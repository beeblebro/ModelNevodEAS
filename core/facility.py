from core.cluster import Cluster
from core.utils import functional, divide_square, g_ratio
from core.amplitude import *


class Facility:
    """Класс для представления установки НЕВОД-ШАЛ"""

    def __init__(self):

        self.clusters = [
            Cluster([-28.4, -7.8, -7.0], 13.3, 12.4),
            Cluster([-28.4, 23.8, -7.0], 13.3, 12.4),
            Cluster([0, 0, 0], 25.1, 13.5),
            Cluster([33.3, 7.8, -14.5]),
            Cluster([35, 46, -14.5]),
            Cluster([-2, -47, -14.5]),
            Cluster([-18, 62, -14.5]),
            Cluster([50, -2, -2]),
            Cluster([50, 26, -2]),
            Cluster([50, 58, -8]),
        ]

        self.average_n = [0, 0, 0]  # Средний из восстановленных векторов
        self.rec_power = None
        self.rec_age = None
        self.rec_x = None
        self.rec_y = None
        self.clust_ok = 0  # Число сработавших кластеров
        self.exp_n = []  # Экспериментальное число частиц
        self.sigma_n = []  # Сигма в функционале
        self.average_x0 = None  # Средневзвешенные x0 и y0
        self.average_y0 = None

    def get_eas(self, eas):
        """Все кластеры получают ШАЛ"""
        for cluster in self.clusters:
            cluster.get_eas(eas)

    def reset(self):
        """Сбрасываем все кластеры"""
        for cluster in self.clusters:
            cluster.reset()
        self.average_n = [0, 0, 0]
        self.clust_ok = 0
        self.exp_n = []
        self.sigma_n = []
        self.average_x0 = None
        self.average_y0 = None
        self.rec_age = None
        self.rec_power = None

    def start(self):
        """Запуск всех кластеров"""
        for cluster in self.clusters:
            if cluster.start():
                self.average_n += cluster.rec_n
                self.clust_ok += 1
        if self.clust_ok == 0:
            print("ERROR: Не сработал ни один кластер")
            return False
        # Получили средний из восстановленных векторов
        self.average_n /= self.clust_ok
        # Средняя амплитуда с поправкой на косинус тета
        fixed_av_ampl = get_av_amplitude() / self.average_n[2]

        particles_sum = 0
        self.average_x0 = 0
        self.average_y0 = 0
        for cl in self.clusters:
            for st in cl.stations:
                st.particles = st.amplitude / fixed_av_ampl

                self.exp_n.append(st.particles)
                self.sigma_n.append(st.sigma_particles)

                particles_sum += st.particles
                self.average_x0 += st.coord[0] * st.particles
                self.average_y0 += st.coord[1] * st.particles

        self.average_x0 /= particles_sum
        self.average_y0 /= particles_sum
        return True

    def rec_params(self):
        """Восстанавливаем точку прихода, мощность и возраст"""
        _x = 0
        _y = 0
        _power = 10 ** 6
        _age = 1.3
        theo_n = self._count_theo(_x, _y, _power, _age)

        _func = functional(self.exp_n, self.sigma_n, theo_n)
        _side = 100

        for steps in range(3):
            step = self._make_step(_side, _x, _y, _power, _age, _func)
            _x = step['x']
            _y = step['y']
            _power = step['power']
            _age = step['age']
            _func = step['func']

            _side /= 3

        self.rec_x = _x
        self.rec_y = _y
        self.rec_power = _power
        self.rec_age = _age

    def _make_step(self, side, start_x, start_y, start_power, start_age, min_func):
        """Делаем шаг по сетке"""
        step_cen = divide_square(start_x, start_y, side)
        step = {'x': [], 'y': [], 'func': [], 'power': [], 'age': []}

        for x in step_cen['x']:
            for y in step_cen['y']:
                theo_n = self._count_theo(x, y, start_power, start_age)
                a = self._power_age_search(x, y, functional(self. exp_n, self.sigma_n,
                                                            theo_n))

                step['x'].append(x)
                step['y'].append(y)

                step['func'].append(a['func'])
                step['power'].append(a['power'])
                step['age'].append(a['age'])

        if min(step['func']) < min_func:

            _index = step['func'].index(min(step['func']))

            return {'x': step['x'][_index],
                    'y': step['y'][_index],
                    'func': step['func'][_index],
                    'power': step['power'][_index],
                    'age': step['age'][_index]}
        else:
            return {'x': start_x,
                    'y': start_y,
                    'func': min_func,
                    'power': start_power,
                    'age': start_age}

    def _power_age_search(self, x, y, min_func):
        """Здесь варьируем мощность и возраст для данной точки"""
        power = 10 ** 6  # Исходное значение мощности
        age = 1.3  # Исходное значение возраста

        for i in range(3):
            s_pow = self._rec_power_age(x, y, power, age, min_func, 'power')
            min_func = s_pow['func']
            power = s_pow['power']

            s_age = self._rec_power_age(x, y, power, age, min_func, 'age')
            min_func = s_age['func']
            age = s_age['age']

        return {'func': min_func, 'power': power, 'age': age}

    def _rec_power_age(self, x, y, power_0, age_0, min_func, flag):
        """Восстановление мощности или возраста"""

        if flag == "power":
            # Восстанавливаем мощность
            step_pow = 5000
            accuracy_pow = 100
            step_age = 0
            accuracy_age = -1
        elif flag == "age":
            # Восстанавливаем возраст
            step_age = 0.1
            accuracy_age = 0.001
            step_pow = 0
            accuracy_pow = -1
        else:
            print("ERROR: Неверный флаг при восстановлении мощност/возраста")
            return False

        # Функционал при меньшем значении восстанавливаемого параметра
        l_func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y,
                                                            power_0 - step_pow,
                                                            age_0 - step_age))

        # Функционал при большем значении восстанавливаемого параметра
        r_func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y,
                                                            power_0 + step_pow,
                                                            age_0 + step_age))

        if r_func > l_func:
            # Нужно уменьшать параметр
            step_pow, step_age = -step_pow, -step_age

        age = age_0 + step_age
        if age >= 2.0 or age <= 0.5:
            age -= step_age

        power = power_0 + step_pow
        if power < 10 ** 4:
            power -= step_pow

        func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y, power, age))

        while abs(step_pow) > accuracy_pow and abs(step_age) > accuracy_age:

            while func <= min_func:

                min_func = func

                power += step_pow
                age += step_age
                if (power < 10 ** 4 or power > 10 ** 8) or (age >= 2.0 or age <= 0.5):
                    break

                func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y,
                                                                             power, age))

            step_pow /= -g_ratio
            step_age /= -g_ratio

        return {'func': func, 'power': power, 'age': age}

    def _count_theo(self, x, y, power, age):
        """Подсчёт теоретическиого числа частиц для каждой станции"""
        theo_n = []

        for cluster in self.clusters:
            cluster.rec_particles(self.average_n, x, y, power, age)
            for station in cluster.stations:
                theo_n.append(station.rec_particles)

        return theo_n

    def draw_func_power(self, x, y, age):
        """Функция для получения зависимость функционала от мощности в данной точке при
        данном возрасте"""
        power = 10**4  # Минимальное значение мощности
        with open('data/power_age_func/func_power.txt', 'w') as file:
            for i in range(10000):
                power += 1000
                func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y,
                                                                             power,  age))

                file.write(str(power) + '\t' + str(func) + '\n')

    def draw_func_age(self, x, y, power):
        """Функция для получения зависимость функционала от возраста в данной точке при
        данной мощности"""
        age = 1.2
        with open('data/power_age_func/func_age.txt', 'w') as file:
            for i in range(10000):
                age += 0.00007
                func = functional(self.exp_n, self.sigma_n, self._count_theo(x, y, power,
                                                                         age))

                file.write(str(age) + '\t' + str(func) + '\n')
