from scipy.optimize import minimize, brute, basinhopping, differential_evolution
from numpy import array
from math import sqrt, acos, atan, pi, sin

from core.cluster import Cluster
from core.utils import divide_square, g_ratio
from core.amplitude import *


class Facility:
    """Класс для представления установки НЕВОД-ШАЛ"""

    def __init__(self, geometry='nevod'):

        if geometry == 'nevod':
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
        elif geometry == 'flat':
            self.clusters = [
                Cluster([-40.0, 40.0, 0.0], 20.0, 20.0),
                Cluster([0.0, 40.0, 0.0], 20.0, 20.0),
                Cluster([40.0, 40.0, 0.0], 20.0, 20.0),
                Cluster([-40.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([0.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([40.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([-40.0, -40.0, 0.0], 20.0, 20.0),
                Cluster([0.0, -40.0, 0.0], 20.0, 20.0),
                Cluster([40.0, -40.0, 0.0], 20.0, 20.0),
            ]

        self.grid_steps = 4  # Число шагов по сетке
        self.power_age_steps = 4  # Число шагов поиска мощности и возраста

        self.average_n = None  # Средний из восстановленных векторов
        self.rec_power = None  # Восстановленная мощность
        self.rec_age = None  # Восстановленный возраст
        self.rec_x = None  # Восстановленне координаты прихода ШАЛ
        self.rec_y = None
        self.rec_theta = None
        self.rec_phi = None
        self.clust_ok = None  # Число сработавших кластеров
        self.exp_n = []  # Экспериментальное число частиц
        self.sigma_n = []  # Сигма в функционале
        self.average_x0 = None  # Средневзвешенные x0 и y0
        self.average_y0 = None

        self.real_n = None  # Настоящий вектор ШАЛ

    def get_eas(self, eas):
        """Все кластеры получают ШАЛ"""
        self.real_n = eas.n
        for cluster in self.clusters:
            cluster.get_eas(eas)

    def reset(self):
        """Сбрасываем все кластеры"""
        for cluster in self.clusters:
            cluster.reset()
        self.average_n = None
        self.exp_n = []
        self.sigma_n = []
        self.clust_ok = None
        self.average_x0 = None
        self.average_y0 = None
        self.rec_age = None
        self.rec_power = None
        self.rec_theta = None
        self.rec_phi = None

    def start(self):
        """Запуск всех кластеров"""
        self.clust_ok = 0
        self.average_n = [0, 0, 0]
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

                # Получаем экспериментальное число частиц в станции
                st.particles = st.amplitude / fixed_av_ampl

                self.exp_n.append(st.particles)
                self.sigma_n.append(st.sigma_particles)

                particles_sum += st.particles
                self.average_x0 += st.coord[0] * st.particles
                self.average_y0 += st.coord[1] * st.particles
        self.average_x0 /= particles_sum
        self.average_y0 /= particles_sum
        return True

    def rec_direction(self):
        """Восстановление направления ШАЛ"""
        cl_ok = 0
        self.average_n = [0, 0, 0]
        for cl in self.clusters:
            if cl.respond:
                if cl.rec_direction():
                    self.average_n += cl.rec_n
                    # self.average_n[0] += cl.rec_n[0]
                    # self.average_n[1] += cl.rec_n[1]
                    # self.average_n[2] += cl.rec_n[2]
                    cl_ok += 1

        self.average_n = array(self.average_n)
        if cl_ok == 0:
            print("ERROR: Не воссталовилось направление")
            return False
        else:
            self.average_n /= cl_ok

            self.rec_theta = acos(self.average_n[2])
            # self.rec_phi = atan(self.average_n[1]/self.average_n[0]) * (180/pi)
            self.rec_phi = acos(self.average_n[0] / sin(self.rec_theta))
            # Переведём в градусы
            self.rec_theta *= (180 / pi)
            self.rec_phi *= (180 / pi)
            return True

    def rec_particles(self):
        """Восстановление числа частиц в станциях"""
        if self.average_n is None or self.average_n[2] == 0:
            print("ERROR: Не восстановлились частицы")
            return False

        fixed_av_ampl = get_av_amplitude() / self.average_n[2]

        particles_sum = 0
        self.average_x0 = 0
        self.average_y0 = 0
        for cl in self.clusters:
            for st in cl.stations:
                # Получаем экспериментальное число частиц в станции
                st.particles = st.amplitude / fixed_av_ampl

                self.exp_n.append(st.particles)
                st.sigma_particles = sqrt(st.particles * get_sqr_sigma())
                if st.sigma_particles == 0:
                    st.sigma_particles = 1.3
                self.sigma_n.append(st.sigma_particles)

                particles_sum += st.particles
                self.average_x0 += st.coord[0] * st.particles
                self.average_y0 += st.coord[1] * st.particles

        # self.exp_n = tuple(self.exp_n)
        # self.sigma_n = tuple(self.sigma_n)
        self.average_x0 /= particles_sum
        self.average_y0 /= particles_sum

        return True

    def rec_params_nelder_mead(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])

        res = minimize(self.func, _args, method='Nelder-Mead',
                       options={'maxiter': 1e6, 'maxfev': 1e6, 'disp': True})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params_powell(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])
        _bnds = ((-50, 50), (-50, 50), (10**5, 10**9), (0.7, 2.0))

        res = minimize(self.func, _args, method='Powell',
                       options={'maxiter': 1e6, 'maxfev': 1e6, 'disp': True,
                                'xtol': 1e-06, 'ftol': 1e-06})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            # self.check_bnds(_bnds)
            return True
        else:
            return False

    def rec_params_bashinhopping(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])
        _bnds = ((-50, 50), (-50, 50), (10**5, 10**9), (0.7, 2.0))
        minimizer_kwargs = {"method": "Powell"}

        res = basinhopping(self.func, _args, disp=True, niter=200,
                           minimizer_kwargs=minimizer_kwargs)

        self.rec_x = res.x[0]
        self.rec_y = res.x[1]
        self.rec_power = res.x[2]
        self.rec_age = res.x[3]
        return True

    def rec_params_diff_evo(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])
        _bnds = ((-50, 50), (-50, 50), (10**5, 10**9), (0.7, 2.0))

        res = differential_evolution(self.func, _bnds, maxiter=1000, disp=False)

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params_bfgs(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])

        res = minimize(self.func, _args, method='BFGS',
                       options={'maxiter': None, 'disp': False, 'gtol': 1e-09})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params_lbfgsb(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])
        _bnds = ((-50, 50), (-50, 50), (10**5, 10**9), (0.7, 2.0))

        res = minimize(self.func, _args, method='L-BFGS-B', bounds=_bnds,
                       options={'maxiter': 1e8, 'disp': False, 'gtol': 1e-09,
                                'maxls': 100, 'ftol': 1e-11, 'maxcor': 20,
                                'maxfun': 20000, 'eps': 1e-10})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params_slsqp(self):

        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.3
        _args = array([_x, _y, _power, _age])
        _bnds = ((-50, 50), (-50, 50), (10**5, 10**9), (0.7, 2.0))

        res = minimize(self.func, _args, method='SLSQP', bounds=_bnds,
                       options={'maxiter': 1e7, 'disp': False, 'ftol': 1e-10})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params(self):
        """Восстанавливаем точку прихода, мощность и возраст"""
        _x = 0
        _y = 0
        _power = 10 ** 5
        _age = 1.3
        theo_n = list(self.count_theo(_x, _y, _power, _age))

        _func = self.functional(theo_n)
        _side = 100

        for steps in range(self.grid_steps):
            step = self.make_step(_side, _x, _y, _power, _age, _func)
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

    def make_step(self, side, start_x, start_y, start_power, start_age, min_func):
        """Делаем шаг по сетке"""
        step_cen = divide_square(start_x, start_y, side)
        step = {'x': [], 'y': [], 'func': [], 'power': [], 'age': []}

        for x in step_cen['x']:
            for y in step_cen['y']:
                theo_n = list(self.count_theo(x, y, start_power, start_age))
                a = self.power_age_search(x, y, self.functional(theo_n))

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

    def power_age_search(self, x, y, min_func):
        """Здесь варьируем мощность и возраст для данной точки"""
        power = 10 ** 5  # Исходное значение мощности
        age = 1.3  # Исходное значение возраста

        for i in range(self.power_age_steps):
            s_pow = self.rec_power_age(x, y, power, age, min_func, 'power')
            min_func = s_pow['func']
            power = s_pow['power']

            s_age = self.rec_power_age(x, y, power, age, min_func, 'age')
            min_func = s_age['func']
            age = s_age['age']

        return {'func': min_func, 'power': power, 'age': age}

    def rec_power_age(self, x, y, power_0, age_0, min_func, flag):
        """Восстановление мощности или возраста"""

        if flag == "power":
            # Восстанавливаем мощность
            step_pow = 5000
            accuracy_pow = 100
            step_age = 0
            accuracy_age = -1
        elif flag == "age":
            # Восстанавливаем возраст
            step_age = 0.5
            accuracy_age = 0.0001
            step_pow = 0
            accuracy_pow = -1
        else:
            print("ERROR: Неверный флаг при восстановлении мощност/возраста")
            return False

        # Функционал при меньшем значении восстанавливаемого параметра
        l_func = self.functional(
                       list(self.count_theo(x, y, power_0 - step_pow, age_0 - step_age)))

        # Функционал при большем значении восстанавливаемого параметра
        r_func = self.functional(
                       list(self.count_theo(x, y, power_0 + step_pow, age_0 + step_age)))

        if r_func > l_func:
            # Нужно уменьшать параметр
            step_pow, step_age = -step_pow, -step_age

        age = age_0 + step_age
        if age >= 2.0 or age <= 0.7:
            age -= step_age

        power = power_0 + step_pow
        if power < 10 ** 5:
            power -= step_pow

        func = self.functional(list(self.count_theo(x, y, power, age)))

        while abs(step_pow) > accuracy_pow and abs(step_age) > accuracy_age:

            while func < min_func:

                min_func = func

                power += step_pow
                age += step_age

                if power < 10 ** 5 or power > 10 ** 8:
                    power -= step_pow
                    break

                if age >= 2.0 or age <= 0.7:
                    age -= step_age
                    break

                func = self.functional(list(self.count_theo(x, y, power, age)))

            # Разворачиваемся и уменьшаем шаг
            step_pow /= -g_ratio
            step_age /= -g_ratio
            # Шагнули обратно с меньшим шагом
            power += step_pow
            age += step_age

            func = self.functional(list(self.count_theo(x, y, power, age)))

        return {'func': func, 'power': power, 'age': age}

    def count_theo(self, x, y, power, age):
        """Подсчёт теоретическиого числа частиц для каждой станции, возвращает
        генератор"""
        for cluster in self.clusters:
            cluster.rec_particles(self.average_n, x, y, power, age)
            for station in cluster.stations:
                yield station.rec_particles

    def functional(self, theo_n):
        """Считает функционал"""
        f = 0

        if len(self.sigma_n) == len(theo_n) == len(self.exp_n):
            for n_e, n_t, sigma in zip(self.exp_n, theo_n, self.sigma_n):
                f += ((n_e - n_t) ** 2) / (sigma ** 2)
        else:
            print("ERROR: Не совпадает число параметров в функционале")
            return False
        return f

    def func(self, params):
        """Функционал одной функцией от параметров ШАЛ в виде ndarray"""
        x = params[0]
        y = params[1]
        power = params[2]
        age = params[3]

        theo_n = list(self.count_theo(x, y, power, age))

        return self.functional(theo_n)

    def draw_func_power(self, x, y, age):
        """Функция для получения зависимость функционала от мощности в данной точке при
        данном возрасте"""
        power = 10**4  # Минимальное значение мощности
        with open('data/power_age_func/func_power.txt', 'w') as file:
            for i in range(10000):
                power += 1000
                func = self.functional(list(self.count_theo(x, y, power, age)))

                file.write(str(power) + '\t' + str(func) + '\n')

    def draw_func_age(self, x, y, power):
        """Функция для получения зависимость функционала от возраста в данной точке при
        данной мощности"""
        age = 1.2
        with open('data/power_age_func/func_age.txt', 'w') as file:
            for i in range(10000):
                age += 0.001
                func = self.functional(list(self.count_theo(x, y, power, age)))

                file.write(str(age) + '\t' + str(func) + '\n')

    def check_bnds(self, bnds):
        """Проверяет восстановленные параметры на попадание в диапазон границ.
            Иначе присваивает им граничные значения"""
        # x0
        if self.rec_x < bnds[0][0]:
            self.rec_x = bnds[0][0]

        if self.rec_x > bnds[0][1]:
            self.rec_x = bnds[0][1]

        # y0
        if self.rec_y < bnds[1][0]:
            self.rec_x = bnds[1][0]

        if self.rec_x > bnds[1][1]:
            self.rec_x = bnds[1][1]

        # Мощность
        if self.rec_power < bnds[2][0]:
            self.rec_x = bnds[2][0]

        if self.rec_power > bnds[2][1]:
            self.rec_x = bnds[2][1]

        # Возраст
        if self.rec_age < bnds[3][0]:
            self.rec_x = bnds[3][0]

        if self.rec_age > bnds[3][1]:
            self.rec_x = bnds[3][1]
