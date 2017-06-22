from scipy.optimize import minimize, basinhopping, differential_evolution
from numpy import array, arctan2, degrees
from numpy.linalg import det, norm
from math import sqrt, acos

from core.cluster import Cluster
from core.utils import divide_square, g_ratio
from core.amplitude import *

light_speed = 0.299792458  # Скорость света [м/нс]


class Facility:
    """Класс для представления установки НЕВОД-ШАЛ"""

    def __init__(self, geometry='nevod', num=None):

        if geometry == 'nevod':
            self.clusters = [
                Cluster([-28.4, -7.8, -7.0], 13.3, 12.4),
                Cluster([-28.4, 23.8, -7.0], 13.3, 12.4),
                Cluster([0.0, 0.0, 0.0], 25.1, 13.5),
                Cluster([33.3, 7.8, -14.5]),
                Cluster([35.0, 46.0, -14.5]),
                Cluster([-2.0, -47.0, -14.5]),
                Cluster([-18.0, 62.0, -14.5]),
                Cluster([60.0, -2.0, -2.0]),
                Cluster([60.0, 26.0, -2.0]),
                Cluster([60.0, 58.0, -8.0]),
            ]
        elif geometry == 'flat':
            self.clusters = [
                Cluster([-60.0, 60.0, 0.0], 20.0, 20.0),
                Cluster([0.0, 60.0, 0.0], 20.0, 20.0),
                Cluster([60.0, 60.0, 0.0], 20.0, 20.0),
                Cluster([-60.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([0.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([60.0, 0.0, 0.0], 20.0, 20.0),
                Cluster([-60.0, -60.0, 0.0], 20.0, 20.0),
                Cluster([0.0, -60.0, 0.0], 20.0, 20.0),
                Cluster([60.0, -60.0, 0.0], 20.0, 20.0),
            ]

        self.num = num  # Номер установки, если создали много штук
        self.eas = None  # ШАЛ, упавший на данную установку
        self.clust_ok = None  # Число сработавших кластеров

        self.average_n = None  # Средний из восстановленных векторов
        self.rec_n = None  # Восстановленный вектор по установке

        self.psi = None  # Угол отклонения среднего
        self.sigma_psi = None  # Среднеквадратичный угол отклонения

        self.rec_power = None  # Восстановленная мощность
        self.rec_age = None  # Восстановленный возраст
        self.rec_x = None  # Восстановленне координаты прихода ШАЛ
        self.rec_y = None
        self.rec_theta = None  # Восстановленные тета и фи
        self.rec_phi = None

        self.exp_n = []  # Экспериментальное число частиц
        self.sigma_n = []  # Сигма в функционале

        self.average_x0 = None  # Средневзвешенные x0 и y0
        self.average_y0 = None

        self.grid_steps = 5  # Число шагов по сетке
        self.power_age_steps = 5  # Число шагов поиска мощности и возраста

    def reset(self):
        """Сбрасываем все кластеры"""
        for cluster in self.clusters:
            cluster.reset()

        self.eas = None
        self.clust_ok = None

        self.average_n = None
        self.rec_n = None

        self.psi = None
        self.sigma_psi = None

        self.rec_power = None
        self.rec_age = None
        self.rec_x = None
        self.rec_y = None
        self.rec_theta = None
        self.rec_phi = None

        self.exp_n = []
        self.sigma_n = []

        self.average_x0 = None
        self.average_y0 = None

    def get_eas(self, eas):
        """Получить данные ШАЛ без запуска"""
        self.eas = eas
        for cluster in self.clusters:
            cluster.get_eas(eas)

    def start(self, eas):
        """Запуск установки"""
        self.eas = eas
        self.clust_ok = 0
        for cluster in self.clusters:
            if cluster.start(eas):
                self.clust_ok += 1

        if self.clust_ok == 0:
            return False
        else:
            return True

    def set_facility_state(self, evt, eas=None):
        """Устанавить состояние установки в соответствии с
        прочитанным событием"""
        self.eas = eas
        self.clust_ok = 0
        for cl_n, cl in enumerate(self.clusters):
            if cl.set_cluster_state(evt['clusters'][cl_n]):
                self.clust_ok += 1

        if self.clust_ok <= 0:
            return False
        else:
            return True

    def rec_direction(self):
        """Восстановление направления ШАЛ"""
        cl_ok = 0  # Считаем кластеры, которые смогли восстановить направление
        self.average_n = array([0.0, 0.0, 0.0])
        average_sqr_n = 0.0  # Средний квадрат
        for cl in self.clusters:
            if cl.respond:
                if cl.rec_direction():
                    self.average_n += cl.rec_n
                    average_sqr_n += norm(cl.rec_n**2, ord=1)
                    cl_ok += 1

        if cl_ok == 0:
            print("ERROR: Не воссталовилось направление")
            return False
        else:
            self.average_n /= cl_ok
            average_sqr_n /= cl_ok

            # Расчитаем зенитный и азимутальный углы
            self.rec_theta = acos(abs(self.average_n[2]))
            self.rec_phi = arctan2(self.average_n[1], self.average_n[0])
            if self.rec_phi < 0:
                self.rec_phi += 2 * pi
            # Переведём в градусы
            self.rec_theta = degrees(self.rec_theta)
            self.rec_phi = degrees(self.rec_phi)

            # Вычислим угол отклонения среднего:
            delta = norm(self.eas.n - self.average_n, ord=2)
            self.psi = degrees(acos(1 - delta**2 / 2))

            # Вычислим среднеквадратичный угол отклонения
            sigma_psi = sqrt(average_sqr_n - norm(self.average_n**2, ord=1))
            self.sigma_psi = degrees(acos(1 - sigma_psi**2 / 2))

            return True

    def rec_particles(self):
        """Восстановление числа частиц в станциях, заполняем экспериментальное
        число частиц для функционала"""
        if self.average_n is None or self.average_n[2] == 0:
            print("ERROR: Не восстановлились частицы")
            return False

        fixed_av_ampl = get_av_amplitude() / abs(self.average_n[2])

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

        self.average_x0 /= particles_sum
        self.average_y0 /= particles_sum

        return True

    def make_times_relative(self):
        """Сделаем времена относительными по установке"""

        # Записали сюда все времена сработавших станций по установке
        min_t = min([st.rndm_time for cl in self.clusters for st in cl.stations
                     if st.respond])

        for cl in self.clusters:
            for st in cl.stations:
                if st.respond:
                    st.rndm_time -= min_t

        return True

    def rec_direction_new(self):
        """Восстанавливает вектор прихода ШАЛ методом наименьших квадратов"""

        # Изменим времена срабатывания станций на относительные
        self.make_times_relative()

        sum0 = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_yy = 0
        sum_xy = 0
        sum_t = 0
        sum_tx = 0
        sum_ty = 0

        for cl in self.clusters:
            for st in cl.stations:
                if st.respond:
                    sqr_sigma_t = pow(st.sigma_t, 2)

                    sum0 += 1 / sqr_sigma_t
                    sum_t += st.rndm_time / sqr_sigma_t
                    sum_x += st.coord[0] / sqr_sigma_t
                    sum_y += st.coord[1] / sqr_sigma_t
                    sum_xx += pow(st.coord[0], 2) / sqr_sigma_t
                    sum_yy += pow(st.coord[1], 2) / sqr_sigma_t
                    sum_xy += (st.coord[0] * st.coord[1]) / sqr_sigma_t
                    sum_tx += (st.rndm_time * st.coord[0]) / sqr_sigma_t
                    sum_ty += (st.rndm_time * st.coord[1]) / sqr_sigma_t

        sqr_light_speed = pow(light_speed, 2)

        sum0 /= sqr_light_speed
        sum_xx /= sqr_light_speed
        sum_xy /= sqr_light_speed
        sum_yy /= sqr_light_speed
        sum_x /= sqr_light_speed
        sum_y /= sqr_light_speed
        sum_tx /= light_speed
        sum_ty /= light_speed
        sum_t /= light_speed

        line_1 = [sum_xx, sum_xy, sum_x]
        line_2 = [sum_xy, sum_yy, sum_y]
        line_3 = [sum_x, sum_y, sum0]
        line_4 = [sum_tx, sum_ty, sum_t]

        det1 = det([line_1, line_2, line_3])
        det2 = det([line_4, line_2, line_3])
        det3 = det([line_1, line_4, line_3])
        # det4 = det([line_1, line_2, line_4])

        a = det2 / det1
        b = det3 / det1

        if (pow(a, 2) + pow(b, 2)) <= 1:
            # Если вектор восстановился успешно
            c = sqrt(1 - pow(a, 2) - pow(b, 2))
            self.rec_n = array([a, b, c])

            # Расчитаем зенитный и азимутальный углы
            self.rec_theta = acos(abs(self.rec_n[2]))
            self.rec_phi = arctan2(self.rec_n[1], self.rec_n[0])
            if self.rec_phi < 0:
                self.rec_phi += 2 * pi
            # Переведём в градусы
            self.rec_theta = degrees(self.rec_theta)
            self.rec_phi = degrees(self.rec_phi)

            # Вычислим угол отклонения среднего:
            delta = norm(self.eas.n - self.rec_n, ord=2)
            self.psi = degrees(acos(1 - delta**2 / 2))

            return True
        else:
            print("ERROR: Не удалось восстановить направление")
            return False

    def rec_params_powell(self):
        """Восстановить параметры ШАЛ методом Пауэлла"""
        _x = self.average_x0
        _y = self.average_y0
        _power = 10 ** 5
        _age = 1.5
        _args = array([_x, _y, _power, _age])

        res = minimize(self.func, _args, method='Powell',
                       options={'maxiter': 1e6, 'maxfev': 1e6, 'disp': True,
                                'xtol': 1e-06, 'ftol': 1e-06})

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return True
        else:
            return False

    def rec_params_diff_evo(self):
        """Восстановление параметров ШАЛ методом дифференциальной эволюции"""
        _bnds = ((-80, 120), (-100, 120), (10**5, 10**8), (0.5, 2.0))

        res = differential_evolution(self.func, _bnds, maxiter=2000, disp=False,
                                     polish=True)

        if res.success:
            self.rec_x = res.x[0]
            self.rec_y = res.x[1]
            self.rec_power = res.x[2]
            self.rec_age = res.x[3]
            return [self.num,
                    self.rec_theta, self.rec_phi,
                    self.rec_x, self.rec_y,
                    self.rec_power, self.rec_age]
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

        # Начальные значения параметров
        x = 0
        y = 0
        power = 10 ** 5
        age = 1.3
        params = [x, y, power, age]

        func = self.func(params)
        side = 100  # Длина стороны квадрата, охватывающего установку

        for steps in range(self.grid_steps):
            step = self.make_step(side, x, y, power, age, func)
            x = step['x']
            y = step['y']
            power = step['power']
            age = step['age']
            func = step['func']
            side /= 3

        self.rec_x = x
        self.rec_y = y
        self.rec_power = power
        self.rec_age = age

        return True

    def make_step(self, side, start_x, start_y, start_power, start_age, min_func):
        """Делаем шаг по сетке"""
        step_cen = divide_square(start_x, start_y, side)
        step = {'x': [], 'y': [], 'func': [], 'power': [], 'age': []}

        for x in step_cen['x']:
            for y in step_cen['y']:

                # Варьируем мощность и возраст для каждой точки
                a = self.power_age_search(x, y, self.func([x, y, start_power, start_age]))

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
        """Общая функция для варьирования мощности или возраста 
        для каждой точки"""

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
        l_func = self.func([x, y, power_0 - step_pow, age_0 - step_age])

        # Функционал при большем значении восстанавливаемого параметра
        r_func = self.func([x, y, power_0 + step_pow, age_0 + step_age])

        if r_func > l_func:
            # Нужно уменьшать параметр
            step_pow, step_age = -step_pow, -step_age

        age = age_0 + step_age
        if age > 2.0 or age < 0.7:
            age -= step_age

        power = power_0 + step_pow
        if power < 10 ** 5 or power > 10**8:
            power -= step_pow

        func = self.func([x, y, power, age])

        while abs(step_pow) > accuracy_pow and abs(step_age) > accuracy_age:

            while func < min_func:

                min_func = func

                power += step_pow
                age += step_age

                if power < 10 ** 5 or power > 10 ** 8:
                    power -= step_pow
                    break

                if age > 2.0 or age < 0.7:
                    age -= step_age
                    break

                func = self.func([x, y, power, age])

            # Разворачиваемся и уменьшаем шаг
            step_pow /= -g_ratio
            step_age /= -g_ratio
            # Шагнули обратно с меньшим шагом
            power += step_pow
            age += step_age

            func = self.func([x, y, power, age])

        return {'func': func, 'power': power, 'age': age}

    def count_theo(self, params):
        """Подсчёт теоретическиого числа частиц для каждой станции, возвращает
        генератор"""
        for cluster in self.clusters:
            cluster.rec_particles(self.average_n, params)
            for station in cluster.stations:
                yield station.rec_particles

    def func(self, params):
        """Функционал одной функцией от параметров ШАЛ в виде ndarray или list"""
        theo_n = list(self.count_theo(params))

        f = 0
        if len(self.sigma_n) == len(theo_n) == len(self.exp_n):
            for n_e, n_t, sigma in zip(self.exp_n, theo_n, self.sigma_n):
                f += ((n_e - n_t) ** 2) / (sigma ** 2)
        else:
            print("ERROR: Не совпадает число параметров в функционале")
            return False
        return f

    def draw_func_power(self, x, y, age):
        """Функция для получения зависимость функционала от мощности в данной точке при
        данном возрасте"""
        power = 10**4  # Минимальное значение мощности
        with open('data/power_age_func/func_power.txt', 'w') as file:
            for i in range(10000):
                power += 1000
                func = self.func([x, y, power, age])

                file.write(str(power) + '\t' + str(func) + '\n')

    def draw_func_age(self, x, y, power):
        """Функция для получения зависимость функционала от возраста в данной точке при
        данной мощности"""
        age = 1.2
        with open('data/power_age_func/func_age.txt', 'w') as file:
            for i in range(10000):
                age += 0.001
                func = self.func([x, y, power, age])

                file.write(str(age) + '\t' + str(func) + '\n')

