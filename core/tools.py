from math import pi, pow, sqrt, exp, cos, sin
from collections import namedtuple
from numpy import array, cross, size
from numpy.linalg import norm
import random as rn

g_ratio = (1 + 5 ** 0.5) / 2  # Золотое сечение

time_dist_func = []
tau = 5 / sqrt(3)
# with open('data/time/dist_func.txt', 'w') as dist_func_file:

for i in range(300):
    # Заполняем функцию распределения
    t = 0.5 * i

    new_element = (2 * pow(tau, 3) - tau * exp(-t / tau) * (pow(t, 2) + 2 * t *
                            tau + 2 * pow(tau, 2))) * (1 / (2 * pow(tau, 3)))

    time_dist_func.append(new_element)

# dist_func_file.write(str(t) + '\t' + str(new_element) + '\n')


def randomize_time():
    """Возвращает случайную поправку ко времени срабатывания станции"""
    gm = rn.random()
    for j in range(300):
        if gm < time_dist_func[j]:
            return j * 0.5
        if j == 299:
            return j * 0.5


def psn(x):
    """Генератор распределения Пуассона"""
    p = exp(-x)
    imax = int(x + 5.0*sqrt(x))
    if imax < 5:
        imax = 5
    gm = rn.random()
    sump = p
    if gm < sump:
        return 0
    for j in range(1, imax + 1):
        p = p * x/j
        sump += p
        if gm < sump:
            return i
    return imax


def get_distance(st_coord, n, x0, y0):
    """Вычисляет расстояние от станции до оси ШАЛ"""
    # Вектор от станции до точки прихода ШАЛ
    a = array([st_coord[0] - x0, st_coord[1] - y0, st_coord[2]])
    b = cross(a, n)  # Векторное произведение вектора a на вектор ШАЛ
    # Возвращаем модуль векторного произведения, т.к длина n равна единице
    return norm(b, ord=2)


def get_theta():
    """Методом Неймана получаем случайный зенитный угол"""
    while True:
        kx1 = rn.random()
        kx2 = rn.random()
        kx11 = 0.0 + ((pi/2) - 0.0) * kx1
        kx22 = 1.0 * kx2

        if kx22 <= pow(cos(kx11), 8.5) * sin(kx11):
            return kx11 * (180/pi)


def get_power():
    """Получаем мощность методом обратных функций"""
    gm = rn.random()
    return ((10**6) / (1 - gm))**(2/3)


def functional(exp_n, sigma_n, theo_n):
    """Считает функционал"""
    f = 0

    if len(sigma_n) == len(theo_n) == len(exp_n):
        for n_e, n_t, sigma in zip(exp_n, theo_n, sigma_n):
            f += ((n_e - n_t)**2) / (sigma**2)
    else:
        print("ERROR: Не совпадает число параметров в функционале")
        return False
    return f


def divide_square(center_x, center_y, side):
    """Делит квадрат на девять равных и возвращает их центры"""
    result_x = []
    result_y = []
    new_side = side / 3
    left_corner_x = center_x - side / 2
    left_corner_y = center_y - side / 2
    for k in [1, 3, 5]:
        # Заполнили первый ряд
        result_y.append(left_corner_y + k * (new_side/2))
        result_x.append(left_corner_x + k * (new_side/2))

    return [result_x, result_y]


def count_theo(clusters, average_n, x, y, energy, age):
    """Подсчёт теоретическиого числа частиц для каждой станции"""
    theo_n = []

    for cluster in clusters:
        cluster.rec_particles(average_n, x, y, energy, age)
        for station in cluster.stations:
            theo_n.append(station.rec_particles)

    return theo_n


def draw_func_power(clusters, average_n, x, y, power, age, exp_n, sigma_n):
    """Функция для получения зависимость функционала от мощности"""
    with open('data/power_age_func/func_power.txt', 'w') as file:
        for i in range(10000):
            power += 1000
            func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x,
                                                         y, power, age))

            file.write(str(power) + '\t' + str(func) + '\n')


def draw_func_age(clusters, average_n, x, y, power, age, exp_n, sigma_n):
    """Функция для получения зависимость функционала от мощности"""
    with open('data/power_age_func/func_age.txt', 'w') as file:
        for i in range(10000):
            age += 0.00007
            func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x,
                                                         y, power, age))

            file.write(str(age) + '\t' + str(func) + '\n')


def search_power(clusters, average_n, x, y, power_0, age_0, exp_n, sigma_n,
                 min_func):
    """Поиск мощности"""
    step1 = 5000
    step2 = 2000
    step3 = 500

    # Функционал при меньшем значении энергии
    l_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   power_0 - step1, age_0))
    # Фкнкционал при большем значении энергии
    r_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   power_0 + step1, age_0))

    if r_func > l_func:
        # Если истинное значение энергии меньше текущего - будем шагать назад
        step1, step2, step3 = -step1, - step2, -step3

    power = power_0 + step1
    if power < 10**4:
        power -= step1

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 power, age_0))

    def mk_step(step):
        nonlocal power
        nonlocal func
        nonlocal min_func

        while func <= min_func:
            min_func = func
            power += step
            if power < 10 ** 4 or power > 10 ** 8:
                break
            func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x,
                                                         y, power, age_0))
        power -= step
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     power, age_0))

    mk_step(step1)
    mk_step(step2)
    mk_step(step3)

    return {'func': func, 'power': power}


def search_age(clusters, average_n, x, y, power_0, age_0, exp_n, sigma_n,
               min_func):
    """Поиск возраста"""
    step1 = 0.1
    step2 = 0.05
    step3 = 0.01

    # Функционал при меньшем возраста возраста
    l_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   power_0, age_0 - step1))
    # Фкнкционал при большем значении возраста
    r_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   power_0, age_0 + step1))

    if r_func > l_func:
        # Истинное значение энергии меньше текущего - будем шагать назад
        step1, step2, step3 = -step1, -step2, -step3

    age = age_0 + step1
    if age >= 2.0 or age <= 0.5:
        age -= step1

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 power_0, age))

    def mk_step(step):
        nonlocal age
        nonlocal func
        nonlocal min_func

        while func <= min_func:
            min_func = func
            age += step
            if age >= 2.0 or age <= 0.5:
                break
            func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x,
                                                         y, power_0, age))
        age -= step
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     power_0, age))

    mk_step(step1)
    mk_step(step2)
    mk_step(step3)

    return {'func': func, 'age': age}


def energy_age_search(clusters, average_n, x, y, exp_n, sigma_n, min_func):
    """Здесь варьируем мощность и возраст для данной точки"""
    power = 10**4  # Исходное значение мощности
    age = 1.3  # Исходное значение возраста

    for i in range(3):

        s_pow = search_power(clusters, average_n, x, y, power, age, exp_n,
                             sigma_n, min_func)
        min_func = s_pow['func']
        power = s_pow['power']

        s_age = search_age(clusters, average_n, x, y, power, age, exp_n, sigma_n
                           , min_func)
        min_func = s_age['func']
        age = s_age['age']

    return {'func': min_func, 'power': power, 'age': age}


def make_step(clusters, average_n, side, start_x, start_y, start_power,
              start_age, exp_n, sigma_n, min_func):

    step_cen = divide_square(start_x, start_y, side)
    step = {'x': [], 'y': [], 'func': [], 'power': [], 'age': []}

    for x in step_cen[0]:
        for y in step_cen[1]:
            theo_n = count_theo(clusters, average_n, x, y, start_power,
                                start_age)
            a = energy_age_search(clusters, average_n, x, y, exp_n, sigma_n,
                                  functional(exp_n, sigma_n, theo_n))

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
