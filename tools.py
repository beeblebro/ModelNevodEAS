from math import pi, pow, sqrt, exp, cos, sin
from numpy import array, cross, size
from numpy.linalg import norm
import random as rn

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


def get_energy():
    """Получаем энергию методом обратных функций"""
    gm = rn.random()
    return ((10**6) / (1 - gm))**(2/3)


def functional(exp_n, sigma_n, theo_n):
    """Считает функционал"""
    f = 0

    if len(sigma_n) == len(theo_n) == len(exp_n):
        for j in range(0, len(sigma_n)):
            f += ((exp_n[j] - theo_n[j])**2) / pow(sigma_n[j], 2)
    else:
        print("Functional error#0")
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
'''
    for k in [1, 3, 5]:
        # Второй ряд
        result_y.append(left_corner_y + 3 * (new_side/2))
        result_x.append(left_corner_x + k * (new_side/2))

    for k in [1, 3, 5]:
        # Третий ряд
        result_y.append(left_corner_y + 5 * (new_side/2))
        result_x.append(left_corner_x + k * (new_side/2))
'''


def count_theo(clusters, average_n, x, y, energy, age):
    """Подсчёт теоретическиого числа частиц для каждой станции"""
    theo_n = []

    for cluster in clusters:
        cluster.rec_particles(x, y, energy, age, average_n)
        for i in range(4):
            theo_n.append(cluster.stations[i].rec_particles)

    return theo_n


def search_energy(clusters, average_n, x, y, energy_0, age_0, exp_n, sigma_n,
                  min_func):
    """Поиск энергии"""
    # search = {'func': [], 'energy': []}
    # Функционал при меньшем значении энергии
    l_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   energy_0 - 1000, age_0))
    # Фкнкционал при большем значении энергии
    r_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   energy_0 + 1000, age_0))

    step1 = 1000
    step2 = 500
    step3 = 100

    if r_func > l_func:
        # Истинное значение энергии меньше текущего - будем шагать назад
        step1 = -step1
        step2 = -step2
        step3 = -step3

    energy = energy_0 + step1

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy, age_0))
    # Шагаем в нужную сторону, пока функционал не станет больше
    while func < min_func:
        min_func = func
        energy += step1
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     energy, age_0))

    # Сделаем шаг в другую сторону
    energy -= step1
    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy, age_0))
    # Продолжим шагать с меньшим шагом
    while func < min_func:
        min_func = func
        energy += step2
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     energy, age_0))

    energy -= step2
    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy, age_0))
    while func < min_func:
        min_func = func
        energy += step3
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     energy, age_0))
    energy -= step3

    # for k in range(10):
    #     energy = energy_0 + 100000 * k
    #
    #     theo_n = count_theo(clusters, average_n, x, y, energy, age_0)
    #
    #     search['func'].append(functional(exp_n, sigma_n, theo_n))
    #     search['energy'].append(energy)
    #
    # if min(search['func']) < min_func:
    #     min_func = min(search['func'])
    #     min_num = search['func'].index(min(search['func']))
    #     new_energy = search['energy'][min_num]
    # else:
    #     new_energy = energy_0

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy, age_0))

    return {'new_func': func, 'new_energy': energy}


def search_age(clusters, average_n, x, y, energy_0, age_0, exp_n, sigma_n,
               min_func):
    """Поиск возраста"""
    # search = {'func': [], 'age': []}
    # Функционал при меньшем возраста возраста
    l_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   energy_0, age_0 - 0.05))
    # Фкнкционал при большем значении возраста
    r_func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                   energy_0, age_0 + 0.05))
    step1 = 0.05
    step2 = 0.01
    step3 = 0.005

    if r_func > l_func:
        # Истинное значение энергии меньше текущего - будем шагать назад
        step1 = -step1
        step2 = -step2
        step3 = -step3

    age = age_0 + step1

    # Шагаем в нужную сторону, пока функционал не станет больше
    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy_0, age))
    while func < min_func:
        min_func = func
        age += step1
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy_0, age))
    # Сделаем шаг в другую сторону
    age -= step1

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy_0, age))
    # Продолжим шагать с меньшим шагом
    while functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                energy_0, age)) < min_func:
        min_func = func
        age += step2
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     energy_0, age))

    age -= step2

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy_0, age))
    while functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                energy_0, age)) < min_func:
        min_func = func
        age += step3
        func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                     energy_0, age))

    age -= step3

    # for k in range(10):
    #     age = age_0 + 0.05 * k
    #
    #     theo_n = count_theo(clusters, average_n, x, y, energy_0, age)
    #
    #     search['func'].append(functional(exp_n, sigma_n, theo_n))
    #     search['age'].append(age)
    #
    # if min(search['func']) < min_func:
    #     min_func = min(search['func'])
    #     min_num = search['func'].index(min(search['func']))
    #     new_age = search['age'][min_num]
    # else:
    #     new_age = age_0

    func = functional(exp_n, sigma_n, count_theo(clusters, average_n, x, y,
                                                 energy_0, age))

    return {'new_func': func, 'new_age': age}


def energy_age_search(clusters, average_n, x, y, exp_n, sigma_n, min_func):
    """Здесь варьируем энергию и возраст для данной точки"""
    energy_0 = 10**6  # Исходное значение энергии
    age_0 = 1.3  # Исходное значение возраста

    a = search_energy(clusters, average_n, x, y, energy_0, age_0, exp_n,
                      sigma_n, min_func)
    min_func = a['new_func']
    energy = a['new_energy']

    b = search_age(clusters, average_n, x, y, energy, age_0, exp_n, sigma_n,
                   min_func)
    min_func = b['new_func']
    age = b['new_age']

    c = search_energy(clusters, average_n, x, y, energy, age, exp_n, sigma_n,
                      min_func)
    min_func = c['new_func']
    energy = c['new_energy']

    d = search_age(clusters, average_n, x, y, energy, age, exp_n, sigma_n,
                   min_func)
    min_func = d['new_func']
    age = d['new_age']

    return {'func': min_func, 'energy': energy, 'age': age}


def make_step(clusters, average_n, side, start_x, start_y, start_energy,
              start_age, exp_n, sigma_n, min_func):
    step_cen = divide_square(start_x, start_y, side)
    step = {'func': [], 'x': [], 'y': [], 'energy': [], 'age': []}
    theo_n = []

    for x in step_cen[0]:
        for y in step_cen[1]:
            for cluster in clusters:
                cluster.rec_particles(start_x, start_y, start_energy, start_age,
                                      average_n)
                for k in range(4):
                    # Получаем теоретические значения числа частиц
                    theo_n.append(cluster.stations[k].rec_particles)
            a = energy_age_search(clusters, average_n, x, y, exp_n, sigma_n,
                                  functional(exp_n, sigma_n, theo_n))
            step['func'].append(a['func'])
            step['x'].append(x)
            step['y'].append(y)
            step['energy'].append(a['energy'])
            step['age'].append(a['age'])

            theo_n = []

    if min(step['func']) < min_func:
        min_func = min(step['func'])
        min_square_num = step['func'].index(min(step['func']))
        new_x = step['x'][min_square_num]
        new_y = step['y'][min_square_num]
        new_energy = step['energy'][min_square_num]
        new_age = step['age'][min_square_num]
        return {'x': new_x, 'y': new_y, 'func': min_func, 'energy': new_energy,
                'age': new_age}
    else:
        return {'x': start_x, 'y': start_y, 'func': min_func,
                'energy': start_energy, 'age': start_age}
