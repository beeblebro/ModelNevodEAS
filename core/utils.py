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
