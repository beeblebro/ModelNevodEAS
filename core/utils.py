# Вспомогательные функции

from math import pi, pow, sqrt, exp, cos, sin, gamma
from numpy import array, cross
from numpy.linalg import norm
from numpy.random import normal, poisson
import random as rn

g_ratio = (1 + 5 ** 0.5) / 2  # Золотое сечение
m_radius = 71  # Радиус Мольера [м]
light_speed = 0.299792458  # Скорость света

time_dist_func = []  # Функция распределения для случайной поправки ко времени
tau = 5 / sqrt(3)  # Параметр тау

for i in range(300):
    # Заполняем функцию распределения
    t = 0.5 * i

    new_element = (2 * pow(tau, 3) - tau * exp(-t / tau) * (pow(t, 2) + 2 * t *
                            tau + 2 * pow(tau, 2))) * (1 / (2 * pow(tau, 3)))

    time_dist_func.append(new_element)


def gen_rndm_time():
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
    imax = int(x + 5.0 * sqrt(x))
    if imax < 5:
        imax = 5
    gm = rn.random()
    sump = p
    if gm < sump:
        return 0
    for j in range(1, imax + 1):
        p *= (x/j)
        sump += p
        if gm < sump:
            return i
    return imax


def get_distance(st_coord, n, x0, y0):
    """Вычисляет расстояние от станции до оси ШАЛ"""
    # Вектор от станции до точки прихода ШАЛ
    a = array([st_coord[0] - x0, st_coord[1] - y0, st_coord[2]])
    b = cross(a, n)  # Векторное произведение вектора a на вектор ШАЛ
    # Возвращаем модуль векторного произведения, т.к длина vector равна единице
    return norm(b, ord=2)


def gen_x_y():
    """Более эффективные координаты ШАЛ"""
    x = rn.uniform(-40, 80)
    y = rn.uniform(-60, 80)

    if y < -20 and not -20 < x < 20:
        return gen_x_y()
    return x, y


def gen_theta():
    """Методом Неймана получаем случайный зенитный угол"""
    while True:
        kx1 = rn.random()
        kx2 = rn.random()
        kx11 = 0.0 + ((pi/2) - 0.0) * kx1
        kx22 = 1.0 * kx2

        if kx22 <= pow(cos(kx11), 8.5) * sin(kx11):
            return kx11 * (180/pi)


def gen_phi():
    """Получаем случайный азимутальный угол"""
    return rn.uniform(0, 360)


def gen_power():
    """Получаем мощность методом обратных функций"""
    beta = 2.5  # Показатель степени
    a = 10**5  # Нижняя граница диапазона мощности
    gm = rn.random()

    power = a * pow((1 - gm), 1/(1 - beta))

    while power < 10**4 or power > 10**9:
        power = a * pow((1 - gm), 1 / (1 - beta))

    return power


def gen_age(power, theta):
    """Генератор значений возраста"""
    theta_rad = theta * (pi/180)  # Перевели тета в радианы

    n_1 = 10**7
    s_1 = 0.5

    n_2 = 10**5
    s_2 = 1.35

    b = (s_2 - (n_2 * s_1)/n_1) / (1 - n_2/n_1)
    k = (s_1 - b)/n_1

    age = k * power + b
    age = normal(age + 0.3 * sin(theta_rad), 0.1)

    while age < 0.5 or age > 2.0:
        age = normal(age + 0.3 * sin(theta_rad), 0.1)

    return age


def divide_square(center_x, center_y, side):
    """Делит квадрат на девять равных и возвращает их центры"""
    result_x = []
    result_y = []
    new_side = side / 3
    left_corner_x = center_x - side / 2
    left_corner_y = center_y - side / 2
    for k in [1, 3, 5]:
        result_y.append(left_corner_y + k * (new_side/2))
        result_x.append(left_corner_x + k * (new_side/2))
    result = {'x': result_x, 'y': result_y}
    return result


def nkg(radius, power, age):
    """Функция пространственного распределения Нишимуры-Каматы-Грейзена"""
    # Разбили формулу на четыре множителя
    m1 = power / pow(m_radius, 2)
    m2 = gamma(4.5 - age) / (2 * pi * gamma(age) * gamma(4.5 - 2 * age))
    m3 = pow((radius + 1e-10) / m_radius, age - 2)
    m4 = pow(1 + (radius + 1e-10) / m_radius, age - 4.5)

    ro = m1 * m2 * m3 * m4
    # Возвращает поверхностную плотность на расстоянии radius от оси ливня
    return ro


def poisson_gauss_gen(n):
    """Генератор Пуассона и Гаусса"""
    if n <= 25:
        return poisson(n)
    else:
        return round(normal(n, sqrt(n)))


def check_effective_area(x, y):
    """Проверяет, попадаю ли x,y в эффективную область"""
    if y > -10:
        if y < 60 and -20 < x < 50:
            return True
    elif y > -40:
        if -10 < x < 5:
            return True
    else:
        return False

