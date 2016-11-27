from math import *
from numpy import *
from numpy.linalg import norm
import random as rn

time_distr_func = []
tau = 5 / sqrt(3)
for i in range(0, 300):
    # Заполняем функцию распределения
    t = 0.5 * i
    new_element = (2 * pow(tau, 3) - tau * exp(-t / tau) * (pow(t, 2) + 2 * t * tau + 2 * pow(tau, 2))) * (1 / (2 * pow(tau, 3)))
    time_distr_func.append(new_element)


def randomize_time():
    """Возвращает случайную поправку ко времени срабатывания станции"""
    gm = rn.random()
    for j in range(0, 300):
        if gm < time_distr_func[j]:
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
    for j in range(1, imax):
        p = p * x/j
        sump += p
        if gm < sump:
            return i
    return imax


def get_distance(st_coord, n, x0, y0):
    """Вычисляет расстояние от станции до оси ШАЛ"""
    a = array([st_coord[0] - x0, st_coord[1] - y0, st_coord[2]])  # Вектор от станции до точки прихода ШАЛ
    b = cross(a, n)  # Векторное произведение вектора a на вектор ШАЛ
    return norm(b, ord=2)  # Возвращаем модуль векторного произведения, т.к длин n равна единице


def get_theta():
    """Методом Неймана получаем случайные зенитный угол с необходимиым распределением"""
    while True:
        kx1 = rn.random()
        kx2 = rn.random()
        kx11 = 0.0 + ((pi/2) - 0.0) * kx1
        kx22 = 1.0 * kx2

        if kx22 <= pow(cos(kx11), 8.5) * sin(kx11):
            return kx11 * (180/pi)


def functional(exp_n, theor_n):
    """Считает функционал"""
    f = 0
    for j in range(0, size(exp_n)):
        f += (exp_n[j] - theor_n[j])**2
    return f

