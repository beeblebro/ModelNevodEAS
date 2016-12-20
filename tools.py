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
    for j in range(1, imax):
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


def functional(exp_n, theo_n, sigma_n):
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
