import random as rn
from math import exp
from numpy import array, histogram, arange, average, power
import matplotlib.pyplot as plt


def ampl_func():
    """Случайная амплитуда методом Неймана"""
    y0 = 0
    A = 63
    xc = 13.5
    W = 4

    while True:
        kx1 = rn.random()
        kx2 = rn.random()
        kx11 = 0 + (80 - 0) * kx1
        kx22 = 63 * kx2
        z = (kx11 - xc)/W

        if kx22 <= y0 + A * exp(-exp(-z) - z + 1):
            return kx11

ampl = []  # Массив для амплитуд от одиночного мюона
with open('a1.txt') as file:
    for line in file:
        ampl.append(float(line))  # Получаем данные из массива

ampl = array(ampl)
sqr_ampl = power(ampl, 2)  # Квадраты амплитуд
average_sqr = average(sqr_ampl)  # Средний квадрат
average_ampl = average(ampl)  # Вычислили среднюю амплитуду
sqr_sigma = average_sqr - average_ampl**2  # Квадрат стандартной ошибки

# print(average_ampl)
# print(sqr_sigma)

# Строим гистограмму
# my_hist = histogram(ampl, density=True)
my_hist = histogram(ampl, bins=arange(80.0), range=(0, 80), density=True)
w = my_hist[0]  # Здесь веса
a = my_hist[1]  # Здесь амплитуды
integral = sum(w)


def get_av_amplitude():
    """Возвращает среднюю амплитуду"""
    return average_ampl


def get_sqr_sigma():
    """Возвращает квадрат стандартной ошибки"""
    return sqr_sigma


def get_amplitude():
    """Функция генерирует амплитуду от одиночного мюона в пКл"""
    gm = rn.random()  # Выбрали случайное число от 0 до 1 (гамма)
    part = integral * gm  # Берём какую-то часть от интеграла

    current_sum = 0  # Текущая сумма
    for k in range(len(a)):
        current_sum += w[k]  # Накапливаем сумму
        if current_sum > part:
            delta = current_sum - part
            correction = delta / w[k]
            return a[k] - correction
    return len(a)

