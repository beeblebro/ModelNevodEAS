import random as rn
from math import exp
from numpy import array, histogram, arange, average, power
import matplotlib.pyplot as plt


# def ampl_func():
#    """Генерирует случайную амплитуду методом Неймана используя аналитическую функцию"""
#    y0 = 0.55306
#    A = 62.68954
#    xc = 13.2051
#    w = 4.01761
#
#    while True:
#        kx1 = rn.random()
#        kx2 = rn.random()
#        kx11 = 0 + (80 - 0) * kx1
#        kx22 = 65 * kx2
#        z = (kx11 - xc)/w
#
#        if kx22 <= y0 + A * exp(-exp(-z) - z + 1):
#            return kx11

ampl = []  # Массив для амплитуд от одиночного мюона
with open('data/analit_func.txt') as file:
    for line in file:
        ampl.append(float(line))  # Получаем данные из массива

apml = array(ampl)
sqr_ampl = power(ampl, 2)  # Квадраты амплитуд
average_sqr = average(sqr_ampl)  # Средний квадрат
average_ampl = average(ampl)  # Вычислили среднюю амплитуду
sqr_sigma = average_sqr - average_ampl**2  # Квадрат стандартной ошибки

# Строим гистограмму
my_hist = histogram(ampl, bins=arange(80.0), range=(0, 80), density=True)
w = my_hist[0]  # Здесь веса
a = my_hist[1]  # Здесь амплитуды


def function(x):
    y0 = 0.55306
    A = 62.68954
    xc = 13.2051
    wn = 4.01761

    z = (x - xc) / wn

    return y0 + A * exp(-exp(-z) - z + 1)

# print(function(13))
# print(w[13])

analit_int = 0
with open('data/ampl_test.txt', 'w') as f:
    for ai in a:
        analit_int += function(ai)

    for i in range(79):
        # print(w[i] - function(a[i]) / analit_int)
        f.write(str(w[i] - function(a[i]) / analit_int))
        f.write('\t')
        f.write(str(a[i]))
        f.write('\n')


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
    for k in range(0, 79):
        current_sum += w[k]  # Накапливаем сумму
        if current_sum > part:
            delta = current_sum - part
            correction = delta / w[k]
            return a[k] - correction
    return a[79]
