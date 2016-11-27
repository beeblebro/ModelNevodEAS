import random as rn
from numpy import *
import matplotlib.pyplot as plt


ampl = []  # Массив для амплитуд от одиночного мюона
with open('data/mu.txt') as file:
    for line in file:
        ampl.append(float(line))  # Получаем данные из массива
my_hist = histogram(ampl, bins=arange(80.0), range=(0, 80), density=True)  # Строим гистограмму
apml = array(ampl)
average_ampl = average(ampl)  # Вычислили среднюю амплитуду


w = my_hist[0]  # Здесь веса
a = my_hist[1]  # Здесь амплитуды


integral = 0
for l in range(0, 79):
    # Считаем интеграл
    integral += w[l]


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


def ampl_func():
    """Генерирует случайную амплитуду используя аналитическую функцию"""
    y0 = 0.55306
    A = 62.68954
    xc = 13.2051
    w = 4.01761

    while True:
        kx1 = rn.random()
        kx2 = rn.random()
        kx11 = 0 + (80 - 0) * kx1
        kx22 = 65 * kx2
        z = (kx11 - xc)/w

        if kx22 <= y0 + A * exp(-exp(-z) - z + 1):
            return kx11


print(ampl_func.__doc__)
'''
f = open('data/amplitude.txt', 'w')
for i in range(0, 750):
    f.write(str(get_amplitude()))
    f.write('\n')
'''