import random as rn
from numpy import *
import matplotlib.pyplot as plt


ampl = []  # Массив для амплитуд от одиночного мюона
with open('data/mu.txt') as file:
    for line in file:
        ampl.append(float(line))  # Получаем данные из массива
my_hist = histogram(ampl, bins=arange(80.0), range=(0, 80), density=True)  # Строим гистограмму

w = my_hist[0]  # Здесь веса
a = my_hist[1]  # Здесь амплитуды
integral = 0
for l in range(0, 79):
    # Считаем интеграл
    integral += w[l] * a[l]


def get_amplitude():
    # Функция генерирует амплитуду от одиночного мюона в пКл
    gm = rn.random()  # Выбрали случайное число от 0 до 1 (гамма)
    part = integral * gm  # Берём какую-то часть от интеграла

    current_sum = 0  # Текущая сумма
    for k in range(0, 79):
        current_sum += w[k] * a[k]  # Накапливаем сумму
        # print(current_sum)
        if current_sum > part:
            return a[k]
    return a[79]


