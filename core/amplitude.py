import random as rn
from math import exp, log, e, pi

# Постоянная Эйлера-Маскерони
euler_masch_const = 0.5772156649015328606
xc = 13.5
wid = 4
# Средняя амплитуда
av_ampl = xc + wid * euler_masch_const
# Дисперсия или квадрат сигма
variance = (pi**2) * (wid**2) / 6


def get_amplitude():
    """Получаем случайную амплитуду методом обратной функции"""
    _enabled_gen = True  # Вкл/Выкл генератор
    if _enabled_gen:
        # Из условия нормировки получаем константу
        const_a = 1 / (wid * (e - exp(1 - exp(xc/wid))))

        gm = rn.random()
        value = xc - wid * log(-log(gm/(const_a * wid * e) + exp(-exp(xc/wid))))
        return value
    else:
        return get_av_amplitude()


def get_av_amplitude():
    """Возвращает среднюю амплитуду"""
    return av_ampl


def get_sqr_sigma():
    """Возвращает квадрат стандартной ошибки"""
    return variance

