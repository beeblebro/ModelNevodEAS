import random as rn
from math import exp, log, e
import matplotlib.pyplot as plt


def get_amplitude():
    """Получаем случайную амплитуду методом обратной функции"""
    _enabled_gen = False  # Вкл/Выкл генератор
    if _enabled_gen:
        xc = 13.5
        wid = 4
        # Из условия нормировки
        const_a = 1 / (wid * (e - exp(1 - exp(xc/wid))))

        gm = rn.random()
        value = xc - wid * log(-log(gm/(const_a * wid * e) + exp(-exp(xc/wid))))
        return value
    else:
        return get_av_amplitude()


def get_av_amplitude():
    """Возвращает среднюю амплитуду"""
    # return average_ampl
    return 15.80949


def get_sqr_sigma():
    """Возвращает квадрат стандартной ошибки"""
    # return sqr_sigma
    return 26.31064**2
