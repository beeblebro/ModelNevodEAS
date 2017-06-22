import random as rn
from math import exp, log, e, pi
from math import modf


def get_ca_amplitude():
    """Возвращает калибровочную амплитуду станции"""
    return 13.5

# Постоянная Эйлера-Маскерони
euler_masch_const = 0.5772156649015328606
xc = get_ca_amplitude()
wid = 4
# Средняя амплитуда
av_ampl = xc + wid * euler_masch_const
# Дисперсия или квадрат сигма
variance = (pi**2) * (wid**2) / 6

det_threshold = xc / 8  # Порог срабатывания детектора


def get_amplitude():
    """Возвращает случайную амплитуду методом обратной функции"""
    _enabled_gen = True  # Вкл/Выкл генератор
    if _enabled_gen:
        # Из условия нормировки получаем константу
        const_a = 1 / (wid * (e - exp(1 - exp(xc/wid))))

        gm = rn.random()
        value = xc - wid * log(-log(gm/(const_a * wid * e) + exp(-exp(xc/wid))))
        return value
    else:
        return get_av_amplitude()


def get_amplitudes(p):
    """Возвращает амплитуду от прохождения p частиц"""
    ampl = 0.0
    if p == 0:
        return ampl
    else:
        for j in range(int(p)):
            # Вычисляем амплитуду в детекторе
            ampl += get_amplitude()
        # Добавим десятичную часть
        ampl += (get_amplitude() * modf(p)[0])

        if ampl < det_threshold:
            ampl = 0.0

        return ampl


def get_av_amplitude():
    """Возвращает среднюю амплитуду"""
    return av_ampl


def get_sqr_sigma():
    """Возвращает квадрат стандартной ошибки"""
    return variance

