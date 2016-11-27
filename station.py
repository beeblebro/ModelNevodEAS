from math import *
from numpy import *


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, coordinates):
        self.coordinates = array(coordinates)
        self.side = 1.6  # Длина стороны станции в метрах
        self.area = 2.56  # Площадь станции в квадратных метрах
        self.failed = False  # Параметр отказа станции
        self.particles = 0  # Экспериментальное число частиц, попавших в стнцию
        self.real_time = 0  # Истинное время срабатывания станции
        self.rndm_time = 0  # Рандомизированное время срабатывания станции
        self.rec_particles = 0  # Восстановленное число частиц
