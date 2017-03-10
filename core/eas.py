from math import pi, sin, cos
from numpy import array
import random as rn


class Eas:
    """Класс для представления широкого атмосферного ливня"""
    light_speed = 0.299792458

    def __init__(self, theta, phi, x0=0, y0=0, energy=10**6, age=1.3):
        self.theta = theta * (pi/180)  # Зенитный угол в радианах
        self.phi = phi * (pi/180)  # Азимутальный в радианах
        self.x0 = x0  # Координаты (x,y) точки
        self.y0 = y0  # прихода ШАЛ
        self.energy = energy  # Энергия
        self.age = age  # Возраст
        self.m_rad = 71  # Радиус Мольера
        self.D = 1000  # Параметр D плоскости ливня
        self.tau = 5 / (3**0.5)  # Параметр тау в формуле временного профиля ШАЛ
        # n - истинный вектор прихода ШАЛ. т.е. нормаль к плоскости ливня
        self.n = array([sin(self.theta) * cos(self.phi),
                        sin(self.theta) * sin(self.phi),
                        cos(self.theta)])

    def generate_energy(self):
        """Сгенерировать энергию"""
        self.energy = ((10 ** 6) / (1 - rn.random())) ** (2 / 3)

