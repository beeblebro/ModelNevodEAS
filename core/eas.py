from math import sin, cos
from numpy import array, radians

from core.utils import get_theta, get_power, get_age, modified_area, get_phi


class Eas:
    """Класс для представления широкого атмосферного ливня"""
    light_speed = 0.299792458

    def __init__(self, params=None):

        if not params:
            # Разыграем параметры ШАЛ, если они не были переданы

            # Углы наклона оси в градусах
            self.theta_deg = get_theta()
            self.phi_deg = get_phi()

            # Углы наклона оси в радианах
            self.theta = radians(self.theta_deg)
            self.phi = radians(self.phi_deg)
            # Параметры ШАЛ
            self.x0, self.y0 = modified_area()
            self.power = get_power()
            self.age = get_age(self.power, self.theta)
        else:
            # Передали параметры в виде словаря

            self.theta_deg = params['theta']
            self.phi_deg = params['phi']

            self.theta = radians(self.theta_deg)
            self.phi = radians(self.phi_deg)

            self.x0 = params['x0']
            self.y0 = params['y0']
            self.power = params['power']
            self.age = params['age']

        self.m_rad = 71  # Радиус Мольера
        self.D = 1000  # Параметр D плоскости ливня
        self.tau = 5 / (3**0.5)  # Параметр тау в формуле временного профиля ШАЛ
        # n - истинный вектор прихода ШАЛ. т.е. нормаль к плоскости ливня
        self.n = array([sin(self.theta) * cos(self.phi),
                        sin(self.theta) * sin(self.phi),
                        cos(self.theta)])

    def get_params_list(self):
        """Получить параметры ШАЛ списком"""
        return [self.theta_deg, self.phi_deg,
                self.x0, self.y0,
                self.power, self.age]

    def get_params_dict(self):
        """Получить параметры ШАЛ словарём"""
        return {'theta': self.theta_deg,
                'phi': self.phi_deg,
                'x0': self.x0,
                'y0': self.y0,
                'power': self.power}
