from math import pow, sqrt, pi, gamma, modf, acos, cos
from numpy import array
from numpy.random import poisson, normal
from numpy.linalg import det

from core.station import Station
from core.amplitude import *
from core.utils import *


class Cluster:
    """Класс для представления кластера установки НЕВОД-ШАЛ"""

    def __init__(self, center, length=15, width=15):
        """Создаём кластер, передаём ему координаты и накрывший его ШАЛ"""
        self.eas = None  # ШАЛ, падающий на кластер
        self.center = array(center)  # Координаты кластера [м]
        self.length = length  # Длина кластреа (вдоль оси Y) [м]
        self.width = width  # Ширина кластера (вдоль оси X) [м]
        self.respond = None  # Отклик кластера
        self.rec_n = None  # Координаты восстановленного вектора
        self.time = None  # Время срабатывания кластера [нс]

        self.stations = (
            # Создаём станции кластера, задаём им координаты
            Station([self.center[0] + self.width / 2,
                     self.center[1] - self.length / 2,
                     self.center[2]]),

            Station([self.center[0] - self.width / 2,
                     self.center[1] - self.length / 2,
                     self.center[2]]),

            Station([self.center[0] - self.width / 2,
                     self.center[1] + self.length / 2,
                     self.center[2]]),

            Station([self.center[0] + self.width / 2,
                     self.center[1] + self.length / 2,
                     self.center[2]]),
        )

    def get_eas(self, eas):
        """Получаем ШАЛ"""
        self.eas = eas

    def start(self):
        """Запуск кластера"""
        if self.model_amplitudes():
            self.model_times()

            self.rec_direction()
            return self.respond
        else:
            return False

    def reset(self):
        """Возвращаем кластер к исходному состоянию"""
        self.eas = None
        self.time = None
        self.rec_n = None
        for station in self.stations:
            station.reset()

    def model_amplitudes(self):
        """Моделирование амплитуд, полученных от станций"""
        _enabled_gen = False  # Включить/выключить генераторы Пуассона и Гаусса

        for st in self.stations:
            dist = get_distance(st.coord, self.eas.n, self.eas.x0, self.eas.y0)
            temp = st.area * self.eas.n[2] * self.nkg(dist, self.eas.power, self.eas.age)

            if _enabled_gen:
                st.particles = self._poisson_gauss_gen(temp)
            else:
                st.particles = temp

            if st.particles == 0:
                # Станция не сработала
                st.respond = False
                # Не сработал и кластер (четырёхкратные совпадения)
                self.respond = False
                st.sigma_particles = 1.3
                st.amplitude = 0
            else:
                # Станция сработала
                st.respond = True
                st.amplitude = 0
                for j in range(int(st.particles)):
                    # Вычисляем амплитуду в станции
                    st.amplitude += get_amplitude()
                # Добавим десятичную часть
                st.amplitude += (get_amplitude() * modf(st.particles)[0])
                st.sigma_particles = sqrt(st.particles * get_sqr_sigma())

        if self.respond is None:
            # Если ни одна станция не сменила отклик кластера на False => кластер сработал
            self.respond = True

        if self.respond:
            return True
        else:
            return False

    def model_times(self):
        """Моделирует времена срабатывания станций"""
        temp = []
        for station in self.stations:
            station.real_time = (sum(self.eas.n * station.coord) +
                                 self.eas.D) / self.eas.light_speed

            station.rndm_time = randomize_time() + station.real_time
            temp.append(station.rndm_time)

        self.time = max(temp)
        min_t = min(temp)

        for station in self.stations:
            station.rndm_time -= min_t

    def rec_direction(self):
        """Восстанавливает вектор прихода ШАЛ методом наименьших квадратов"""

        sum0 = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_yy = 0
        sum_xy = 0
        sum_t = 0
        sum_tx = 0
        sum_ty = 0

        for i in range(4):
            sqr_sigma_t = pow(self.stations[i].sigma_t, 2)

            sum0 += 1 / sqr_sigma_t
            sum_t += self.stations[i].rndm_time / sqr_sigma_t
            sum_x += self.stations[i].coord[0] / sqr_sigma_t
            sum_y += self.stations[i].coord[1] / sqr_sigma_t
            sum_xx += pow(self.stations[i].coord[0], 2) / sqr_sigma_t
            sum_yy += pow(self.stations[i].coord[1], 2) / sqr_sigma_t
            sum_xy += (self.stations[i].coord[0] * self.stations[i].coord[1]) / sqr_sigma_t
            sum_tx += (self.stations[i].rndm_time * self.stations[i].coord[0]) / sqr_sigma_t
            sum_ty += (self.stations[i].rndm_time * self.stations[i].coord[1]) / sqr_sigma_t

        sqr_light_speed = pow(self.eas.light_speed, 2)

        sum0 /= sqr_light_speed
        sum_xx /= sqr_light_speed
        sum_xy /= sqr_light_speed
        sum_yy /= sqr_light_speed
        sum_x /= sqr_light_speed
        sum_y /= sqr_light_speed
        sum_tx /= self.eas.light_speed
        sum_ty /= self.eas.light_speed
        sum_t /= self.eas.light_speed

        line_1 = [sum_xx, sum_xy, sum_x]
        line_2 = [sum_xy, sum_yy, sum_y]
        line_3 = [sum_x, sum_y, sum0]
        line_4 = [sum_tx, sum_ty, sum_t]

        det1 = det([line_1, line_2, line_3])
        det2 = det([line_4, line_2, line_3])
        det3 = det([line_1, line_4, line_3])
        # det4 = det([line_1, line_2, line_4])

        a = det2 / det1
        b = det3 / det1

        if (pow(a, 2) + pow(b, 2)) <= 1:
            # Если вектор восстановился успешно
            c = sqrt(1 - pow(a, 2) - pow(b, 2))
            self.rec_n = array([a, b, c])
            self.respond = True
            return True
        else:
            self.respond = False
            print("ERROR: Не удалось восстановить направление")
            return False

    def rec_particles(self, n, x, y, power, age):
        """Восстанавливаем число частиц в каждой станции, предполагая мощность
        координаты прихода ШАЛ, мощность и возраст. А вычисление расстояние до станций 
        от оси ШАЛ просиходит с помощью восстанолвенного вектора"""
        for station in self.stations:
            dist = get_distance(station.coord, n, x, y)
            station.rec_particles = station.area * n[2] * self.nkg(dist, power, age)

    def nkg(self, radius, power, age):
        """Функция пространственного распределения Нишимуры-Каматы-Грейзена"""
        # Разбили формулу на четыре множителя
        m1 = power / pow(self.eas.m_rad, 2)
        m2 = gamma(4.5 - age) / (2 * pi * gamma(age) * gamma(4.5 - 2 * age))
        m3 = pow(radius / self.eas.m_rad, age - 2)
        m4 = pow(1 + radius / self.eas.m_rad, age - 4.5)

        ro = m1 * m2 * m3 * m4
        # Возвращает поверхностную плотность на расстоянии radius от оси ливня
        return ro

    @staticmethod
    def _poisson_gauss_gen(n):
        """Генератор Пуассона и Гаусса"""
        if n <= 25:
            return poisson(n)
        else:
            return round(normal(n, sqrt(n)))
