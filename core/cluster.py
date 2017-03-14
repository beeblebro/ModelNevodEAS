from math import pow, sqrt, pi, gamma, acos, cos
from numpy import array
from numpy.random import poisson, normal
from numpy.linalg import det
import json

from core.tools import randomize_time, get_distance
from core.amplitude import get_amplitude, get_av_amplitude, get_sqr_sigma
from core.station import Station


class Cluster:
    """Класс для представления кластера установки НЕВОД-ШАЛ"""

    def __init__(self, center, eas, length=15, width=15):
        """Создаём кластер, передаём ему координаты и накрывший его ШАЛ"""
        self.center = array(center)  # Координаты кластера [м]
        self.length = length  # Длина кластреа (вдоль оси Y) [м]
        self.width = width  # Ширина кластера (вдоль оси X) [м]
        self.eas = eas  # ШАЛ, накрывший кластр

        self.failed = False  # Параметр отказа кластера
        self.rec_n = array([0, 0, 0])  # Координаты восстановленного вектора
        self.time = 0  # Меньшее из времён срабатываний станций [нс]

        self.times = []  # Времена срабатывания

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

    def time_func(self):
        """Вычисляет времена времена срабатывания станций"""
        for station in self.stations:
            station.real_time = (sum(self.eas.n * station.coordinates) +
                                 self.eas.D) / self.eas.light_speed

            station.rndm_time = randomize_time() + station.real_time
            self.times.append(station.rndm_time)

        self.time = min(self.times)

        for station in self.stations:
            station.rndm_time -= self.time

    def least_squares(self):
        """Восстанавливает вектор прихода ШАЛ методом наименьших квадратов"""
        self.time_func()  # Получили времена срабатывания станций

        # Если хотя бы одно станции не хватило частиц - кластер не сработал
        if not self.model_amplitudes():
            return False

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
            sum_x += self.stations[i].coordinates[0] / sqr_sigma_t
            sum_y += self.stations[i].coordinates[1] / sqr_sigma_t
            sum_xx += pow(self.stations[i].coordinates[0], 2) / sqr_sigma_t
            sum_yy += pow(self.stations[i].coordinates[1], 2) / sqr_sigma_t
            sum_xy += (self.stations[i].coordinates[0] * self.stations[i].coordinates[1]) / sqr_sigma_t
            sum_tx += (self.stations[i].rndm_time * self.stations[i].coordinates[0]) / sqr_sigma_t
            sum_ty += (self.stations[i].rndm_time * self.stations[i].coordinates[1]) / sqr_sigma_t

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
            return True
        else:
            self.failed = True
            print("ERROR: Не удалось восстановить направление")
            return False

    def nkg(self, radius):
        """Функция пространственного распределения Нишимуры-Каматы-Грейзена"""
        # Разбили формулу на четыре множителя
        m1 = self.eas.energy / pow(self.eas.m_rad, 2)
        m2 = gamma(4.5 - self.eas.age) / (2 * pi * gamma(self.eas.age) *
                                          gamma(4.5 - 2 * self.eas.age))

        m3 = pow(radius / self.eas.m_rad, self.eas.age - 2)
        m4 = pow(1 + radius / self.eas.m_rad, self.eas.age - 4.5)

        ro = m1 * m2 * m3 * m4
        # Возвращает поверхностную плотность на расстоянии radius от оси ливня
        return ro

    def NKG(self, radius, energy, age):
        """НКГ, в которой варьируются мощность и возраст"""
        # Разбили формулу на четыре множителя
        m1 = energy / pow(self.eas.m_rad, 2)
        m2 = gamma(4.5 - age) / (2 * pi * gamma(age) * gamma(4.5 - 2 * age))
        m3 = pow(radius / self.eas.m_rad, age - 2)
        m4 = pow(1 + radius / self.eas.m_rad, age - 4.5)

        ro = m1 * m2 * m3 * m4
        # Возвращает поверхностную плотность на расстоянии radius от оси ливня
        return ro

    def model_amplitudes(self):
        """Моделирование амплитуд, полученных от станций"""
        for i in range(4):
            dist = get_distance(self.stations[i].coordinates, self.eas.n,
                                self.eas.x0, self.eas.y0)
            # Не радномизированное число частиц в станции
            temp = self.stations[i].area * self.eas.n[2] * self.nkg(dist)

            if temp <= 25:
                # Используем Пуассон
                self.stations[i].particles = poisson(temp)
            else:
                # Используем Гаусс
                self.stations[i].particles = round(normal(temp, sqrt(temp)))

            if self.stations[i].particles < 0:
                print("ERROR: Отрицательное число частиц в Cluster")

            # Станция не сработала
            if self.stations[i].particles == 0:
                self.stations[i].failed = True
                # Не сработал и кластер (четырёхкратные совпадения)
                self.failed = True
                self.stations[i].sigma_particles = 1.3
                self.stations[i].amplitude = 0.0
            # Станция сработала
            else:
                for j in range(self.stations[i].particles):
                    # Вычислили амплитуду в станции
                    self.stations[i].amplitude += get_amplitude()

        if not self.failed:
            return True
        else:
            return False

    def model_amplitudes_simple(self):
        """Упрощенная версия моделирования амплитуд.
        Пригодится, если интересует только факт срабатывания"""
        for i in range(4):
            dist = get_distance(self.stations[i].coordinates, self.eas.n,
                                self.eas.x0, self.eas.y0)
            # Число частиц в станции
            self.stations[i].particles = poisson(self.stations[i].area *
                                                 self.nkg(dist))

            # Станция не сработала
            if self.stations[i].particles == 0:
                self.stations[i].failed = True
                # Не сработал и кластер (четырёхкратные совпадения)
                self.failed = True
                return False
            # Станция сработала
            else:
                return True

    def recons_particles(self, sug_x, sug_y, average_n):
        """Восстанавливаем число частиц в станциях. Получаем точку, в которой
        предполагаем центр ШАЛ и средний восстановленный вектор ШАЛ"""
        for station in self.stations:
            dist = get_distance(station.coordinates, average_n, sug_x, sug_y)
            station.rec_particles = station.area * self.eas.n[2] * self.nkg(dist)

    def rec_particles(self, average_n, sug_x, sug_y, sug_energy, sug_age):
        """Восстанавливаем число частиц, варьируя возраст и мощность"""
        for station in self.stations:
            dist = get_distance(station.coordinates, average_n, sug_x, sug_y)
            station.rec_particles = station.area * self.eas.n[2] * \
                                    self.NKG(dist, sug_energy, sug_age)

    def record_event(self):
        """Записывает событие в формате json"""
        event = {
            'energy': 'self.eas.energy',
            'age': self.eas.age,
            'theta': self.eas.theta,
            'phi': self.eas.phi,
            'real vector': str(self.eas.n),
        }

        filename = 'data/EAS_events.json'
        with open(filename, 'w') as f_obj:
            json.dump(event, f_obj, sort_keys=True, indent=4, separators=(',', ': '))
