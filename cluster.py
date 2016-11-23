from math import *
from numpy import *
from numpy.linalg import norm, det
import random as rn

from tools import randomize_time, psn, get_distance
from amplitude import  get_amplitude


class Cluster:
    # Класс для представления кластера установки НЕВОД-ШАЛ
    area = 2.56  # Площадь станции
    width = 20  # Ширина кластера
    length = 20  # Длина кластера
    sigma_t = (5, 5, 5, 5)  # Ошибки в определении времени
    failed = False  # Параметр отказа кластера
    rec_n = array([0, 0, 0])  # Координаты восстановленного вектора
    time = 0.0  # Меньшее из времён срабатываний станций
    st_failed = [0, 0, 0, 0]  # Параметры отказа станций
    st_particles = [0, 0, 0, 0]  # Число частиц, попавших в каждую станцию
    st_amplitudes = [0, 0, 0, 0]  # Амплитуда в каждой станции
    real_time = [0, 0, 0, 0]  # Истинное время срабатывания станций
    rndm_time = [0, 0, 0, 0]  # Рандомизированное время срабатывания станций

    def __init__(self, center, eas, length=20, width=20):
        self.center = array(center)
        self.length = length
        self.width = width
        self.eas = eas

        self.st_coord = [
            # Вычисляем координаты станций кластера
            array([self.center[0] + self.width / 2, self.center[1] - self.length / 2, self.center[2]]),
            array([self.center[0] - self.width / 2, self.center[1] - self.length / 2, self.center[2]]),
            array([self.center[0] - self.width / 2, self.center[1] + self.length / 2, self.center[2]]),
            array([self.center[0] + self.width / 2, self.center[1] + self.length / 2, self.center[2]])
        ]

    def time_func(self):
        # Вычисляет времена истинные и рандомизированные времена срабатывания станций
        for i in range(0, 4):
            self.real_time[i] = (norm(self.eas.n * self.st_coord[i], ord=1) + self.eas.D) / self.eas.light_speed
            self.rndm_time[i] = randomize_time() + self.real_time[i]

    def least_squares(self):
        # Восстанавливает вектор прихода ШАЛ методом наименьших квадратов
        self.time_func()  # Получили времена срабатывания станций
        if not self.coutn_particles():  # Если хотя бы одно станции не хватило частиц - кластер не сработал
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

        for i in range(0, 4):
            sum0 += 1 / pow(self.sigma_t[i], 2)
            sum_t += self.rndm_time[i] / pow(self.sigma_t[i], 2)
            sum_x += self.st_coord[i][0] / pow(self.sigma_t[i], 2)
            sum_y += self.st_coord[i][1] / pow(self.sigma_t[i], 2)
            sum_xx += pow(self.st_coord[i][0], 2) / pow(self.sigma_t[i], 2)
            sum_yy += pow(self.st_coord[i][1], 2) / pow(self.sigma_t[i], 2)
            sum_xy += (self.st_coord[i][0] * self.st_coord[i][1]) / pow(self.sigma_t[i], 2)
            sum_tx += (self.rndm_time[i] * self.st_coord[i][0]) / pow(self.sigma_t[i], 2)
            sum_ty += (self.rndm_time[i] * self.st_coord[i][1]) / pow(self.sigma_t[i], 2)

        sum0 /= pow(self.eas.light_speed, 2)
        sum_xx /= pow(self.eas.light_speed, 2)
        sum_xy /= pow(self.eas.light_speed, 2)
        sum_yy /= pow(self.eas.light_speed, 2)
        sum_x /= pow(self.eas.light_speed, 2)
        sum_y /= pow(self.eas.light_speed, 2)
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
            return False

    def nkg(self, radius):
        # Функция пространственного распределения Нишимуры-Каматы-Грейзена
        m1 = self.eas.energy / pow(self.eas.m_rad, 2)  # Разбили формулу на четыре множителя
        m2 = gamma(4.5 - self.eas.age) / (2 * pi * gamma(self.eas.age) * gamma(4.5 - 2 * self.eas.age))
        m3 = pow(radius / self.eas.m_rad, self.eas.age - 2)
        m4 = pow(1 + radius / self.eas.m_rad, self.eas.age - 4.5)

        ro = m1 * m2 * m3 * m4
        return ro  # Возвращает поверхностную плотность на расстоянии radius от оси ливня

    def coutn_particles(self):
        # Считает число частиц, попавших в каждую станцию кластера
        for i in range(0, 4):
            dist = get_distance(self.st_coord[i], self.eas.n, self.eas.x0, self.eas.y0)
            self.st_particles[i] = psn(self.area * self.nkg(dist))
            self.st_amplitudes[i] = self.st_particles[i] * get_amplitude()
            if self.st_particles[i] == 0:  # Станция не сработала
                self.st_particles = [0, 0, 0, 0]
                self.st_failed[i] = 1
                self.failed = True
                return False
        return True

    def rec_particles(self):
        #  Восстанавливаем число частиц в станциях. Получаем точку, которую
        pass
