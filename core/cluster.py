from math import pow, sqrt
from numpy import array
from numpy.random import poisson, normal
from numpy.linalg import det

from core.station import Station
from core.utils import get_distance, nkg, light_speed


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
        self.st_ok = None  # Число сработавших станций
        self.matches_required = 4  # Крастность совпадений

        self.stations = (
            # Создаём станции кластера, задаём им координаты и номера
            Station(1,
                    [self.center[0] - self.width / 2,
                     self.center[1] + self.length / 2,
                     self.center[2]]),

            Station(2,
                    [self.center[0] + self.width / 2,
                     self.center[1] + self.length / 2,
                     self.center[2]]),

            Station(3,
                    [self.center[0] - self.width / 2,
                     self.center[1] - self.length / 2,
                     self.center[2]]),

            Station(4,
                    [self.center[0] + self.width / 2,
                     self.center[1] - self.length / 2,
                     self.center[2]]),
        )

    def get_eas(self, eas):
        """Получаем ШАЛ без запуска"""
        self.eas = eas

    def start(self, eas):
        """Запуск кластера"""
        # Получаем ШАЛ
        self.eas = eas
        # Счётчик сработавших станций
        self.st_ok = 0
        for st in self.stations:
            # Запускаем станции
            if st.start(self.eas):
                self.st_ok += 1

        if self.st_ok >= self.matches_required:
            # Кластер сработал (четырёхкратные совпадения)
            self.respond = True
        else:
            # Кластер не сработал
            self.respond = False

        return self.respond

    def reset(self):
        """Возвращаем кластер к исходному состоянию"""
        self.respond = None
        self.eas = None
        self.time = None
        self.rec_n = None
        self.st_ok = None
        for station in self.stations:
            station.reset()

    def set_cluster_state(self, evt_cluster):
        """Устанавливаем состояние кластера  в соответствии
        с прочитанным событием"""
        self.st_ok = 0
        for st_n, st in enumerate(self.stations):
            if st.set_station_state(evt_cluster['st'][st_n]):
                self.st_ok += 1

        if self.st_ok == 4:
            self.respond = True
        else:
            self.respond = False

        return self.respond

    def make_times_relative(self):
        """Делает времена срабатывания станций относительными"""
        temp = [st.rndm_time for st in self.stations if st.respond]
        min_t = min(temp)
        for st in self.stations:
            if st.respond:
                st.rndm_time -= min_t

    def rec_direction(self):
        """Восстанавливает вектор прихода ШАЛ методом наименьших квадратов"""

        # Изменим времена срабатывания станций на относительные
        self.make_times_relative()

        sum0 = 0
        sum_x = 0
        sum_y = 0
        sum_xx = 0
        sum_yy = 0
        sum_xy = 0
        sum_t = 0
        sum_tx = 0
        sum_ty = 0

        for st in self.stations:
            sqr_sigma_t = pow(st.sigma_t, 2)

            sum0 += 1 / sqr_sigma_t
            sum_t += st.rndm_time / sqr_sigma_t
            sum_x += st.coord[0] / sqr_sigma_t
            sum_y += st.coord[1] / sqr_sigma_t
            sum_xx += pow(st.coord[0], 2) / sqr_sigma_t
            sum_yy += pow(st.coord[1], 2) / sqr_sigma_t
            sum_xy += (st.coord[0] * st.coord[1]) / sqr_sigma_t
            sum_tx += (st.rndm_time * st.coord[0]) / sqr_sigma_t
            sum_ty += (st.rndm_time * st.coord[1]) / sqr_sigma_t

        sqr_light_speed = pow(light_speed, 2)

        sum0 /= sqr_light_speed
        sum_xx /= sqr_light_speed
        sum_xy /= sqr_light_speed
        sum_yy /= sqr_light_speed
        sum_x /= sqr_light_speed
        sum_y /= sqr_light_speed
        sum_tx /= light_speed
        sum_ty /= light_speed
        sum_t /= light_speed

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

    def rec_particles(self, n, params):
        """Восстанавливаем число частиц в каждой станции, предполагая мощность
        координаты прихода ШАЛ, мощность и возраст. А вычисление расстояния до станций 
        от оси ШАЛ просиходит с помощью восстанолвенного вектора"""
        for station in self.stations:
            dist = get_distance(station.coord, n, params[0], params[1])
            station.rec_particles = station.area * n[2] * nkg(dist, params[2],
                                                              params[3])

    @staticmethod
    def _poisson_gauss_gen(n):
        """Генератор Пуассона и Гаусса"""
        if n <= 25:
            return poisson(n)
        else:
            return round(normal(n, sqrt(n)))
