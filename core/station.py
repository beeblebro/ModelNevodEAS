from numpy import array, sign
from core.utils import get_distance, nkg, poisson_gauss_gen, randomize_time, light_speed
from core.amplitude import *
from math import sqrt

h_side = 0.8  # Половина стороны станции [м]
side = 1.6  # Длина стороны станции [м]
det_area = 0.64  # Площадь детектора [м2]
det_threshold = get_ca_amplitude() / 8  # Порог срабатывания детектора
st_threshold = get_ca_amplitude() / 2  # Порог срабатывания станции


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, number, coordinates):

        self.num = number  # Нормер станции в кластере
        self.coord = array(coordinates)  # Координаты станции
        self.area = 2.56  # Площадь станции [м2]
        self.sigma_t = 3.7  # Точность временной привязки [нс]
        self.respond = None  # Отклик станции
        self.particles = None  # Экспериментальное число частиц
        self.real_time = None  # Истинное время срабатывания станции [нс]
        self.rndm_time = None  # Рандомизированное время станции [нс]
        self.rec_particles = None  # Восстановленное число частиц
        self.amplitude = None  # Амплитуда сигнала от станции [пКл]
        self.sigma_particles = None  # Ошибка для функционала по частицам
        self.add_ampl = None  # Амплитуда сигнала от доп.ФЭУ

        # Детекторы станции
        self.detectors = (

            {
                'coord': (self.coord[0] + (-1)**(self.num + 1) * h_side,
                          self.coord[1] + sign(2.5 - self.num) * h_side,
                          self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] + (-1)**self.num * h_side,
                          self.coord[1] + sign(2.5 - self.num) * h_side,
                          self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] + (-1)**self.num * h_side,
                          self.coord[1] + sign(self.num - 2.5) * h_side,
                          self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] + (-1)**(self.num + 1) * h_side,
                          self.coord[1] + sign(self.num - 2.5) * h_side,
                          self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},
        )

        self.add_det = {
            # Дополнительный соседствует с 4ым ФЭУ
            'coord': self.detectors[3]['coord'],
            'time': None,
            'ampl': None,
            'particles': None,
            'respond': None
        }

    def start(self, eas):
        """Запуск станции"""
        if self.model_ampl(eas):
            self.model_times(eas)

        return self.respond

    def reset(self):
        """Возвращаемся к исходному состоянию"""
        self.respond = None
        self.particles = None
        self.real_time = None
        self.rndm_time = None
        self.rec_particles = None
        self.amplitude = None
        self.sigma_particles = None

        for i in range(len(self.detectors)):
            # Обнуляем словари детекторов
            self.detectors[i]['time'] = None
            self.detectors[i]['ampl'] = None
            self.detectors[i]['respond'] = None
            self.detectors[i]['particles'] = None

        # Обнуляем данные дополнительного ФЭУ
        self.add_det['time'] = None
        self.add_det['ampl'] = None
        self.add_det['respond'] = None
        self.add_det['particles'] = None

    def set_station_state(self, evt_station):
        """Устанавливаем состояние станции в соответствии с
        прчитанным событием"""
        self.amplitude = evt_station['ampl']
        self.rndm_time = evt_station['time']

        if self.amplitude > st_threshold and self.rndm_time is not None:
            self.respond = True
        else:
            self.respond = False
            self.amplitude = 0.0
            self.rndm_time = None

        return self.respond

    def model_ampl(self, eas):
        """Моделирование амплитуд в счётчиках станций"""
        enabled_gen = True  # Включить/выключить генераторы Пуассона и Гаусса
        self.amplitude = 0.0
        self.sigma_particles = 0

        for det in self.detectors:
            dist = get_distance(det['coord'], eas.n, eas.x0, eas.y0)
            temp = det_area * abs(eas.n[2]) * nkg(dist, eas.power, eas.age)

            if enabled_gen:
                p = poisson_gauss_gen(temp)
            else:
                p = temp

            det['particles'] = p
            det['ampl'] = get_amplitudes(p)

            if det['ampl'] < det_threshold:
                # Порог не превышен - детектор не сработал
                det['respond'] = False
            else:
                # Детектор сработал
                det['respond'] = True

                # Накапливаем сигму по частицам и амплитуду станции
                self.sigma_particles += sqrt(p * get_sqr_sigma())
                self.amplitude += det['ampl']

        # Моделируем отклик дополнительного на основании отклика его соседа
        self.add_det['particles'] = self.detectors[3]['particles']
        self.add_det['respond'] = self.detectors[3]['respond']
        self.add_det['ampl'] = get_amplitudes(self.add_det['particles']) / 100

        if self.sigma_particles == 0:
            self.sigma_particles = 1.3

        if self.amplitude > 0:
            self.respond = True
        else:
            self.respond = False

        return self.respond

    def model_times(self, eas):
        """Моделирует времена срабатывания станций"""
        temp = []
        for det in self.detectors:
            # Моделируем времена срабатывания детекторов
            if det['respond']:
                det['time'] = (sum(eas.n * det['coord']) + eas.D) / light_speed
                det['time'] += randomize_time()
                temp.append(det['time'])
            else:
                det['time'] = None

        self.rndm_time = min(temp)

        # Смоделируем время срабатывания доп. ФЭУ
        if self.add_det['respond']:
            self.add_det['time'] = (sum(eas.n * self.add_det['coord']) + eas.D) / light_speed
            self.add_det['time'] += randomize_time()
        else:
            self.add_det['time'] = None

        return self.respond
