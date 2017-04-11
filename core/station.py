from collections import namedtuple
from numpy import array
from core.utils import get_distance, nkg, poisson_gauss_gen, randomize_time
from core.amplitude import get_amplitude, get_sqr_sigma
from math import modf, sqrt

h_side = 0.8  # Половина стороны станции [м]
max_ampl = 1500.0  # Максимальная амплитуда пКл
side = 1.6  # Длина стороны станции [м]
det_area = 0.64  # Площадь детектора [м2]
# multiplicity_of_matches = 3  # Кратность совпадений


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, number, coordinates):

        self.num = number  # Нормер станции в кластере
        self.coord = array(coordinates)  # Координаты станции
        self.area = 2.56  # Площадь станции [м2]
        self.respond = None  # Отклик станции
        self.particles = None  # Экспериментальное число частиц
        self.real_time = None  # Истинное время срабатывания станции [нс]
        self.rndm_time = None  # Рандомизированное время станции [нс]
        self.rec_particles = None  # Восстановленное число частиц
        self.sigma_t = 5  # Ошибка определения времени [нс]
        self.amplitude = None  # Амплитуда сигнала от станции [пКл]
        self.sigma_particles = None  # Ошибка для функционала по частицам
        self.add_ampl = None  # Амплитуда сигнала от доп.ФЭУ

        # self.Det = namedtuple('Det', 'coord time ampl particles respond')

        # Детекторы станции
        self.detectors = (

            {
                'coord': (self.coord[0] - h_side, self.coord[1] + h_side, self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] + h_side, self.coord[1] + h_side, self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] + h_side, self.coord[1] - h_side, self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},

            {
                'coord': (self.coord[0] - h_side, self.coord[1] - h_side, self.coord[2]),
                'time': None,
                'ampl': None,
                'particles': None,
                'respond': None},
        )

        self.add_det = {
            'coord': self.detectors[self.num - 1]['coord'],
            'time': None,
            'ampl': None,
            'particles': None,
            'respond': None
        }

    def reset(self):
        """Возвращаемся к исходному состоянию"""
        self.respond = None
        self.particles = None
        self.real_time = None
        self.rndm_time = None
        self.rec_particles = None
        self.amplitude = None
        self.sigma_particles = None

        for i in range(4):
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

        # Обнуляем именованные кортежи детекторов
        # self.det1 = None
        # self.det2 = None
        # self.det3 = None
        # self.det4 = None
        # self.det5 = None

    def start(self, eas):
        """Запуск станции"""
        if self.model_ampl(eas):
            self.model_times(eas)

        # self.pack_namedtuples()
        return self.respond

    def model_ampl(self, eas):
        """Моделирование амплитуд в счётчиках станций"""
        _enabled_gen = True  # Включить/выключить генераторы Пуассона и Гаусса
        self.amplitude = 0.0
        self.sigma_particles = 0.0

        for det in self.detectors:
            dist = get_distance(det['coord'], eas.n, eas.x0, eas.y0)
            temp = det_area * eas.n[2] * nkg(dist, eas.power, eas.age)
            # p - временная переменная для хранения числа частиц
            # a - для амплитуды
            if _enabled_gen:
                p = poisson_gauss_gen(temp)
            else:
                p = temp

            if p == 0:
                # Детектор не сработал
                det['ampl'] = 0.0
                det['respond'] = False
                # self.respond = False
            else:
                # Детектор сработал
                det['respond'] = True
                a = 0.0
                for j in range(int(p)):
                    # Вычисляем амплитуду в детекторе
                    a += get_amplitude()
                # Добавим десятичную часть
                a += (get_amplitude() * modf(p)[0])

                if a > max_ampl:
                    a = max_ampl

                self.sigma_particles += sqrt(p * get_sqr_sigma())
                self.amplitude += a

                det['particles'] = p
                det['ampl'] = a

        # Дополнительный крадёт данные у ФЭУ-соседа
        self.add_det['ampl'] = self.detectors[self.num - 1]['ampl'] / 100
        self.add_det['respond'] = self.detectors[self.num - 1]['respond']

        if self.sigma_particles == 0:
            self.sigma_particles = 1.3

        if self.amplitude > 0.0:
            self.respond = True
        else:
            self.respond = False

        return self.respond

    def model_times(self, eas):
        """Моделирует времена срабатывания станций"""
        temp = []
        for det in self.detectors:
            if det['ampl'] != 0.0:
                det['time'] = (sum(eas.n * det['coord']) + eas.D) / eas.light_speed
                det['time'] += randomize_time()
                temp.append(det['time'])
            else:
                det['time'] = None

        self.rndm_time = min(temp)

        # Дополнительный крадёт время у соседа
        self.add_det['time'] = self.detectors[self.num - 1]['time']
        return self.respond

    # def pack_namedtuples(self):
    #     """Представим детекторы в виде именованных кортежей"""
    #     self.det1 = self.Det(**self.detectors[0])
    #     self.det2 = self.Det(**self.detectors[1])
    #     self.det3 = self.Det(**self.detectors[2])
    #     self.det4 = self.Det(**self.detectors[3])
    #
    #     self.det5 = self.Det(**self.add_det)
