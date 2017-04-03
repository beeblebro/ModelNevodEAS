from numpy import array


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, coordinates):

        # h_side = 0.8  # Половина стороны станции [м]

        self.coord = array(coordinates)
        self.side = 1.6  # Длина стороны станции [м]
        self.area = 2.56  # Площадь станции [м2]
        self.respond = None  # Отклик станции
        self.particles = None  # Экспериментальное число частиц
        self.real_time = None  # Истинное время срабатывания станции [нс]
        self.rndm_time = None  # Рандомизированное время станции [нс]
        self.rec_particles = None  # Восстановленное число частиц
        self.sigma_t = 5  # Ошибка определения времени [нс]
        self.amplitude = None  # Амплитуда сигнала от станции [пКл]
        self.sigma_particles = None  # Ошибка для функционала по частицам
        self.dyn_range = 800  # Динамический диапазон

        # Детекторы станции TODO
        # self.detectors = (
        #
        #     {
        #         'coord': (self.coord[0] - h_side, self.coord[1] - h_side, self.coord[2]),
        #         'time': None,
        #         'ampl': None},
        #
        #     {
        #         'coord': (self.coord[0] - h_side, self.coord[1] + h_side, self.coord[2]),
        #         'time': None,
        #         'ampl': None},
        #
        #     {
        #         'coord': (self.coord[0] + h_side, self.coord[1] - h_side, self.coord[2]),
        #         'time': None,
        #         'ampl': None},
        #
        #     {
        #         'coord': (self.coord[0] + h_side, self.coord[1] + h_side, self.coord[2]),
        #         'time': None,
        #         'ampl': None},
        #
        #     {   # Дополнительный ФЭУ
        #         'coord': (self.coord[0] + h_side, self.coord[1] + h_side, self.coord[2]),
        #         'time': None,
        #         'ampl': None,
        #         'additional': True}
        # )

    def reset(self):
        """Возвращаемся к исходному состоянию"""
        self.respond = None
        self.particles = None
        self.real_time = None
        self.rndm_time = None
        self.rec_particles = None
        self.amplitude = None
        self.sigma_particles = None
