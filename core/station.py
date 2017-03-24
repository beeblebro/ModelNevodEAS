from numpy import array


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, coord):
        self.coord = array(coord)
        self.side = 1.6  # Длина стороны станции [м]
        self.area = 2.56  # Площадь станции [м2]
        self.respond = None  # Отклик станции
        self.particles = None  # Экспериментальное число частиц, попавших в стнцию
        self.real_time = None  # Истинное время срабатывания станции [нс]
        self.rndm_time = None  # Рандомизированное время срабатывания станции [нс]
        self.rec_particles = None  # Восстановленное число частиц
        self.sigma_t = 5  # Ошибка определения времени [нс]
        self.amplitude = None  # Амплитуда сигнала от станции [пКл]
        self.sigma_particles = None  # Ошибка для функционала по частицам
        self.dyn_range = 800  # Динамический диапазон

    def reset(self):
        """Возвращаемся к исходному состоянию"""
        self.respond = None
        self.particles = None
        self.real_time = None
        self.rndm_time = None
        self.rec_particles = None
        self.amplitude = None
        self.sigma_particles = None