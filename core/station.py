from numpy import array


class Station:
    """Класс для представления станции установки НЕВОД-ШАЛ"""

    def __init__(self, coordinates):
        self.coordinates = array(coordinates)
        self.side = 1.6  # Длина стороны станции [м]
        self.area = 2.56  # Площадь станции [м2]
        self.failed = False  # Параметр отказа станции
        self.particles = 0  # Экспериментальное число частиц, попавших в стнцию
        self.real_time = 0  # Истинное время срабатывания станции [нс]
        self.rndm_time = 0  # Рандомизированное время срабатывания станции [нс]
        self.rec_particles = 0  # Восстановленное число частиц
        self.sigma_t = 5  # Ошибка определения времени [нс]
        self.amplitude = 0  # Амплитуда сигнала от станции [пКл]
        self.sigma_particles = 0  # Ошибка для функционала по частицам
        self.dyn_range = 800  # Динамический диапазон