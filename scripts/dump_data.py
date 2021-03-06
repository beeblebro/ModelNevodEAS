import struct
import json
import random as rn

from core import Facility
from core import Eas
from core import gen_age, gen_theta, gen_power, gen_x_y
from core import get_ca_amplitude


def create_det():
    """Структура для записи данных детектора"""
    return {
        'ampl': None,
        'time': None
    }


def create_station(dets):
    """Структура для записи данных станции"""
    return {
        'x': None,
        'y': None,
        'z': None,
        'ca_ampl': None,
        'ampl': None,
        'time': None,
        'add_det_ampl': None,
        'add_det_time': None,
        'det': dets
    }


def create_cluster(stations):
    """Структура для записи данных кластера"""
    return {
        'st': stations
    }


def create_params():
    """Задаём параметры ливня"""
    theta = gen_theta()  # Тета в градусах
    phi = rn.uniform(0, 360)  # и фи в градусах
    x0, y0 = gen_x_y()
    power = gen_power()
    age = gen_age(power, theta)
    return {
        'theta': theta,
        'phi': phi,
        'x0': x0,
        'y0': y0,
        'power': power,
        'age': age
    }


def run_facility(nevod_eas, params):
    """Запускаем установку и создаём ШАЛ"""
    NevodEAS = nevod_eas
    # Создали ШАЛ
    eas = Eas(params)
    # Передаём ШАЛ и запускаем установку
    NevodEAS.start(eas)
    return NevodEAS.clusters


def save_to_json(file, event_json):
    """Сохраняем событие в файл формата jsonl"""
    json.dump(event_json, file)
    file.write('\n')


def save_to_bin(file, event_json):
    """Сохраняем в бианрный файл"""
    # Счётчик байтов
    b_count = 0
    # Записали номер события
    b_count += file.write(struct.pack('q', event_json['num']))
    # Записали параметры ШАЛ
    b_count += file.write(struct.pack('dddddd',
                                      event_json['params']['theta'],
                                      event_json['params']['phi'],
                                      event_json['params']['x0'],
                                      event_json['params']['y0'],
                                      event_json['params']['power'],
                                      event_json['params']['age']))

    # Меняем None на -1.0 для времени, чтобы записать в бинарник
    for cl in event_json['clusters']:

        for st in cl['st']:

            b_count += file.write(struct.pack('ddd', st['x'], st['y'],
                                  st['z']))

            b_count += file.write(struct.pack('d', st['ca_ampl']))
            b_count += file.write(struct.pack('d', st['ampl']))
            if st['time'] is None:
                b_count += file.write(struct.pack('d', -1.0))
            else:
                b_count += file.write(struct.pack('d', st['time']))

            for det in st['det']:
                b_count += file.write(struct.pack('d', det['ampl']))
                if det['time'] is None:
                    b_count += file.write(struct.pack('d', -1.0))
                else:
                    b_count += file.write(struct.pack('d', det['time']))

            b_count += file.write(struct.pack('d', st['add_det_ampl']))
            if st['add_det_time'] is None:
                b_count += file.write(struct.pack('d', -1.0))
            else:
                b_count += file.write(struct.pack('d', st['add_det_time']))

    # print(b_count)
    return b_count


def create_event(clusters, count_event, params):
    """Создаём и заполняем стурктуру события"""
    clusters_json = []
    for i in range(len(clusters)):
        stations = []
        for j in range(4):
            dets = []
            for k in range(4):
                dets.append(create_det())
            stations.append(create_station(dets))
        clusters_json.append(create_cluster(stations))

    for cl_n, cl in enumerate(clusters):

        for st_n, st in enumerate(cl.stations):
            clusters_json[cl_n]['st'][st_n]['x'] = st.coord[0]
            clusters_json[cl_n]['st'][st_n]['y'] = st.coord[1]
            clusters_json[cl_n]['st'][st_n]['z'] = st.coord[2]

            clusters_json[cl_n]['st'][st_n]['ca_ampl'] = get_ca_amplitude()
            clusters_json[cl_n]['st'][st_n]['ampl'] = st.amplitude
            clusters_json[cl_n]['st'][st_n]['time'] = st.rndm_time

            clusters_json[cl_n]['st'][st_n]['add_det_ampl'] = st.add_det['ampl']
            clusters_json[cl_n]['st'][st_n]['add_det_time'] = st.add_det['time']

            for det_n, det in enumerate(st.detectors):
                clusters_json[cl_n]['st'][st_n]['det'][det_n]['ampl'] = det['ampl']
                clusters_json[cl_n]['st'][st_n]['det'][det_n]['time'] = det['time']

    return {
        # Номер события
        'num': count_event,
        # Параметры ШАЛ
        'params': params,
        # Данные срабатывания счётчиков
        'clusters': clusters_json
    }


if __name__ == '__main__':
    nevod_eas = Facility(configuration='nevod')  # Создали НЕВОД
    flat_eas = Facility(configuration='flat')  # Создали плоскую установку
    count_event = 10000

    # Файлы для установки НЕВОД
    file_bin = open('../input_data/model_10k.bin', 'wb')
    file_json = open('../input_data/model_10k.jsonl', 'w')
    # Файлы для плоской установки
    file_json_flat = open('../input_data/model_10k_flat.jsonl', 'w')
    file_bin_flat = open('../input_data/model_10k_flat.bin', 'wb')

    for event_number in range(count_event):
        params = create_params()
        # Запускаем установки
        clusters = run_facility(nevod_eas, params)
        clusters_flat = run_facility(flat_eas, params)

        event = create_event(clusters, event_number, params)
        event_flat = create_event(clusters_flat, event_number, params)

        save_to_bin(file_bin, event)
        save_to_json(file_json, event)

        save_to_bin(file_bin_flat, event_flat)
        save_to_json(file_json_flat, event_flat)

        print(event_number)
        nevod_eas.reset()
        flat_eas.reset()

    file_bin.close()
    file_json.close()
    file_bin_flat.close()
    file_json_flat.close()
