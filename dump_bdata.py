import struct
import json
import random as rn

from core import Facility
from core import Eas
from core import get_age, get_theta, get_power

fbin = open('modelEAS_flat_10k_11.04.bin', 'wb')
fjson = open('modelEAS_flat_10k_11.04.jsonl', 'w')
fp = open('distr.txt', 'w')

NevodEAS = Facility()
for tries in range(10000):
    theta = get_theta()  # Тета в градусах
    phi = rn.uniform(0, 360)  # и фи в градусах
    x0 = rn.uniform(-50, 50)
    y0 = rn.uniform(-50, 50)
    power = get_power()
    age = get_age(power, theta)
    # Создали ШАЛ
    eas = Eas(theta, phi, x0, y0, power, age)

    # Запишем параметры ШАЛ, чтобы построить распределения
    fp.write(str(theta) + '\t' + str(phi) + '\t'
             + str(x0) + '\t' + str(y0) + '\t'
             + str(power) + '\t' + str(age) + '\n')

    # Установка получает ливень
    NevodEAS.get_eas(eas)
    # Запускаем установку
    NevodEAS.start()

    dets = {
        'ampl': None,
        'time': None
    }

    station = {
        'ampl': None,
        'time': None,
        'add_det_ampl': None,
        'add_det_time': None,
        'det': [dets for i in range(4)]
    }

    cluster = {
        'st': [station for j in range(4)]
    }

    clusters = [cluster for k in range(len(NevodEAS.clusters))]

    for cl_n, cl in enumerate(NevodEAS.clusters):

        for st_n, st in enumerate(cl.stations):
            clusters[cl_n]['st'][st_n]['ampl'] = st.amplitude
            clusters[cl_n]['st'][st_n]['time'] = st.rndm_time

            clusters[cl_n]['st'][st_n]['add_det_ampl'] = st.add_det['ampl']
            clusters[cl_n]['st'][st_n]['add_det_time'] = st.add_det['time']

            for det_n, det in enumerate(st.detectors):
                clusters[cl_n]['st'][st_n]['det'][det_n]['ampl'] = det['ampl']
                clusters[cl_n]['st'][st_n]['det'][det_n]['time'] = det['time']

    data = {
        # Номер события
        'num': tries,
        # Параметры ШАЛ
        'params': {
            'theta': theta,
            'phi': phi,
            'x0': x0,
            'y0': y0,
            'power': power,
            'age': age
        },
        # Данные срабатывания счётчиков
        'cluster': clusters
    }

    json.dump(data, fjson)
    fjson.write('\n')

    # Заменим None на -1, чтобы записать в бинарник
    for cl in NevodEAS.clusters:

        for st in cl.stations:
            if st.rndm_time is None:
                st.rndm_time = -1.0
            if st.amplitude is None:
                st.amplitude = 0.0

            if st.add_det['time'] is None:
                st.add_det['time'] = -1.0

            if st.add_det['ampl'] is None:
                st.add_det['ampl'] = 0.0

            for det in st.detectors:
                if det['time'] is None:
                    det['time'] = -1.0

                if det['ampl'] is None:
                    det['ampl'] = 0.0

    # Счётчик байтов
    b_count = 0
    # Записали номер события
    b_count += fbin.write(struct.pack('L', tries))
    # Записали параметры ШАЛ
    b_count += fbin.write(struct.pack('dddddd', theta, phi, x0, y0, power, age))
    for cl in NevodEAS.clusters:
        for st in cl.stations:
            # Время и амплитуда срабатывания станции
            b_count += fbin.write(struct.pack('dd', st.rndm_time, st.amplitude))

            for det in st.detectors:
                # Время и амплитуда детекторов
                b_count += fbin.write(struct.pack('dd', det['time'], det['ampl']))

            # Время и амплитуда дополнительного ФЭУ
            b_count += fbin.write(struct.pack('dd', st.add_det['time'], st.add_det['ampl']))

    print(tries)
    NevodEAS.reset()

fp.close()
fbin.close()
fjson.close()
