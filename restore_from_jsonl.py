from multiprocessing import Pool
import csv
import json
import time

from core import Facility, Eas

input_file = open('input_data/model_10k.jsonl', 'r')


def record_params(params_writer):
    """Запишем параметры ШАЛ"""
    for line in input_file:
        event = json.loads(line)
        params_writer.writerow(list(event['params'].values()))


def main_process(event_num):
    Nevod = Facility()
    for line_num, line in enumerate(input_file):
        if line_num == event_num:
            event = json.loads(line)
            eas = Eas(event['params'])
            Nevod.set_facility_state(event, eas)
            Nevod.rec_direction()
            Nevod.rec_particles()
            result = Nevod.rec_params_diff_evo()
            print(event_num)
            return [event_num] + result
        else:
            continue


def record(res, w_r):
    w_r.writerows(res)


if __name__ == '__main__':
    start_time = time.time()
    # Число считываемых событий
    event_num = 20
    # Файлы для записи
    params_file = open('output_data/params.txt', 'w')
    params_writer = csv.writer(params_file, dialect="excel-tab")
    restore_file = open('output_data/restore.txt', 'w')
    restore_writer = csv.writer(restore_file, dialect="excel-tab")

    record_params(params_writer)
    zoo = [i for i in range(event_num)]
    p = Pool(8)
    res = p.map(main_process, zoo)
    # Запись результоатов
    record(res, restore_writer)

    params_file.close()
    restore_file.close()
    input_file.close()

    print("--- %s seconds ---" % (time.time() - start_time))
