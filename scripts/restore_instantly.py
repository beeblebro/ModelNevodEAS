# Смоделировать и восстановить ШАЛ сразу же, без сохранения события
from multiprocessing import Pool
import csv
import time

from core import Facility, Eas


def main_process(event_num):
    eas = Eas()
    params = eas.get_params_list()
    nevod = Facility()
    nevod.start(eas)
    nevod.rec_direction()
    nevod.rec_particles()
    result = nevod.rec_params_diff_evo()
    # print(event_num)
    return [event_num] + params, [event_num] + result


def record(res, w_p, w_r):
    """Запись заданных и восстановленных параметров"""
    for evt in res:
        w_p.writerow(evt[0])
        w_r.writerow(evt[1])


if __name__ == '__main__':
    start_time = time.time()
    # Число моделируемых событий
    event_num = 20
    # Файлы для записи
    params_file = open('../output_data/params.txt', 'w')
    params_writer = csv.writer(params_file, dialect="excel-tab")
    restore_file = open('../output_data/restore.txt', 'w')
    restore_writer = csv.writer(restore_file, dialect="excel-tab")

    zoo = [i for i in range(event_num)]
    p = Pool(8)
    res = p.map(main_process, zoo)
    # Запись результатов
    record(res, params_writer, restore_writer)
    params_file.close()
    restore_file.close()
    print("--- %s seconds ---" % (time.time() - start_time))
