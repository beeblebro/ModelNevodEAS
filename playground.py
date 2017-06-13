from multiprocessing import Pool
import csv
import json
import time

from core import Facility


def set_state(zoo, f_data, w_p):
    """Установим отклики для всех установок и запишем параметры"""

    for line, facility in zip(f_data, zoo):
        event = json.loads(line)

        theta = event['params']['theta']
        phi = event['params']['phi']
        x0 = event['params']['x0']
        y0 = event['params']['y0']
        power = event['params']['power']
        age = event['params']['age']
        params = [facility.num, theta, phi, x0, y0, power, age]
        w_p.writerow(params)
        facility.set_facility_state(event)


def run(facility):
    facility.rec_direction()
    facility.rec_particles()
    result = facility.rec_params_diff_evo()
    return result


def record(res, w_r):
    w_r.writerows(res)


if __name__ == '__main__':
    start_time = time.time()

    f_data = open('data/model_10k.jsonl', 'r')
    f_params = open('params.txt', 'w')
    w_params = csv.writer(f_params, dialect="excel-tab")
    f_restore = open('restore.txt', 'w')
    w_restore = csv.writer(f_restore, dialect="excel-tab")

    zoo = [Facility(geometry='nevod', num=i) for i in range(10)]
    set_state(zoo, f_data, w_params)
    p = Pool(8)
    res = p.map(run, zoo)
    record(res, w_restore)

    f_data.close()
    f_params.close()
    f_restore.close()

    print("--- %s seconds ---" % (time.time() - start_time))
