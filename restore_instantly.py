# Смоделировать и восстановить ШАЛ сразу же, без сохранения события
from multiprocessing import Pool
import csv
import random as rn
import time

from core import Facility
from core import Eas
from core.utils import get_theta, get_power, get_age, modified_area


def run(facility):
    theta = get_theta()  # Тета
    phi = rn.uniform(0, 360)  # и фи в градусах
    x0, y0 = modified_area()
    power = get_power()
    age = get_age(power, theta)
    params = [theta, phi, x0, y0, power, age]
    eas = Eas(theta, phi, x0, y0, power, age)
    facility.start(eas)
    facility.rec_direction()
    facility.rec_particles()
    result = facility.rec_params_diff_evo()
    return params, result


def record(res, w_p, w_r):
    for evt in res:
        w_p.writerow(evt[0])
        w_r.writerow(evt[1])


if __name__ == '__main__':
    start_time = time.time()
    f_params = open('params.txt', 'w')
    w_params = csv.writer(f_params, dialect="excel-tab")
    f_restore = open('restore.txt', 'w')
    w_restore = csv.writer(f_restore, dialect="excel-tab")
    zoo = [Facility() for i in range(100)]
    p = Pool(8)
    res = p.map(run, zoo)
    record(res, w_params, w_restore)
    f_params.close()
    f_restore.close()
    print("--- %s seconds ---" % (time.time() - start_time))
