from cluster import Cluster
from eas import Eas
from tools import get_theta

import random as rn

with open('data/efficiency/eff.txt', 'w') as eff_file:
    for i in range(100):
        en_pow = 2 + 0.05 * i
        energy = 10 ** en_pow

        eff_3 = 0  # Эффективность срабатывания третьего кластера
        eff_10 = 0  # десятого кластера
        eff_s = 0  # Установки

        exp_num = 5000  # Число испытаний для каждог значения энергии

        print(str(i) + '\t' + str(en_pow))

        for experiments in range(exp_num):
            theta = get_theta()
            phi = rn.uniform(0, 360)
            x0 = rn.uniform(-50, 50)
            y0 = rn.uniform(-50, 50)
            eas = Eas(theta, phi, x0, y0, energy)
            # print('\t' + str(experiments))
            clusters = [
                Cluster([-28.4, -7.8,  -7.0], eas, 13.3, 12.4),
                Cluster([-28.4, 23.8,  -7.0], eas, 13.3, 12.4),
                Cluster([    0,    0,     0], eas, 25.1, 13.5),
                Cluster([ 33.3,  7.8, -14.5], eas),
                Cluster([   35,   46, -14.5], eas),
                Cluster([   -2,  -47, -14.5], eas),
                Cluster([  -18,   62, -14.5], eas),
                Cluster([   50,   -2,    -2], eas),
                Cluster([   50,   26,    -2], eas),
                Cluster([   50,   58,    -8], eas),
            ]

            clust_ok = 0  # Считаем сработавшие кластеры
            for cluster in clusters:
                if cluster.model_amplitudes_simple():
                    clust_ok += 1

            eff_s += clust_ok

            if not clusters[2].failed:
                eff_3 += 1

            if not clusters[9].failed:
                eff_10 += 1

        eff_3 /= exp_num
        eff_10 /= exp_num
        eff_s /= exp_num

        eff_file.write(str(en_pow) + '\t' + str(eff_s) + '\t' + str(eff_3) +
                       '\t' + str(eff_10) + '\n')
