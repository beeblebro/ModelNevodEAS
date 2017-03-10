from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta

import random as rn

with open('data/efficiency/eff.txt', 'w') as eff_file:
    for i in range(100):
        en_pow = 3 + 0.05 * i
        energy = 10 ** en_pow

        eff_3 = 0  # Эффективность срабатывания третьего кластера
        eff_10 = 0  # десятого кластера
        eff_s = 0  # Установки

        exp_num = 500  # Число испытаний для каждог значения энергии

        print(str(i) + '\t' + str(en_pow))

        for experiments in range(exp_num):
            theta = get_theta()
            phi = rn.uniform(0, 360)
            x0 = rn.uniform(-50, 50)
            y0 = rn.uniform(-50, 50)
            eas = Eas(theta, phi, x0, y0, energy)
            # print('\t' + str(experiments))
            clusters = [
                Cluster([25, 25, 15], eas),
                Cluster([125, 25, 10], eas),
                Cluster([275, 25, 0], eas),
                Cluster([25, 125, 0], eas),
                Cluster([125, 125, 15], eas),
                Cluster([275, 125, 10], eas),
                Cluster([25, 275, 10], eas),
                Cluster([125, 275, 10], eas),
                Cluster([275, 275, 15], eas),
            ]

            clust_ok = 0  # Считаем сработавшие кластеры
            for cluster in clusters:
                if cluster.model_amplitudes_simple():
                    clust_ok += 1

            eff_s += clust_ok

            if not clusters[2].failed:
                eff_3 += 1

            if not clusters[8].failed:
                eff_10 += 1

        eff_3 /= exp_num
        eff_10 /= exp_num
        eff_s /= exp_num

        eff_file.write(str(en_pow) + '\t' + str(eff_s) + '\t' + str(eff_3) +
                       '\t' + str(eff_10) + '\n')
