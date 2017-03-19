from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta
from numpy import array

import random as rn


with open('data/time/delta_time.txt', 'w') as time_file:
    for experiments in range(1):
        theta = get_theta()
        phi = rn.uniform(0, 360)
        x0 = rn.uniform(-50, 50)
        y0 = rn.uniform(-50, 50)
        eas = Eas(theta, phi, x0, y0)

        clust_times = []  # Времена срабатывания кластеров

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

        clust_ok = 0
        for cluster in clusters:
            if cluster.rec_direction():
                clust_times.append(cluster.time)
                clust_ok += 1

        if clust_ok == 0:
            continue

        min_time = min(clust_times)

        clust_times = array(clust_times)
        clust_times -= min_time

        delta_t = max(clust_times) - min(clust_times)
        print(clust_times)
        print(experiments)
        # print(delta_t)

        time_file.write(str(theta) + '\t' + str(phi) + '\t' + str(delta_t)
                        + '\n')
