"""Хотим посмотреть отношение экспериментального числа частиц к расчётному, но с правильной осью ШАЛ"""

from core.cluster import Cluster
from core.eas import Eas
from numpy import *
import random as rn
from core.tools import get_theta, functional, divide_square


f = open('data/amplitude_test.txt', 'w')
for experiments in range(1000):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(0, 275)
    y0 = rn.uniform(0, 275)

    eas = Eas(theta, phi, x0, y0)

    clusters = [
        Cluster([-28.4, -7.8, -7.0], eas, 13.3, 12.4),
        Cluster([-28.4, 23.8, -7.0], eas, 13.3, 12.4),
        Cluster([0, 0, 0], eas, 25.1, 13.5),
        Cluster([33.3, 7.8, -14.5], eas),
        Cluster([35, 46, -14.5], eas),
        Cluster([-2, -47, -14.5], eas),
        Cluster([-18, 62, -14.5], eas),
        Cluster([50, -2, -2], eas),
        Cluster([50, 26, -2], eas),
        Cluster([50, 58, -8], eas),
    ]

    for cluster in clusters:
        if cluster.model_amplitudes():
            cluster.recons_particles(x0, y0, eas.n)
            for station in cluster.stations:
                if station.particles != 0:
                    test = station.particles / station.rec_particles
                    print(test)
                    f.write(str(test))
                    f.write('\n')

f.close()
