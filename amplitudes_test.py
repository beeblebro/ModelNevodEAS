"""Хотим посмотреть отношение экспериментального числа частиц к расчётному, но с правильной осью ШАЛ"""

from cluster import Cluster
from eas import Eas
from numpy import *
import random as rn
from tools import get_theta, functional, divide_square


f = open('data/amplitude_test.txt', 'w')
for experiments in range(1000):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(0, 275)
    y0 = rn.uniform(0, 275)

    eas = Eas(theta, phi, x0, y0)

    clusters = [
        Cluster([ 25,  25, 15], eas),
        Cluster([125,  25, 10], eas),
        Cluster([275,  25,  0], eas),
        Cluster([ 25, 125,  0], eas),
        Cluster([125, 125, 15], eas),
        Cluster([275, 125, 10], eas),
        Cluster([ 25, 275, 10], eas),
        Cluster([125, 275, 10], eas),
        Cluster([275, 275, 15], eas),
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
