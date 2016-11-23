from cluster import Cluster
from eas import Eas
from numpy import array
import random as rn
from tools import get_theta

for experiments in range(0, 1000):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(0, 300)
    y0 = rn.uniform(0, 300)
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
        Cluster([275, 275, 15], eas)
    ]

    for cluster in clusters:
        if cluster.least_squares():
            print(cluster.rec_n)

