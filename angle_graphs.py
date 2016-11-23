from cluster import Cluster
from eas import Eas
from numpy import array

eas = Eas(5, 45, 100, 100)
print(eas.n)
print("")

average = array([0.0, 0.0, 0.0])

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

clust_ok = 0
for cluster in clusters:
    if cluster.least_squares():
        print(cluster.st_amplitudes)
        average += cluster.rec_n
        clust_ok += 1
average /= clust_ok
print("")
print(average)
