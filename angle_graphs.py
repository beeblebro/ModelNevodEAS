from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta

from math import sqrt, pi, acos
from numpy import array, square
from numpy.linalg import norm
import random as rn

psi_list = [[], []]
sigma_psi_list = [[], []]

with open('data/angle_parameters/angles.txt', 'w') as angle_file:
    for experiments in range(100000):
        theta = get_theta()
        phi = rn.uniform(0, 360)
        x0 = rn.uniform(-50, 50)
        y0 = rn.uniform(-50, 50)
        eas = Eas(theta, phi, x0, y0)

        average = array([0.0, 0.0, 0.0])  # Средний из восстановленных векторов
        average_sqr = 0.0  # Среднее квадратов

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
            if cluster.least_squares():
                average += cluster.rec_n
                average_sqr += sum(square(cluster.rec_n))
                clust_ok += 1

        if clust_ok == 0:
            continue

        average /= clust_ok
        average_sqr /= clust_ok  # Среднее квадратов
        sqr_average = sum(square(average))  # Квадрат среднего вектора

        delta = norm(average - eas.n)  # Модуль разности среднего и истинного
        # Среднеквадратичный разброс векторов
        sigma_n = sqrt(average_sqr - sqr_average)

        # Переведём в градусную меру
        psi = (180 / pi) * acos(1 - 0.5 * pow(delta, 2))
        sigma_psi = (180 / pi) * acos(1 - 0.5 * pow(sigma_n, 2))

        psi_list[0].append(theta)
        psi_list[1].append(psi)

        sigma_psi_list[0].append(theta)
        sigma_psi_list[1].append(sigma_psi)

        print(experiments)
        # print(psi)
        # print(sigma_psi)

        angle_file.write(str(theta) + '\t' + str(phi) + '\t'
                         + str(psi) + '\t' + str(sigma_psi) + '\n')

with open('data/angle_parameters/psi_theta.txt', 'w') as psi_theta_file:
    for th in range(20):
        current_psi = []
        av = 0
        sig = 0
        count = 0
        for i in range(len(psi_list[1])):
            if psi_list[0][i] <= th * 5:
                psi_list[0][i] = 999
                av += psi_list[1][i]
                count += 1
                current_psi.append(psi_list[1][i])

        if count == 0:
            continue

        av /= count
        for j in range(len(current_psi)):
            sig += (av - current_psi[j])**2
        sig /= count
        psi_theta_file.write(str(th * 5) + '\t' + str(av) + '\t' + str(sig) + '\n')

with open('data/angle_parameters/sigma_psi_theta.txt', 'w') as sigma_psi_theta_file:
    for th in range(20):
        current_sig_psi = []
        av = 0
        sig = 0
        count = 0
        for i in range(len(sigma_psi_list[1])):
            if sigma_psi_list[0][i] <= th * 5:
                sigma_psi_list[0][i] = 999
                av += sigma_psi_list[1][i]
                count += 1
                current_sig_psi.append(sigma_psi_list[1][i])

        if count == 0:
            continue

        av /= count
        for j in range(len(current_sig_psi)):
            sig += (av - current_sig_psi[j])**2
        sig /= count
        sigma_psi_theta_file.write(str(th * 5) + '\t' + str(av) + '\t' + str(sig) + '\n')
