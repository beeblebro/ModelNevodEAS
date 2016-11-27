from cluster import Cluster
from eas import Eas
from numpy import *
import random as rn
from tools import get_theta, functional

for experiments in range(0, 1):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(0, 275)
    y0 = rn.uniform(0, 275)
    print(theta, phi, x0, y0)
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

    experimental_n = []  # Экспериментальное число частиц
    theoretical_n = []  # Теоритическое число частиц

    average_n = [0, 0, 0]  # средний из восстановленных векторов
    average_x = 0  # Средневзвешанные x и y
    average_y = 0

    particles_sum = 0  # Сумма частиц по всем станциям

    clust_ok = 0  # Число сработавших кластеров
    for cluster in clusters:
        if cluster.least_squares():
            clust_ok += 1
            average_n += cluster.rec_n
            # cluster.coutn_particles()
            particles_sum += sum(cluster.st_particles)
            for i in range(0, 4):
                experimental_n.append(cluster.st_particles[i])
                average_x += cluster.st_coord[i][0] * cluster.st_particles[i]
                average_y += cluster.st_coord[i][1] * cluster.st_particles[i]
    average_n /= clust_ok
    average_x /= particles_sum
    average_y /= particles_sum

    print(average_x, average_y)

    for cluster in clusters:
        if not cluster.failed:
            cluster.recons_particles(average_x, average_y, average_n)  # Получаем теоретические значения числа частиц
            for i in range(0, 4):
                theoretical_n.append(cluster.rec_particles[i])

    func = functional(experimental_n, theoretical_n)
    theoretical_n = []

    func_step1 = [[] for i in range(3)]  # Функционалы первого шага

    for x in [50, 150, 250]:
        for y in [50, 150, 250]:
            for cluster in clusters:
                if not cluster.failed:
                    cluster.recons_particles(x, y, average_n)  # Получаем теоретические значения числа частиц
                    for k in range(0, 4):
                        theoretical_n.append(cluster.rec_particles[k])
            func_step1[0].append(functional(experimental_n, theoretical_n))
            func_step1[1].append(x)
            func_step1[2].append(y)
            theoretical_n = []

    if min(func_step1[0]) < func:
        # print(min(func_step1), func_step1.index(min(func_step1)))
        min_sqare_num = func_step1[0].index(min(func_step1[0]))
        new_x = func_step1[1][min_sqare_num]
        new_y = func_step1[2][min_sqare_num]
        print(new_x, new_y)

