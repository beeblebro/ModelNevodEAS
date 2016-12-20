from math import sqrt
import random as rn

from cluster import Cluster
from eas import Eas
from tools import get_theta, functional, divide_square


f = open('data/arrival_point.txt', 'w')
for experiments in range(1):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(0, 275)
    y0 = rn.uniform(0, 275)
    print(x0, y0)
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

    experimental_n = []  # Экспериментальное число частиц
    theoretical_n = []  # Теоритическое число частиц
    sigma_n = []  # Сигма в функционале

    average_n = [0, 0, 0]  # средний из восстановленных векторов
    average_x = 0  # Средневзвешанные x и y
    average_y = 0

    particles_sum = 0  # Сумма частиц по всем станциям

    clust_ok = 0  # Число сработавших кластеров
    for cluster in clusters:
        if cluster.least_squares():
            clust_ok += 1
            average_n += cluster.rec_n
            for i in range(4):

                particles_sum += cluster.stations[i].particles
                experimental_n.append(cluster.stations[i].particles)
                sigma_n.append(cluster.stations[i].sigma_particles)

                average_x += cluster.stations[i].coordinates[0] * \
                             cluster.stations[i].particles

                average_y += cluster.stations[i].coordinates[1] * \
                             cluster.stations[i].particles

    average_n /= clust_ok
    average_x /= particles_sum
    average_y /= particles_sum

    print(average_x, average_y)

    for cluster in clusters:
        if not cluster.failed:
            cluster.recons_particles(average_x, average_y, eas.n)
            for i in range(4):
                # Получаем теоретические значения числа частиц
                theoretical_n.append(cluster.stations[i].rec_particles)

    func = functional(experimental_n, theoretical_n, sigma_n)
    theoretical_n = []

    func_step1 = [[] for i in range(3)]  # Функционалы первого шага

    for x in [50.0, 150.0, 250.0]:
        for y in [50.0, 150.0, 250.0]:
            for cluster in clusters:
                if not cluster.failed:
                    cluster.recons_particles(x, y, average_n)
                    for k in range(4):
                        # Получаем теоретические значения числа частиц
                        theoretical_n.append(cluster.stations[k].rec_particles)
            func_step1[0].append(functional(experimental_n, theoretical_n, sigma_n))
            func_step1[1].append(x)
            func_step1[2].append(y)
            theoretical_n = []

    if min(func_step1[0]) < func:
        func = min(func_step1[0])
        # print(min(func_step1), func_step1.index(min(func_step1)))
        min_sqare_num = func_step1[0].index(min(func_step1[0]))
        new_x = func_step1[1][min_sqare_num]
        new_y = func_step1[2][min_sqare_num]
        # print(new_x, new_y)
    else:
        new_x = average_x
        new_y = average_y

    step2_cen = divide_square(new_x, new_y, 100)
    func_step2 = [[] for i in range(3)]  # Функционалы второго шага

    for x in step2_cen[0]:
        for y in step2_cen[1]:
            for cluster in clusters:
                if not cluster.failed:
                    # Получаем теоретические значения числа частиц
                    cluster.recons_particles(x, y, average_n)
                    for k in range(4):
                        theoretical_n.append(cluster.stations[k].rec_particles)
            func_step2[0].append(functional(experimental_n, theoretical_n, sigma_n))
            func_step2[1].append(x)
            func_step2[2].append(y)
            theoretical_n = []

    if min(func_step2[0]) < func:
        func = min(func_step2[0])
        # print(min(func_step1), func_step1.index(min(func_step1)))
        min_sqare_num = func_step2[0].index(min(func_step2[0]))
        new_x = func_step2[1][min_sqare_num]
        new_y = func_step2[2][min_sqare_num]
    print(new_x, new_y)

    delta = sqrt((x0 - average_x)**2 + (y0 - average_y)**2)
    f.write(str(delta))
    f.write('\n')

f.close()
