from math import sqrt
import random as rn

from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta, functional, make_step
from core.amplitude import get_av_amplitude, get_sqr_sigma


f = open('data/arrival_point/dist_from_average.txt', 'w')
for experiments in range(5):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(-50, 50)
    y0 = rn.uniform(-50, 50)
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

    exp_n = []  # Экспериментальное число частиц
    theo_n = []  # Теоритическое число частиц
    sigma_n = []  # Сигма в функционале

    average_n = [0, 0, 0]  # средний из восстановленных векторов
    average_x = 0  # Средневзвешанные x и y
    average_y = 0

    particles_sum = 0  # Сумма частиц по всем станциям

    clust_ok = 0  # Число сработавших кластеров
    # Запускаем кластеры, считаем средний вектор
    for cluster in clusters:
        if cluster.rec_direction():
            clust_ok += 1
            average_n += cluster.rec_n

    if clust_ok == 0:
        # Не сработал ни один кластер
        print("GLOBAL FAIL!")
        continue

    # Восстановили вектор прихода ШАЛ
    average_n /= clust_ok
    # Косиунус восстановленного зенитного угла
    cos_rec_theta = average_n[2]
    # Среднняя амплитуда, скорректрованная на толщину
    fixed_av_amplitude = get_av_amplitude() / average_n[2]

    for cluster in clusters:
        for i in range(4):
            # Считаем экспериментальные частицы в каждой станции
            cluster.stations[i].particles = cluster.stations[i].amplitude /\
                                            fixed_av_amplitude
            #  Считаем сигмы
            cluster.stations[i].sigma_particles = sqrt(cluster.stations[i].particles * get_sqr_sigma())
            if cluster.stations[i].sigma_particles == 0:
                cluster.stations[i].sigma_particles = 1.3

            exp_n.append(cluster.stations[i].particles)
            sigma_n.append(cluster.stations[i].sigma_particles)

            particles_sum += cluster.stations[i].particles
            average_x += cluster.stations[i].coordinates[0] * \
                         cluster.stations[i].particles
            average_y += cluster.stations[i].coordinates[1] * \
                         cluster.stations[i].particles

    average_x /= particles_sum
    average_y /= particles_sum

    # Защитим экспериментальные данные от изменений
    exp_n = tuple(exp_n)
    sigma_n = tuple(sigma_n)

    for cluster in clusters:
        cluster.rec_particles(average_n, average_x, average_y, pow(10, 4), 1.0)
        for i in range(4):
            # Получаем теоретические значения числа частиц
            theo_n.append(cluster.stations[i].rec_particles)

    func = functional(exp_n, sigma_n, theo_n)
    # print(func)

    step_1 = make_step(clusters, average_n, 100, 0, 0, exp_n, sigma_n, func)

    if step_1[0] == 0 and step_1[1] == 0:
        step_1[0] = average_x
        step_1[1] = average_y

    step_2 = make_step(clusters, average_n, 100 / 3, step_1[0], step_1[1],
                       exp_n, sigma_n, step_1[2])

    step_3 = make_step(clusters, average_n, 100 / 9, step_2[0], step_2[1],
                       exp_n, sigma_n, step_2[2])

    new_x = step_3[0]
    new_y = step_3[1]

    result_energy = step_3[3]
    result_age = step_3[4]

    # print(x0, y0)
    # print(average_x, average_y)
    # print(new_x, new_y)

    print(experiments)
    #  delta = sqrt((x0 - new_x)**2 + (y0 - new_y)**2)
    #  delta_av = sqrt((x0 - average_x)**2 + (y0 - average_y)**2)
    f.write(str(result_energy) + '\t' + str(result_age))
    f.write('\n')

f.close()

