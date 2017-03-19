from math import sqrt
import random as rn

from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta, functional, make_step, count_theo, \
    draw_func_power, draw_func_age
from core.amplitude import get_av_amplitude, get_sqr_sigma

f = open('data/power_age_func/power_age.txt', 'w')
for experiments in range(1):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    x0 = rn.uniform(-50, 50)
    y0 = rn.uniform(-50, 50)
    power = 5 * 10 ** 6
    age = 1.4

    eas = Eas(theta, phi, x0, y0, power, age)

    # Предполагаемые изначально мощность и возраст
    start_power = 10 ** 4
    start_age = 1.3

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
    # theo_n = []  # Теоритическое число частиц
    sigma_n = []  # Сигма в функционале

    average_n = [0, 0, 0]  # средний из восстановленных векторов
    average_x = 0  # Средневзвешанные x и y
    average_y = 0

    particles_sum = 0  # Сумма частиц по всем станциям

    clust_ok = 0  # Число сработавших кластеров
    # Запускаем кластеры, считаем средний вектор
    for cl in clusters:
        if cl.start():
            clust_ok += 1
            average_n += cl.rec_n

    if clust_ok == 0:
        # Не сработал ни один кластер
        print("ERROR: Не сработал ни один кластер")
        continue

    # Восстановили вектор прихода ШАЛ
    average_n /= clust_ok
    # Среднняя амплитуда, скорректрованная на толщину
    fixed_av_amplitude = get_av_amplitude() / average_n[2]

    for cl in clusters:
        for st in cl.stations:
            # Считаем экспериментальные частицы в каждой станции
            st.particles = st.amplitude / fixed_av_amplitude

            if st.sigma_particles < 0:
                print("ERROR: Отрицательное число частиц в основном цикле")

            #  Считаем сигмы
            st.sigma_particles = sqrt(st.particles * get_sqr_sigma())
            if st.sigma_particles == 0:
                st.sigma_particles = 1.3

            exp_n.append(st.particles)
            sigma_n.append(st.sigma_particles)

            particles_sum += st.particles
            average_x += st.coord[0] * st.particles
            average_y += st.coord[1] * st.particles

    average_x /= particles_sum
    average_y /= particles_sum

    # Защитим экспериментальные данные от изменений
    exp_n = tuple(exp_n)
    sigma_n = tuple(sigma_n)

    draw_func_power(clusters, eas.n, x0, y0, start_power, age, exp_n, sigma_n)
    draw_func_age(clusters, eas.n, x0, y0, power, start_age - 0.1, exp_n, sigma_n)

    theo_n = count_theo(clusters, average_n, average_x, average_y, start_power,
                        start_age)

    func = functional(exp_n, sigma_n, theo_n)

    step_1 = make_step(clusters, average_n, 100, average_x, average_y, start_power, start_age,
                       exp_n, sigma_n, func)

    # if step_1['x'] == 0 and step_1['y'] == 0:
    #     step_1['x'] = average_x
    #     step_1['y'] = average_y

    step_2 = make_step(clusters, average_n, 100 / 3, step_1['x'], step_1['y'],
                       step_1['power'], step_1['age'], exp_n, sigma_n,
                       step_1['func'])

    step_3 = make_step(clusters, average_n, 100 / 9, step_2['x'], step_2['y'],
                       step_2['power'], step_2['age'], exp_n, sigma_n,
                       step_2['func'])

    delta = sqrt((step_3['x'] - x0)**2 + (step_3['y'] - y0)**2)
    delta_age = step_3['age'] - age
    delta_power = step_3['power'] - power

    print(experiments)
    f.write(str(delta) + '\t' + str(delta_power) + '\t' + str(delta_age))
    f.write('\n')

f.close()

