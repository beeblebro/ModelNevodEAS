from math import sqrt
import random as rn

from core.cluster import Cluster
from core.eas import Eas
from core.tools import get_theta, functional, make_step, count_theo, \
    draw_func_power, draw_func_age
from core.amplitude import get_av_amplitude, get_sqr_sigma

f = open('data/power_age_func/power_age.txt', 'w')
for experiments in range(5):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    # x0 = rn.uniform(-50, 50)
    # y0 = rn.uniform(-50, 50)
    x0 = 25
    y0 = 25
    power = 10 ** 6
    age = 1.45

    eas = Eas(theta, phi, x0, y0, power, age)

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
    average_n = eas.n
    # Среднняя амплитуда, скорректрованная на толщину
    fixed_av_amplitude = get_av_amplitude() / average_n[2]

    for cl in clusters:
        for st in cl.stations:
            # Считаем экспериментальные частицы в каждой станции
            # st.particles = st.amplitude / fixed_av_amplitude

            if st.sigma_particles < 0:
                print("ERROR: Отрицательное число частиц в основном цикле")

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

    # draw_func_power(clusters, eas.n, x0, y0, start_power, age, exp_n, sigma_n)
    # draw_func_age(clusters, eas.n, x0, y0, power, start_age - 0.1, exp_n, sigma_n)

    _x = 0
    _y = 0
    _power = 10 ** 6
    _age = 1.3
    theo_n = count_theo(clusters, average_n, _x, _y, _power, _age)

    _func = functional(exp_n, sigma_n, theo_n)
    _side = 100

    for steps in range(3):
        step = make_step(clusters, average_n, _side, _x, _y, _power, _age, exp_n
                         , sigma_n, _func)
        _x = step['x']
        _y = step['y']
        _power = step['power']
        _age = step['age']
        _func = step['func']

        _side /= 3

    delta = sqrt((_x - x0)**2 + (_y - y0)**2)
    delta_age = _age - age
    delta_power = _power - power

    print(experiments)
    f.write(str(_power) + '\t' + str(_age))
    f.write('\n')

f.close()

