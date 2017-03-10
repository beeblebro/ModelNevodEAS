from core.amplitude import get_amplitude
import random as rn

with open('data/one_det.txt', 'w') as f:
    for experiments in range(1000):

        particles = rn.randint(1, 20)  # Теоретическое число частиц
        amplitude = 0  # Амплитуда

        for i in range(particles):  # За каждую частичку прибавляем амплитуду
            amplitude += get_amplitude()
        exp_particles = amplitude / 17.0147121406

        f.write(str(exp_particles / particles) + '\t' + str(particles) + '\n')

        print(exp_particles / particles)
