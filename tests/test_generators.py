from numpy.random import poisson, normal
from math import sqrt
from core.amplitude import get_amplitude, get_av_amplitude, ampl_func


def test_poisson(lam, num):

    with open("poisson_test.txt", 'w') as fp:
        for i in range(num):
            if lam <= 25:
                fp.write(str(poisson(lam)) + "\n")
            else:
                fp.write(str(normal(lam, sqrt(lam))) + "\n")


def test_am(t, num):

    with open("am_test.txt", 'w') as fa:
        for i in range(num):
            res = 0
            for j in range(t):
                res += get_amplitude()
            res /= get_av_amplitude()
            fa.write(str(res) + "\n")


def test_both(t, num):

    with open("both.txt", "w") as f:
        for i in range(num):
            if t <= 25:
                t1 = poisson(t)
            else:
                t1 = normal(t, sqrt(t))

            res = 0
            for j in range(t1):
                res += get_amplitude()
            res /= get_av_amplitude()
            f.write(str(res) + "\n")

test_poisson(2, 10**5)
# test_am(2, 10**5)
# test_both(2, 10**5)
