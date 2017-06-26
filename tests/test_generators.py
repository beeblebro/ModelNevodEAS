from numpy.random import poisson, normal
from math import sqrt
from core.amplitude import gen_amplitude, get_av_amplitude


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
                res += gen_amplitude()
            res /= get_av_amplitude()
            fa.write(str(res) + "\n")


def amplitudes(num):

    with open("am_test.txt", 'w') as fs:
        for i in range(num):
            res = gen_amplitude()
            fs.write(str(res) + "\n")


def test_both(t, num):

    with open("both.txt", "w") as f:
        for i in range(num):
            if t <= 25:
                t1 = poisson(t)
            else:
                t1 = round(normal(t, sqrt(t)))

            res = 0
            for j in range(t1):
                res += gen_amplitude()
            res /= get_av_amplitude()
            f.write(str(res) + "\n")


def func():
    with open("am_test.txt", 'w') as f:
        x = 0
        for i in range(10**7):
            x += i * 0.000008
            f.write(str(x) + "\t" + str(function(x)) + "\n")


# func()
# test_poisson(2, 10**6)
# test_am(5, 10**6)
test_both(30, 10**6)
# amplitudes(10**7)
