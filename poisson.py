from numpy.random import poisson


with open('data/poisson/psn.txt', 'w') as file:
    for i in range(100000):
        file.write(str(poisson(1)) + '\t'
                   + str(poisson(4)) + '\t' + str(poisson(10)) + '\n')
