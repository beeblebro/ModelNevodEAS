# Тестируем генератор значений энергии

from core.tools import get_energy
from math import log10

with open('data/energy/energy.txt', 'w') as file:
    for i in range(100000):
        file.write(str(get_energy()) + '\n')
