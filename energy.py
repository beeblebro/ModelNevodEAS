# Тестируем генератор значений энергии

from core.tools import get_power
from math import log10

with open('data/power/power.txt', 'w') as file:
    for i in range(100000):
        file.write(str(get_power()) + '\n')
