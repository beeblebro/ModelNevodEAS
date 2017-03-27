import random as rn

from core import Facility
from core import Eas
from core.utils import get_theta

f = open('data/power_age_func/power_age.txt', 'w')

# Создали устаовку
NevodEAS = Facility()


for experiments in range(5):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    # x0 = rn.uniform(-50, 50)
    # y0 = rn.uniform(-50, 50)
    x0 = 30
    y0 = -50
    power = 5 * 10 ** 6
    age = 1.5
    # Создали ШАЛ
    eas = Eas(theta, phi, x0, y0, power, age)

    # Установка получает ливень
    NevodEAS.get_eas(eas)
    # Запускаем установку
    if not NevodEAS.start():
        # Пропустим итерацию цикла, если ничего не сработало
        print("ERROR: Не сработал ни один кластер")
        continue

    # Восстанавливаем точку прихода, мощность и возраст
    NevodEAS.new_rec_params_bfgs()

    # NevodEAS.draw_func_age(x0, y0, power)

    print(experiments)
    f.write(str(NevodEAS.rec_power) + '\t' + str(NevodEAS.rec_age))
    f.write('\n')
    NevodEAS.reset()
f.close()

