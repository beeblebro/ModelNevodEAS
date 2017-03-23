import random as rn

from core import Facility
from core import Eas
from core.utils import get_theta

f = open('data/power_age_func/power_age.txt', 'w')

# Создали устаовку
NevodEAS = Facility()


for experiments in range(100):
    theta = get_theta()
    phi = rn.uniform(0, 360)
    # x0 = rn.uniform(-50, 50)
    # y0 = rn.uniform(-50, 50)
    x0 = 25
    y0 = 25
    power = 10 ** 6
    age = 1.45
    # Создали ШАЛ
    eas = Eas(theta, phi, x0, y0, power, age)

    # Установка получает ливень
    NevodEAS.get_eas(eas)
    # Запускаем установку
    if not NevodEAS.start():
        # Пропустим итерацию цикла, если ничего не сработало
        continue
    # Восстанавливаем точку прихода, мощность и возраст
    NevodEAS.rec_params()

    # draw_func_power(NevodEAS.clusters, eas.n, eas.x0, eas.y0, start_power, eas,age,
    # NevodEAS.exp_n, NevodEAS.sigma_n)
    # draw_func_age(NevodEAS.clusters, eas.n, eas.x0, eas.y0, start_power, eas,age,
    # NevodEAS.exp_n, NevodEAS.sigma_n)

    print(experiments)
    f.write(str(NevodEAS.rec_power) + '\t' + str(NevodEAS.rec_age))
    f.write('\n')
    NevodEAS.reset()
f.close()

