import random as rn

from core import Facility
from core import Eas
from core.utils import get_theta, get_power, get_age

f = open('data/power_age_func/power_age.txt', 'w')

# Создали устаовку
NevodEAS = Facility()

for tries in range(20):
    theta = get_theta()  # Тета
    phi = rn.uniform(0, 360)  # и фи в градусах
    x0 = rn.uniform(-50, 50)
    y0 = rn.uniform(-50, 50)
    power = get_power()
    age = get_age(power, theta)
    # Создали ШАЛ
    eas = Eas(theta, phi, x0, y0, power, age)

    # Установка получает ливень
    NevodEAS.get_eas(eas)
    # Запускаем установку
    if not NevodEAS.start():
        # Пропустим итерацию цикла, если ничего не сработало
        NevodEAS.reset()
        continue

    # Восстанавливаем точку прихода, мощность и возраст
    if not NevodEAS.new_rec_params_bfgs():
        print("ERROR: Не удалось восстановить параметры ШАЛ")
        NevodEAS.reset()
        continue

    f.write(str(x0) + '\t' + str(NevodEAS.rec_x) + '\t'
            + str(y0) + '\t' + str(NevodEAS.rec_y) + '\t'
            + str(power) + '\t' + str(NevodEAS.rec_power) + '\t'
            + str(age) + '\t' + str(NevodEAS.rec_age) + '\n')
    NevodEAS.reset()

f.close()

