import json

from core import Facility
from core import Eas


f_data = open('model_10k.jsonl', 'r')
f_restore = open('restore.txt', 'w')
f_params = open('params.txt', 'w')
f_delta = open('delta.txt', 'w')
f_func = open('func.txt', 'w')

# Создали устаовку
NevodEAS = Facility(geometry='nevod')

for line in f_data:
    event = json.loads(line)

    num = event['num']
    theta = event['params']['theta']
    phi = event['params']['phi']
    x0 = event['params']['x0']
    y0 = event['params']['y0']
    power = event['params']['power']
    age = event['params']['age']

    f_params.write(str(theta) + '\t' +
                   str(phi) + '\t' +
                   str(x0) + '\t' +
                   str(y0) + '\t' +
                   str(power) + '\t' +
                   str(age) + '\n')

    eas = Eas(theta, phi, x0, y0, power, age)

    # f_params.write(str(eas.n[0]) + '\t' + str(eas.n[1]) + '\t' + str(eas.n[2]) + '\n')
    # f_params.write(str(theta) + '\t' + str(phi) + '\n')

    if not NevodEAS.set_facility_state(event):
        continue

    # f_restore.write(str(NevodEAS.average_n[0]) + '\t' + str(NevodEAS.average_n[1]) +
    #                 '\t' + str(NevodEAS.average_n[2]) + '\n')

    # f_restore.write(str(NevodEAS.rec_theta) + '\t' + str(NevodEAS.rec_phi) + '\n')

    if NevodEAS.rec_params_diff_evo():
        # Восстановление
        params = [x0, y0, power, age]
        rec_params = [NevodEAS.rec_x, NevodEAS.rec_y, NevodEAS.rec_power, NevodEAS.rec_age]
        delta = [x0 - x for x0, x in zip(params, rec_params)]

        # Восстановленные данные
        f_restore.write(str(NevodEAS.rec_theta) + '\t' + str(NevodEAS.rec_phi) + '\t' +
                        str(rec_params[0]) + '\t' + str(rec_params[1]) + '\t' +
                        str(rec_params[2]) + '\t' + str(rec_params[3]) + '\n')

        # Дельты
        f_delta.write(str(delta[0]) + '\t' + str(delta[1]) + '\t' +
                      str(delta[2]) + '\t' + str(delta[3]) + '\n')

        # Остаточная сумма
        f_func.write(str(NevodEAS.func(rec_params)) + '\n')

    NevodEAS.reset()

    print(num)
    if num == 999:
        break

f_data.close()
f_restore.close()
f_params.close()
f_func.close()
f_delta.close()
