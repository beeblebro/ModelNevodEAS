import json

from core import Facility
from core import Eas


f_data = open('model_10k.jsonl', 'r')
f_dump = open('rec_n.txt', 'w')
f_dump_params = open('true_n.txt', 'w')
f_dump_func = open('func.txt', 'w')

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

    # f_dump_params.write(str(theta) + '\t' +
    #                     str(phi) + '\t' +
    #                     str(x0) + '\t' +
    #                     str(y0) + '\t' +
    #                     str(power) + '\t' +
    #                     str(age) + '\n')

    eas = Eas(theta, phi, x0, y0, power, age)
    NevodEAS.get_eas(eas)

    # f_dump_params.write(str(round(theta, 4)) + '\t' +
    #                     str(round(phi, 4)) + '\n')

    f_dump_params.write(str(round(eas.n[0], 4)) + '\t' +
                        str(round(eas.n[1], 4)) + '\t' +
                        str(round(eas.n[2], 4)) + '\n')

    NevodEAS.set_facility_state(event)

    f_dump.write(str(round(NevodEAS.average_n[0], 4)) + '\t' +
                 str(round(NevodEAS.average_n[1], 4)) + '\t' +
                 str(round(NevodEAS.average_n[2], 4)) + '\n')

    # f_dump.write(str(round(NevodEAS.rec_theta, 4)) + '\t' +
    #              str(round(NevodEAS.rec_phi, 4)) + '\n')

    # if NevodEAS.rec_params_powell():
    #
    #     params = [x0, y0, power, age]
    #     rec_params = [NevodEAS.rec_x, NevodEAS.rec_y, NevodEAS.rec_power, NevodEAS.rec_age]
    #     # delta = [x0 - x for x0, x in zip(params, rec_params)]
    #     # print(delta)
    #
    #     f_dump.write(str(NevodEAS.rec_theta) + '\t' + str(NevodEAS.rec_phi) + '\t' +
    #                  str(rec_params[0]) + '\t' + str(rec_params[1]) + '\t' +
    #                  str(rec_params[2]) + '\t' + str(rec_params[3]) + '\n')
    #
    #     # func = NevodEAS.func(params)
    #     func_rec = NevodEAS.func(rec_params)
    #     # delta_func = func - func_rec
    #
    #     f_dump_func.write(str(func_rec) + '\n')

    NevodEAS.reset()

    print(num)
    if num == 999:
        break

f_data.close()
f_dump.close()
f_dump_params.close()
f_dump_func.close()