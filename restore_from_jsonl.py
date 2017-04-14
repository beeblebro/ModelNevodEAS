import random as rn
import json

from core import Facility
from core import Eas
from core.utils import get_theta, get_power, get_age

f_data = open('model_10k_flat.jsonl', 'r')
f_dump = open('restore_flat.txt', 'w')
f_dump_params = open('params_flat.txt', 'w')

# Создали устаовку
NevodEAS = Facility(geometry='flat')

for line in f_data:
    event = json.loads(line)

    num = event['num']
    theta = event['params']['theta']
    phi = event['params']['phi']
    x0 = event['params']['x0']
    y0 = event['params']['y0']
    power = event['params']['power']
    age = event['params']['age']

    f_dump_params.write(str(theta) + '\t' +
                        str(phi) + '\t' +
                        str(x0) + '\t' +
                        str(y0) + '\t' +
                        str(power) + '\t' +
                        str(age) + '\n')

    eas = Eas(theta, phi, x0, y0, power, age)
    NevodEAS.get_eas(eas)

    for cl_n, cl in enumerate(NevodEAS.clusters):

        st_ok = 0
        for st_n, st in enumerate(cl.stations):
            st.amplitude = event['clusters'][cl_n]['st'][st_n]['ampl']
            st.rndm_time = event['clusters'][cl_n]['st'][st_n]['time']

            # Расставим правильные отклики, т.к. их не сохраняли
            if st.amplitude is not None and st.rndm_time is not None:
                st.respond = True
                st_ok += 1
            else:
                st.respond = False

        if st_ok == 4:
            cl.respond = True
        else:
            cl.respond = False

    NevodEAS.rec_direction()
    NevodEAS.rec_particles()

    if NevodEAS.rec_params_diff_evo():
        f_dump.write(str(NevodEAS.rec_theta) + '\t' + str(NevodEAS.rec_phi) + '\t' +
                     str(NevodEAS.rec_x) + '\t' + str(NevodEAS.rec_y) + '\t' +
                     str(NevodEAS.rec_power) + '\t' + str(NevodEAS.rec_age) + '\n')

    NevodEAS.reset()

    print(num)
    if num == 99:
        break

f_data.close()
f_dump.close()
f_dump_params.close()

