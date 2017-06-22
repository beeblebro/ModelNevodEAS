# Отборать восстановленные события из params.txt и restore.txt и сохранить
# в nice_restore.txt
import csv
from core import check_effective_area

params_file = open("output_data/params.txt", 'r')
restore_file = open("output_data/restore.txt", 'r')
file = open("output_data/nice_restore.txt", 'w')

params_data = []
restore_data = []
func_data = []

for column in zip(*[line for line in csv.reader(params_file, dialect="excel-tab")]):
    params_data.append(column)

for column in zip(*[line for line in csv.reader(restore_file, dialect="excel-tab")]):
    restore_data.append(column)

true_x = map(float, params_data[3])
true_y = map(float, params_data[4])
true_power = map(float, params_data[5])
true_age = map(float, params_data[6])

rec_x = map(float, restore_data[3])
rec_y = map(float, restore_data[4])
rec_power = map(float, restore_data[5])
rec_age = map(float, restore_data[6])

for r_x, r_y, r_age, r_power, t_x, t_y, t_age, t_power in zip(rec_x, rec_y, rec_age, rec_power, true_x, true_y, true_age, true_power):
    if check_effective_area(r_x, r_y):

        file.write(str(t_x) + '\t' + str(t_y) + '\t' +
                   str(t_age) + '\t' + str(t_power) + '\t' +
                   str(r_x) + '\t' + str(r_y) + '\t' +
                   str(r_age) + '\t' + str(r_power) + '\n')

params_file.close()
restore_file.close()
file.close()
