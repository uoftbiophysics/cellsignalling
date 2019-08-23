import os
import re
fname = os.getcwd() + os.sep + 'MLELookUpTable_dense_NewT.m'
new_fname = os.getcwd() + os.sep + 'MLELookUpTable_dense_NewT_Fixed.m'

def get_numbers(line):
    line = line.strip()
    line = re.sub('{','',line)
    line = re.sub('}', '', line)
    line = re.sub("\*\^", 'E', line)
    line = line.split(',')
    line[0] = int(line[0])
    line[1] = float(line[1])
    line[2] = float(line[2])
    return line[:3]

fixed_file = ''
with open(fname, 'r') as f:
    for m in range(10, 9990, 10):
        chunk = ' {'
        line = f.readline()
        numbers = get_numbers(line)
        numbers = [n / 10. for n in numbers]
        nmck = [m, int(numbers[0]), numbers[1], numbers[2]]
        newline = '{' + str(nmck[0]) + ', {' + str(nmck[1]) + ', {' + str(nmck[2]) + ', ' + str(nmck[3]) + '}}},\n'
        chunk += newline
        for n in range(20, 9990, 10):
            line = f.readline()
            numbers = get_numbers(line)
            numbers = [n/10. for n in numbers]
            nmck = [m, int(numbers[0]), numbers[1], numbers[2]]
            newline = '  {' + str(nmck[0]) + ', {' + str(nmck[1]) + ', {' + str(nmck[2]) + ', ' + str(nmck[3]) + '}}},\n'
            chunk += newline
        chunk = chunk[:-2] + '},\n'
        fixed_file += chunk
fixed_file = '{' + fixed_file[1:-2] + '}'

with open(new_fname, 'w+') as f:
    f.write(fixed_file)



