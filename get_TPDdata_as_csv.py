"""Takes a timestamp as argument and writes the corresponding TPD experiment
to a CSV file.
"""
import sys
import numpy
import TPD

timestamp = sys.argv[1]
data = TPD.Experiment(timestamp, caching=False)
exp = data.isolate_experiments()

filename = '{} {} - %s_%s'.format(timestamp, data.name)

easy = False # If true: provide "as function of temperature"

# Header
header = 'time (s),value (arb)'
if easy is True:
    header += ',temperature (K)'
header += '\r\n'

# Body
for number in exp.keys():
    print('Writing experiment {}'.format(number))
    new_filename = filename.split('--')[0] + '_Region_{}.csv'.format(number)
    f = open(new_filename %('M30', number), 'a')
    g = open(new_filename %('Sample temperature', number), 'a')

    x = exp[number]['M30'][0]
    y = exp[number]['M30'][1]
    if easy is True:
        z = exp[number]['M30'][2]

    f.write(header)
    for i in range(len(y)):
        line = '{},{}\r\n'.format(x[i]-x[0], y[i])
        f.write(line)
    f.close()

    x = exp[number]['Sample temperature'][0]
    y = exp[number]['Sample temperature'][1]
    if easy is True:
        z = exp[number]['Sample temperature'][2]

    g.write(header)
    for i in range(len(y)):
        line = '{},{}\r\n'.format(x[i]-x[0], y[i])
        g.write(line)
    g.close()

print('Done!')
