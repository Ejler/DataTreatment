"""Remove data points from volvo ISS resulting from single spikes.
"""
import numpy as np
import matplotlib.pyplot as plt

import common_toolbox as ct
from cinfdata import Cinfdata

#mode = 'test'
mode = 'plot'
#mode = 'save'
settings = {
    #9489: 220,
    #9494: 400,
    #9504: 1500,
    #9507: 800,
    #9509: 1000,
    #9512: 600,
    #9516: 1000,
    #9518: 500,
    #9527: 600,
    9545: 350,
}

db = Cinfdata('dummy', use_caching=False)
for ID, LEVEL in settings.items():
    data = db.get_data(ID)
    time, signal = data[:, 0], data[:, 1]

    smooth = ct.smooth(signal)
    diff = np.abs(np.diff(signal-smooth))

    index = np.where(diff > LEVEL)[0] + 1
    i = np.where(time > 0)[0]
    
    if mode in ['test', 'plot']:
        plt.title(str(ID))
        plt.plot(data[:, 0], data[:, 1], 'b-', label='Raw data')
        plt.plot(time[index], signal[index], 'mo', markerfacecolor='w')

    index = [x for x in i if x not in index]
    if mode == 'plot':
        plt.plot(time[index], signal[index], 'm-')
        plt.plot(time[index], ct.smooth(signal[index]), 'k-')
    elif mode == 'test':
        plt.plot(time[1:], diff)
        plt.axhline(y=LEVEL)
    if mode in ['test', 'plot']:
        plt.legend()
        plt.show()

    if mode == 'save':
        name = f'iss_smoothed_id{ID}.csv'
        with open(name, 'w') as f:
            x, y = time[index], ct.smooth(signal[index])
            for i, j in list(zip(x, y)):
                f.write(f'{i},{j}\r\n')
