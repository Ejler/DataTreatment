"""
A common toolbox containing useful functions for data analysis and plotting
"""
import numpy as np
import matplotlib.pyplot as plt

colors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']*5

def flip_x(ax=None):
    """Flip the axis of `Axes` object ax"""
    if ax is None:
        ax = plt.gca()
    x1, x2 = ax.get_xlim()
    ax.set_xlim([x2, x1])

def get_range(x, lim1, lim2):
    """Return the index range of x between lim1 and lim2"""
    index1 = np.where(x <= lim2)[0]
    index2 = np.where(x >= lim1)[0]
    index = np.intersect1d(index1, index2)
    if len(index) == 0:
        print('"get_range" didn\'t find any data within the limits!')
    return index

def copy_range(x, lim1, lim2):
    """Return a copy of x between lim1 and lim2"""
    x_copy = x[np.where(x >= lim1)[0]]
    x_copy = x_copy[np.where(x_copy <= lim2)[0]]
    return x_copy

def smooth(data, width=1):
    """Average `data` with `width` neighbors"""
    if len(data) == 0:
        print('Empty data! ', len(data))
        return data
    smoothed_data = np.zeros(len(data))
    smoothed_data[width:-width] = data[2*width:]
    for i in range(2*width):
        smoothed_data[width:-width] += data[i:-2*width+i]
        if i < width:
            smoothed_data[i] = sum(data[0:i+width+1])/len(data[0:i+width+1])
            smoothed_data[-1-i] = sum(data[-1-i-width:])/len(data[-1-i-width:])
    smoothed_data[width:-width] = smoothed_data[width:-width]/(2*width+1)
    return smoothed_data

def weighted_smooth(data, repeat=1):
    """I[n] = 1/4(I[n-1] + 2I[n] + I[n+1]"""
    for iteration in range(repeat):
        copy = np.zeros(data.shape)
        copy += data*2
        copy[1:] += data[:-1]
        copy[:-1] += data[1:]
        copy[1:-1] /= 4
        copy[0] /= 3
        copy[-1] /= 3
        data = copy.copy()
    return data
