#encoding:utf-8
from cinfdata import Cinfdata
import common_toolbox as ct
import time, datetime
import pickle
from os import path
import numpy as np
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.stats import linregress
from scipy.integrate import simps

"""
Author: Jakob Ejler Sørensen
"""

CACHING = None
CACHE_DIR = '/home/jejsor/.tpd_cache'

class color():
    """
Instance to supply colors for plots. See:
https://matplotlib.org/3.1.0/gallery/color/named_colors.html
"""
    def __init__(self):
        # Red: reducing
        # Green: inert
        # Blue: oxidizing
        self.mass2color_dict = {
            '2': 'firebrick',
            '4': 'olivedrab',
            '12': 'black',
            '15': 'dimgray',
            '16': 'lightskyblue',
            '17': 'deepskyblue',
            '18': 'royalblue',
            '28': 'darkgreen',
            '29': 'olive',
            '30': 'darkorange',
            '32': 'mediumblue',
            '40': 'limegreen',
            '44': 'darkorchid', # matplotlib.colors.CSS4_COLORS[
        }

    def __getitem__(self, key):
        if isinstance(key, str):
            if key.lower().startswith('m'):
                new_key = key.lower().lstrip('m')
            else:
                new_key = key
        elif isinstance(key, int):
            new_key = str(key)
        else:
            raise TypeError('Mass must be a string or an integer.')
        if not new_key in self.mass2color_dict.keys():
            print(key, 'not found in color dict. Using standard sequential colors!')
            return
        else:
            name = self.mass2color_dict[new_key]
            color = mcolors.CSS4_COLORS[name]
        return color

def show_color_scheme():
    plt.figure()
    x = np.array([0, 1])
    y = np.array([1, 1])
    global color
    counter = 0
    mass_colors = color()
    for label, name in mass_colors.mass2color_dict.items():
        counter += 1
        if isinstance(label, str):
            label = label.lower().lstrip('m')
        plt.plot(x, y*int(label), color=mass_colors[label], linewidth=4, label='-'.join([label, name]))
    plt.plot(x, y*(counter+1), color=mass_colors['m1'], linewidth=4, label='m1')
    plt.legend()
    plt.title('Chosen color scheme for TPD masses')
    plt.show()

def get_range(X, time):
    A = np.where(X > time[0])[0][0]
    B = np.where(X < time[1])[0][-1]
    return np.arange(A, B)

def time2temp(data):
    time = data[:,0]
    temp = data[:,1]
    # Cycle through points and average groups
    times = np.array([])
    temps = np.array([])
    stds = np.array([])
    for i in range(len(time)):
        if i == 0:
            first = i
            value = temp[i]
            continue
        if temp[i] - value != 0:
            # Collect oversampling of temperature in a single point:
            # Extract region
            region = np.arange(first, i)
            times = np.append(times, time[region].mean())
            stds = np.append(stds, time[region].std())
            temps = np.append(temps, value)
            # Reset counters
            first = i
            value = temp[i]
    region = np.arange(first, i)
    times = np.append(times, time[region].mean())
    stds = np.append(stds, time[region].std())
    temps = np.append(temps, value)
    # Create interpolate function
    #interpolate = interp1d(times, temps, kind='linear')
    interpolate = interp1d(times, temps, kind='linear', bounds_error=False, fill_value='extrapolate')
    return interpolate


def import_data(ID, setup='omicron'):
    db = Cinfdata(setup, use_caching=True)
    data = db.get_data(ID)
    try:
        data[:,0] = data[:,0]/1000
    except:
        meta = db.get_metadata(ID)
        print('"{}" data for ID {} is empty and set to "None"'.format(meta['mass_label'], ID))
        data = None
    return data

#ALT_INFO = ['Analysis pressure', 'Preparation pressure',
#            'Sample temperature', 'Setpoint temperature',
#            'Voltage setpoint', 'Voltage monitor',
#            'Current setpoint', 'Current monitor']

class Experiment(object):
    """Imports a TPD data set from Omicron into a TPD class to implement
different practical functions for data treatment"""

    def __init__(self, timestamp, caching=False, set_label=None, temp_label='Sample temperature', recalculate=False):
        """Save all data in 'DATA' dict object"""

        # Various attributes
        self.pickle_format = 'tpd__{}__{}.pickle'
        self.exps = {}
        #self.caching = True if CACHING is True or caching is True else False
        self.caching = caching
        if CACHING is not None:
            self.caching = CACHING
        self.timestamp = {}
        self.timestamp['date'] = timestamp
        self.timestamp['unix'] = time.mktime(datetime.datetime.strptime(
                                             timestamp, "%Y-%m-%d %H:%M:%S").timetuple())
        # Try to load from file
        cache_file = []
        load_from_pickle = False
        if self.caching and not recalculate:
            list_of_files = path.os.listdir(CACHE_DIR)
            for name in list_of_files:
                if name.startswith('tpd__{}__'.format(timestamp)):
                    cache_file.append(name)
            if len(cache_file) == 1:
                with open(path.join(CACHE_DIR, cache_file[0]), 'rb') as f:
                    copy = pickle.load(f)
                    self.name = copy.name
                    self.labels = copy.labels
                    self.data = copy.data
                    self.exps = copy.exps
                load_from_pickle = True
                print(' *** Loaded TPD data from pickle: {} ***'.format(self.name))
            elif len(cache_file) > 1:
                print('Found multiple pickles matching timestamp:\n{}'.format(cache_file))

        # Get data from database
        if load_from_pickle is False:
            # Connect to database
            db = Cinfdata('omicron',
                        use_caching=caching,
                        grouping_column='time',
                        )
            # Get data in unit seconds
            group_data = db.get_data_group(timestamp, scaling_factors=(1E-3, None))
            group_meta = db.get_metadata_group(timestamp)

            # Get a list of labels and group data
            print(group_meta)
            self.name = group_meta[list(group_meta.keys())[0]]['Comment'].replace('/', '')
            self.labels = [group_meta[key]['mass_label'] for key in group_data.keys()]
            self.data = {group_meta[key]['mass_label']: group_data[key] for key in group_data.keys()}
            print('Loaded data from Experiment: "{}"'.format(self.name))
            self.isolate_experiments(set_label=set_label, temp_label=temp_label)


    def isolate_experiments(self, set_label=None, temp_label='Sample temperature'):
        """
1) Isolate regions of heating ramps
2) Organize into returnable data
3) Convert time to temperature
    - intentionally oversample the temperature
    - collect groups into points of average/std
    - interpolate points to correlate arbitrary time with a temperature
"""

        # 1) Find regions of heating ramps
        if set_label is None:
            for i in self.labels:
                if 'setpoint' in i.lower().split() and 'temperature' in i.lower().split():
                    set_label = i
                    break
        if set_label is None:
            print('Ramp label is not found automatically. Specify "set_label"!')
            raise ValueError()
        if len(self.data[set_label]) == 0:
            raise ValueError('No heat ramps detected')
        marker_start, marker_end = 0, 0
        counter = 0
        regions = []
        for i in self.data[set_label][:,0]:
            if counter == 0:
                last = i
                marker_start = counter
                counter += 1
                continue

            # Main check
            if i - last > 15: # difference of more than 15 seconds
                marker_end = counter
                region = np.arange(marker_start, marker_end)
                regions.append(region)
                marker_start = counter

            # End of loop
            last = i
            counter += 1
        region = np.arange(marker_start, counter)
        regions.append(region)
        number_of_regions = len(regions)
        #print(number_of_regions)

        # 2) Organize into returnable data
        # Get temperature data
        if not temp_label in self.labels:
            for i in self.labels:
                if 'sample' in i.lower() and 'temperature' in i.lower():
                    temp_label = i
                    print('Temperature label deviating from default found: "{}"'.format(temp_label))
                    break
            else:
                print('"Sample temperature" not found in loaded dataset. Please specify "temp_label"!')
                return None
        interpolate_temp = time2temp(self.data[temp_label])
        # Local heating rate
        number = 0
        t_window = 10

        # For fast temperature logging, this seems to work very good. Hasn't been tested on
        # sparse datasets *** OBSOLETE ***
        while number < 3:
            t_window += 1
            number = len(ct.get_range(self.data[temp_label][:, 0], self.data[temp_label][-1, 0]-t_window, self.data[temp_label][-1, 0]))
        self.slope = self.data[temp_label].copy()
        self.slope[:number+1, 1] = np.nan
        for i in range(len(self.slope[:, 0]) -1, number, -1):
            coeff = np.polyfit(self.data[temp_label][i-number:i+1, 0], self.data[temp_label][i-number:i+1, 1], 1)
            self.slope[i, 1] = coeff[0]
        interpolate_slope = time2temp(self.slope)

        # Set up data structure for experiments
        for i in range(number_of_regions):
            self.exps[i] = {}
            region = regions[i]
            time_ref = self.data[set_label][:,0][region]
            for label in self.labels:
                # skip if empty
                #print(self.data[label], label)
                #if self.data[label] == None:
                #if len(self.data[label]) == 0:
                #    continue
                if len(self.data[label]) == 0:
                    continue
                self.exps[i][label] = {}
                local_range = get_range(self.data[label][:, 0], [time_ref[0]-3, time_ref[-1]+3])
                for x in [0, 1]:
                    self.exps[i][label][x] = self.data[label][:, x][local_range]

            # 3) Add temperature axis to each data set
            for label in self.labels:
                # skip if empty
                #if self.data[label] == None:
                #if len(self.data[label]) == 0:
                #    continue
                if len(self.data[label]) == 0:
                    continue
                print(i, label)
                self.exps[i][label][2] = interpolate_temp(self.exps[i][label][0])
                # Extract information about local heating rate
                dat = self.exps[i][label]
                dat[3] = np.empty(len(dat[0]))
                dat[3][0] = np.nan
                dat[3][1:] = ct.smooth(np.diff(dat[2])/np.diff(dat[0]), width=1)

        # Save results in pickle:
        self.save_pickle()

    def save_pickle(self):
        """Save class instance to pickle file"""
        if self.caching:
            if not path.exists(CACHE_DIR):
                print('TPD.CACHE_DIR doesn\'t exist. Create it or change to existing path.')
                raise OSError('Create path for me!')
            pickle_name = self.pickle_format.format(self.timestamp['date'], self.name)
            file_string = path.join(CACHE_DIR, pickle_name)
            with open(file_string, 'wb') as f:
                pickle.dump(self, f, pickle.HIGHEST_PROTOCOL)
            print('Data saved as {}'.format(pickle_name))

        else:
            print('Set TPD.CACHING = True to save instance as pickle. Nothing saved')

    def load_pickle(self):
        """Reload self from pickle file"""
        if CACHING is True:
            # Check if data is loadable
            path

    def get_timesnippet(self, label, limit):
        """Return (time, data) of data set 'label' cropped to fit into time 'limit'"""
        index = ct.get_range(self.data[label][:, 0], limit[0], limit[1])
        return self.data[label][:, 0][index], self.data[label][:, 1][index]

def get_total_signal(experiment, mass_label='M30'):
    """ Integrate time based QMS signal from 'isolate_experiments()' 
Resembles a lot"""

    # Check type of data input
    if isinstance(experiment, Experiment):
        # Redefine data
        experiment = experiment.exps
    elif isinstance(experiment, dict):
        pass # Probably good
    else:
        raise TypeError('"data" should be a TPD.Experiment instance for default behaviour' \
'or data set like TPD.Experiment.data["Prep Pressure"] for specific behaviour.')
    if not mass_label in experiment[0]:
        raise KeyError('mass_label {} not found in experiment dictionary {}'.format(mass_label, experiment.keys()))

    # For every experiment, add integrated signal as 'COVERAGE'
    coverage = np.zeros(len(experiment.keys()))
    for i in experiment.keys():
        ranges = []
        x = experiment[i][mass_label][0]
        y = experiment[i][mass_label][1]
        N = len(x)
        indice = np.arange(10)
        indice = np.append(indice, range(N-10, N))
        slope, intercept, rvalue, pvalue, std_err = linregress(x[indice], y[indice])
        baseline = intercept + slope*x
        coverage[i] = simps(y - baseline, x=x) 
    return coverage

def Langmuir(dose):
    """Take a dose in mbar*s and return it in Langmuir"""
    return dose*100/1.33e-4

def get_doses(data, times, unit='Langmuir', subtract_background=False):
    """Integrate doses.
Data can be a TPD.Experiment object, in which case it defaults to self.data['Prep Pressure']
Otherwise, data can be the numpy array you want to integrate. Pass it for instance:
    self.data['M30'] to integrate the mass 30 signal.
Times is a list of time sequences to integrate. Ex: times=[[12, 25], [1200, 1310]]
Signal in data is integrated between these times with a linear background subtraction.

subtract_background (bool): if False (default), the raw signal is integrated with simps method,
    else a linear background is fitted and subtracted from between the endpoints times[i]

Returns: a list of doses with each index corresponding to the index in times."""

    # Check type of data input
    if isinstance(data, Experiment):
        # Redefine data
        data = data.data['Prep Pressure']
    elif isinstance(data, np.ndarray):
        pass # Good
    else:
        raise TypeError('"data" should be a TPD.Experiment instance for default behaviour' \
'or data set like TPD.Experiment.data["Prep Pressure"] for specific behaviour.')
    
    # Initialize variables and ranges
    ranges = []
    for deltat in times:
        ranges.append(get_range(data[:,0], deltat))
    doses = np.zeros(len(ranges))
    
    # Loop to subtract background and integrate doses
    i = 0
    for r in ranges:
        x = data[:,0][r]
        y = data[:,1][r]
        if subtract_background:
            N = len(x)
            indice = np.arange(10)
            indice = np.append(indice, range(N-10, N))
            slope, intercept, rvalue, pvalue, std_err = linregress(x[indice], y[indice])
            baseline = intercept + slope*x
            doses[i] = simps(y - baseline, x=x)
        else:
            doses[i] = simps(y, x=x)
        i += 1

    # Convet to return unit
    if unit.lower().startswith('l'): # Langmuir
        doses = Langmuir(doses)
    return doses


def get_doseage(data, times):
    """Changed to get_doses """
    
    # Initialize variables and ranges
    ranges = []
    for deltat in times:
        ranges.append(get_range(data[:,0], deltat))
    doses = np.zeros(len(ranges))
    
    # Loop to subtract background and integrate doses
    i = 0
    for r in ranges:
        x = data[:,0][r]
        y = data[:,1][r]
        N = len(x)
        indice = np.arange(10)
        indice = np.append(indice, range(N-10, N))
        slope, intercept, rvalue, pvalue, std_err = linregress(x[indice], y[indice])
        baseline = intercept + slope*x
        doses[i] = simps(y - baseline, x=x)
        i += 1
    return doses



def tick_function(T):
    """ Dummy function for 'ADD_SCALE_KELVIN' """

    K = T + 273.15
    return ["%d" % z for z in K]



def add_scale_Kelvin(ax, axis='x'):
    """ Convert a Celcius axis to Kelvin and coplot on opposite axis """
    
    if axis not in ['x', 'y']:
        raise ValueError('Choose \'x\' or \'y\' as AXIS')
    if axis == 'x':
        ax_twin = ax.twiny()
        ax_twin.set_xlim(ax.get_xlim())
        new_ticks = ax.get_xticks()
    elif axis == 'y':
        ax_twin = ax.twinx()
        ax_twin.set_ylim(ax.get_ylim())
        new_ticks = ax.get_yticks()
    new_start = 100 - 273.15
    new_ticklocations = np.arange(len(new_ticks)-1)*np.diff(new_ticks) + new_start
    if axis == 'x':
        ax_twin.set_xticks(new_ticklocations)
        ax_twin.set_xticklabels(tick_function(new_ticklocations))
        ax_twin.set_xlabel('Temperature (K)')
    elif axis == 'y':
        ax_twin.set_yticks(new_ticklocations)
        ax_twin.set_yticklabels(tick_function(new_ticklocations))
        ax_twin.set_ylabel('Temperature (K)')



def Veff(H, r = 0.5):
    """Effort to estimate the volume between sniffer and sample"""
    return r * H**2 + r**2 * H + 1./3 * H**3



def generate_standard_plots(exp, coverage=None, doses=None, selection=[]):
    """ Generate standard plots for a single TPD experiment. """

    # Check if selection has been given
    if len(selection) < 1:
        selection = exp.keys()

    GUIDE_LINES = (np.arange(10)+1) * 100
    fig, ax = {}, {}    
    # Plot general data overview plot for each sub-TPD
    for i in selection:
        fig_id = 'overview TPD {}'.format(i)
        fig[fig_id] = plt.figure(fig_id)

        # Plot actual TPD result
        ax[i] = {}
        ax[i][1] = fig[fig_id].add_subplot(221)
        data = exp[i]['M30']
        ax[i][1].plot(data[2], data[1], 'r-')
        ax[i][1].set_xlabel('Temperature (C)')
        ax[i][1].set_ylabel('SEM current (A)')
        add_scale_Kelvin(ax[i][1])

        # Plot raw signal vs time
        ax[i][2] = fig[fig_id].add_subplot(222)
        ax[i][2].plot(data[0] - data[0][0], data[1], 'ro-')
        ax[i][2].set_xlabel('Time (s)')
        ax[i][2].set_ylabel('SEM current (A)')
        if not coverage is None and not doses is None:
            ax[i][2].text(0.5, 0.85, 'Dose:        {:.3} mbar s\nCoverage: {:.3} C'.format(doses[i], coverage[i]), transform=ax[i][2].transAxes)

        # Plot temperature vs time
        ax_id = 'temperature {}'.format(i)
        ax[i][3] = fig[fig_id].add_subplot(223)
        data = exp[i]['Sample temperature']
        slope = exp[i]['Slope']
        ax[i][3].plot(data[0] - data[0][0], data[1], 'b-')
        ax[i][3].plot(data[0] - data[0][0], slope, 'k-')
        value_of_slope = (slope[-1]- slope[0])/(data[0][-1]-data[0][0])
        ax[i][3].text(.05, 0.7, 'Heat rate: {:.3} K/s'.format(value_of_slope), transform=ax[i][3].transAxes)
        ax[i][3].set_xlabel('Time (s)')
        ax[i][3].set_ylabel('Temperature (C)')
        add_scale_Kelvin(ax[i][3], axis='y')

        # Plot difference in temperature and slope vs time
        ax[i][4] = fig[fig_id].add_subplot(224)
        ax[i][4].axhline(y=0, linestyle='solid', color='k')
        for j in [-1, 1]:
            ax[i][4].axhline(y=j, linestyle='dashed', color='k')
        ax[i][4].plot(data[0] - data[0][0], data[1] - slope, 'b-')
        ax[i][4].set_xlabel('Time (s)')
        ax[i][4].set_ylabel('Temp$_{meas}$ - Temp$_{avg}$')

        
        # Plot power supply data
        fig_id = 'Filament data {}'.format(i)
        fig[fig_id] = plt.figure(fig_id)
        ax[i][5] = fig[fig_id].add_subplot(121)
        data = exp[i]['Voltage monitor']
        ax[i][5].plot(data[0] - data[0][0], data[1], 'k-')
        ax[i][5].set_xlabel('Time (s)')
        ax[i][5].set_ylabel('Voltage (V)')
        ax[i][5].set_ylim([0, 6.5])
        ax[i][5].arrow(.7, .4, 0.05, 0, color='m', transform=ax[i][5].transAxes)
        ax[i][5].text(.7, .4, 'Current ', color='m', verticalalignment='center', horizontalalignment='right', transform=ax[i][5].transAxes)
        ax[i][5].arrow(.55, .3, -0.05, 0, color='k', transform=ax[i][5].transAxes)
        ax[i][5].text(.55, .3, ' Voltage', color='k', verticalalignment='center', horizontalalignment='left', transform=ax[i][5].transAxes)
        ax_twinx = ax[i][5].twinx()
        
        data = exp[i]['Current monitor']
        ax_twinx.plot(data[0] - data[0][0], data[1], 'm-')
        ax_twinx.set_ylim([0, 6.0])
        ax_twinx.set_ylabel('Current (A)')

        # Plot PID data
        ax[i][6] = fig[fig_id].add_subplot(122)
        data = exp[i]['P error']
        x = data[0] - data[0][0]
        y = data[1]
        ax[i][6].plot(x[25:], y[25:], 'k-') # Hardcoded index range due to earlier extended experiment region
        ax[i][6].set_xlabel('Time (s)')
        ax[i][6].set_ylabel('P error')
        ax_twinx = ax[i][6].twinx()
        
        data = exp[i]['I error']
        x = data[0] - data[0][0]
        y = data[1]
        ax_twinx.plot(x[25:], y[25:], 'm-') # Hardcoded index range due to earlier extended experiment region
        ax_twinx.set_ylabel('I error')

        # Plot guiding lines
        XLIM = [0, data[0][-1]-data[0][0]]
        for j in [2, 3, 4, 5, 6]:
            print(i,j)
            axis = ax[i][j]
            for X in GUIDE_LINES:
                axis.axvline(x=X, linestyle='-.', color='k', alpha=0.3)
            axis.set_xlim(XLIM)
    
    # Plot main result
    fig_id = 'Result'
    fig[fig_id] = plt.figure(fig_id)
    ax_id = 'Result'
    ax[ax_id] = fig[fig_id].add_subplot(111)
    for i in selection:
        data = exp[i]['M30']
        ax[ax_id].plot(data[2], data[1])
    ax[ax_id].set_xlabel('Temperature (C)')
    ax[ax_id].set_ylabel('SEM current (A)')
    ax[ax_id].text(0.05, 0.90, 'm/z = 30', color='k')
    add_scale_Kelvin(ax[ax_id])

    if not coverage is None and not doses is None:
        # Plot uptake curve
        fig_id = 'Uptake curve'
        fig[fig_id] = plt.figure(fig_id)
        ax_id = 'Uptake curve'
        ax[ax_id] = fig[fig_id].add_subplot(111)
        ax[ax_id].plot(doses/1.33e-6, coverage, 'ko')
        ax[ax_id].set_xlabel('Doseage (1.33e-6 mbar s)')
        ax[ax_id].set_ylabel('Coverage (C)')
        x1, x2 = ax[ax_id].get_xlim()
        ax[ax_id].set_xlim([0, x2])

    # Show plots
    plt.show()

if __name__ == '__main__':
    data = Experiment('2020-02-26 12:49:44', caching=False)
    exp = data.isolate_experiments()

    colors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']

    # Plot experiment(s)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    for i in exp.keys():
        dat = exp[i]['M30']
        ax.plot(dat[2], dat[1]/1e-12, marker='o', color=colors[i])
    ax.set_xlabel('Temperature ($^\circ$C)')
    ax.set_ylabel('SEM current (pA)')
    add_scale_Kelvin(ax)

    plt.show()
