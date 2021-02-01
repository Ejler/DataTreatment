# pylint: disable=invalid-name,too-many-arguments
"""FIXME"""
import matplotlib.pyplot as plt
import numpy as np
import functools
from scipy.optimize import curve_fit
from scipy.integrate import simps
import common_toolbox as ct

class Experiment():
    """Load an ISS experiment exported as text or VAMAS file.

Author: Jakob Ejler Sorensen
Version: 5.0
Date: 2020 Oct 21
    """

    def __init__(self, filename, mass=4, theta=146.7, E0=1000, default_scan=0):
        """Initialize the class"""
        # Constants
        self.settings = dict()
        self.settings['mass'] = mass
        self.settings['theta'] = theta
        self.settings['E0'] = E0
        self.default_scan = default_scan

        # Initialize variables
        self.energy = dict()
        self.cps = dict()
        self.dwell = dict()
        self.mode = dict()
        self.mode_value = dict()
        self.note = dict()
        filename = str(filename)
        self.filename = filename

        # Convenience function variables
        self.peak_positions = None
        self.peak_heights_raw = None
        self.peak_heights_bg = None
        self._background = None
        self.background_settings = {
            'type': None,
            'ranges': None,
            'on': False,
            }

        #----------------------------------------------------------------------
        # Read data from textfile:
        if filename.endswith('.txt'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            self.format = 'Text file'

            start_points = [i for i, line in enumerate(lines) if line == 'Energy\tCounts\r\n']
            self.scans = len(start_points)
            if self.scans == 0:
                raise ImportError('File apparently empty!')

            if lines[0].lower().startswith('note'):
                self.note = lines[0].split('=')[1].lstrip(' ')
            # Copy data points
            counter = 0
            for start in start_points:
                line = lines[start-4].split('\t')
                for i, word in enumerate(line):
                    if word.lower() == 'dwell':
                        self.dwell[counter] = float(lines[start-3].split('\t')[i])
                if not start == start_points[-1]:
                    interval = range(start+1, start_points[counter+1]-4)
                else:
                    interval = range(start+1, len(lines))
                self.energy[counter] = np.zeros(len(interval))
                self.cps[counter] = np.zeros(len(interval))
                counter_inner = 0
                for index in interval:
                    line = lines[index].rstrip().split('\t')
                    self.energy[counter][counter_inner] = float(line[0])
                    self.cps[counter][counter_inner] = float(line[1])
                    counter_inner += 1
                self.cps[counter] = self.cps[counter]/self.dwell[counter]
                counter += 1
        #----------------------------------------------------------------------
        # Read data from old VAMAS block file
        elif filename.endswith('.vms'):
            # Open filename with ISS data
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()
            # Old format:
            if lines[6].lower().startswith('experiment type'):
                self.setup = 'omicron'
                self.format = 'Old VAMAS'
                print('Loading file: ' + filename)
                blocks_4 = [i for i, line in enumerate(lines) if (line.strip() == '-1') \
and (lines[i+1].lower().strip() == 'kinetic energy')]
                blocks_2_ISS = [i for i, line in enumerate(lines) if (line.strip() == 'ISS') \
and (lines[i+1].strip() == '')]
                print(lines[9].rstrip())
                self.scans = len(blocks_4)
                if len(blocks_4) == int(lines[9].rstrip()) \
and len(blocks_4) == len(blocks_2_ISS):
                    self.scans = len(blocks_4)
                else:
                    msg = 'Error: Identified {} "Block 4", {} "Block 2", but "Block 1" says: {}'
                    msg = msg.format(len(blocks_4), len(blocks_2_ISS), int(lines[9].rstrip()))
                    raise ImportError(msg)

                # Copy data points
                self.note = dict()
                for counter, block in enumerate(blocks_4):
                    if not len(lines[blocks_2_ISS[counter] - 1]) == 5:
                        self.note[counter] = lines[blocks_2_ISS[counter] - 1].rstrip()
                    else:
                        self.note[counter] = ''
                    self.mode[counter] = lines[block-11].rstrip()
                    self.mode_value[counter] = float(lines[block-10].rstrip())
                    self.dwell[counter] = float(lines[block+9].rstrip())
                    data_points = int(lines[block+16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[block+4].rstrip())
                    E_start = float(lines[block+3].rstrip())
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float(lines[block+19+counter_inner]) \
/self.dwell[counter]
                self.note[counter] = ''
                print(self.energy.keys())
                print('Comments: {}'.format(self.note))
                print('Dwell time: {}'.format(self.dwell))
                print('Modes: {}'.format(self.mode))
                print('Mode values: {}'.format(self.mode_value))
            #----------------------------------------------------------------------
            # New format
            if lines[6].lower().startswith('created with'):
                self.setup = 'omicron'
                self.format = 'New VAMAS'
                ENDING = '_1-Detector_Region.vms'
                filen = filename.rstrip(ENDING)
                counter = 0
                while True:
                    try:
                        f = open(filen + '--' + str(counter+1) + ENDING)
                        counter += 1
                    except:
                        #print('{} files detected of series:'.format(counter))
                        #print('* ' + filen + '--' + str(1) + ENDING + ' *')
                        COUNTER = counter
                        break
                # Open filename with ISS data
                self.scans = COUNTER
                for counter in range(COUNTER):
                    new_filename = filen + '--' + str(counter+1) + ENDING
                    f = open(new_filename, 'r')
                    lines = f.readlines()
                    f.close()
                    print('Loading file: ' + new_filename)
                    blocks_4 = [i for i, line in enumerate(lines) if (line.rstrip() == '-1') \
and (lines[i+1].lower().rstrip() == 'kinetic energy')]
                    #print(lines[9].rstrip())
                    if len(blocks_4) > 1:
                        print('*** Interesting! More than 1 scan has been detected in above file!')
                    # Copy data points
                    i = blocks_4[0]
                    ###
                    if counter == 0:
                        _counter = 0
                        while True:
                            if lines[_counter].startswith('CREATION COMMENT START'):
                                comment_start = _counter
                                break
                            else:
                                _counter += 1
                                if _counter > len(lines):
                                    break
                        _counter = 0
                        while True:
                            if lines[_counter].startswith('CREATION COMMENT END'):
                                comment_end = _counter
                                break
                            else:
                                _counter += 1
                                if _counter > len(lines):
                                    break
                    self.note = lines[comment_start+1:comment_end]
                    ###
                    self.mode[counter] = lines[i-11].rstrip()
                    self.mode_value[counter] = float(lines[i-10].rstrip())
                    self.dwell[counter] = float(lines[i+9].rstrip())
                    data_points = int(lines[i+16])
                    self.cps[counter] = np.zeros(data_points)
                    E_step = float(lines[i+4].rstrip())
                    E_start = float(lines[i+3].rstrip())
                    self.energy[counter] = np.arange(data_points)*E_step + E_start
                    for counter_inner in range(data_points):
                        self.cps[counter][counter_inner] = float(lines[i+19+counter_inner]) \
/self.dwell[counter]
        #----------------------------------------------------------------------
        # Import Thetaprobe .avg data
        elif filename.endswith('.avg'):
            self.setup = 'thetaprobe'
            with open(filename, 'r', encoding='latin-1') as f:
                lines = f.readlines()

            # Check for ISS
            info = {line.split(':')[0].strip(): line.split('=')[1].strip() for line in lines if line.startswith('DS_')}
            if info['DS_ANPROPID_LENS_MODE_NAME'] != "'ISS'":
                print('{} does not appear to be an ISS experiment!'.format(self.filename))
                print('Expected \'ISS\', but encountered: {}'.format(info['DS_ANPROPID_LENS_MODE_NAME']))
                raise ImportError('File not an ISS experiment!')
            if info['DS_EXT_SUPROPID_CREATED'] == info['DS_EXT_SUPROPID_SAVED']:
                #print('Created and saved dates are identical - checking for empty dataset...')
                check_empty = True
            else:
                check_empty = False

            # Metadata
            self.note[0] = info['DS_EXT_SUPROPID_SUBJECT']
            self.date = info['DS_EXT_SUPROPID_CREATED']
            self.dwell[0] = float(info['DS_ACPROPID_ACQ_TIME'])
            self.mode[0] = int(info['DS_ANPROPID_MODE'])
            self.mode_value[0] = float(info['DS_ANPROPID_PASS'])
            if info['DS_GEPROPID_VALUE_LABEL'] == "'Counts'":
                normalize = True # normalize to "counts per second"
            else:
                normalize = False

            # Data
            #data_info = {}
            line_number = [i for i, line in enumerate(lines) if line.startswith('$DATAAXES')]
            if len(line_number) > 1:
                print('Reading file: {}'.format(self.filename))
                raise ImportError('Import of multiple dataaxes not implemented yet!')
            else:
                line_number = line_number[0]
            keys = [key.strip() for key in lines[line_number-1].split('=')[1].split(',')]
            values = [key.strip() for key in lines[line_number+1].split('=')[1].split(',')]
            data_info = {key: value for key, value in list(zip(keys, values))}

            start, end = float(data_info['start']), float(data_info['end'])

            #space_info = {}
            line_number = [i for i, line in enumerate(lines) if line.startswith('$SPACEAXES')]
            if len(line_number) > 1:
                print('Reading file: {}'.format(self.filename))
                raise ImportError('Import of multiple dataaxes not implemented yet!')
            else:
                line_number = line_number[0]
            keys = [key.strip() for key in lines[line_number-1].split('=')[1].split(',')]
            values = [key.strip() for key in lines[line_number+1].split('=')[1].split(',')]
            space_info = {key: value for key, value in list(zip(keys, values))}

            num = int(space_info['numPoints'])
            if space_info['linear'] != 'LINEAR':
                print('Reading file: {}'.format(self.filename))
                raise ImportError('Check .avg file if energy axis is linear!')

            # Generate xy-data
            self.energy[0] = np.linspace(start, end, num)
            self.cps[0] = self.energy[0]*np.nan

            line_number = [i for i, line in enumerate(lines) if line.startswith('$DATA=')]
            if len(line_number) > 1:
                msg = 'Reading file: {}'.format(self.filename)
                raise ImportError('Import of multiple dataaxes not implemented yet!')
            else:
                line_number = line_number[0]

            for j in range(num):
                if j%4 == 0: # values are grouped in chunks of 4
                    line_number += 1
                    line = lines[line_number].split('=')[1].split(',')
                try:
                    self.cps[0][j] = float(line[j%4])
                except ValueError:
                    pass # #empty# values
            if check_empty:
                if not np.any(np.isfinite(self.cps[0])):
                    raise ImportError('Dataset from {} is empty!'.format(self.filename))
                else:
                    print('Dataset appeared to be empty from the saved timestamps, but is not empty.')
            if normalize:
                self.cps[0] /= self.dwell[0]
        else:
            raise IOError('File: "{}" not found or fileending not accepted.'.format(self.filename))

        # Print loaded settings
        print('Successfully loaded file: {}'.format(filename))
        string = 'Used settings:\nProbing mass: {} amu\nScatter angle: {}\nPrimary energy: {} eV'
        #print(string.format(*[self.settings[key] for key in ['mass', 'theta', 'E0']]))

    @property
    def x(self):
        return self.energy[self.default_scan]

    @x.setter
    def x(self, var):
        if not var in self.energy.keys():
            print('"{}" not an available key! {}'.format(var, self.energy.keys()))
        self.default_scan = var

    @property
    def y(self):
        return self.cps[self.default_scan]

    @y.setter
    def y(self, var):
        if not var in self.energy.keys():
            print('"{}" not an available key! {}'.format(var, self.energy.keys()))
        self.default_scan = var

    @property
    def xy(self):
        return np.vstack((self.x, self.y)).T

    def get_xy(self, index):
        return np.vstack((self.energy[index], self.cps[index])).T

    @property
    def background(self):
        if self._background is not None:
            return self._background[self.default_scan]
        else:
            return None


    def ConvertEnergy(self, mass):
        """Converts a measured energy to mass of surface atom
corresponding the settings stored in the experiment.
        """
        angle = self.settings['theta'] * np.pi/180
        return self.settings['E0'] * ((self.settings['mass']*np.cos(angle) + \
np.sqrt(mass**2 - self.settings['mass']**2*np.sin(angle)**2))/(mass + self.settings['mass']))**2


    def PlotAllScans(self, exclude=[None], color=None):
        """Plot all elements in file in single figure."""
        selection = [i for i in self.energy.keys() if not i in exclude]
        if not color:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i])
        else:
            for i in selection:
                plt.plot(self.energy[i], self.cps[i], color=color)
        plt.xlabel('Kinetic energy (eV)')
        plt.ylabel('Counts per second')


    def Normalize(self, interval='Max', exclude=[None], unit='Mass', delta_e=10):
        """Normalize to highest value in interval=[value1, value2]"""
        self.delta_e = delta_e
        if isinstance(interval, int):
            self.normalization_criteria = interval
        elif isinstance(interval, str):
            if interval == 'Total':
                self.normalization_criteria = 'all'
            elif interval.lower().startswith('max'):
                self.normalization_criteria = 'max'
            elif interval == 'Au':
                self.normalization_criteria = 196.
        if not isinstance(interval, list):
            if self.normalization_criteria == 'all':
                selection = [i for i in range(self.scans) if not i in exclude]
                for __counter in selection:
                    total = simps(self.cps[__counter], self.energy[__counter])
                    self.cps[__counter] /= total
            elif self.normalization_criteria == 'max':
                selection = [i for i in range(self.scans) if not i in exclude]
                for __counter in selection:
                    ydata = ct.smooth(self.cps[__counter], width=2)
                    norm_value = max(ydata)
                    self.cps[__counter] /= norm_value
            else:
                interval = [0, 0]
                if unit.lower() == 'mass':
                    interval[0] = self.ConvertEnergy(self.normalization_criteria) - self.delta_e
                    interval[1] = self.ConvertEnergy(self.normalization_criteria) + self.delta_e
                elif unit.lower() == 'energy':
                    interval[0] = self.normalization_criteria - self.delta_e
                    interval[1] = self.normalization_criteria + self.delta_e
                selection = [i for i in range(self.scans) if (not i in exclude) and \
                             (not interval[0] > max(self.energy[i])) and (not interval[1] < min(self.energy[i]))]
                for __counter in selection:
                    range_1 = np.where(self.energy[__counter] < interval[1])[0]
                    range_2 = np.where(self.energy[__counter] > interval[0])[0]
                    energy_range = np.intersect1d(range_1, range_2)
                    value = max(self.cps[__counter][energy_range])
                    self.cps[__counter] = self.cps[__counter]/value


    def AddMassLines(self, masses, ax=None, offset=0, color='k', labels=True, linestyle='dotted', **kwargs):
        """Add vertical lines for mass references."""
        energies = self.ConvertEnergy(np.array(masses))
        if ax is None:
            ax = plt.gca()
        [x1, x2, y1, y2] = ax.axis()
        for energy, mass in zip(energies, masses):
            ax.axvline(x=energy-offset, ymin=0, ymax=1, linestyle=linestyle, color=color, **kwargs)
            if labels:
                ax.text(float(energy)/x2, 0.95, 'm-{}'.format(mass), transform=ax.transAxes)


    def AddRegions(self):
        """Add regions indicating the whereabouts of 3d, 4d, 5d metals and the
lanthanides and actinides."""
        ax = plt.gca()
        d3 = [45, 65]
        d4 = [89, 112]
        d5 = [178, 201]
        lant = [139, 175]
        act = [227, 260]
        for i in [d3, d4, d5]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]),
                       color='k', alpha=0.2)
        for i in [lant, act]:
            ax.axvspan(xmin=self.ConvertEnergy(i[0]), xmax=self.ConvertEnergy(i[1]),
                       color='y', alpha=0.2)


# Convenience functions
def iss_to_time_axis(iss_data=[], ax=None):
    """Fold a series of iss measurements unto a time axis to separate them graphically"""

    colors = ['k', 'r', 'g', 'b', 'm', 'y', 'c']

    # Allow to not pass single argument as a list object
    if isinstance(iss_data, np.ndarray):
        iss_data = [iss_data]

    # Create axes object
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    # Plot iss
    separator = 0
    print('Time-plot of ISS made:\n-----')
    for j, iss in enumerate(iss_data):
        print(j, iss.filename)
        for i in iss.energy.keys():
            x, y = iss.energy[i], iss.cps[i]
            ax.plot(x + separator, y, color=colors[j])
            separator += max(x) + 10
    ax.set_xlabel('Energy (eV) -> Time')
    ax.set_ylabel('Signal (cts/s)')

    # return handles to figure
    fig = plt.gcf()
    return fig, ax

def get_range(X, limit):
    """Return index of list X where limit[0] < x < limit[1] """
    try:
        A = np.where(X > limit[0])[0][0]
        B = np.where(X < limit[1])[0][-1]
        return np.arange(A, B)
    except IndexError:
        return None

def get_peaks(iss_data=[], peaks=None):
    """FIXME"""
    AVG = 3.5
    for iss in iss_data:
        iss.peak_positions = peaks
        iss.peaks_raw = dict()
        if iss.background_settings['on']:
            iss.peaks_bg = dict()
        for key in iss.energy.keys():
            iss.peaks_raw[key] = dict()
            if iss.background_settings['on']:
                iss.peaks_bg[key] = dict()
            for element, pos in peaks.items():
                indice = get_range(iss.energy[key], [pos-AVG, pos+AVG])
                if indice is None:
                    #print('Peaks', iss.filename, key, (element, pos), 'None')
                    average = -500
                    average = np.nan
                else:
                    average = np.average(iss.cps[key][indice])
                    #print('Peaks', iss.filename, key, (element, pos), average)
                iss.peaks_raw[key][element] = average
                if iss.background_settings['on']:
                    if indice is None:
                        iss.peaks_bg[key][element] = np.nan
                    else:
                        iss.peaks_bg[key][element] = iss.peaks_raw[key][element] - \
iss._background[key][int(np.average(indice))]


class LoadSet():
    """Load a set of ISS data """
    def __init__(self, input_files=[]):
        """input_files (list of dicts/dict): 'path' to filenames"""
        if isinstance(input_files, dict):
            self.input_files = [input_files]
        elif isinstance(input_files, list):
            #if not isinstance(input_files[0], dict):
            #    raise InputError('input_files must be a list of dictionaries!')
            self.input_files = input_files

        # Data variables
        self.data = {}

    def add_data(self, input_files):
        """Add data to load."""
        if isinstance(input_files, dict):
            self.input_files.append(input_files)
        elif isinstance(input_files, list):
            if not isinstance(input_files[0], dict):
                raise ValueError('input_files must be a list of dictionaries!')
            for item in input_files:
                self.input_files.append(item)

    def load(self, normalize=None, exclude=[], unit='Mass'):
        """Compile and return dictionary of data"""

        # Load data
        loaded_files = [x.filename for x in self.data.values()]
        for input_file in self.input_files:
            path = input_file['path']
            if len(path) > 0:
                if path[-1] != '/':
                    path += '/'
            for key, filename in input_file.items():
                if key == 'path':
                    continue
                if path+filename in loaded_files:
                    continue
                else:
                    print('Input file: {}'.format(path+filename))
                self.data[key] = Experiment(path + filename)

        # Normalize data if chosen
        if not normalize is None:
            for key in self.data:
                self.data[key].Normalize(normalize, exclude, unit)

        # Return data
        return self.data

def plot(data, index, ax=None, args=[], kwargs={}):
    """Plot index from Experiment 'data'"""
    if ax is None:
        ax = plt.gca()
    x, y = data.energy[index], data.cps[index]
    ax.plot(x, y, *args, **kwargs)

def subtract_single_background(xy, ranges=[], avg=3):
    """Subtract the background from a single spectrum"""
    x, y = xy[:, 0], xy[:, 1]
    background = np.copy(y)
    for limit in ranges:
        indice = get_range(x, limit)
        # if first index is chosen
        # OR
        # if last ten indice are included
        if indice[0] == 0 or indice[-1] > len(x) - 10:
            print('Uhh', indice[0], indice[-1], limit)
            background[indice] = 0
        elif len(indice) == 0:
            print('Did not find data within limit: {}'.format(limit))
        else:
            y1 = np.average(y[indice[0]-avg:indice[0]+avg])
            y2 = np.average(y[indice[-1]-avg:indice[-1]+avg])
            a_coeff = (y2-y1)/(limit[1]-limit[0])
            b_coeff = y1 - a_coeff*limit[0]
            #print(y1, y2, a_coeff, b_coeff)
            background[indice] = x[indice]*a_coeff + b_coeff
    return background


def subtract_backgrounds(iss_data=[], ranges=[], btype='linear', avg=3):
    """Subtract a linear background from defined 'ranges'.
Return data above backgrounds."""

    AVG = avg + 0.5
    for iss in iss_data:
        iss._background = dict()
        for key in iss.energy.keys():
            iss._background[key] = subtract_single_background(iss.get_xy(key), ranges, avg)

        # Tell settings that function is done
        iss.background_settings['type'] = btype
        iss.background_settings['ranges'] = ranges
        iss.background_settings['on'] = True

def reorganize_peaks(iss_data=[], bg=True):
    """Take the list of dicts() of peaks and gather into a dict of elements
TODO: add possibility of only selecting specific sets (iss.energy[i])"""
    peaks = dict()
    for element in iss_data[0].peak_positions.keys():
        if bg:
            peaks[element] = [iss.peaks_bg[i][element] for iss in iss_data \
for i in iss.peaks_bg.keys()]
            peaks[element] = np.array(peaks[element])
        else:
            peaks[element] = [iss.peaks_raw[i][element] for iss in iss_data \
for i in iss.peaks_raw.keys()]
            peaks[element] = np.array(peaks[element])
    total = np.zeros(len(peaks[element]))
    for element in peaks.keys():
        total += peaks[element]
    return peaks, total

def convolute_peaks(iss_data=[], peaks=[], sigma=[30, 30], amp=[1000, 1000]):
    """
Take a set of iss data with background subtraction and fit a dictionary of peak positions.
Return: None. Peaks are saved as attributes of ISS object under 'peak_fit' as a dict.
"""
    def gauss(x, x0, sigma, amp):
        """Gaussian function"""
        return amp * np.exp(-((x - x0)**2)/(2*sigma**2))

    def fit_gauss(x, y, x01, x02, sigma1, sigma2, amp1, amp2):

        def double_gauss(x, sigma1, sigma2, amp1, amp2):
            return gauss(x, x01, sigma1, amp1) + gauss(x, x02, sigma2, amp2)

        param = [sigma1, sigma2, amp1, amp2]
        ret, _null = curve_fit(double_gauss, x, y, p0=param)
        return ret

    for iss in iss_data:
        iss.peak_fit = dict()
        for i in range(iss.scans):
            signal = iss.cps[i] - iss._background[i]
            energy = iss.energy[i][:]
            iss.peak_fit[i] = fit_gauss(energy, signal, peaks[0], peaks[1], sigma[0], sigma[1], \
amp[0], amp[1])
            print(iss.peak_fit[i])

### NEW ###
def fit_references(data, references, e_ranges, labels=None):
    """Fit 'references' to 'data' in the energy window (x) 'e_range'"""
    # Restrict data to within e_range
    index = []
    for e_range in e_ranges:
        index.append(ct.get_range(data[:, 0], *e_range))
    index = functools.reduce(np.union1d, index)
    #return index
    data = data[index]
    mod_ref = {label: ref[index] for label, ref in references.items()}

    # Define fitting relationship
    def func(x, *args):
        signal = x*0
        for arg, ref in list(zip(args, mod_ref.values())):
            signal += arg*ref[:, 1]
        return signal

    # Calculate best fit
    ret, _ = curve_fit(func, data[:, 0], data[:, 1], p0=[0.5]*len(references), bounds=(0, np.inf))
    dict_ = {label: ret[i] for i, label in enumerate(references.keys())}
    pretty_dict = dict_as_percent(dict_)
    for label, percent in pretty_dict.items(): print(label, percent)
    return dict_

def fit_references_v2(iss, references, e_ranges, labels=None):
    """Fit 'references' to 'data' in the energy window (x) 'e_range'

Subtract background from references to only use peak intensities for fit.
Subtract beckground in the regions you want to fit.
"""
    #plt.show()
    # Restrict data to within e_range
    index = []
    for e_range in e_ranges:
        index.append(ct.get_range(iss.x, *e_range))
    index = functools.reduce(np.union1d, index)

    # Subtract background from references
    modified_references = {}
    for label, (ref, limit) in references.items():

        
        modified_references[label] = np.copy(ref)
        #modified_references[label][:, 1] -= peak[index]

        background = subtract_single_background(ref, ranges=[limit])
        modified_references[label][:, 1] -= background
        #modified_references[label] = modified_references[label][index]
        #plt.plot(ref[:, 0], ref[:, 1] - background, label=label)

    # Subtract background from data
    subtract_backgrounds(iss_data=[iss], ranges=e_ranges)
    data = iss.xy[index]
    data[:, 1] -= iss.background[index]
    #plt.plot(data[:, 0], data[:, 1], 'ko', label='Data')
    #plt.plot(iss.x, iss.background, label='Data (raw)')


    #plt.legend()
    #plt.show()
    #modified_references = {label: ref[index] for label, ref in references.items()}

    # Define fitting relationship
    def func(x, *args):
        signal = x*0
        for arg, ref in list(zip(args, modified_references.values())):
            signal += arg*ref[:, 1][index]
        return signal

    # Calculate best fit
    ret, _ = curve_fit(func, data[:, 0], data[:, 1], p0=[0.5]*len(references), bounds=(0, np.inf))
    dict_ = {label: ret[i] for i, label in enumerate(references.keys())}
    #for label, percent in pretty_dict.items(): print(label, percent)
    return dict_, modified_references

def align_spectra(iss_data, limits=[350, 520], masses=[16, 18], key='oxygen'):
    """Shift the iss data within 'limits' region to snap maximum signal
unto nearest mass in list 'masses'."""

    plt.figure('aligned')
    for data in iss_data:
        # Initialize attributes
        try:
            data.shifted
        except AttributeError:
            data.shifted = {}
        try:
            data.smoothed
        except AttributeError:
            data.smoothed = {}

        # Get index of region of interest
        index = ct.get_range(data.x, *limits)
        # Find maximum in region
        ys = ct.smooth(data.y, width=4)
        data.smoothed[key] = ys
        maximum = max(ys[index])
        if not np.isfinite(maximum):
            data.good = False
            continue # "Broken" dataset
        data.good = True
        i_max = np.where(ys == maximum)[0]
        x_max = data.x[i_max]

        # Find difference from reference
        energies = data.ConvertEnergy(np.array(masses))
        try:
            distance = np.abs(x_max - energies)
        except ValueError:
            print('ValueError')
            print(data.filename)
            print(data.x)
            print(data.y[np.isfinite(data.y)])
            print(data.smoothed[key])
            print(maximum)
            print(x_max)
            print(energies)
            plt.figure()
            plt.plot(data.x, data.y)
            plt.show()
            raise

        # Snap to nearest line
        data.shifted[key] = data.x - (x_max - energies)[np.where(distance == min(distance))[0]]
        plt.plot(data.shifted[key], data.y)
        plt.plot(data.shifted[key], ys)
    data.AddMassLines(masses)

    # Return a nd.array of modified xy data
    return np.vstack((data.shifted[key], data.smoothed[key])).T

def from_fit(references, fit):
    y = None
    for label, value in references.items():
        if y is None:
            y = value[:, 1] * fit[label]
        else:
            y += fit[label] * value[:, 1]
    return y

def interpret_fit(dictionary):
    total = sum([val for val in dictionary.values()])
    new_dict = {label: value/total*100 for label, value in dictionary.items()}
    results = {}
    results['info'] = """Interpretation of keys:
    'isotope_ratio': Percentage of O-18 of total oxygen present
    'oxide_ratio': Percentage of total oxygen signal compared to the Ru (oxide) signal
    'metallic_ratio': Percentage metallic Ru signal compared to total Ru signal
"""
    print(new_dict)

    O18 = dictionary['O-18']/(dictionary['O-18'] + dictionary['O-16'])
    print('Percent 18-O: {} %'.format(O18*100))
    results['isotope_ratio'] = O18*100

    #oxide = (dictionary['O-16'] + dictionary['O-18']) / dictionary['Metallic2']
    #print('Degree of oxide: {} %'.format(oxide*100))
    #results['oxide_ratio'] = oxide*100

    #metallic = dictionary['Metallic1'] / (dictionary['Metallic1'] + dictionary['Metallic2'])
    #print('Degree metallic: {} %'.format(metallic*100))
    #results['metallic_ratio'] = metallic*100
    #print()
    return new_dict, results
    

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 2:
        print(sys.argv[1])
        experiment = Experiment(sys.argv[1])
        experiment.PlotAllScans()
        plt.show()
    else:
        print('Pass a single ISS file as extra argument to plot all the data scans')
