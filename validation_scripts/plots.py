'''
Making bar plots from columns in a tab delimited file
Call this as :
>> from validation_scripts import plots
>> bp = plots.BarPlotter("FATCAT_STATS")
>> bp.bar_plots(0, [2,5,7]) # where 0 is the index to the x value (in this case the first column in the file)
# and [2,5,7] are the y indices to this plot (the respective columns)

To dos:
Currently the y bars are on top of each other.
Allowing different colors for different bars.
The option to save is included, but not completed
'''
import os, sys
import numpy as np
import matplotlib                       # for plotting files
import matplotlib.pyplot as plt         # shortcut for pyplot
import matplotlib.mlab as mlab

class BarPlotter():
    '''
    Base class for bar plots.
    '''
    def __init__(self, filename):
        '''
        Input argument to the constructor is the filename.
        Make sure that the filename is valid, that the first row contains title labels,
        and that it contains tab delimited values of equal length.
        '''
        self.file_name = filename
        # Get columns into single lists from the file
        self.file_list = [line.strip() for line in open(self.file_name).readlines()]
        self.file_list = [line.split("\t") for line in self.file_list]
        self.file_list = np.transpose(np.array(self.file_list))
        self.file_list = [[listvals[0], list(listvals[1:])] for listvals in self.file_list]
        # print self.file_list # to see what self.file_list looks like here

    def bar_plots(self, x_index, y_indices_list, save=True):
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        x_vals = self.file_list[x_index][1]
        #print x_vals
        ind = np.arange(len(x_vals)) # Gives the start for each bar in a group
        orig_ind = ind[:] # Keep a copy of this location for later
        width = 0.35 # width of each bar
        y_bars = []
        for y_index in y_indices_list:
            y_vals = [float(val) for val in self.file_list[y_index][1]]
            rects = ax.bar(ind, y_vals, width)
            y_bars.append(rects)
            ind += width
        ax.set_xticks(orig_ind)
        ax.set_xticklabels(x_vals, rotation = 30)
        plt.show()


if __name__ == "__main__":
    bp = BarPlotter ('/home/ajith/Documents/Write/Programming/Python/SjolanderLab/FAT_CAT/Previous_version/FATCAT_STATS')
    bp.bar_plots(0, [1, 2, 1])
    
