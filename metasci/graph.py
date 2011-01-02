"""Functions to aid in graphing."""

#Requires MatPlotLib
#from pylab import *
import pylab as plt
from matplotlib import rc

rc('text', usetex=True)
rc('font', family='roman')

### Other MetaSci functions ###
from . import StairStep 
from . import Nines2Log
from . import LogLabelConvert

from . import LogUniformBins
from . import LinearUniformBins

import mathematics as msm


BWLines = ['k-', 'r-', 'k--', 'r--', 'k-.', 'r-.', 'k:', 'r:', 'k.', 'r.', 'k+', 'r+', 'kx', 'rx']
BWLinesLen = len(BWLines)

###################################
### Figure out the right values ###
###################################

def set_defaults(**kwgraph):
    """Sets undefined defaults.

    Keyword Args:
        * `name` (str): Figure name.
        * `xlabel` (str): x-axis label.
        * `ylabel` (str): y-axis label.
        * `datalabel` (str or list):  data set label or list of them (['data set 1 label', 'data set 2 label', ...]).
        * `colorline` (str): line type and color, ie 'k-'.
        * `axis` (4-list): A list defining the axis limits, [xmin, xmax, ymin, ymax].
        * `xaxis` (2-list): A list defining the x-axis limits, [xmin, xmax].
        * `yaxis` (2-list): A list defining the y-axis limits, [ymin, ymax].
        * `xmin` (numeric): The x-axis minimum, xmin.
        * `xmax` (numeric): The x-axis maximum, xmin.
        * `ymin` (numeric): The y-axis minimum, xmin.
        * `ymax` (numeric): The y-axis maximum, xmin.
        * `scale` (str): axis scale, ie 'linear' or 'log' or 'xlog' or 'ylog'.
        * `base` (int or float): Base of `scale`, ie 10.0.
        * `unit` (str): String for data's units, ie '[unitless]'. 
        * `write` (True or False): Boolean on whether to write figure to disk.
        * `show` (True or False): Boolean to display the figure.

    Returns:
        * `kwgraph` (dict): An updated keyword graph dictionary.    
    """

    if 'name' not in kwgraph.keys():
        kwgraph['name'] = 'Plot'

    if 'xlabel' not in kwgraph.keys():
        kwgraph['xlabel'] = ''

    if 'ylabel' not in kwgraph.keys():
        kwgraph['ylabel'] = ''

    if not 'datalabel' in kwgraph.keys():
        kwgraph['datalabel'] =  ""

    if not 'colorline' in kwgraph.keys():
        kwgraph['colorline'] =  "k-"

    if 'axis' not in kwgraph.keys():
        kwgraph['axis'] = []

    if 'xaxis' not in kwgraph.keys():
        kwgraph['xaxis'] = []

    if 'yaxis' not in kwgraph.keys():
        kwgraph['yaxis'] = []

    if 'xmin' not in kwgraph.keys():
        kwgraph['xmin'] = None

    if 'xmax' not in kwgraph.keys():
        kwgraph['xmax'] = None

    if 'ymin' not in kwgraph.keys():
        kwgraph['ymin'] = None

    if 'ymax' not in kwgraph.keys():
        kwgraph['ymax'] = None

    if 'scale' not in kwgraph.keys():
        kwgraph['scale'] = 'linear'

    if 'base' not in kwgraph.keys():
        kwgraph['base'] = 10.0

    if 'unit' not in kwgraph.keys():
        kwgraph['unit'] = ''

    if 'write' not in kwgraph.keys():
        kwgraph['write'] = True

    if 'show' not in kwgraph.keys():
        kwgraph['show'] = False

    return kwgraph

def set_axis(axis=[], xaxis=[], yaxis=[], xmin=None, xmax=None, ymin=None, ymax=None, **kwgraph):
    """Sets the axes of the current plot.

    Keyword Args:
        * `axis` (4-list): A list defining the axis limits, [xmin, xmax, ymin, ymax].
        * `xaxis` (2-list): A list defining the x-axis limits, [xmin, xmax].
        * `yaxis` (2-list): A list defining the y-axis limits, [ymin, ymax].
    """

    if axis:
        plt.axis(axis)
    else:
        ax = list(plt.axis())

        if xaxis:
            ax[0] = xaxis[0]
            ax[1] = xaxis[1]
        elif xmin != None:
            ax[0] = xmin
        elif xmax != None:
            ax[1] = xmax

        if yaxis:
            ax[2] = yaxis[0]
            ax[3] = yaxis[1]
        elif ymin != None:
            ax[2] = ymin
        elif ymax != None:
            ax[3] = ymax

        plt.axis(ax)

    return

def set_legend_label(xlabel='', ylabel='', datalabel='', **kwgraph):
    """Sets the legend, x-axis, and y-axis labels when available.

    Keyword Args:
        * `xlabel` (str): x-axis label, ignores if unset.
        * `ylabel` (str): y-axis label, ignores if unset.
        * `datalabel` (str or list):  data set label or list of them (['data set 1 label', 'data set 2 label', ...]).
          Legend does not appear if not set.
    """

    if xlabel:
        plt.xlabel(xlabel)

    if ylabel:
        plt.ylabel(ylabel)

    if datalabel:
        plt.legend(loc=0)

    return

def set_scale(scale='linear', base=10.0, **kwgraph):
    """Sets the scale on the figure.

    Keyword Args:
        * `scale` (str): axis scale, ie 'linear', 'log', 'xlog', 'ylog', 'logx', 'logy' or 'nines'.
        * `base` (int or float): Base of scale, ie 10.0.
    """
    if 'log' == scale:
        plt.xscale('log', basex=base)
        plt.yscale('log', basey=base)
    elif scale in ['xlog', 'logx']:
        plt.xscale('log', basex=base)
        plt.yscale('linear')
    elif scale in ['ylog', 'logy']:
        plt.xscale('linear')
        plt.yscale('log', basey=base)
    else:
        plt.xscale('linear')
        plt.yscale('linear')
    return

##########################################
### Function to clear and write graphs ###
##########################################
def shortname(figstr):
    """Returns the short name of a figure.

    Args:
        * `figstr` (str): Figure string, like 'project/figs/an_image.png'

    Returns:
        * `shortname` (str): Figure name without the path, ie 'an_image.png'.
    """
    return figstr.split("/")[-1]

def writeGraph(name='Plot', write=True, show=False, **kwgraph):
    """Saves and clears the current figure.

    Keyword Args:
        * `name` (str): Figure name.

    Returns:
        * `name` (str): The location name of this figure.
    """

    if write:
        plt.savefig(name + ".eps")
        plt.savefig(name + ".png")

    if show:
        plt.show()

    if (write or show):
        plt.clf()

    return name


##################################################
### Functions that graph to the current figure ###
##################################################

def SimpleGraph(x, y, **kwgraph):
    """Graphs a single x-y data plot.

    Args:
        * `x` (numeric sequence): Independnet variable data set.
        * `y` (numeric sequence): Dependnet variable data set.

    Keyword Args:
        * `kwgraph` (dict): A dictionary of keywords that are used to specify the figure.
          Please see the set_defaults() function.
    """
    kwgraph = set_defaults(**kwgraph)

    plt.plot(x, y, kwgraph['colorline'], label=kwgraph['datalabel'])
    set_scale(**kwgraph)
    set_axis(**kwgraph)
    set_legend_label(**kwgraph)

    writeGraph(**kwgraph)

    return 


###################
### Nuke Grpahs ###
###################

def StairStepEnergy(data, energy_bins=[], G=3, elower=(10.0**-9), eupper=10.0, **kwgraph):
    """Plots nuclear data which is separated into G energy bins.

    Args:
       * `data` (numeric sequence) Data Set; length G [units].

    Keyword Args:
       * `energy_bins` (numeric sequence or None): Energy Bins; length G+1 [MeV].
       * `G` (int): Number of energy groups.
       * `elower` (float): Lower Energy Bound [MeV], if energy_bins not given.
       * `eupper` (float): Upper Energy Bound [MeV], if energy_bins not given.
       * `kwgraph` (dict): Other graphing keywords.
    """
    if 'scale' not in kwgraph.keys():
        kwgraph['scale'] = 'log'
    if 'xlabel' not in kwgraph.keys():
        kwgraph['xlabel'] = 'Energy [MeV]'

    kwgraph = set_defaults(**kwgraph)

    if not energy_bins:
        if kwgraph['scale'] == 'linear':
            energy_bins = LinearUniformBins(elower, eupper, G)
        elif kwgraph['scale'] == 'log':
            energy_bins = LogUniformBins(elower, eupper, G)
        else:
            print "The energy bins could not be calculated because the scale variable is not valid:", kwgraph['scale']

    #Determine axes bounds
    elow  = kwgraph['base']**msm.orderfloor(energy_bins[0], kwgraph['base'])
    ehigh = kwgraph['base']**msm.orderceil(energy_bins[-1], kwgraph['base'])
    dlow  = kwgraph['base']**msm.orderfloor(msm.min_above(data, 0.0), kwgraph['base'])        
    dhigh = kwgraph['base']**msm.orderceil(max(data), kwgraph['base'])        

    if None != kwgraph['xmin'] < elow:
        elow = kwgraph['xmin']
    if ehigh < kwgraph['xmax'] != None:
        ehigh = kwgraph['xmax']

    if None != kwgraph['ymin'] < dlow:
        dlow = kwgraph['ymin']
    if dhigh < kwgraph['ymax'] != None:
        dhigh = kwgraph['ymax']

    if (not kwgraph['axis']) and (not kwgraph['xaxis']):
        kwgraph['xaxis'] = [elow, ehigh]
    if (not kwgraph['axis']) and (not kwgraph['yaxis']):
        kwgraph['yaxis'] = [dlow, dhigh]

    set_scale(**kwgraph)

    #Make dataset plotable if there are zeros...
    for n in range(len(data)):
        if data[n] == 0.0:
            data[n] = dlow * 0.1
            if dlow == 0.0:
                data[n] = 10.0**-300

    #plot data
    x, y = StairStep(energy_bins, data, G)
    plt.plot(x, y, kwgraph['colorline'], label=kwgraph['datalabel'])

    set_axis(**kwgraph)
    set_legend_label(**kwgraph)

    writeGraph(**kwgraph)

    return 
