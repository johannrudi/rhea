#!/usr/bin/env python

###############################################################################
# Creates plots of the parallel distribution of work (#elements, #cnodes, etc.)
# and memory usage.
#
# Author:             Johann Rudi <johann@ices.utexas.edu>
###############################################################################

#
# Load Modules
#

import sys
import time
import argparse
import math
import numpy
import matplotlib.pyplot

#
# Set Constant Parameters
#

# set script name
THIS_SCRIPT_NAME = 'workload_plot'

# set regex patterns for reading petsc files
REGEX_INT = r'([\-\+]?\d+)'
REGEX_REAL = r'([\-\+]?[\d\.]+[e]?[\-\+]?\d*)'
REGEX_PETSC_MALLOC = (r'\[' + REGEX_INT  +
                      r'\]\s*(Current space PetscMalloc\(\)ed)\s+' +
                      REGEX_REAL + r',\s+(...)')
REGEX_PETSC_MEM = (r'\[' + REGEX_INT  + r'\]\s*(Current process memory)\s+' +
                   REGEX_REAL + r'\s+(...)')

# set figure width and height (inch)
FIGURE_WIDTH = 8.0
FIGURE_HEIGHT_PER_PROC_MIN = 0.1
FIGURE_HEIGHT_PER_PROC_MAX = 0.5

# set margin for axes (inch)
FIGURE_AXES_MARGIN = 1.2

# set figure resolution (dots per inch)
FIGURE_DPI = 10.0
OUTPUT_DPI_DEFAULT = 100.0
OUTPUT_MAX_PIXELS = 32768.0

# set width of bars
BAR_WIDTH_DEFAULT = 0.8
BAR_AXES_MARGIN = 0.4

# qualitative colorset (printer friendly, colorblind save) (by colorbrewer2.org)
COLOR_LIST = {
    'red'   : [228, 26, 28],
    'blue'  : [ 55,126,184],
    'green' : [ 77,175, 74],
    'violet': [152, 78,163],
    'orange': [255,127,  0],
    'yellow': [255,255, 51],
    'brown' : [166, 86, 40],
    'pink'  : [247,129,191],
    'gray'  : [153,153,153]
}
for k, c in COLOR_LIST.items():
    COLOR_LIST[k] = tuple(numpy.array(c, dtype='f') / 255.0)

# set colors for plots
COLOR_ELEM_DISTRIB     = COLOR_LIST['red']
COLOR_CNODE_DISTRIB    = COLOR_LIST['orange']
COLOR_MEM_TOTAL        = COLOR_LIST['blue']
COLOR_MEM_OWN          = COLOR_LIST['yellow']
COLOR_MEM_PETSC_TOTAL  = COLOR_LIST['violet']
COLOR_MEM_PETSC_MALLOC = COLOR_LIST['pink']

# set list of plot names
PLOT_NAMES = {
    'elements': 'Distribution of elements',
    'cnodes'  : 'Distribution of continuous nodes',
    'memory'  : 'Memory usage (10^6 Bytes)'
}

###############################################################################

# prints messages in stdout
def print_info(message):
    print '[' + THIS_SCRIPT_NAME + '] ' + message

# prints error messages
def print_error(message):
    print_info('Error: ' + message)

# prints warning messages
def print_warning(message):
    print_info('Warning: ' + message)

###############################################################################

# status output
print_info('Start [' + time.strftime('%Y-%m-%d %H:%M:%S') + ']')

#
# Process Arguments
#

parser = argparse.ArgumentParser()
parser.add_argument(
        'file',
        help='path to the workload file'
)
parser.add_argument(
        '-p', '--petscfile',
        help='path to the petsc file'
)
parser.add_argument(
        '-b', '--barwidth', type=float,
        help='width of a bar in a bar chart, range (0,1]'
)
parser.add_argument(
        '-o', '--output',
        help='path for output'
)
parser.add_argument(
        '-d', '--outdpi',
        help='resolution of output in DPI'
)
args = parser.parse_args()

# get input path
in_filepath = args.file

# process arguments
if args.barwidth:
    barwidth = float(args.barwidth)
else:
    barwidth = BAR_WIDTH_DEFAULT

if args.petscfile:
    in_petsc_filepath = args.petscfile
else:
    in_petsc_filepath = None

if args.output:
    out_filepath = args.output
else:
    out_filepath = in_filepath[0:in_filepath.rfind('.')]

if args.outdpi:
    out_dpi = float(args.outdpi)
else:
    out_dpi = OUTPUT_DPI_DEFAULT

#
# Read Input File(s)
#

# read main workload file
print_info('Read file ' + in_filepath)
W = numpy.loadtxt(in_filepath, comments='#')

# read petsc file
if in_petsc_filepath is not None:  # if petsc file given
    print_info('Read PETSc file ' + in_petsc_filepath)

    # read petsc mallocs
    W_petsc_malloc = numpy.fromregex(
            in_petsc_filepath, REGEX_PETSC_MALLOC,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('malloc', numpy.float64), ('t2', numpy.bool_)])

    # read petsc memory usage
    W_petsc_mem = numpy.fromregex(
            in_petsc_filepath, REGEX_PETSC_MEM,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('mem', numpy.float64), ('t2', numpy.bool_)])

    petsc_data_given = True
elif 9 < W.shape[1]:
    petsc_data_given = True
else:
    petsc_data_given = False

###DEV###
#print(W)

#
# Prepare Data
#

print_info('Prepare data')

# prepare common data
mpirank = W[:,0]
n_procs = mpirank[-1] + 1

# prepare data for distribution of elements
elem_distrib = W[:,2]
max_elements = max(elem_distrib)

# prepare data for distribution of continuous nodes
cnode_distrib = W[:,3]
max_cnodes = max(cnode_distrib)

# prepare data for memory usage
mem_usage_own = 1.0e-6 * W[:,7]
mem_usage_total = 1.0e-6 * W[:,8]
max_mem_total = max(mem_usage_total)
if petsc_data_given:
    if in_petsc_filepath is not None:  # if petsc file given
        petsc_malloc = 1.0e-6 * W_petsc_malloc['malloc']
        mem_usage_petsc = 1.0e-6 * W_petsc_mem['mem']
    else:
        petsc_malloc = 1.0e-6 * W[:,9]
        mem_usage_petsc = 1.0e-6 * W[:,10]
    max_mem_petsc = max(mem_usage_petsc)

#
# Prepare Plotting
#

# calculate height per process in figure
fig_height_per_proc = max(
        FIGURE_HEIGHT_PER_PROC_MIN,
        FIGURE_HEIGHT_PER_PROC_MAX - 0.1 * math.log(n_procs)/math.log(16))
print_info('Set space per rank in figure to %g in' % (fig_height_per_proc))

# calculate figure size (inch)
fig_width = FIGURE_WIDTH
fig_height = n_procs * fig_height_per_proc + 2.0 * FIGURE_AXES_MARGIN
print_info('Set figure size to %.0f x %.0f in' % (fig_width,  fig_height))

# calculate axes size and position (relative)
axes_x = FIGURE_AXES_MARGIN/fig_width
axes_y = FIGURE_AXES_MARGIN/fig_height
axes_width = (fig_width - 2.0 * FIGURE_AXES_MARGIN) / fig_width
axes_height = (n_procs * fig_height_per_proc) / fig_height
print_info('Set axes size to %g x %g' % (axes_width,  axes_height))

# set default bar positions
bar_pos = -mpirank - 0.5 * barwidth

# set bar positions for memory usage plots
if petsc_data_given:
    barwidth_mem = 0.5 * barwidth
    bar_pos_mem = -mpirank
    bar_pos_petsc = -mpirank - barwidth_mem
else:  # if no petsc file given
    barwidth_mem = barwidth
    bar_pos_mem = bar_pos

# initialize figures and axes
figure = {}
axes = {}
for k, t in PLOT_NAMES.items():
    figure[k] = matplotlib.pyplot.figure(
            num=None, figsize=(fig_width, fig_height), dpi=FIGURE_DPI)
    axes[k] = figure[k].add_axes(
            [axes_x, axes_y, axes_width, axes_height])

#
# Plot Distribution of Elements
#

print_info('Create "elements" plot')

# plot distribution of elements
line_max = axes['elements'].plot(
        [max_elements, max_elements],
        [bar_pos[0] + barwidth, bar_pos[-1]],
        color=COLOR_ELEM_DISTRIB, zorder=1)
bars_elem_distrib = axes['elements'].barh(
        bar_pos, elem_distrib, barwidth,
        color=COLOR_ELEM_DISTRIB, zorder=10)

#
# Plot Distribution of Continuous Nodes
#

print_info('Create "cnodes" plot')

# plot distribution of continuous nodes
line_max = axes['cnodes'].plot(
        [max_cnodes, max_cnodes],
        [bar_pos[0] + barwidth, bar_pos[-1]],
        color=COLOR_CNODE_DISTRIB, zorder=1)
bars_cnodes_distrib = axes['cnodes'].barh(
        bar_pos, cnode_distrib, barwidth,
        color=COLOR_CNODE_DISTRIB, zorder=10)

#
# Plot Memory Usage
#

print_info('Create "memory" plot')

# plot memory usage
line_max = axes['memory'].plot(
        [max_mem_total, max_mem_total],
        [bar_pos_mem[0] + barwidth_mem, bar_pos_mem[-1]],
        color=COLOR_MEM_TOTAL, zorder=1)
bars_mem_total = axes['memory'].barh(
        bar_pos_mem, mem_usage_total, barwidth_mem,
        color=COLOR_MEM_TOTAL, zorder=10)
bars_mem_own = axes['memory'].barh(
        bar_pos_mem, mem_usage_own, barwidth_mem,
        color=COLOR_MEM_OWN, zorder=11)

if petsc_data_given:
    line_max = axes['memory'].plot(
            [max_mem_petsc, max_mem_petsc],
            [bar_pos_petsc[0] + barwidth_mem, bar_pos_petsc[-1]],
            color=COLOR_MEM_PETSC_TOTAL, zorder=2)
    bars_mem_petsc = axes['memory'].barh(
            bar_pos_petsc, mem_usage_petsc, barwidth_mem,
            color=COLOR_MEM_PETSC_TOTAL, zorder=12)
    if 0 < len(petsc_malloc):
        bars_petsc_malloc = axes['memory'].barh(
                bar_pos_petsc, petsc_malloc, barwidth_mem,
                color=COLOR_MEM_PETSC_MALLOC, zorder=13)
    else:
        bars_petsc_malloc = []

#
# Format Plots
#

# set common formatting
ytick = -mpirank
ytick_labels = ['rank {:.0f}'.format(rank) for rank in mpirank]

# format all axes
for k, a in axes.items():
    # add title
    a_title = a.set_title(PLOT_NAMES[k] + '\n\n')

    # set y-ticks
    a.set_yticks(ytick)
    a.set_yticklabels(ytick_labels)
    a.set_ylim(1.0 - (n_procs + 0.5 * barwidth + BAR_AXES_MARGIN),
               0.5 * barwidth + BAR_AXES_MARGIN)

    # add x-ticks on top of axes
    a_twin_y = a.twiny()
    a_twin_y.set_frame_on(False)
    a_twin_y.set_xticks(a.get_xticks())

    # add legend to memory usage plot
    if 'memory' == k:
        legend = [bars_mem_own, bars_mem_total]
        legend_label = ['mem own', 'mem total']
        if petsc_data_given:
            legend.extend([bars_petsc_malloc, bars_mem_petsc])
            legend_label.extend(['PETSc malloc', 'PETSc mem'])

        matplotlib.pyplot.legend(legend, legend_label, loc='upper left',
                                 bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0,
                                 prop={'size':8})

#
# Output of Plots
#

# set output size of pictures
out_width = fig_width * out_dpi
out_height = fig_height * out_dpi
if OUTPUT_MAX_PIXELS < out_width or OUTPUT_MAX_PIXELS < out_height:
    # reduce resolution if picture is too large
    if out_width < out_height:  # if height is too large
        out_dpi = math.floor(OUTPUT_MAX_PIXELS / fig_height)
    else:  # if width is too large
        out_dpi = math.floor(OUTPUT_MAX_PIXELS / fig_width)
    print_warning('Picture size is too large, reduce resolution to %.0f DPI' %
                  (out_dpi))
    out_width = fig_width * out_dpi
    out_height = fig_height * out_dpi
print_info('Set file output picture size to %.0f x %.0f px, %.0f DPI' %
           (out_width, out_height, out_dpi))

# save all plots to disk
for k, f in figure.items():
    print_info('Save "' + k + '" plot')
    f.savefig(out_filepath + '.' + k + '.png', dpi=int(out_dpi))

# show plot
#matplotlib.pyplot.show()

# status output
print_info('Finished [' + time.strftime('%Y-%m-%d %H:%M:%S') + ']')

