#!/usr/bin/python

import os, pathlib, re, sys
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Setup
###############################################################################

# set environment
script_path = os.path.abspath(sys.argv[0])
script_dir  = os.path.dirname(script_path)
script_base = os.path.splitext(os.path.basename(script_path))[0] # w/o ext.
argn = len(sys.argv) - 1

# set input arguments
assert 3 <= argn
param_dim_path  = pathlib.Path(sys.argv[1])
rhea_out_path   = pathlib.Path(sys.argv[2])
iter_files_base = pathlib.Path(sys.argv[3])

# set output path
plot_base        = rhea_out_path.parent / script_base
plot_summary     = plot_base.stem+"_summary.png"
plot_parameters  = plot_base.stem+"_parameters.png"

# define strings of inversion solver summary
ITER_SUMMARY_HEADER='[rhea] Newton solve summary: Inversion'
ITER_SUMMARY_FOOTER='[rhea] ========================================'
ITER_SUMMARY_DELIMITERS=[',', ';', '(', ')', '[', ']']

# define column indices of inversion solver summary
ITER_SUMMARY_COLUMNS = {
    'itn':               0,
    'obj':               1,
    'obj_data1':         2,
    'obj_data2':         3,
    'obj_prior':         4,
    'grad':              5,
    'grad_data1':        6,
    'grad_prior':        7,
    'step_converged':    8,
    'step_n_reductions': 9,
    'step_length':       10
}

# set plotting parameters
PLOT_TITLES_GLO = [
    "Scaling in upper mantle",
    "Activation energy in upper mantle",
    "Scaling in lower mantle",
   #"Activation energy in lower mantle",
    "Stress exponent",
    "Yield stress"
]
PLOT_TITLE_WEAK = "Weak zone factor"

# set physical parameters for nondimensialization
phys_glo_dim = np.array([
    1.0,        # scaling in upper mantle
    8.314*1400, # activation energy in upper mantle
    1.0,        # scaling in lower mantle
   #8.314*1400, # activation energy in lower mantle
    1.0,        # stress exponent
    1.0e20*1.0e-6/6371.0e3**2  # yield stress
])

###############################################################################
# Functions
###############################################################################

def print_info(message: str):
    '''Prints messages in stdout.'''
    print("["+script_base+"] "+message)

def parse_solver_filename(base: pathlib.Path, itn: int, name: str, suffix='.txt'):
    return '{}'.format(base) + ('_itn%02d_' % itn) + name + suffix

def extract_iterations_summary(rhea_out_path: pathlib.Path, header=ITER_SUMMARY_HEADER, footer=ITER_SUMMARY_FOOTER,
                               delimiters=ITER_SUMMARY_DELIMITERS):
    # extract lines
    lines = None
    with open(rhea_out_path, 'r') as f:
        for line in f:
            if lines is None and header in line: # if header found
                lines = []
            elif lines is not None:
                if footer in line: # if footer found
                    break
                elif '----' in line: # if skip line
                    continue
                else: # otherwise add line to output
                    lines.append(line.lstrip('[rhea] ').rstrip())
    # convert lines to numbers
    summary = None
    if lines is None:
        return summary
    for line in lines:
        for delim in delimiters:
            line = line.replace(delim, '')
        array = np.fromstring(line, sep=' ')
        if summary is None:
            summary = array[np.newaxis,:]
        elif summary.shape[1] == len(array):
            summary = np.append(summary, array[np.newaxis,:], axis=0)
        else:
            print_info('extract_iterations_summary, exclude line: {}'.format(line))
    return summary

def load_parameters(iter_files_base: pathlib.Path, n_itn: int):
    parameters = None
    for itn in range(n_itn):
        filename = parse_solver_filename(iter_files_base, itn, 'parameters')
        p = np.loadtxt(filename)
        if parameters is None:
            parameters = p[np.newaxis,:]
        else:
            parameters = np.append(parameters, p[np.newaxis,:], axis=0)
    return parameters

def load_steps(iter_files_base: pathlib.Path, n_itn: int):
    steps = None
    for itn in range(1, n_itn):
        filename = parse_solver_filename(iter_files_base, itn, 'step')
        s = np.loadtxt(filename)
        if steps is None:
            steps = np.empty_like(s[np.newaxis,:])
            steps[:] = np.nan
        steps = np.append(steps, s[np.newaxis,:], axis=0)
    return steps

###############################################################################
# Main
###############################################################################

print_info("INPUT")
print_info("- Parameter dim:        {}".format(param_dim_path))
print_info("- Rhea output:          {}".format(rhea_out_path))
print_info("- Iteration files base: {}".format(iter_files_base))
print_info("OUTPUT")
print_info("- Plot iterations:      {}".format(plot_summary))
print_info("- Plot parameters:      {}".format(plot_parameters))

# load files
param_dim    = np.loadtxt(param_dim_path)
iter_summary = extract_iterations_summary(rhea_out_path)
assert iter_summary is not None
n_itn        = iter_summary.shape[0]
parameters   = load_parameters(iter_files_base, n_itn)
steps        = load_steps(iter_files_base, n_itn)

# set sizes
n_parameters      = len(param_dim)
n_parameters_glo  = len(phys_glo_dim)
n_parameters_weak = n_parameters - n_parameters_glo
assert parameters.shape[1] == n_parameters
assert parameters.shape == steps.shape

# extract columns from summary
iter_itn               = iter_summary[:,ITER_SUMMARY_COLUMNS['itn']].astype(int)
iter_obj               = iter_summary[:,ITER_SUMMARY_COLUMNS['obj']]
iter_obj_data1         = iter_summary[:,ITER_SUMMARY_COLUMNS['obj_data1']]
iter_obj_prior         = iter_summary[:,ITER_SUMMARY_COLUMNS['obj_prior']]
iter_grad              = iter_summary[:,ITER_SUMMARY_COLUMNS['grad']]
if ITER_SUMMARY_COLUMNS['step_length'] <= iter_summary.shape[1]:
    iter_step_length_eff   = iter_summary[:,ITER_SUMMARY_COLUMNS['step_length']]
    iter_step_n_reductions = iter_summary[:,ITER_SUMMARY_COLUMNS['step_n_reductions']].astype(int)
    # compute step metrics
    iter_step_backtrack    = np.ones_like(iter_step_length_eff) * (0.5 ** iter_step_n_reductions)
    iter_step_length_init  = iter_step_length_eff / iter_step_backtrack
    iter_step_rel          = np.mean(np.abs(steps), axis=1)
else:
    iter_step_length_eff  = None
    iter_step_length_init = None
    iter_step_rel         = None

########################################
# Plot Objective, Gradient, Step Length
########################################

# create figure
fig, ax = plt.subplots(4, 1, figsize=(6, 8))

# plot objective
ax[0].plot(iter_itn, iter_obj, label='objective')
ax[0].plot(iter_itn, iter_obj_data1, label='data')
ax[0].plot(iter_itn, iter_obj_prior, label='prior')
ax[0].set_ylabel('Objective')
ax[0].legend()

# plot gradient
ax[1].plot(iter_itn, iter_grad)
ax[1].set_yscale('log')
ax[1].set_ylabel('Gradient')

# plot step length
if iter_step_length_init is not None:
    ax[2].plot(iter_itn[1:], iter_step_length_init[1:])
if iter_step_length_eff is not None:
    ax[2].plot(iter_itn[1:], iter_step_length_eff[1:])
ax[2].set_yscale('log')
ax[2].set_ylabel('Step length')

# plot relative step size
if iter_step_rel is not None:
    ax[3].plot(iter_itn[1:], iter_step_rel[1:])
ax[3].set_yscale('log')
ax[3].set_ylabel('Relative step size')

# format axes
for a in ax:
    a.set_xlim(iter_itn[0], iter_itn[-1])
    a.grid(True)

# set spacing between subplots
fig.set_tight_layout({'pad': 0.5})

# save plots
fig.savefig(plot_summary, dpi=360)

########################################
# Plot Parameters, Steps
########################################

# create figure
fig, ax = plt.subplots(n_parameters, 2, figsize=(8, 16))

for i in range(n_parameters):
    ax[i,0].plot(iter_itn, parameters[:,i])
    ax[i,1].plot(iter_itn[1:], np.abs(steps[1:,i]))
    ax[i,1].set_yscale('log')
    if i < n_parameters_glo:
        ax[i,0].set_ylabel(PLOT_TITLES_GLO[i])
    else:
        ax[i,0].set_ylabel(PLOT_TITLE_WEAK + ' %d' % (i - n_parameters_glo))

# format axes
for i in range(n_parameters):
    for j in range(2):
        ax[i,j].set_xlim(iter_itn[0], iter_itn[-1])
        ax[i,j].grid(True)

# set spacing between subplots
fig.set_tight_layout({'pad': 0.5})

# save plots
fig.savefig(plot_parameters, dpi=360)

# show all plots
plt.show()
