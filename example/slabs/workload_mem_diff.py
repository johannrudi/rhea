#!/usr/bin/env python

###############################################################################
# Takes the difference of the memory usage columns of two files.
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

#
# Set Constant Parameters
#

# set script name
THIS_SCRIPT_NAME = 'workload_diff_mem'

# set regex patterns for reading petsc files
REGEX_INT = r'([\-\+]?\d+)'
REGEX_REAL = r'([\-\+]?[\d\.]+[e]?[\-\+]?\d*)'
REGEX_PETSC_MALLOC = (r'\[' + REGEX_INT  +
                      r'\]\s*(Current space PetscMalloc\(\)ed)\s+' +
                      REGEX_REAL + r',\s+(...)')
REGEX_PETSC_MEM = (r'\[' + REGEX_INT  + r'\]\s*(Current process memory)\s+' +
                   REGEX_REAL + r'\s+(...)')

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
        'prevfile',
        help='path to the previous file that will be subtracted'
)
parser.add_argument(
        'file',
        help='path to the main file'
)
parser.add_argument(
        'outfile',
        help='path for output'
)
parser.add_argument(
        '-p', '--petsc',
        help="include memory diff from petsc files",
        action="store_true")
args = parser.parse_args()

# process arguments
in_filepath = args.file
in_prev_filepath = args.prevfile
out_filepath = args.outfile

#
# Read Input File(s)
#

# read main workload file
print_info('Read main file ' + in_filepath)
W = numpy.loadtxt(in_filepath, comments='#')

# read previous workload file
print_info('Read prev file ' + in_prev_filepath)
W_prev = numpy.loadtxt(in_prev_filepath, comments='#')

# read petsc files
if args.petsc:
    in_petsc_filepath = in_filepath.replace('.txt', '.petsc.txt')
    in_petsc_prev_filepath = in_prev_filepath.replace('.txt', '.petsc.txt')

    print_info('Read PETSc file ' + in_petsc_filepath)
    W_petsc_malloc = numpy.fromregex(
            in_petsc_filepath, REGEX_PETSC_MALLOC,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('malloc', numpy.float64), ('t2', numpy.bool_)])
    W_petsc_mem = numpy.fromregex(
            in_petsc_filepath, REGEX_PETSC_MEM,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('mem', numpy.float64), ('t2', numpy.bool_)])

    print_info('Read prev PETSc file ' + in_petsc_prev_filepath)
    W_petsc_prev_malloc = numpy.fromregex(
            in_petsc_prev_filepath, REGEX_PETSC_MALLOC,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('malloc', numpy.float64), ('t2', numpy.bool_)])
    W_petsc_prev_mem = numpy.fromregex(
            in_petsc_prev_filepath, REGEX_PETSC_MEM,
            [('mpirank', numpy.uint32), ('t1', numpy.bool_),
             ('mem', numpy.float64), ('t2', numpy.bool_)])

#
# Create Diff
#

# initialize empty matrix
w_diff_shape = list(W.shape)
if args.petsc:
    w_diff_shape[1] += 2
W_diff = numpy.zeros(shape=w_diff_shape)

# set diff from main files
W_diff[:,0:7] = W[:,0:7]
W_diff[:,7] = 0
W_diff[:,8] = W[:,8] - W_prev[:,8]

# set diff from petsc files
if args.petsc:
    W_diff[:,9] = W_petsc_malloc['malloc'] - W_petsc_prev_malloc['malloc']
    W_diff[:,10] = W_petsc_mem['mem'] - W_petsc_prev_mem['mem']

#
# Write Output
#

head = ' mpirank, maxlevel, num elements, num owned cnodes, ' + \
       'num shared cnodes, num velocity unknowns, ' + \
       'num pressure unknowns, [unused], memory usage diff (Bytes)'

format = '%4lu  %2lu  %6lu  %9lu  %9lu  %9lu  %9lu  %10lu  %10lu'

if args.petsc:
    head += ', petsc malloc diff, petsc mem diff'
    format += '  %10lu  %10lu'

print_info('Write diff to file "' + out_filepath + '"')
numpy.savetxt(out_filepath, W_diff, fmt=format, header=head, comments='#')

# status output
print_info('Finished [' + time.strftime('%Y-%m-%d %H:%M:%S') + ']')

