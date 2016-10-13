#!/usr/bin/env python

###############################################################################
# Computes the error between the input weak zone data and the output
# of weak zone points after being processed by ymir program.
#
# Author:             Johann Rudi <johann@ices.utexas.edu>
###############################################################################

#
# Load Modules
#

import argparse
import sys
import numpy

#
# Set Constant Parameters
#

# set script name
THIS_SCRIPT_NAME = 'check_weakzone_data'

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

#
# Process Arguments
#

parser = argparse.ArgumentParser()
parser.add_argument(
        'in_weak_data',
        help='path to the input weak zone file'
)
parser.add_argument(
        'out_weak_data',
        help='path to the ouput weak zone file'
)
args = parser.parse_args()

# get paths
in_filepath = args.in_weak_data
out_filepath = args.out_weak_data

#
# Read Files
#

# read input temperature
print_info('Read file ' + in_filepath)
in_weak = numpy.loadtxt(in_filepath, comments='#')

# read ouput temperature
print_info('Read file ' + out_filepath)
out_weak = numpy.loadtxt(out_filepath, comments='#')

#
# Check Temperature Data
#

# check number of entries
if not numpy.all(numpy.equal(in_weak.shape, out_weak.shape)):
    print_error('Dimensions not matching')
    sys.exit(1)

# compute error
error = abs(in_weak - out_weak)
error_min = error.min()
error_max = error.max()
error_rms = numpy.sqrt(numpy.sum(numpy.square(error)) / len(error))

# print error
print_info('Weak zone data error min %1.3e, max %1.3e' % (error_min, error_max))
print_info('Weak zone data error RMS %1.3e' % (error_rms))

