#!/usr/bin/env python

###############################################################################
# Computes the error between the input temperature data and the output
# of the temperature after being processed by ymir program.
#
# Note: A single file with output temperature data is required.
#   Create it by merging the multiple files, e.g.,
#     $ cat earth_temp_data.*.txt > earth_temp_data.txt
#   and then run this script.
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
THIS_SCRIPT_NAME = 'check_temperature_data'

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
        'in_temp_data',
        help='path to the input temperature file'
)
parser.add_argument(
        'out_temp_data',
        help='path to the ouput temperature file'
)
args = parser.parse_args()

# get paths
in_filepath = args.in_temp_data
out_filepath = args.out_temp_data

#
# Read Files
#

# read input temperature
print_info('Read file ' + in_filepath)
in_mat = numpy.loadtxt(in_filepath, comments='#')
in_temp = in_mat[:,1]

# read ouput temperature
print_info('Read file ' + out_filepath)
out_mat = numpy.loadtxt(out_filepath, comments='#')
out_temp = out_mat[:,4]

#
# Check Temperature Data
#

# check number of entries
if len(in_temp) != len(out_temp):
    print_error('Dimensions not matching')
    sys.exit(1)

# compute error
error = abs(in_temp - out_temp)
error_min = error.min()
error_max = error.max()
error_rms = numpy.sqrt(numpy.sum(numpy.square(error)) / len(error))

# print error
print_info('Temperature data error min %1.3e, max %1.3e' %
           (error_min, error_max))
print_info('Temperature data error RMS %1.3e' % (error_rms))

