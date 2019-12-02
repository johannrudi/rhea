#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Setup
###############################################################################

# set plotting parameters
XLIM = [0, 240]
XTICKS_INCREMENT = 30

# set environment
argn = len(sys.argv) - 1
script_path = os.path.abspath(sys.argv[0])
script_dir  = os.path.dirname(script_path)
script_base = os.path.splitext(os.path.basename(script_path))[0] # w/o ext.

# set input arguments
assert 1 <= argn
data_in  = sys.argv[1]
data_dir = os.path.dirname(os.path.abspath(data_in))

# set output path
if 2 <= argn:
    plot_out = sys.argv[2]
else:
    plot_out = data_dir+"/"+script_base+".png"

###############################################################################
# Functions
###############################################################################

def print_info(message):
    '''Prints messages in stdout.'''
    print("["+script_base+"] "+message)

def rotate(x, y, angle):
    '''Rotates coordinates in the (x,y)-plane by a given angle.'''
    assert x.size == y.size
    rot = np.zeros((x.size, 2))
    rot[:,0] = np.cos(angle) * x - np.sin(angle) * y
    rot[:,1] = np.sin(angle) * x + np.cos(angle) * y
    return rot

def get_longitude(x, y):
    '''Computes longitude in [0,360) from Cartesian coordinates (x,y).'''
    longitude = np.arctan2(y, x)
    idx_negative = longitude < -1.0e-12
    longitude[idx_negative] = 2.0*np.pi + longitude[idx_negative]
    return np.maximum(0.0, longitude * 180.0/np.pi)

def get_tangential_velocity(x, y, vx, vy):
    '''Computes tangential component of the velocity (vx,vy) at (x,y).'''
    a = np.zeros((x.size, 2))
    a[:,0] = x
    a[:,1] = y
    b = np.zeros((vx.size, 2))
    b[:,0] = vx
    b[:,1] = vy
    return np.sqrt(vx**2 + vy**2) * np.sign(np.cross(a, b))

###############################################################################
# Main
###############################################################################

print_info("Data input: "+data_in)
print_info("Plot output: "+plot_out)

# load file
M = np.loadtxt(data_in, skiprows=1, delimiter=",")
obs_weight = M[:, 0]
vel_fwd    = M[:, 1: 4]
vel_adj    = M[:, 4: 7]
vel_obs    = M[:, 7:10]
vel_misfit = M[:,10:13]
point      = M[:,13:16]
point_idx  = M[:,16]

# set coordinates
x = -point[:,0]
z = point[:,2]
radius = np.sqrt(x**2 + z**2)
angle  = np.arctan2(z, x)
print_info("rotate by %g" % angle[0])
rot       = rotate(x, z, -angle[0])
longitude = get_longitude(rot[:,0], rot[:,1])

# set velocities
fwd    = get_tangential_velocity(x, z, -vel_fwd[:,0], +vel_fwd[:,2])
obs    = get_tangential_velocity(x, z, -vel_obs[:,0], +vel_obs[:,2])
misfit = get_tangential_velocity(x, z, -vel_misfit[:,0], +vel_misfit[:,2])

# mask non existing observational data
idx_nan = np.logical_and(-1.0e-16 < obs, obs < +1.0e-16)
obs[idx_nan] = np.nan
misfit[idx_nan] = np.nan

# set vector field variables
idx_vec = np.arange(0, x.size, 50)
qx = x[idx_vec]
qy = z[idx_vec]
qu = -vel_fwd[idx_vec,0]
qv = +vel_fwd[idx_vec,2]
qc = np.hypot(qu, qv)

# create figure
fig, ax = plt.subplots(3, 1, figsize=(6, 9))

# plot spherical domain and forward velocity
ax[0].plot(x, z, color='0.5')
ax[0].quiver(qx, qy, qu, qv, qc)
ax[0].axis("equal")
ax[0].set_ylabel("model velocity (spherical)")

# plot forward and observational velocities
ax[1].plot([XLIM[0], XLIM[1]], [0, 0], color='0.5', linestyle=':', linewidth=0.5)
ax[1].plot(longitude, fwd, 'b-')
ax[1].plot(longitude, obs, 'g-')
ax[1].set_xlim(XLIM[1], XLIM[0]) # decreasing
ax[1].set_xticks(np.arange(XLIM[0], XLIM[1]+1, XTICKS_INCREMENT))
ax[1].set_ylabel("model vel. vs. data [mm/yr]")

# plot data misfit
ax[2].plot([XLIM[0], XLIM[1]], [0, 0], color='0.5', linestyle=':', linewidth=0.5)
ax[2].plot(longitude, misfit, 'r-')
ax[2].set_xlim(XLIM[1], XLIM[0]) # decreasing
ax[2].set_xticks(np.arange(XLIM[0], XLIM[1]+1, XTICKS_INCREMENT))
ax[2].set_ylabel("weighted data misfit [mm/yr]")

# set figure annotations
ax[0].set_title(data_in)
ax[2].set_xlabel("longitude [degree]")

# save & show plots
fig.savefig(plot_out)
#plt.show()
