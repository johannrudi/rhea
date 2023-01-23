#!/usr/bin/python

import os, pathlib, re, sys
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
data_in       = pathlib.Path(sys.argv[1])
data_dir      = data_in.parent
rhea_out_path = pathlib.Path(sys.argv[2])

# set output path
if 3 <= argn:
    plot_out = pathlib.Path(sys.argv[3])
else:
    plot_out = data_dir / script_base

# define strings for rhea output
VEL_NOISE_HEADER='[rhea] rhea_inversion: obs_velocity_weight'
VEL_NOISE_FOOTER='[rhea] ========================================'
VEL_NOISE_PREFIX='[rhea] plate_idx'
VEL_NOISE_DELIMITERS=[',', ';', ':', '(', ')', '[', ']']

# set plot file suffix
PLOT_SUFFIX = ['png', 'eps']

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

def extract_velocity_noise(rhea_out_path: pathlib.Path, header=VEL_NOISE_HEADER, footer=VEL_NOISE_FOOTER,
                           linePrefix=VEL_NOISE_PREFIX, delimiters=VEL_NOISE_DELIMITERS):
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
                    lines.append(line.lstrip(linePrefix).rstrip())
    # convert lines to numbers
    mat = None
    if lines is None:
        return mat
    for line in lines:
        for delim in delimiters:
            line = line.replace(delim, '')
        row = np.fromstring(line, sep=' ')
        if mat is None:
            mat = row[np.newaxis,:]
        else:
            mat = np.append(mat, row[np.newaxis,:], axis=0)
    return mat

def interpolate_noise(vel_noise, obs):
    noise     = np.ones_like(obs) * np.nan
    noise_len = vel_noise.shape[0]
    n         = 0
    v_prev    = np.nan
    for i, v in enumerate(obs):
        if not np.isnan(v):  # if inside a plate
            assert n < noise_len, (n, noise_len)
            noise[i] = vel_noise[n]
        elif not np.isnan(v_prev) and np.isnan(v):  # if end of plate
            n += 1
        v_prev = v
    return noise

###############################################################################
# Main
###############################################################################

print_info("INPUT")
print_info("- Velocity data: {}".format(data_in))
print_info("- Rhea output:   {}".format(rhea_out_path))
print_info("OUTPUT")
print_info("- Plot velocities: {}".format(plot_out))

# load files
M = np.loadtxt(data_in, skiprows=1, delimiter=",")
obs_weight = M[:, 0]
vel_fwd    = M[:, 1: 4]
vel_adj    = M[:, 4: 7]
vel_obs    = M[:, 7:10]
vel_misfit = M[:,10:13]
point      = M[:,13:16]
point_idx  = M[:,16]
vel_noise  = extract_velocity_noise(rhea_out_path)[1:,2]  # remove first row because first plate has not obs data

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
obs[idx_nan]    = np.nan
misfit[idx_nan] = np.nan

# interpolate noise to match observational data
noise = interpolate_noise(vel_noise, obs)
noise_95ci_l = obs - 2*noise
noise_95ci_h = obs + 2*noise

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
ax[0].set_ylabel("Model velocity (spherical)")
ax[0].set_title(data_in)

# plot forward and observational velocities
ax[1].plot([XLIM[0], XLIM[1]], [0, 0], color='0.5', linestyle=':', linewidth=0.5)
ax[1].fill_between(longitude, noise_95ci_l, noise_95ci_h, label='noise, 95% confidence', color='limegreen', alpha=0.3)
ax[1].plot(longitude, fwd, label='velocity, model', color='cornflowerblue', linestyle='solid', linewidth=2)
ax[1].plot(longitude, obs, label='velocity, data', color='darkgreen', linestyle='dashed', linewidth=2, alpha=0.9)
ax[1].set_xlim(XLIM[1], XLIM[0]) # decreasing
ax[1].set_xticks(np.arange(XLIM[0], XLIM[1]+1, XTICKS_INCREMENT))
ax[1].set_xlabel("Longitude [degree]")
ax[1].set_ylabel("Model vs. data [mm/yr]")
ax[1].legend(loc='upper right', fontsize='small')

# plot data misfit
ax[2].plot([XLIM[0], XLIM[1]], [0, 0], color='0.5', linestyle=':', linewidth=0.5)
ax[2].plot(longitude, misfit, 'r-')
ax[2].set_xlim(XLIM[1], XLIM[0]) # decreasing
ax[2].set_xticks(np.arange(XLIM[0], XLIM[1]+1, XTICKS_INCREMENT))
ax[2].set_xlabel("Longitude [degree]")
ax[2].set_ylabel("Weighted data misfit [mm/yr]")

# set spacing between subplots
fig.tight_layout()

# save & show plots
for suffix in PLOT_SUFFIX:
    fig.savefig('{}.{}'.format(plot_out, suffix), dpi=360)

# plot only forward and observatioanl velocities
fig, ax = plt.subplots(1, 1, figsize=(6, 2.5))
ax.plot([XLIM[0], XLIM[1]], [0, 0], color='0.5', linestyle=':', linewidth=0.5)
ax.fill_between(longitude, noise_95ci_l, noise_95ci_h, label='noise, 95% confidence', color='limegreen', alpha=0.3)
ax.plot(longitude, fwd, label='velocity, model', color='cornflowerblue', linestyle='solid', linewidth=2)
ax.plot(longitude, obs, label='velocity, data', color='darkgreen', linestyle='dashed', linewidth=2, alpha=0.9)
ax.set_xlim(XLIM[1], XLIM[0]) # decreasing
ax.set_xticks(np.arange(XLIM[0], XLIM[1]+1, XTICKS_INCREMENT))
ax.set_xlabel("Longitude [degree]")
ax.set_ylabel("Model vs. data [mm/yr]")
ax.legend(loc='upper right', fontsize='small')
fig.tight_layout()
for suffix in PLOT_SUFFIX:
    fig.savefig('{}_model_vs_data.{}'.format(plot_out, suffix), dpi=360)

#plt.show()
