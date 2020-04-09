#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Setup
###############################################################################

# set environment
argn = len(sys.argv) - 1
script_path = os.path.abspath(sys.argv[0])
script_dir  = os.path.dirname(script_path)
script_base = os.path.splitext(os.path.basename(script_path))[0] # w/o ext.

# set input arguments
assert 5 <= argn
param_dim_txt = sys.argv[1]
prior_max_txt = sys.argv[2]
prior_cov_txt = sys.argv[3]
post_max_txt  = sys.argv[4]
post_cov_txt  = sys.argv[5]

# set output path
if 6 <= argn:
    plot_base = sys.argv[6]
else:
    data_dir = os.path.dirname(os.path.abspath(post_max_txt))
    plot_base = data_dir+"/"+script_base

# set plotting parameters
PLOT_N_NODES = 256+1
PLOT_TITLES_GLO = [
    "Scaling in upper mantle",
    "Activation energy in upper mantle",
    "Scaling in lower mantle",
    "Activation energy in lower mantle",
    "Stress exponent",
    "Yield stress"
]
PLOT_TITLE_WEAK = "Weak zone factor"

# set physical parameters for nondimensialization
phys_glo_dim = np.array([
    1.0,        # scaling in upper mantle
    8.314*1400, # activation energy in upper mantle
    1.0,        # scaling in lower mantle
    8.314*1400, # activation energy in lower mantle
    1.0,        # stress exponent
    1.0e20*1.0e-6/6371.0e3**2  # yield stress
])

# set indices of densities that are lognormal
idx_lognormal_glo = [0, 2]

###############################################################################
# Functions
###############################################################################

def print_info(message):
    '''Prints messages in stdout.'''
    print("["+script_base+"] "+message)

def adjust_mean_of_lognormal(mean, variance, param_dim):
    # shift by mode:
    #   mode = (mean*param_dim - variance*param_dim^2)/param_dim
    return mean + variance*param_dim

def transfer_param_to_physical(m, param_dim, phys_dim, lognormal=False):
    if lognormal:
        return phys_dim * np.exp(param_dim * m)
    else:
        return phys_dim * param_dim * m

def transfer_param_to_physical_multi(m, param_dim, phys_glo_dim, idx_lognormal):
    # apply dimensional scaling of inversion parameters
    p = param_dim[:,np.newaxis] * m
    # apply transfer function
    p[idx_lognormal,:] = np.exp(p[idx_lognormal,:])
    # apply physical dimensional scaling
    n_glo = len(phys_glo_dim)
    p[0:n_glo,:] = phys_glo_dim[:,np.newaxis] * p[0:n_glo,:]
    # return parameters in physical dimensions
    return p

def transfer_physical_to_param_multi(p, param_dim, phys_glo_dim, idx_lognormal):
    m = np.copy(p)
    # remove physical dimensional scaling
    n_glo = len(phys_glo_dim)
    m[0:n_glo,:] = m[0:n_glo,:]/phys_glo_dim[:,np.newaxis]
    # invert transfer function
    m[idx_lognormal,:] = np.log(m[idx_lognormal,:])
    # remove dimensional scaling of inversion parameters
    m = m/param_dim[:,np.newaxis]
    # return inversion parameters
    return m

def transfer_physical_to_param_lin_deriv(p, param_dim, phys_glo_dim):
    return np.ones_like(p)/phys_glo_dim/param_dim

def transfer_physical_to_param_log_deriv(p, param_dim, phys_glo_dim=None):
    return 1.0/p/param_dim

def transfer_physical_to_param_deriv_multi(p, param_dim, phys_glo_dim, idx_lognormal):
    dm = np.ones_like(p)
    # remove physical dimensional scaling
    n_glo = len(phys_glo_dim)
    dm[0:n_glo,:] = dm[0:n_glo,:]/phys_glo_dim[:,np.newaxis]
    # invert nonlinear transfer function
    dm[idx_lognormal,:] = 1.0/p[idx_lognormal,:]
    # remove dimensional scaling of inversion parameters
    dm = dm/param_dim[:,np.newaxis]
    # return derivative of the inverse transfer function
    return dm

def gaussian_1d(m, mean, variance):
    '''Computes the univariate Gaussian distribution.'''
    return 1.0/np.sqrt(2.0*np.pi) / np.sqrt(variance[:,np.newaxis]) * \
           np.exp(-0.5*(m - mean[:,np.newaxis])**2/variance[:,np.newaxis])

def pdf_1d(m, mean, variance, p, param_dim, phys_glo_dim, idx_lognormal):
    dm = np.abs(transfer_physical_to_param_deriv_multi(
            p, param_dim, phys_glo_dim, idx_lognormal))
    return dm*gaussian_1d(m, mean, variance)

def gaussian_2d(X, Y, Mean, Cov, dX=None, dY=None):
    '''Computes the multivariate Gaussian distribution.

    From M0 and M1, we construct a single array by packing the meshed arrays of
    variables x_1, x_2, x_3, ..., x_k into its _last_ dimension.'''
    # pack M0 and M1 into a single 3-dimensional array
    Pos = np.empty(X.shape + (2,))
    Pos[:,:,0] = X
    Pos[:,:,1] = Y
    # compute values pertaining to the covariance matrix
    Cov_inv = np.linalg.inv(Cov)
    cov_det = np.abs(np.linalg.det(Cov))
    # calculate: (Pos-Mean)T.Cov-1.(Pos-Mean)
    # Note: `einsum` calculates in a vectorized way across all the input variables
    #power = np.einsum('...k,kl,...l->...', Pos-Mean, Cov_inv, Pos-Mean)
    power = np.einsum('ijk,kl->ijl', Pos-Mean, Cov_inv)
    power = np.einsum('ijk,ijk->ij', power, Pos-Mean)
    return 1.0/(2.0*np.pi*np.sqrt(cov_det)) * np.exp(-0.5*power)

def pdf_2d(indices, m, mean, cov,
           p, param_dim, phys_glo_dim, idx_lognormal):
    i = indices[0]
    j = indices[1]
    #TODO applying derivatives does not work with off-diagonal entries in covariance
    # set derivatives of transfer function for index `i`
#   dim = param_dim[i]
#   if i < len(phys_glo_dim):
#       phys_dim = phys_glo_dim[i]
#   else:
#       phys_dim = 1.0
#   if i in idx_lognormal:
#       dmi = np.abs(transfer_physical_to_param_log_deriv(p[i,:], dim, phys_dim))
#   else:
#       dmi = np.abs(transfer_physical_to_param_lin_deriv(p[i,:], dim, phys_dim))
    # set derivatives of transfer function for index `j`
#   dim = param_dim[j]
#   if j < len(phys_glo_dim):
#       phys_dim = phys_glo_dim[j]
#   else:
#       phys_dim = 1.0
#   if j in idx_lognormal:
#       dmj = np.abs(transfer_physical_to_param_log_deriv(p[j,:], dim, phys_dim))
#       print(j,"is log")
#   else:
#       dmj = np.abs(transfer_physical_to_param_lin_deriv(p[j,:], dim, phys_dim))
    # compute pdf
#   dMi, dMj = np.meshgrid(dmi, dmj)
    Mi, Mj = np.meshgrid(m[i,:], m[j,:])
    Mean = np.array([mean[i], mean[j]])
    Cov = np.array([ [cov[i,i], cov[i,j]],
                     [cov[j,i], cov[j,j]] ])
#   return dMi*dMj*gaussian_2d(Mi, Mj, Mean, Cov)
    return gaussian_2d(Mi, Mj, Mean, Cov)

def compute_correlation(cov):
    inv_stddev = np.diag(1.0/np.sqrt(np.diag(cov)))
    return np.dot(inv_stddev, np.dot(cov, inv_stddev))

###############################################################################
# Main
###############################################################################

print_info("INPUT")
print_info("- Parameter dim: "+param_dim_txt)
print_info("- Prior max: "+prior_max_txt)
print_info("- Prior covariance: "+prior_cov_txt)
print_info("- Posterior max: "+post_max_txt)
print_info("- Posterior covariance: "+post_cov_txt)
print_info("OUTPUT")
print_info("- Plot: "+plot_base)

# load files
param_dim = np.loadtxt(param_dim_txt)
prior_max = np.loadtxt(prior_max_txt)
prior_cov = np.loadtxt(prior_cov_txt)
post_max  = np.loadtxt(post_max_txt)
post_cov  = np.loadtxt(post_cov_txt)

# set sizes and indices
n_glo    = len(phys_glo_dim)
n_weak   = len(param_dim) - n_glo
idx_glo  = range(n_glo)
idx_weak = range(n_glo, n_glo+n_weak)
idx_lognormal = idx_lognormal_glo + idx_weak

# set mean and covariance of prior
prior_cov_diag = np.diag(prior_cov)
prior_mean = np.copy(prior_max)
for i in idx_lognormal:
    prior_mean[i] = adjust_mean_of_lognormal(prior_mean[i], prior_cov_diag[i], param_dim[i])
# set mean and covariance of posterior
post_cov_diag = np.diag(post_cov)
post_mean = np.copy(post_max)
for i in idx_lognormal:
    post_mean[i] = adjust_mean_of_lognormal(post_mean[i], post_cov_diag[i], param_dim[i])

# set physical dimensions
phys_dim = np.concatenate([phys_glo_dim, np.ones(n_weak)])
# convert variances to physical dimensions
prior_cov_diag_dim = np.zeros_like(prior_cov_diag)
post_cov_diag_dim = np.zeros_like(post_cov_diag)
for i in range(n_glo+n_weak):
    prior_cov_diag_dim[i] = transfer_param_to_physical(
            np.sqrt(prior_cov_diag[i]), param_dim[i], phys_dim[i],
            i in idx_lognormal)
    post_cov_diag_dim[i] = transfer_param_to_physical(
            np.sqrt(post_cov_diag[i]), param_dim[i], phys_dim[i],
            i in idx_lognormal)

np.set_printoptions(precision=3, linewidth=150)
print_info("Compare prior max vs. posterior max")
print(np.array([prior_max, post_max]))
print_info("Compare prior mean vs. posterior mean")
print(np.array([prior_mean, post_mean]))
print_info("Compare prior stddev vs. posterior stddev")
print(np.array([np.sqrt(prior_cov_diag), np.sqrt(post_cov_diag)]))

# create range of parameters for plotting
m_glo  = np.linspace(0.2, 1.8, num=PLOT_N_NODES)
m_weak = np.linspace(-2.0, 0.0, num=PLOT_N_NODES)
m = np.concatenate([np.tile(m_glo, (n_glo, 1)),
                    np.tile(m_weak, (n_weak, 1))], axis=0)

# transfer inversion parameters to physical values
p = transfer_param_to_physical_multi(m, param_dim, phys_glo_dim, idx_lognormal)
print_info(
    "Error of transfer fn. vs. its inverse: max abs err %g" % \
    np.max(np.abs(m - transfer_physical_to_param_multi(p, param_dim, phys_glo_dim, idx_lognormal)))
)

# compute max coordinates
m_prior_max = prior_max[:,np.newaxis]
m_post_max  = post_max[:,np.newaxis]
p_prior_max = transfer_param_to_physical_multi(m_prior_max, param_dim, phys_glo_dim, idx_lognormal)
p_post_max  = transfer_param_to_physical_multi(m_post_max, param_dim, phys_glo_dim, idx_lognormal)
pdf_prior_1d_max = \
    np.abs(transfer_physical_to_param_deriv_multi(
            p_prior_max, param_dim, phys_glo_dim, idx_lognormal)) * \
    gaussian_1d(m_prior_max, prior_mean, prior_cov_diag)
pdf_post_1d_max  = \
    np.abs(transfer_physical_to_param_deriv_multi(
            p_post_max, param_dim, phys_glo_dim, idx_lognormal)) * \
    gaussian_1d(m_post_max, post_mean, post_cov_diag)

# compute 1-dim probability densities
pdf_prior_1d = pdf_1d(m, prior_mean, prior_cov_diag,
                      p, param_dim, phys_glo_dim, idx_lognormal)
pdf_post_1d  = pdf_1d(m, post_mean, post_cov_diag,
                      p, param_dim, phys_glo_dim, idx_lognormal)

# integrate densities
int_pdf_prior = np.empty(n_glo+n_weak)
int_pdf_post  = np.empty(n_glo+n_weak)
for i in range(n_glo+n_weak):
    int_pdf_prior[i] = np.trapz(pdf_prior_1d[i,:], p[i,:])
    int_pdf_post[i]  = np.trapz(pdf_post_1d[i,:], p[i,:])
print_info("Compare integrals of prior pdf's vs. posterior pdf's")
print(np.array([int_pdf_prior, int_pdf_post]))

########################################
# 1D Plots
########################################

# create title strings
title = np.empty(n_glo+n_weak, dtype=object)
title[0:n_glo] = PLOT_TITLES_GLO
for i in range(n_weak):
    title[n_glo+i] = PLOT_TITLE_WEAK+(" %i" % i)

def plot_pdfs_1d(ax, p,
                 pdf_prior, p_prior_max, pdf_prior_max,
                 pdf_post, p_post_max, pdf_post_max,
                 rescale_pdf_max=False):
    if rescale_pdf_max:
        pdf_prior = pdf_prior/pdf_prior_max
        pdf_post  = pdf_post/pdf_post_max
        pdf_prior_max = 1.0
        pdf_post_max  = 1.0
    # plot prior max point
    x = np.repeat(p_prior_max, 2)
    y = [0.0, pdf_prior_max]
    ax.plot(x, y, color='0.3', linewidth=1, linestyle='--')
    # plot posterior max point
    x = np.repeat(p_post_max, 2)
    y = [0.0, pdf_post_max]
    ax.plot(x, y, color='blue', linewidth=1, linestyle='--')
    # plot prior and posterior densities
    ax.plot(p, pdf_prior, label='prior', color='0.3', linewidth=2)
    ax.plot(p, pdf_post, label='posterior', color='blue', linewidth=2)

# create figure
plot_n_rows = np.max([n_glo, n_weak])
fig, ax = plt.subplots(plot_n_rows, 2, figsize=(10, 10))

# create plots
for i in range(n_glo):
    this_title = PLOT_TITLES_GLO[i]
    ax[i,0].set_title(this_title)
    plot_pdfs_1d(ax[i,0], p[i,:],
                 pdf_prior_1d[i,:], p_prior_max[i], pdf_prior_1d_max[i],
                 pdf_post_1d[i,:], p_post_max[i], pdf_post_1d_max[i],
                 rescale_pdf_max=True)
    handles, labels = ax[i,0].get_legend_handles_labels()
    labels = (
        "m=%.2e s=%.2e" % (p_prior_max[i], prior_cov_diag_dim[i]),
        "m=%.2e s=%.2e" % (p_post_max[i], post_cov_diag_dim[i])
    )
    ax[i,0].legend(handles, labels, loc='upper left', fontsize='xx-small')
    ax[i,0].ticklabel_format(axis='x', style='sci', scilimits=(1,3))
    ax[i,0].set_xlim(p[i,0], p[i,-1])
    ax[i,0].set_ylim(0.0, 1.1)
for ii in range(n_weak):
    i = n_glo + ii
    this_title = PLOT_TITLE_WEAK+(" %i" % ii)
    ax[ii,1].set_title(this_title)
    plot_pdfs_1d(ax[ii,1], p[i,:],
                 pdf_prior_1d[i,:], p_prior_max[i], pdf_prior_1d_max[i],
                 pdf_post_1d[i,:], p_post_max[i], pdf_post_1d_max[i],
                 rescale_pdf_max=True)
    handles, labels = ax[ii,1].get_legend_handles_labels()
    labels = (
        "m=%.2e s=%.2e" % (p_prior_max[i], prior_cov_diag_dim[i]),
        "m=%.2e s=%.2e" % (p_post_max[i], post_cov_diag_dim[i])
    )
    ax[ii,1].legend(handles, labels, loc='upper left', fontsize='xx-small')
    ax[ii,1].set_ylim(0.0, 1.1)

# set log scale for x-axis
for i in idx_lognormal_glo:
    ax[i,0].set_xscale('log')
for ii in range(n_weak):
    ax[ii,1].set_xscale('log')

# set grid lines
for i in range(plot_n_rows):
    for j in [0,1]:
        ax[i,j].grid(True)

# turn off unused axis
for i in range(n_glo, plot_n_rows):
    ax[i,0].axis('off')
for i in range(n_weak, plot_n_rows):
    ax[i,1].axis('off')

# set spacing between subplots
fig.set_tight_layout({'pad': 0.5})

# save plots
fig.savefig(plot_base+"_pdf_1d.png", dpi=360)

########################################
# 2D Plots
########################################

from matplotlib.lines import Line2D

def plot_add_confidence_ellipse(ax, mean, cov, n_std=3.0,
                                param_dim=(1.0, 1.0), phys_glo_dim=(1.0, 1.0),
                                lognormal=(False, False),
                                edgecolor='k', **kwargs):
    """
    Plot a covariance confidence ellipse.
    https://matplotlib.org/gallery/statistics/confidence_ellipse.html
    https://carstenschelp.github.io/2018/09/14/Plot_Confidence_Ellipse_001.html
    """
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # create reference ellipse
    theta = np.linspace(0, 2.0*np.pi, 256+1)
    radius_x = np.sqrt(1 + pearson)
    radius_y = np.sqrt(1 - pearson)
    ellipse = np.array([radius_x * np.cos(theta),
                        radius_y * np.sin(theta)])
    # rotate ellipse counterclockwise by 45 degrees
    rot = np.array([[np.cos(0.25*np.pi), -np.sin(0.25*np.pi)],
                    [np.sin(0.25*np.pi),  np.cos(0.25*np.pi)]])
    ellipse = np.dot(rot, ellipse)
    # scale
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    ellipse[0,:] = scale_x * ellipse[0,:]
    ellipse[1,:] = scale_y * ellipse[1,:]
    # shift by mean
    ellipse = ellipse + mean[:,np.newaxis]
    # apply transfer to physical dimensions
    ellipse[0,:] = transfer_param_to_physical(
            ellipse[0,:], param_dim[0], phys_glo_dim[0], lognormal[0])
    ellipse[1,:] = transfer_param_to_physical(
            ellipse[1,:], param_dim[1], phys_glo_dim[1], lognormal[1])
    # plot ellipse
    return ax.plot(ellipse[0,:], ellipse[1,:], color=edgecolor, **kwargs)

def plot_confidence_2d(ax, indices, m,
                       prior_max, prior_cov, post_max, post_cov,
                       param_dim, phys_glo_dim, idx_lognormal):
    i = indices[0]
    j = indices[1]
    param_dim = (param_dim[i], param_dim[j])
    phys_glo_dim = (phys_glo_dim[i], phys_glo_dim[j])
    lognormal = (i in idx_lognormal, j in idx_lognormal)
    # plot prior confidences
    shift = np.array([prior_max[i], prior_max[j]])
    cov = np.array([ [prior_cov[i,i], prior_cov[i,j]],
                     [prior_cov[j,i], prior_cov[j,j]] ])
    for n_stddev in [1.0]: #[1.0, 2.0]:
        plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
                                    param_dim, phys_glo_dim, lognormal,
                                    edgecolor='0.3', linewidth=2)
    # plot prior max
    p_prior_max = (
        transfer_param_to_physical(prior_max[i], param_dim[0],
                                   phys_glo_dim[0], lognormal[0]),
        transfer_param_to_physical(prior_max[j], param_dim[1],
                                   phys_glo_dim[1], lognormal[1])
    )
    ax.plot(p_prior_max[0], p_prior_max[1], color='0.3', marker='o')
    # plot posterior confidences
    #   (1*sigma, 2*sigma, 3*sigma) = (0.682, 0.954, 0.997) confidence
    shift = np.array([post_max[i], post_max[j]])
    cov = np.array([ [post_cov[i,i], post_cov[i,j]],
                     [post_cov[j,i], post_cov[j,j]] ])
    for n_stddev in [1.0, 2.0, 3.0]:
        plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
                                    param_dim, phys_glo_dim, lognormal,
                                    edgecolor='blue', linewidth=2)
    # plot posterior max
    p_post_max = (
        transfer_param_to_physical(post_max[i], param_dim[0],
                                   phys_glo_dim[0], lognormal[0]),
        transfer_param_to_physical(post_max[j], param_dim[1],
                                   phys_glo_dim[1], lognormal[1])
    )
    ax.plot(p_post_max[0], p_post_max[1], color='blue', marker='.')
    # set axis to logscale
    if lognormal[0]:
        ax.set_xscale('log')
    if lognormal[1]:
        ax.set_yscale('log')

#def generate_plot_data_2d(indices, m, prior_mean, prior_cov, post_mean, post_cov,
#                          p, param_dim, phys_glo_dim, idx_lognormal):
#    i = indices[0]
#    j = indices[1]
#    Pi, Pj = np.meshgrid(p[i,:], p[j,:])
#    Pdf_prior = pdf_2d([i,j], m, prior_mean, prior_cov,
#                       p, param_dim, phys_glo_dim, idx_lognormal)
#    Pdf_post  = pdf_2d([i,j], m, post_mean, post_cov,
#                       p, param_dim, phys_glo_dim, idx_lognormal)
#    return (Pi, Pj, Pdf_prior, Pdf_post)
#
#def plot_pdfs_2d(ax, indices, m,
#                 prior_mean, prior_cov, p_prior_max,
#                 post_mean, post_cov, p_post_max,
#                 p, param_dim, phys_glo_dim, idx_lognormal):
#    i = indices[0]
#    j = indices[1]
#    # generate data
#    (Pi, Pj, Pdf_prior, Pdf_post) = generate_plot_data_2d(
#            indices, m, prior_mean, prior_cov, post_mean, post_cov,
#            p, param_dim, phys_glo_dim, idx_lognormal)
#    # plot prior
#    levels = np.array([0.3, 0.6, 0.9]) * np.max(Pdf_prior)
#    CS = ax.contour(Pi, Pj, Pdf_prior, levels, colors='0.3', linewidths=2)
#    ax.plot(p_prior_max[i], p_prior_max[j], color='0.3', marker='o')
#    # plot posterior
#    levels = np.array([0.3, 0.6, 0.9]) * np.max(Pdf_post)
#    CS = ax.contour(Pi, Pj, Pdf_post, levels, colors='blue', linewidths=2)
#    ax.plot(p_post_max[i], p_post_max[j], color='blue', marker='o')

# create custom lines for legend
custom_lines = [Line2D([0], [0], color='0.3', linewidth=2),
                Line2D([0], [0], color='blue', linewidth=2)]

# create figure
fig, ax = plt.subplots(2, 2, figsize=(8, 8))

# plot (stress exponent vs. yield stress)
plot_confidence_2d(ax[0,0], [4,5], m,
                   prior_max, prior_cov, post_max, post_cov,
                   param_dim, phys_glo_dim, idx_lognormal)
ax[0,0].set_xlabel(PLOT_TITLES_GLO[4])
ax[0,0].set_ylabel(PLOT_TITLES_GLO[5])
ax[0,0].legend(custom_lines, ('prior', 'posterior'),
               loc='upper left', fontsize='small')

# plot (stress exponent vs. activation energy)
plot_confidence_2d(ax[0,1], [4,1], m,
                   prior_max, prior_cov, post_max, post_cov,
                   param_dim, phys_glo_dim, idx_lognormal)
ax[0,1].set_xlabel(PLOT_TITLES_GLO[4])
ax[0,1].set_ylabel(PLOT_TITLES_GLO[1])
ax[0,1].legend(custom_lines, ('prior', 'posterior'),
               loc='upper left', fontsize='small')

# plot (stress exponent vs. scaling in upper mantle)
plot_confidence_2d(ax[1,0], [4,0], m,
                   prior_max, prior_cov, post_max, post_cov,
                   param_dim, phys_glo_dim, idx_lognormal)
ax[1,0].set_xlabel(PLOT_TITLES_GLO[4])
ax[1,0].set_ylabel(PLOT_TITLES_GLO[0])
ax[1,0].legend(custom_lines, ('prior', 'posterior'),
               loc='upper left', fontsize='small')

# plot (stress exponent vs. scaling in lower mantle)
plot_confidence_2d(ax[1,1], [4,2], m,
                   prior_max, prior_cov, post_max, post_cov,
                   param_dim, phys_glo_dim, idx_lognormal)
ax[1,1].set_xlabel(PLOT_TITLES_GLO[4])
ax[1,1].set_ylabel(PLOT_TITLES_GLO[2])
ax[1,1].legend(custom_lines, ('prior', 'posterior'),
               loc='upper left', fontsize='small')

# set grid lines
for i in [0,1]:
    for j in [0,1]:
        ax[i,j].grid(True)

# set spacing between subplots
fig.set_tight_layout({'pad': 0.5})

# save & show plots
fig.savefig(plot_base+"_confidence_2d.png", dpi=360)
#plt.show()

########################################
# Correlation Plots
########################################

# compute correlation matrix
post_corr = compute_correlation(post_cov)

# create figure
fig, ax = plt.subplots(1, 1)

# plot singular values and vectors
c = ax.pcolor(np.abs(post_corr), cmap='hot')

# annotate plots
ax.set_title("Correlation matrix of posterior")
for row in range(len(post_corr)):
    for col in range(len(post_corr)):
        if row == col:
            continue
        x = col + 0.5
        y = row + 0.5
        if 0.5 < np.abs(post_corr[row,col]):
            textcolor = "k"
        else:
            textcolor = "w"
        v = "%+.2f" % post_corr[row,col]
        ax.text(x, y, v, color=textcolor,
                ha="center", va="center", fontsize='x-small')
cb1 = fig.colorbar(c, ax=ax, label="magnitude of correlation")
ax.set_ylim(len(post_corr), 0) # invert y-axis

# save & show plots
fig.savefig(plot_base+"_correlation.png")
#plt.show()

########################################
# SVD Plots
########################################

# compute SVD of posterior covarance
U, s, Vh = np.linalg.svd(post_cov)

# want to plot SVD of inverse covariance
s_plot = 1.0/np.fliplr(s[np.newaxis,:])
U_plot = np.fliplr(np.abs(U))

# create figure
fig, ax = plt.subplots(2, 1)

# plot singular values and vectors
c0 = ax[0].pcolor(s_plot, cmap='hot')
c1 = ax[1].pcolor(U_plot, cmap='hot')

# annotate plots
ax[0].set_title("SVD of inverse covariance matrix (Hessian)")
ax[0].tick_params(
    axis='y',          # changes apply to the y-axis
    which='both',      # both major and minor ticks are affected
    left=False,        # ticks along the left edge are off
    right=False,       # ticks along the right edge are off
    labelleft=False)   # labels along the left edge are off
for k in range(len(s)):
    if k == 0:
        textcolor = "k"
    else:
        textcolor = "w"
    ax[0].text((k+0.5), 0.5, "%.2e"%s_plot[0,k], color=textcolor,
               ha="center", va="center", fontsize='x-small')
cb0 = fig.colorbar(c0, ax=ax[0], label="singular values")
cb1 = fig.colorbar(c1, ax=ax[1], label="magnitude of singular vectors")
ax[1].set_ylim(len(s), 0) # invert y-axis

# set spacing between subplots
fig.set_tight_layout({'pad': 0.5})

# save & show plots
fig.savefig(plot_base+"_svd_inv_covariance.png")
plt.show()
