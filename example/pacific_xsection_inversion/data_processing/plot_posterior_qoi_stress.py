#!/usr/bin/python

from posterior_util import (transfer_param_to_physical,
                            transfer_param_to_physical_multi,
                            transfer_physical_to_param_multi,
                            transfer_physical_to_param_deriv_multi,
                            gaussian_1d, pdf_1d, plot_pdfs_1d,
                            plot_marginal_2d, plot_conditional_2d,
                            compute_correlation,
                            OOMFormatter)
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Environment
###############################################################################

# set environment
argn = len(sys.argv) - 1
script_path = os.path.abspath(sys.argv[0])
script_dir  = os.path.dirname(script_path)
script_base = os.path.splitext(os.path.basename(script_path))[0] # w/o ext.

def print_info(message):
    '''Prints messages in stdout.'''
    print("["+script_base+"] "+message)

###############################################################################
# Setup
###############################################################################

# set input arguments
assert 3 <= argn
post_cov_txt  = sys.argv[1]
qoi_max_txt   = sys.argv[2]
qoi_jac_txt   = sys.argv[3]
invert_post_cov = True

# set output path
if 4 <= argn:
    plot_base = sys.argv[4]
else:
    data_dir = os.path.dirname(os.path.abspath(post_cov_txt))
    plot_base = data_dir+"/"+script_base

# set plotting parameters
PLOT_N_NODES = 1024+1
PLOT_ORDER_MARGINAL_1D = 6
PLOT_TITLE_QOI = 'QOI'
PLOT_TITLE_MARGINAL_1D = [
    'Average normal component of stress along weak zone',
    'Average tangential component of stress along weak zone',
]
PLOT_XLABEL = 'Component of stress [MPa]'
PLOT_XLABEL_MARGINAL_1D = [
    'Normal stress [MPa]',
    'Tangential stress [MPa]',
]
PLOT_LEGEND_MARGINAL_1D = [
    'Chile',
    'Mariana',
    'Ryuku',
]
PLOT_MARGINALS_ALL    = True
PLOT_CONDITIONALS_ALL = True
PLOT_SUFFIX = ['png', 'eps']

# reduce stress tangential components to norm
REDUCE_TANG_TO_NORM = True

# set physical parameters for nondimensialization
stress_dim = [
    1.0e20*1.0e-6/6371.0e3**2
]

# set axis limits in physical dimensions
stress_lim = np.array([
    [0.0, 1.5e8],
    [0.0, 1.5e8],
    [0.0, 1.5e8],
    [0.0, 0.4e8],
    [0.0, 0.4e8],
    [0.0, 0.4e8],
])

# set counts
n_qoi_all = 9
n_weak    = 3

# set indices
idx_stressNorm  = 0  # until 2
idx_stressTangX = 3  # until 5
idx_stressTangZ = 6  # until 8

# set colors
COLOR_PRIOR='dimgray'
COLOR_POST=np.array([117,112,179])/255.0 #'darkviolet'
COLOR_POST_MARGINAL_1D = [
    np.array([ 27,158,119])/255.0,
    np.array([217, 95,  2])/255.0,
    np.array([117,112,179])/255.0,
]
COLOR_POST_COND=np.array([217, 95,  2])/255.0 #'orange'

# set line styles
LINESTYLE_POST_MARGINAL_1D = [
    (0, (1,1)), # densly dotted
    'dashed',
    'dashdot',
]

###############################################################################
# Functions
###############################################################################

def reduce_stress_tangential_to_norm(qoi_max_all, qoi_cov_all, n_weak, n_samples=100000):
    assert qoi_max_all.size == 3*n_weak, qoi_max_all.size
    assert qoi_cov_all.shape == (3*n_weak, 3*n_weak), qoi_cov_all.size
    # sample from gaussian
    rng     = np.random.default_rng()
    samples = rng.multivariate_normal(qoi_max_all, qoi_cov_all, n_samples)
    # compute magnitude of sampled tangentials
    magn = np.empty((n_samples,n_weak))
    x    = samples[:,1*n_weak:2*n_weak]
    z    = samples[:,2*n_weak:3*n_weak]
    for i in range(n_weak):
        magn[:,i] = np.sqrt(x[:,i]**2 + z[:,i]**2)
    # init output
    n       = 2*n_weak
    samples = np.concatenate((samples[:,:n_weak], magn), axis=1)
    qoi_cov = np.cov(samples, rowvar=False)
    qoi_max = np.zeros(n)
    # copy mean of stress normal
    qoi_max[:n_weak] = qoi_max_all[:n_weak]
    #qoi_cov[:n_weak,:n_weak] = qoi_cov_all[:n_weak,:n_weak]
    # compute mean of magnitude
    x    = qoi_max_all[1*n_weak:2*n_weak]
    z    = qoi_max_all[2*n_weak:3*n_weak]
    for i in range(n_weak):
        qoi_max[n_weak+i] = np.sqrt(x[i]**2 + z[i]**2)
    # return reduced max and covariance
    return qoi_max, qoi_cov, samples

###############################################################################
# Main
###############################################################################

print_info("INPUT")
print_info("- Posterior covariance: "+post_cov_txt)
print_info("- QOI max:              "+qoi_max_txt)
print_info("- QOI covariance:       "+qoi_jac_txt)
print_info("OUTPUT")
print_info("- Plot: "+plot_base)

# load files
post_cov  = np.loadtxt(post_cov_txt)
qoi_max   = np.loadtxt(qoi_max_txt)
qoi_jac   = np.loadtxt(qoi_jac_txt)

# set size
n_qoi = len(qoi_max)

# invert matrices
if invert_post_cov:
    post_cov = np.linalg.inv(0.5*(post_cov + post_cov.T))

# we assume that the jacobian was computed with a finite difference unit
# perturbation of [0...0, 0.1, 0...0], thus scale by 10
qoi_jac = 10.0*qoi_jac

# compute QOI covariance matrix: Cov_qoi = Jac * Cov_post * Jac^T
qoi_cov = np.dot(qoi_jac, np.dot(post_cov, qoi_jac.T))

# reduce stress tangentials to their norm
qoi_samples = None
if REDUCE_TANG_TO_NORM and n_qoi == n_qoi_all:
    qoi_max, qoi_cov, qoi_samples = reduce_stress_tangential_to_norm(qoi_max, qoi_cov, n_weak)
    n_qoi = len(qoi_max)

# set mean and covariance
qoi_cov_diag = np.diag(qoi_cov)
qoi_mean     = np.copy(qoi_max)

# set physical dimensions
param_dim = np.array(n_qoi*[1.0])
phys_dim  = np.array(n_qoi*stress_dim)
phys_lim  = stress_lim  #np.repeat(stress_lim[np.newaxis,:], n_qoi, axis=0)

# convert variances to standard deviations in physical dimensions
qoi_std_diag_dim = np.zeros_like(qoi_cov_diag)
for i in range(n_qoi):
    qoi_std_diag_dim[i] = transfer_param_to_physical(
            np.sqrt(qoi_cov_diag[i]), param_dim[i], phys_dim[i], False)

# convert samples to physical dimensions
qoi_samples_dim = None
if qoi_samples is not None:
    qoi_samples_dim = transfer_param_to_physical_multi(qoi_samples.T, param_dim, phys_dim, []).T

print_info("PARAMETERS")
print_info("- QOI count: %i" % n_qoi)
print_info("DISTRIBUTIONS")
np.set_printoptions(precision=3, linewidth=150)
print_info("QOI max")
print(qoi_max)
print_info("QOI mean")
print(qoi_mean)
print_info("QOI stddev")
print(np.sqrt(qoi_cov_diag))

# create range of parameters for plotting
m_qoi = np.linspace(0.5, 4.0*np.mean(qoi_mean), num=PLOT_N_NODES)
m = np.tile(m_qoi, (n_qoi, 1))

# transfer inversion parameters to physical values
p = transfer_param_to_physical_multi(m, param_dim, phys_dim, [])
print_info(
    "Error of transfer fn. vs. its inverse: max abs err %g" % \
    np.max(np.abs(m - transfer_physical_to_param_multi(p, param_dim, phys_dim, [])))
)

# compute max coordinates
m_max = qoi_max[:,np.newaxis]
p_max = transfer_param_to_physical_multi(m_max, param_dim, phys_dim, [])
pdf_1d_max  = \
    np.abs(transfer_physical_to_param_deriv_multi(
            p_max, param_dim, phys_dim, [])) * \
    gaussian_1d(m_max, qoi_mean, qoi_cov_diag)

# compute 1-dim probability densities
pdf_1d = pdf_1d(m, qoi_mean, qoi_cov_diag,
                p, param_dim, phys_dim, [])

# integrate densities
int_pdf = np.empty(n_qoi)
for i in range(n_qoi):
    int_pdf[i] = np.trapz(pdf_1d[i,:], p[i,:])
print_info("Integrals of QOI pdf's")
print(int_pdf)

########################################
# Text Output
########################################

with open(plot_base+'_data.txt', 'w') as f:
    f.write('QOI: index, max, standard deviation\n')
    for i in range(n_qoi):
        f.write('%d, %.6e, %.6e\n' % (i, p_max[i], qoi_std_diag_dim[i]))

########################################
# 1D Plots
########################################

# create figure
plot_n_combine = 3
plot_n_rows = n_qoi//plot_n_combine
fig, ax = plt.subplots(plot_n_rows, 1, figsize=(8, 4))

# create plots
for i in range(n_qoi):
    plot_row   = i // plot_n_combine
    plot_layer = i %  plot_n_combine
    #ax[plot_row].set_title(PLOT_TITLE_MARGINAL_1D[plot_row])
    x        = p[i,:]/10.0**PLOT_ORDER_MARGINAL_1D
    x1_max   = p_max[i]/10.0**PLOT_ORDER_MARGINAL_1D
    pdf1     = pdf_1d[i,:] # /pdf_1d_max[i]
    pdf1_max = np.amax(pdf_1d[i,:]) # 1.0
    plot_pdfs_1d(ax[plot_row], x,
                 None, None, None, None, # no prior
                 pdf1, x1_max, pdf1_max, COLOR_POST_MARGINAL_1D[plot_layer],
                 pdf1_linestyle=LINESTYLE_POST_MARGINAL_1D[plot_layer])
    ax[plot_row].set_xlabel(PLOT_XLABEL_MARGINAL_1D[plot_row])
    ax[plot_row].set_xlim(phys_lim[i,0]/10.0**PLOT_ORDER_MARGINAL_1D,
                          phys_lim[i,1]/10.0**PLOT_ORDER_MARGINAL_1D)
    #ax[plot_row].xaxis.set_major_formatter(OOMFormatter(PLOT_ORDER_MARGINAL_1D, '%1.1f'))
    ax[plot_row].get_yaxis().set_visible(False)
    if False:  #qoi_samples_dim is not None:
        ax[plot_row].hist(qoi_samples_dim[:,i], bins=32, density=True, stacked=True)
    else:
        ylim = ax[plot_row].get_ylim()
        if ylim[1] < 1.1*pdf1_max:
            ax[plot_row].set_ylim(0.0, 1.1*pdf1_max)

# set grid lines
for plot_row in range(plot_n_rows):
    handles, labels = ax[plot_row].get_legend_handles_labels()
    labels = PLOT_LEGEND_MARGINAL_1D
    ax[plot_row].legend(handles, labels, loc='upper right', fontsize='small')
    ax[plot_row].grid(True)

# set spacing between subplots
fig.set_tight_layout({'pad': 2.0})

# save plots
for suffix in PLOT_SUFFIX:
    fig.savefig('{}_marginal_1d.{}'.format(plot_base, suffix), dpi=360)
#plt.show()

########################################
# 2D Plots
########################################

# create figure with all 2D maginals of the posterior
if PLOT_MARGINALS_ALL:
    fig, ax = plt.subplots(n_qoi, n_qoi, figsize=(16, 16))
    for axcol in range(0, n_qoi):
        ax[0,axcol].set_title("%s %i" % (PLOT_TITLE_QOI, axcol))
        for axrow in range(0, n_qoi):
            axcurr = ax[axrow,axcol]
            if axcol == 0:
                axcurr.set_ylabel("%s %i" % (PLOT_TITLE_QOI, axrow))
            if axrow == axcol:
                axcurr.get_xaxis().set_ticks([])
                axcurr.get_yaxis().set_ticks([])
                continue
            plot_marginal_2d(axcurr, [axcol,axrow], m,
                             None, None, None, # no prior
                             qoi_max, qoi_cov, COLOR_POST,
                             param_dim, phys_dim, [], [])
            axcurr.set_xlim(phys_lim[axcol,0], phys_lim[axcol,1])
            axcurr.set_ylim(phys_lim[axrow,0], phys_lim[axrow,1])
            axcurr.grid(True)
    # set spacing between subplots
    fig.tight_layout()
    # save & show plots
    for suffix in PLOT_SUFFIX:
        fig.savefig('{}_marginal_2d_all.{}'.format(plot_base, suffix), dpi=360)
    #plt.show()

# create figure with all 2D conditionals of the posterior
if PLOT_CONDITIONALS_ALL:
    fig, ax = plt.subplots(n_qoi, n_qoi, figsize=(16, 16))
    for axcol in range(0, n_qoi):
        ax[0,axcol].set_title("%s %i" % (PLOT_TITLE_QOI, axcol))
        for axrow in range(0, n_qoi):
            axcurr = ax[axrow,axcol]
            if axcol == 0:
                axcurr.set_ylabel("%s %i" % (PLOT_TITLE_QOI, axrow))
            if axrow == axcol:
                axcurr.get_xaxis().set_ticks([])
                axcurr.get_yaxis().set_ticks([])
                continue
            plot_conditional_2d(axcurr, [axcol,axrow], m,
                                None, None, None, # no prior
                                qoi_max, qoi_cov, COLOR_POST_COND,
                                param_dim, phys_dim, [], [])
            axcurr.set_xlim(phys_lim[axcol,0], phys_lim[axcol,1])
            axcurr.set_ylim(phys_lim[axrow,0], phys_lim[axrow,1])
            axcurr.grid(True)
    # set spacing between subplots
    fig.tight_layout()
    # save & show plots
    for suffix in PLOT_SUFFIX:
        fig.savefig('{}_conditional_2d_all.{}'.format(plot_base, suffix), dpi=360)
    #plt.show()

########################################
# Correlation Plots
########################################

# compute correlation matrix
corr = compute_correlation(qoi_cov)

# create figure
fig, ax = plt.subplots(1, 1)

# plot singular values and vectors
c = ax.pcolor(np.abs(corr), cmap='hot')

# annotate plots
ax.set_title("Correlation matrix of %s" % PLOT_TITLE_QOI)
for row in range(len(corr)):
    for col in range(len(corr)):
        if row == col:
            continue
        x = col + 0.5
        y = row + 0.5
        if 0.5 < np.abs(corr[row,col]):
            textcolor = "k"
        else:
            textcolor = "w"
        v = "%+.2f" % corr[row,col]
        ax.text(x, y, v, color=textcolor,
                ha="center", va="center", fontsize='x-small')
cb1 = fig.colorbar(c, ax=ax, label="magnitude of correlation")
ax.set_ylim(len(corr), 0) # invert y-axis

# save & show plots
fig.savefig(plot_base+"_correlation.png")
plt.show()
