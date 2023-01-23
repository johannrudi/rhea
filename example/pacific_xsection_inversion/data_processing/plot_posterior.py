#!/usr/bin/python

from posterior_util import (adjust_mean_of_lognormal,
                            transfer_param_to_physical,
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
import scipy.linalg
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
assert 5 <= argn
param_dim_txt = sys.argv[1]
prior_max_txt = sys.argv[2]
prior_cov_txt = sys.argv[3]
post_max_txt  = sys.argv[4]
post_cov_txt  = sys.argv[5]
invert_prior_cov = False
invert_post_cov = True

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
    "Activation energy [J/mol]", #in upper mantle [J/mol]",
    "Scaling in lower mantle",
   #"Activation energy in lower mantle [J/mol]",
    "Stress exponent",
    "Yield stress [MPa]"
]
PLOT_TITLE_WEAK = "Weak zone factor"
PLOT_MARGINALS_ALL    = False
PLOT_CONDITIONALS_ALL = False
PLOT_SUFFIX = ['png', 'eps']

# set colors
COLOR_PRIOR=np.array([150,150,150])/255.0 #'dimgray'
COLOR_POST=np.array([152,78,163])/255.0 #'darkviolet'
COLOR_POST_COND=np.array([27,158,119])/255.0 #'darkcyan'

# set physical parameters for nondimensialization
phys_glo_dim = np.array([
    1.0,        # scaling in upper mantle
    8.314*1400, # activation energy in upper mantle
    1.0,        # scaling in lower mantle
   #8.314*1400, # activation energy in lower mantle
    1.0,        # stress exponent
    1.0e20*1.0e-6/6371.0e3**2 / 1.0e6  # yield stress in MPa
])

# set axis limits in physical dimensions
phys_glo_lim = np.array([
    [1.0e6, 1.0e10], # scaling in upper mantle
    [2.0e5, 8.0e5],  # activation energy in upper mantle
    [1.0e0, 1.0e4],  # scaling in lower mantle
   #[2.0e5, 8.0e5],  # activation energy in lower mantle
    [1.0  , 6.0],    # stress exponent
    [1.0e1, 2.0e2]   # yield stress in MPa
])

# set axis formatter
axis_formatter = [
    None,                    # scaling in upper mantle
    OOMFormatter(5, "%.0f"), # activation energy in upper mantle
    None,                    # scaling in lower mantle
   #None,                    # activation energy in lower mantle
    None,                    # stress exponent
    None,                    # yield stress
]

# set indices
idx_scaling_um    = 0
idx_activation_um = 1
idx_scaling_lm    = 2
idx_activation_lm = None
idx_stress_exp    = 3
idx_yield_stress  = 4

# set indices of densities that are lognormal
idx_lognormal_glo = [idx_scaling_um, idx_scaling_lm]

# set indices of densities that are inverted: (.)^-1
idx_inverted = [idx_stress_exp]

###############################################################################
# Main
###############################################################################

print_info("INPUT")
print_info("- Parameter dim:        "+param_dim_txt)
print_info("- Prior max:            "+prior_max_txt)
print_info("- Prior covariance:     "+prior_cov_txt)
print_info("- Posterior max:        "+post_max_txt)
print_info("- Posterior covariance: "+post_cov_txt)
print_info("OUTPUT")
print_info("- Plot: "+plot_base)

# load files
param_dim = np.loadtxt(param_dim_txt)
prior_max = np.loadtxt(prior_max_txt)
prior_cov = np.loadtxt(prior_cov_txt)
post_max  = np.loadtxt(post_max_txt)
post_cov  = np.loadtxt(post_cov_txt)

# compute posterior covariance without prior
if invert_prior_cov:
    prior_cov_sqrt = scipy.linalg.sqrtm(np.linalg.inv(prior_cov))
else:
    prior_cov_sqrt = scipy.linalg.sqrtm(prior_cov)
if invert_post_cov:
    post_cov_withPriorPC = np.dot(np.dot(prior_cov_sqrt, np.copy(post_cov)), prior_cov_sqrt)
else:
    post_cov_withPriorPC = np.dot(np.dot(prior_cov_sqrt, np.linalg.inv(post_cov)), prior_cov_sqrt)
post_cov_withPriorPC = np.linalg.inv(0.5*(post_cov_withPriorPC + post_cov_withPriorPC.T))
post_cov_withPriorPC_diag = np.diag(post_cov_withPriorPC)
print(post_cov_withPriorPC_diag)

# invert matrices
if invert_prior_cov:
    prior_cov = np.linalg.inv(prior_cov)
if invert_post_cov:
    post_cov = np.linalg.inv(0.5*(post_cov + post_cov.T))

# set sizes and indices
n_glo    = len(phys_glo_dim)
n_weak   = len(param_dim) - n_glo
idx_glo  = np.arange(n_glo)
idx_weak = np.arange(n_glo, n_glo+n_weak)
idx_lognormal = np.hstack((idx_lognormal_glo, idx_weak))

# set mean and covariance of prior
prior_cov_diag = np.diag(prior_cov)
prior_mean     = np.copy(prior_max)
for i in idx_lognormal:
    prior_mean[i] = adjust_mean_of_lognormal(prior_mean[i], prior_cov_diag[i], param_dim[i])
# set mean and covariance of posterior
post_cov_diag = np.diag(post_cov)
post_mean     = np.copy(post_max)
for i in idx_lognormal:
    post_mean[i] = adjust_mean_of_lognormal(post_mean[i], post_cov_diag[i], param_dim[i])

# set physical dimensions
phys_dim = np.concatenate([phys_glo_dim, np.ones(n_weak)])
phys_lim = np.concatenate(
        [phys_glo_lim, np.repeat([[1.0e-10, 1.0]], n_weak, axis=0)],
        axis=0)
# convert variances to physical dimensions
prior_std_diag_dim = np.zeros_like(prior_cov_diag)
post_std_diag_dim  = np.zeros_like(post_cov_diag)
for i in range(n_glo+n_weak):
    prior_std_diag_dim[i] = transfer_param_to_physical(
            np.sqrt(prior_cov_diag[i]), param_dim[i], phys_dim[i],
            i in idx_lognormal)
    post_std_diag_dim[i] = transfer_param_to_physical(
            np.sqrt(post_cov_diag[i]), param_dim[i], phys_dim[i],
            i in idx_lognormal)

print_info("PARAMETERS")
print_info("- global parameters count:   %i" % n_glo)
print_info("- weakzone parameters count: %i" % n_weak)
print_info("DISTRIBUTIONS")
np.set_printoptions(precision=3, linewidth=150)
print_info("Compare prior max vs. posterior max")
print(np.array([prior_max, post_max]))
print_info("Compare prior mean vs. posterior mean")
print(np.array([prior_mean, post_mean]))
print_info("Compare prior stddev vs. posterior stddev")
print(np.array([np.sqrt(prior_cov_diag), np.sqrt(post_cov_diag)]))
print_info("Posterior stddev with prior preconditioning")
print(post_cov_withPriorPC_diag)

# create range of parameters for plotting
m_glo  = np.linspace(0.01, 4.0, num=PLOT_N_NODES)
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

# apply inverted transfer function
p[idx_inverted,:]            = np.flipud(1.0/p[idx_inverted,:])
p_prior_max[idx_inverted]    = 1.0/p_prior_max[idx_inverted]
p_post_max[idx_inverted]     = 1.0/p_post_max[idx_inverted]
pdf_prior_1d[idx_inverted,:] = np.flipud(pdf_prior_1d[idx_inverted,:])
pdf_post_1d[idx_inverted,:]  = np.flipud(pdf_post_1d[idx_inverted,:])

# integrate densities
int_pdf_prior = np.empty(n_glo+n_weak)
int_pdf_post  = np.empty(n_glo+n_weak)
for i in range(n_glo+n_weak):
    int_pdf_prior[i] = np.trapz(pdf_prior_1d[i,:], p[i,:])
    int_pdf_post[i]  = np.trapz(pdf_post_1d[i,:], p[i,:])
print_info("Compare integrals of prior pdf's vs. posterior pdf's")
print(np.array([int_pdf_prior, int_pdf_post]))

########################################
# Text Output
########################################

with open(plot_base+'_data.txt', 'w') as f:
    f.write('Info: index, axis min, axis max, lognormal, inverted, name\n')
    for i in range(n_glo):
        f.write( '%d, %.6e, %.6e, %d, %d, "%s"\n' %
                 (i, phys_lim[i,0], phys_lim[i,1],
                  i in idx_lognormal_glo, i in idx_inverted,
                  PLOT_TITLES_GLO[i]) )
    for i in range(n_weak):
        f.write( '%d, %.6e, %.6e, %d, %d, "%s"\n' %
                 (n_glo+i, phys_lim[i,0], phys_lim[i,1], 1, 0,
                 PLOT_TITLE_WEAK+(" %i" % i)) )
    f.write('\n')

    f.write('Prior: index, max, standard deviation\n')
    for i in range(n_glo + n_weak):
        f.write('%d, %.6e, %.6e\n' % (i, p_prior_max[i], prior_std_diag_dim[i]))
    f.write('\n')

    f.write('Posterior: index, max, standard deviation\n')
    for i in range(n_glo + n_weak):
        f.write('%d, %.6e, %.6e\n' % (i, p_post_max[i], post_std_diag_dim[i]))

########################################
# 1D Plots
########################################

# create title strings
title = np.empty(n_glo+n_weak, dtype=object)
if 0 < n_glo:
    title[0:n_glo] = PLOT_TITLES_GLO
for i in range(n_weak):
    title[n_glo+i] = PLOT_TITLE_WEAK+(" %i" % i)

# create figure
plot_n_rows = np.max([n_glo, n_weak])
fig, ax = plt.subplots(plot_n_rows, 2, figsize=(10, 8))

# create plots
for i in range(n_glo):
    this_title = PLOT_TITLES_GLO[i]
    ax[i,0].set_title(this_title)
    x = p[i,:]
    x0_max = p_prior_max[i]
    x1_max = p_post_max[i]
    pdf0 = pdf_prior_1d[i,:]/pdf_prior_1d_max[i]
    pdf1 = pdf_post_1d[i,:]/pdf_post_1d_max[i]
    pdf0_max = 1.0
    pdf1_max = 1.0
    plot_pdfs_1d(ax[i,0], x,
                 pdf0, x0_max, pdf0_max, COLOR_PRIOR,
                 pdf1, x1_max, pdf1_max, COLOR_POST)
    handles, labels = ax[i,0].get_legend_handles_labels()
    labels = ["m=%.2e s=%.2e" % (x0_max, prior_std_diag_dim[i]),
              "m=%.2e s=%.2e" % (x1_max, post_std_diag_dim[i])]
    ax[i,0].legend(handles, labels, loc='upper right', fontsize='xx-small')
    ax[i,0].ticklabel_format(axis='x', style='sci', scilimits=(1,3))
    ax[i,0].set_xlim(phys_glo_lim[i,0], phys_glo_lim[i,1])
    ax[i,0].set_ylim(0.0, 1.1)
    if axis_formatter[i] is not None:
        ax[i,0].xaxis.set_major_formatter(axis_formatter[i])
for ii in range(n_weak):
    i = n_glo + ii
    this_title = PLOT_TITLE_WEAK+(" %i" % ii)
    ax[ii,1].set_title(this_title)
    x = p[i,:]
    x0_max = p_prior_max[i]
    x1_max = p_post_max[i]
    pdf0 = pdf_prior_1d[i,:]/pdf_prior_1d_max[i]
    pdf1 = pdf_post_1d[i,:]/pdf_post_1d_max[i]
    pdf0_max = 1.0
    pdf1_max = 1.0
    plot_pdfs_1d(ax[ii,1], x,
                 pdf0, x0_max, pdf0_max, COLOR_PRIOR,
                 pdf1, x1_max, pdf1_max, COLOR_POST)
    handles, labels = ax[ii,1].get_legend_handles_labels()
    labels = ["m=%.2e s=%.2e" % (x0_max, prior_std_diag_dim[i]),
              "m=%.2e s=%.2e" % (x1_max, post_std_diag_dim[i])]
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
for suffix in PLOT_SUFFIX:
    fig.savefig('{}_marginal_1d.{}'.format(plot_base, suffix), dpi=360)

########################################
# 2D Plots
########################################

from matplotlib.lines import Line2D

# create custom lines for legend
custom_lines = [Line2D([0], [0], color=COLOR_PRIOR, linewidth=2, linestyle='dashed'),
                Line2D([0], [0], color=COLOR_POST, linewidth=2, linestyle='solid')]

# create figure
fig, ax = plt.subplots(2, 2, figsize=(8, 8))

# plot (stress exponent vs. yield stress)
if idx_stress_exp is not None and idx_yield_stress is not None:
    plot_marginal_2d(ax[0,0], [idx_stress_exp,idx_yield_stress], m,
                     prior_max, prior_cov, COLOR_PRIOR,
                     post_max, post_cov, COLOR_POST,
                     param_dim, phys_glo_dim, idx_lognormal, idx_inverted)
    ax[0,0].set_xlabel(PLOT_TITLES_GLO[idx_stress_exp])
    ax[0,0].set_ylabel(PLOT_TITLES_GLO[idx_yield_stress])
    ax[0,0].legend(custom_lines, ('prior', 'posterior'),
                   loc='upper left', fontsize='small')
    ax[0,0].set_xlim(phys_glo_lim[idx_stress_exp,0],
                     phys_glo_lim[idx_stress_exp,1])
    ax[0,0].set_ylim(phys_glo_lim[idx_yield_stress,0],
                     phys_glo_lim[idx_yield_stress,1])
    if axis_formatter[idx_yield_stress] is not None:
        ax[0,0].yaxis.set_major_formatter(axis_formatter[idx_yield_stress])

# plot (stress exponent vs. activation energy)
if idx_stress_exp is not None and idx_activation_um is not None:
    plot_marginal_2d(ax[0,1], [idx_stress_exp,idx_activation_um], m,
                     prior_max, prior_cov, COLOR_PRIOR,
                     post_max, post_cov, COLOR_POST,
                     param_dim, phys_glo_dim, idx_lognormal, idx_inverted)
    ax[0,1].set_xlabel(PLOT_TITLES_GLO[idx_stress_exp])
    ax[0,1].set_ylabel(PLOT_TITLES_GLO[idx_activation_um])
    ax[0,1].legend(custom_lines, ('prior', 'posterior'),
                   loc='upper left', fontsize='small')
    ax[0,1].set_xlim(phys_glo_lim[idx_stress_exp,0],
                     phys_glo_lim[idx_stress_exp,1])
    ax[0,1].set_ylim(phys_glo_lim[idx_activation_um,0],
                     phys_glo_lim[idx_activation_um,1])
    if axis_formatter[idx_activation_um] is not None:
        ax[0,1].yaxis.set_major_formatter(axis_formatter[idx_activation_um])

# plot (stress exponent vs. scaling in upper mantle)
if idx_stress_exp is not None and idx_scaling_um is not None:
    plot_marginal_2d(ax[1,0], [idx_stress_exp,idx_scaling_um], m,
                     prior_max, prior_cov, COLOR_PRIOR,
                     post_max, post_cov, COLOR_POST,
                     param_dim, phys_glo_dim, idx_lognormal, idx_inverted)
    ax[1,0].set_xlabel(PLOT_TITLES_GLO[idx_stress_exp])
    ax[1,0].set_ylabel(PLOT_TITLES_GLO[idx_scaling_um])
    ax[1,0].legend(custom_lines, ('prior', 'posterior'),
                   loc='upper left', fontsize='small')
    ax[1,0].set_xlim(phys_glo_lim[idx_stress_exp,0],
                     phys_glo_lim[idx_stress_exp,1])
    ax[1,0].set_ylim(phys_glo_lim[idx_scaling_um,0],
                     phys_glo_lim[idx_scaling_um,1])

# plot (stress exponent vs. scaling in lower mantle)
if idx_stress_exp is not None and idx_scaling_lm is not None:
    plot_marginal_2d(ax[1,1], [idx_stress_exp,idx_scaling_lm], m,
                     prior_max, prior_cov, COLOR_PRIOR,
                     post_max, post_cov, COLOR_POST,
                     param_dim, phys_glo_dim, idx_lognormal, idx_inverted)
    ax[1,1].set_xlabel(PLOT_TITLES_GLO[idx_stress_exp])
    ax[1,1].set_ylabel(PLOT_TITLES_GLO[idx_scaling_lm])
    ax[1,1].legend(custom_lines, ('prior', 'posterior'),
                   loc='upper left', fontsize='small')
    ax[1,1].set_xlim(phys_glo_lim[idx_stress_exp,0],
                     phys_glo_lim[idx_stress_exp,1])
    ax[1,1].set_ylim(phys_glo_lim[idx_scaling_lm,0],
                     phys_glo_lim[idx_scaling_lm,1])

# set grid lines
for i in [0,1]:
    for j in [0,1]:
        ax[i,j].grid(True)
# set spacing between subplots
fig.set_tight_layout({'pad': 2.0})
# save & show plots
for suffix in PLOT_SUFFIX:
    fig.savefig('{}_marginal_2d.{}'.format(plot_base, suffix), dpi=360)
#plt.show()

# create figure with all 2D maginals of the posterior
if PLOT_MARGINALS_ALL:
    fig, ax = plt.subplots(n_glo+n_weak, n_glo+n_weak, figsize=(26, 22))
    for axcol in range(0, n_glo+n_weak):
        if axcol < n_glo:
            ax[0,axcol].set_title(PLOT_TITLES_GLO[axcol])
        else:
            ax[0,axcol].set_title(PLOT_TITLE_WEAK+(" %i" % (axcol-n_glo)))
        for axrow in range(0, n_glo+n_weak):
            axcurr = ax[axrow,axcol]
            if axcol == 0:
                if axrow < n_glo:
                    axcurr.set_ylabel(PLOT_TITLES_GLO[axrow])
                else:
                    axcurr.set_ylabel(PLOT_TITLE_WEAK+(" %i" % (axrow-n_glo)))
            if axrow == axcol:
                axcurr.get_xaxis().set_ticks([])
                axcurr.get_yaxis().set_ticks([])
                continue
            plot_marginal_2d(axcurr, [axcol,axrow], m,
                             prior_max, prior_cov, COLOR_PRIOR,
                             post_max, post_cov, COLOR_POST,
                             param_dim, phys_dim, idx_lognormal, idx_inverted)
            axcurr.set_xlim(phys_lim[axcol,0], phys_lim[axcol,1])
            axcurr.set_ylim(phys_lim[axrow,0], phys_lim[axrow,1])
            axcurr.grid(True)
    # set spacing between subplots
    fig.set_tight_layout({'pad': 0.5})
    # save & show plots
    for suffix in PLOT_SUFFIX:
        fig.savefig('{}_marginal_2d_all.{}'.format(plot_base, suffix), dpi=360)
    #plt.show()

# create figure with all 2D conditionals of the posterior
if PLOT_CONDITIONALS_ALL:
    fig, ax = plt.subplots(n_glo+n_weak, n_glo+n_weak, figsize=(26, 22))
    for axcol in range(0, n_glo+n_weak):
        if axcol < n_glo:
            ax[0,axcol].set_title(PLOT_TITLES_GLO[axcol])
        else:
            ax[0,axcol].set_title(PLOT_TITLE_WEAK+(" %i" % (axcol-n_glo)))
        for axrow in range(0, n_glo+n_weak):
            axcurr = ax[axrow,axcol]
            if axcol == 0:
                if axrow < n_glo:
                    axcurr.set_ylabel(PLOT_TITLES_GLO[axrow])
                else:
                    axcurr.set_ylabel(PLOT_TITLE_WEAK+(" %i" % (axrow-n_glo)))
            if axrow == axcol:
                axcurr.get_xaxis().set_ticks([])
                axcurr.get_yaxis().set_ticks([])
                continue
            plot_conditional_2d(axcurr, [axcol,axrow], m,
                                prior_max, prior_cov, COLOR_PRIOR,
                                post_max, post_cov, COLOR_POST_COND,
                                param_dim, phys_dim, idx_lognormal, idx_inverted)
           #axcurr.set_xlim(phys_lim[axcol,0], phys_lim[axcol,1])
           #axcurr.set_ylim(phys_lim[axrow,0], phys_lim[axrow,1])
            axcurr.grid(True)
    # set spacing between subplots
    fig.set_tight_layout({'pad': 0.5})
    # save & show plots
    for suffix in PLOT_SUFFIX:
        fig.savefig('{}_conditional_2d_all.{}'.format(plot_base, suffix), dpi=360)
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
