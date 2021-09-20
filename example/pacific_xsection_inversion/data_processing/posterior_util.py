import numpy as np

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

def _transfer_physical_to_param_lin_deriv(p, param_dim, phys_glo_dim):
    return np.ones_like(p)/phys_glo_dim/param_dim

def _transfer_physical_to_param_log_deriv(p, param_dim, phys_glo_dim=None):
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
#       dmi = np.abs(_transfer_physical_to_param_log_deriv(p[i,:], dim, phys_dim))
#   else:
#       dmi = np.abs(_transfer_physical_to_param_lin_deriv(p[i,:], dim, phys_dim))
    # set derivatives of transfer function for index `j`
#   dim = param_dim[j]
#   if j < len(phys_glo_dim):
#       phys_dim = phys_glo_dim[j]
#   else:
#       phys_dim = 1.0
#   if j in idx_lognormal:
#       dmj = np.abs(_transfer_physical_to_param_log_deriv(p[j,:], dim, phys_dim))
#       print(j,"is log")
#   else:
#       dmj = np.abs(_transfer_physical_to_param_lin_deriv(p[j,:], dim, phys_dim))
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

########################################
# 1D Plots
########################################

def plot_pdfs_1d(ax, x,
                 pdf0, x0_max, pdf0_max, pdf0_color,
                 pdf1, x1_max, pdf1_max, pdf1_color):
    assert ax is not None
    assert x is not None
    # plot pdf 0 (prior)
    if not (pdf0 is None or x0_max is None or pdf0_max is None or pdf0_color is None):
        # plot max point
        xx = np.repeat(x0_max, 2)
        yy = [0.0, pdf0_max]
        ax.plot(xx, yy, color=pdf0_color, linewidth=1, linestyle='--')
        # plot density
        ax.plot(x, pdf0, label='prior', color=pdf0_color, linewidth=2)
    # plot pdf 1 (posterior)
    if not (pdf1 is None or x1_max is None or pdf1_max is None or pdf1_color is None):
        # plot max point
        xx = np.repeat(x1_max, 2)
        yy = [0.0, pdf1_max]
        ax.plot(xx, yy, color=pdf1_color, linewidth=1, linestyle='--')
        # plot density
        ax.plot(x, pdf1, label='posterior', color=pdf1_color, linewidth=2)

########################################
# 2D Plots
########################################

def plot_add_confidence_ellipse(ax, mean, cov, n_std=3.0,
                                param_dim=(1.0, 1.0), phys_dim=(1.0, 1.0),
                                lognormal=(False, False),
                                inverted=(False, False),
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
            ellipse[0,:], param_dim[0], phys_dim[0], lognormal[0])
    if inverted[0]:
        ellipse[0,:] = 1.0/np.abs(ellipse[0,:])
    ellipse[1,:] = transfer_param_to_physical(
            ellipse[1,:], param_dim[1], phys_dim[1], lognormal[1])
    if inverted[1]:
        ellipse[1,:] = 1.0/np.abs(ellipse[1,:])
    # plot ellipse
    return ax.plot(ellipse[0,:], ellipse[1,:], color=edgecolor, **kwargs)

def cov_marginal_2d(indices, cov):
    i = indices[0]
    j = indices[1]
    return np.array([ [cov[i,i], cov[i,j]],
                      [cov[j,i], cov[j,j]] ])

def plot_marginal_2d(ax, indices, m,
                     prior_max, prior_cov, prior_color,
                     post_max, post_cov, post_color,
                     param_dim, phys_dim, idx_lognormal, idx_inverted):
    assert ax is not None
    assert indices is not None
    assert m is not None
    i = indices[0]
    j = indices[1]
    param_dim = [param_dim[i], param_dim[j]]
    phys_dim  = [phys_dim[i], phys_dim[j]]
    lognormal = [i in idx_lognormal, j in idx_lognormal]
    inverted  = [i in idx_inverted, j in idx_inverted]
    # plot prior marginal
    #   (1*sigma, 2*sigma, 3*sigma) = (0.682, 0.954, 0.997) confidence
    if not (prior_max is None or prior_cov is None or prior_color is None):
        shift = np.array([prior_max[i], prior_max[j]])
        cov   = cov_marginal_2d(indices, prior_cov)
        for n_stddev in [1.0, 2.0, 3.0]:
            plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
                                        param_dim, phys_dim, lognormal, inverted,
                                        edgecolor=prior_color, linewidth=2)
        # plot prior max
        p_prior_max = np.array([
            transfer_param_to_physical(prior_max[i], param_dim[0],
                                       phys_dim[0], lognormal[0]),
            transfer_param_to_physical(prior_max[j], param_dim[1],
                                       phys_dim[1], lognormal[1])
        ])
        p_prior_max[np.where(inverted)] = 1.0/p_prior_max[np.where(inverted)]
        ax.plot(p_prior_max[0], p_prior_max[1], color=prior_color, marker='.')
    # plot posterior marginal
    #   (1*sigma, 2*sigma, 3*sigma) = (0.682, 0.954, 0.997) confidence
    if not (post_max is None or post_cov is None or post_color is None):
        shift = np.array([post_max[i], post_max[j]])
        cov   = cov_marginal_2d(indices, post_cov)
        for n_stddev in [1.0, 2.0, 3.0]:
            plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
                                        param_dim, phys_dim, lognormal, inverted,
                                        edgecolor=post_color, linewidth=2)
        # plot posterior max
        p_post_max = np.array([
            transfer_param_to_physical(post_max[i], param_dim[0],
                                       phys_dim[0], lognormal[0]),
            transfer_param_to_physical(post_max[j], param_dim[1],
                                       phys_dim[1], lognormal[1])
        ])
        p_post_max[np.where(inverted)] = 1.0/p_post_max[np.where(inverted)]
        ax.plot(p_post_max[0], p_post_max[1], color=post_color, marker='.')
    # set axis to logscale
    if lognormal[0]:
        ax.set_xscale('log')
    if lognormal[1]:
        ax.set_yscale('log')

def cov_conditional_2d(indices, cov):
    i = indices[0]
    j = indices[1]
    cov_d1 = np.array([ [cov[i,i], cov[i,j]],
                        [cov[j,i], cov[j,j]] ])
    cov_u  = np.delete(np.array([cov[i,:], cov[j,:]]), (i,j), axis=1)
    cov_d2 = np.delete(np.delete(cov, (i,j), axis=0), (i,j), axis=1)
    return cov_d1 - np.dot(np.dot(cov_u, np.linalg.inv(cov_d2)), cov_u.T)

def plot_conditional_2d(ax, indices, m,
                        prior_max, prior_cov, prior_color,
                        post_max, post_cov, post_color,
                        param_dim, phys_dim, idx_lognormal, idx_inverted):
    i = indices[0]
    j = indices[1]
    param_dim = [param_dim[i], param_dim[j]]
    phys_dim  = [phys_dim[i], phys_dim[j]]
    lognormal = [i in idx_lognormal, j in idx_lognormal]
    inverted  = [i in idx_inverted, j in idx_inverted]
#   # plot posterior marginal
#   shift = np.array([post_max[i], post_max[j]])
#   cov   = cov_marginal_2d(indices, post_cov)
#  #for n_stddev in [1.0, 2.0, 3.0]:
#  #    plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
#  #                                param_dim, phys_dim, lognormal, inverted,
#  #                                edgecolor=COLOR_POST, linewidth=1)
#   # plot prior max
#   p_post_max = np.array([
#       transfer_param_to_physical(post_max[i], param_dim[0],
#                                  phys_dim[0], lognormal[0]),
#       transfer_param_to_physical(post_max[j], param_dim[1],
#                                  phys_dim[1], lognormal[1])
#   ])
#   p_post_max[np.where(inverted)] = 1.0/p_post_max[np.where(inverted)]
#   ax.plot(p_post_max[0], p_post_max[1], color=COLOR_POST, marker='.')
    # plot posterior conditional
    #   (1*sigma, 2*sigma, 3*sigma) = (0.682, 0.954, 0.997) confidence
    shift = np.array([post_max[i], post_max[j]])
    cov   = cov_conditional_2d(indices, post_cov)
    for n_stddev in [1.0, 2.0, 3.0]:
        plot_add_confidence_ellipse(ax, shift, cov, n_stddev,
                                    param_dim, phys_dim, lognormal, inverted,
                                    edgecolor=post_color, linewidth=1)
    # plot posterior max
    p_post_max = np.array([
        transfer_param_to_physical(post_max[i], param_dim[0],
                                   phys_dim[0], lognormal[0]),
        transfer_param_to_physical(post_max[j], param_dim[1],
                                   phys_dim[1], lognormal[1])
    ])
    p_post_max[np.where(inverted)] = 1.0/p_post_max[np.where(inverted)]
    ax.plot(p_post_max[0], p_post_max[1], color=post_color, marker='.')
    # set axis to logscale
    if lognormal[0]:
        ax.set_xscale('log')
    if lognormal[1]:
        ax.set_yscale('log')

#TODO deprecated code
#
#def generate_plot_data_2d(indices, m, prior_mean, prior_cov, post_mean, post_cov,
#                          p, param_dim, phys_dim, idx_lognormal):
#    i = indices[0]
#    j = indices[1]
#    Pi, Pj = np.meshgrid(p[i,:], p[j,:])
#    Pdf_prior = pdf_2d([i,j], m, prior_mean, prior_cov,
#                       p, param_dim, phys_dim, idx_lognormal)
#    Pdf_post  = pdf_2d([i,j], m, post_mean, post_cov,
#                       p, param_dim, phys_dim, idx_lognormal)
#    return (Pi, Pj, Pdf_prior, Pdf_post)
#
#def plot_pdfs_2d(ax, indices, m,
#                 prior_mean, prior_cov, p_prior_max,
#                 post_mean, post_cov, p_post_max,
#                 p, param_dim, phys_dim, idx_lognormal):
#    i = indices[0]
#    j = indices[1]
#    # generate data
#    (Pi, Pj, Pdf_prior, Pdf_post) = generate_plot_data_2d(
#            indices, m, prior_mean, prior_cov, post_mean, post_cov,
#            p, param_dim, phys_dim, idx_lognormal)
#    # plot prior
#    levels = np.array([0.3, 0.6, 0.9]) * np.max(Pdf_prior)
#    CS = ax.contour(Pi, Pj, Pdf_prior, levels, colors=COLOR_PRIOR, linewidths=2)
#    ax.plot(p_prior_max[i], p_prior_max[j], color=COLOR_PRIOR, marker='o')
#    # plot posterior
#    levels = np.array([0.3, 0.6, 0.9]) * np.max(Pdf_post)
#    CS = ax.contour(Pi, Pj, Pdf_post, levels, colors=COLOR_POST, linewidths=2)
#    ax.plot(p_post_max[i], p_post_max[j], color=COLOR_POST, marker='o')

