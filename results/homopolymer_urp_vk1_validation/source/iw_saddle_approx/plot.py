#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
sys.path.insert(0, '/Users/lyx/Develop/cl_ma_rho/scripts')

import numpy as np
from scipy.io import loadmat

import matplotlib as mpl
# To avoid launching interactive plot, such as wxAgg.
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import mpltex

from DFM import get_phi, convert_3d_to_1d
from DFM import parse_input
from DFM import multiline_plot
from DFM import one_basis_fit, get_basis_set
from DFM import get_iw_K1, get_iw_rpa, get_iw_fga, get_iw_sga, get_iw_gsd
from DFM import fit_error


def multiline_plot(ax, x, data, xlabel, ylabel, labels=None,
                   loc='best', linestyles=[], ms=4, lw=1,
                   xscale=None, yscale=None, xlim=None, ylim=None):
    linestyle1 = mpltex.linestyle_generator(lines=['-'], markers=[],
                                            hollow_styles=[])
    linestyle2 = mpltex.linestyle_generator(lines=[], markers=['o'],
                                            hollow_styles=[])

    n = len(data)  # number of curves
    if type(labels) is str:
        labels = [labels] * n
    has_label = True
    if not labels:
        has_label = False
        labels = [''] * n
    if type(linestyles) is int:
        linestyles = [linestyles] * n
    elif not linestyles:
        linestyles = [1] * n

    i = 0
    for d in data:
        if linestyles[i] == 2:
            ax.plot(x[i], d, label=labels[i], ms=ms, **linestyle2.next())
        else:
            line, = ax.plot(x[i], d, label=labels[i], linewidth=lw, ms=1.5,
                            **linestyle1.next())
            if labels[i] == 'RPA':
                line.set_dashes([3, 1.2])  # dash
                # plt.setp(line, linestyle='--')
            elif labels[i] == 'SGA':
                line.set_dashes([1, 1.2])  # dot
                # plt.setp(line, linestyle=':')
            elif labels[i] == 'FGA':
                line.set_dashes([3, 1.2, 1, 1.2])  # dash-dot
                # plt.setp(line, linestyle='-.')
            elif labels[i] == 'GSD':
                line.set_dashes([3, 1.2, 1, 1.2, 1, 1.2])  # dash-dot-dot
            elif labels[i] == 'URP':
                line.set_linewidth(0.5)
            elif labels[i] == 'VK1':
                line.set_linewidth(1.0)
            elif labels[i] == 'MAP':
                pass
        i += 1

    ax.locator_params(nbins=5)
    ax.xaxis.labelpad = 1.5
    ax.yaxis.labelpad = 1.5
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if n > 1 and has_label:
        ax.legend(loc=loc, fontsize=6.5)

    return ax


@mpltex.aps_decorator
def fit_saddle(data_dir='./'):
    show_scft = True
    show_rpa2 = True
    show_sga = True
    show_fga = True
    show_gsd = True
    show_rpa1 = True
    show_K1 = True
    show_match = True
    show_K2 = False
    show_K3 = False

    data_file = 'phitype5_L3_a0.5_b1_c1'
    inputfile = os.path.join(data_dir, 'input_phitype5.txt')
    # data_file = 'phitype11_L3.5_a0.4_b1'
    # inputfile = os.path.join(data_dir, 'input_phitype11.txt')
    params = parse_input(inputfile)
    phi_type = params['phi_type']
    Nx = params['Nx']
    Ny = params['Ny']
    Nz = params['Nz']
    Lx = params['Lx']
    a = params['a']
    b = params['b']
    c = params['c']
    d = params['d']

    dim = 3
    if Ny == 1:
        dim = dim - 1
    if Nz == 1:
        dim = dim - 1

    x, phi = get_phi(phi_type, [Nx, Lx, a, b, c, d])

    matfile = os.path.join(data_dir, data_file)
    mat = loadmat(matfile)
    if not np.allclose(mat['x'][0, :], x):
        print 'x array in input file does not match that in the result file!'
        exit(1)
    iw = mat['w_snap']
    if dim == 3:
        iw = convert_3d_to_1d(iw)
    elif dim == 1:
        iw = iw[0, :]
    mu = mat['mu'][0, -1]

    x_array = []
    y_array = []
    labels = []
    linestyles = []
    if show_scft:
        x_array.append(x)
        y_array.append(iw-mu)
        labels.append('SCFT')
        linestyles.append(2)

    # Approximation schemes
    Nx_approx = 4 * Nx
    x_approx, phi_approx = get_phi(phi_type, [Nx_approx, Lx*1.0, a, b, c, d])
    # RPA
    iw_approx = get_iw_rpa(phi_approx, Lx, t=2)
    iw_approx2 = get_iw_rpa(phi, Lx, t=2)
    print 'RPA error:', fit_error(iw_approx2, iw-mu)
    if show_rpa2:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('RPA')
        linestyles.append(1)
    # SGA
    iw_approx = get_iw_sga(phi_approx, Lx)
    iw_approx2 = get_iw_sga(phi, Lx)
    print 'SGA error:', fit_error(iw_approx2, iw-mu)
    if show_sga:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('SGA')
        linestyles.append(1)
    # FGA
    iw_approx = get_iw_fga(phi_approx, Lx)
    iw_approx2 = get_iw_fga(phi, Lx)
    print 'FGA error:', fit_error(iw_approx2, iw-mu)
    if show_fga:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('FGA')
        linestyles.append(1)
    # GSD
    iw_approx = get_iw_gsd(phi_approx, Lx)
    iw_approx2 = get_iw_gsd(phi, Lx)
    print 'GSD error:', fit_error(iw_approx2, iw-mu)
    if show_gsd:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('GSD')
        linestyles.append(1)
    # URPA
    iw_approx = get_iw_rpa(phi_approx, Lx, t=1)
    iw_approx2 = get_iw_rpa(phi, Lx, t=1)
    print 'URPA error:', fit_error(iw_approx2, iw-mu)
    if show_rpa1:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('URP')
        linestyles.append(1)
    # K1 - 1
    iw_approx = get_iw_K1(phi_approx, Lx, t=1)
    iw_approx2 = get_iw_K1(phi, Lx, t=1)
    print 'K1-1 error:', fit_error(iw_approx2, iw-mu)
    if show_K1:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('VK1')
        linestyles.append(1)
    # K1 - 2
    iw_approx = get_iw_K1(phi_approx, Lx, t=2)
    iw_approx2 = get_iw_K1(phi, Lx, t=2)
    print 'K1-2 error:', fit_error(iw_approx2, iw-mu)
    if show_K2:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('K1-2')
        linestyles.append(1)
    # K1 - 3
    iw_approx = get_iw_K1(phi_approx, Lx, t=3)
    iw_approx2 = get_iw_K1(phi, Lx, t=3)
    print 'K1-3 error:', fit_error(iw_approx2, iw-mu)
    if show_K3:
        x_array.append(x_approx)
        y_array.append(iw_approx)
        labels.append('K1-3')
        linestyles.append(1)

    # 1-basis fit
    basis_set = get_basis_set('Tang-Freed1', phi, Lx)
    basis1 = basis_set[0]
    K1 = one_basis_fit(basis1, iw - mu + np.log(phi))
    iw_fit = K1 * basis1 - np.log(phi)
    print 'fit K:', K1
    print 'fit error:', fit_error(iw_fit, iw-mu)
    if show_match:
        x_array.append(x)
        y_array.append(iw_fit)
        labels.append('MAP')
        linestyles.append(1)

    fig, ax = plt.subplots()
    multiline_plot(ax, x_array, y_array, '$x$', '$iw^*$',
                   linestyles=linestyles, ms=2, lw=0.75,
                   labels=labels, loc='upper right',
                   xlim=[0, 3])
    minorLocator_x = AutoMinorLocator(5)
    minorLocator_y = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorLocator_x)
    ax.yaxis.set_minor_locator(minorLocator_y)
    ax.text(0.05, 0.9, '$k=2\pi/3$', transform=ax.transAxes)
    ax.text(0.05, 0.8, '$a=0.5$', transform=ax.transAxes)
    ax.text(0.05, 0.7, '$\sigma=1$', transform=ax.transAxes)

    # inset plots for prescribed density
    ax_in = plt.axes([.42, .25, .3, .25], axisbg='w')
    multiline_plot(ax_in, [x], [phi], None, None,
                   linestyles=1, lw=1,
                   xlim=[0, 3])
    # ax_in.plot(x, iw_fit)
    # ax_in.set_ylabel('$\phi$')

    fig.tight_layout(pad=0.15)
    fig.savefig('iw_saddle_approx')


if __name__ == '__main__':
    fit_saddle()
