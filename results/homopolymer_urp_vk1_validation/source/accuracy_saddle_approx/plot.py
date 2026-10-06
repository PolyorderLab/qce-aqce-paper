#!/Users/lyx/Library/Enthought/Canopy_64bit/User/bin/python
# -*- coding: utf-8 -*-
import numpy as np

import matplotlib as mpl
# To avoid launching interactive plot, such as wxAgg.
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import mpltex


def multiline_plot(ax, x, data, xlabel, ylabel, labels=None,
                   loc='best', linestyles=[], ms=4, lw=1,
                   xscale=None, yscale=None, xlim=None, ylim=None):
    linestyle1 = mpltex.linestyle_generator(lines=['-'], markers=[],
                                            hollow_styles=[])
    linestyle2 = mpltex.linestyle_generator(lines=[], markers=['o'],
                                            hollow_styles=[])
    linestyle3 = mpltex.linestyle_generator(lines=[],
                                            markers=['v', '^', '<', '>', 'o', 's', 'D'],
                                            hollow_styles=[True])
    linestyle4 = mpltex.linestyle_generator(lines=['-'], markers=['o'],
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
            # linestyle1.next()
            # linestyle3.next()
            # linestyle4.next()
        elif linestyles[i] == 3:
            extra_styles = linestyle3.next()
            extra_styles['mew'] = 0.4
            ax.plot(x[i], d, label=labels[i], linewidth=lw, ms=ms,
                    **extra_styles)
            # linestyle1.next()
            # linestyle2.next()
            # linestyle4.next()
        elif linestyles[i] == 4:
            ax.plot(x[i], d, label=labels[i], linewidth=lw, ms=ms,
                    **linestyle4.next())
            # linestyle1.next()
            # linestyle2.next()
            # linestyle3.next()
        else:
            ax.plot(x[i], d, label=labels[i], linewidth=lw, ms=ms,
                    **linestyle1.next())
            # linestyle2.next()
            # linestyle3.next()
            # linestyle4.next()
        i += 1

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
        ax.legend(loc=loc)

    return ax


def plot_ax1(ax):
    datafile = './error_K_L_a0.01.dat'
    rawdata = np.loadtxt(datafile)
    L = rawdata[:, 0]
    err_rpa = rawdata[:, 1]
    err_sga = rawdata[:, 2]
    err_fga = rawdata[:, 3]
    err_gsd = rawdata[:, 4]
    err_urpa = rawdata[:, 5]
    err_k1 = rawdata[:, 6]
    err_k2 = rawdata[:, 7]
    err_k3 = rawdata[:, 8]
    err_fit = rawdata[:, 9]
    K1_fit = rawdata[:, 10]
    y_array = [err_rpa, err_sga, err_fga, err_gsd,
               err_urpa, err_k1,
               err_fit]
    y_array = [np.log10(y) for y in y_array]
    k = 2 * np.pi / L
    x_array = [1./k] * len(y_array)

    multiline_plot(ax, x_array, y_array, None, '$\log\\varepsilon$',
                   xlim=[0, 1.4],
                   linestyles=3, ms=2.5, lw=0.5)
    ax.locator_params(nbins=5)
    minorLocator_x = AutoMinorLocator(5)
    minorLocator_y = AutoMinorLocator(5)
    # minorLocator_x = MultipleLocator(5)
    # minorLocator_y = MultipleLocator(5)
    ax.xaxis.set_minor_locator(minorLocator_x)
    ax.yaxis.set_minor_locator(minorLocator_y)
    ax.text(0.15, 0.1, '$a=0.01$', transform=ax.transAxes)
    ax.text(0.85, 0.1, '(a)', transform=ax.transAxes)
    ax.yaxis.labelpad = 1.5


def plot_ax2(ax):
    datafile = './error_K_a_L8.dat'
    rawdata = np.loadtxt(datafile)
    a = rawdata[:, 0]
    err_rpa = rawdata[:, 1]
    err_sga = rawdata[:, 2]
    err_fga = rawdata[:, 3]
    err_gsd = rawdata[:, 4]
    err_urpa = rawdata[:, 5]
    err_k1 = rawdata[:, 6]
    err_k2 = rawdata[:, 7]
    err_k3 = rawdata[:, 8]
    err_fit = rawdata[:, 9]
    K1_fit = rawdata[:, 10]
    y_array = [err_rpa, err_sga, err_fga, err_gsd,
               err_urpa, err_k1,
               err_fit]
    y_array = [np.log10(y) for y in y_array]
    x_array = [a] * len(y_array)

    multiline_plot(ax, x_array, y_array, None, None,
                   xlim=[0, 0.95],
                   linestyles=3, ms=2.5, lw=0.5)
    ax.text(0.85, 0.1, '(c)', transform=ax.transAxes)
    ax.text(0.15, 0.1, '$k=\pi/4$', transform=ax.transAxes)


def plot_ax3(ax):
    datafile = './error_K_L_a0.8.dat'
    rawdata = np.loadtxt(datafile)
    L = rawdata[:, 0]
    err_rpa = rawdata[:, 1]
    err_sga = rawdata[:, 2]
    err_fga = rawdata[:, 3]
    err_gsd = rawdata[:, 4]
    err_urpa = rawdata[:, 5]
    err_k1 = rawdata[:, 6]
    err_k2 = rawdata[:, 7]
    err_k3 = rawdata[:, 8]
    err_fit = rawdata[:, 9]
    K1_fit = rawdata[:, 10]
    y_array = [err_rpa, err_sga, err_fga, err_gsd,
               err_urpa, err_k1,
               err_fit]
    y_array = [np.log10(y) for y in y_array]
    k = 2 * np.pi / L
    x_array = [1./k] * len(y_array)

    multiline_plot(ax, x_array, y_array, '$k^{-1}$', '$\log\\varepsilon$',
                   linestyles=3, ms=2.5, lw=0.5)
    ax.text(0.85, 0.1, '(b)', transform=ax.transAxes)
    ax.text(0.15, 0.1, '$a=0.8$', transform=ax.transAxes)
    ax.xaxis.labelpad = 1.5
    ax.yaxis.labelpad = 1.5


def plot_ax4(ax):
    datafile = './error_K_a_L0.5.dat'
    rawdata = np.loadtxt(datafile)
    a = rawdata[:, 0]
    err_rpa = rawdata[:, 1]
    err_sga = rawdata[:, 2]
    err_fga = rawdata[:, 3]
    err_gsd = rawdata[:, 4]
    err_urpa = rawdata[:, 5]
    err_k1 = rawdata[:, 6]
    err_k2 = rawdata[:, 7]
    err_k3 = rawdata[:, 8]
    err_fit = rawdata[:, 9]
    K1_fit = rawdata[:, 10]
    y_array = [err_rpa, err_sga, err_fga, err_gsd,
               err_urpa, err_k1,
               err_fit]
    y_array = [np.log10(y) for y in y_array]
    x_array = [a] * len(y_array)

    multiline_plot(ax, x_array, y_array, '$a$', None,
                   linestyles=3, ms=2.5, lw=0.5)
    ax.locator_params(nbins=5)
    minorLocator_x = AutoMinorLocator(5)
    minorLocator_y = AutoMinorLocator(5)
    # minorLocator_x = MultipleLocator(5)
    # minorLocator_y = MultipleLocator(5)
    ax.xaxis.set_minor_locator(minorLocator_x)
    ax.yaxis.set_minor_locator(minorLocator_y)
    ax.text(0.15, 0.1, '$k=4\pi$', transform=ax.transAxes)
    ax.text(0.85, 0.1, '(d)', transform=ax.transAxes)
    ax.xaxis.labelpad = 1.5


@mpltex.aps_decorator
def main():
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,
                                                 figsize=(3.375, 2.54),
                                                 sharex='col',
                                                 sharey='row')
    plot_ax1(ax1)
    plot_ax2(ax2)
    plot_ax3(ax3)
    plot_ax4(ax4)
    fig.subplots_adjust(wspace=0, hspace=0)
    # fig.set_tight_layout(True)
    fig.tight_layout(pad=0.35, w_pad=0, h_pad=0)
    fig.savefig('accuracy_saddle_approx')


if __name__ == '__main__':
    main()
