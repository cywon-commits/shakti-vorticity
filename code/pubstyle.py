"""Common style for the PRB figures: final print widths, readable fonts, compact colorbars, panel labels."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

COL = 3.375      # PRB single column (inches)
TWO = 7.0        # PRB double column
FIGDIR = '../paper/figs/'

plt.rcParams.update({
    'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'legend.frameon': False, 'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5, 'lines.linewidth': 1.0, 'lines.markersize': 3.5,
    'savefig.dpi': 400, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02, 'pdf.fonttype': 42,
    'axes.grid': False, 'figure.dpi': 150,
})

# vorticity colours (diverging: blue -1, grey 0, orange +1), unchanged from the working figures
OR, BL, ORl, BLl = '#e8743a', '#2f6db5', '#f6c6a8', '#aac6ea'
CMAPV = LinearSegmentedColormap.from_list('ob', [BL, '#e6e6e6', OR])
# categorical order (validated: CVD-safe with markers as secondary encoding)
CAT = ['#1b6ca8', '#d9822b', '#7b3294', '#2e8b57']
MK = ['o', 's', '^', 'D']
SEQ = {24: '#c6dbef', 32: '#6baed6', 48: '#2171b5', 64: '#08306b', 96: '#2171b5', 128: '#08306b'}

def label(ax, s, x=-0.02, y=1.0, **kw):
    """panel label (a), (b) ... in the upper-left corner outside the axes"""
    ax.text(x, y, s, transform=ax.transAxes, ha='right', va='bottom', fontsize=9, fontweight='bold', **kw)

def label_in(ax, s, x=0.02, y=0.98, **kw):
    ax.text(x, y, s, transform=ax.transAxes, ha='left', va='top', fontsize=9, fontweight='bold', **kw)

def small_cbar(fig, mappable, ax, label=None, loc='right', size=('4%', '45%'), pad=0.02, ticks=None, orientation='vertical'):
    """compact colour bar attached next to an axes"""
    if orientation == 'vertical':
        cax = ax.inset_axes([1 + pad, 0.5 - float(size[1][:-1]) / 200, float(size[0][:-1]) / 100, float(size[1][:-1]) / 100])
    else:
        w = float(size[1][:-1]) / 100; h = float(size[0][:-1]) / 100
        cax = ax.inset_axes([0.5 - w / 2, -pad - h, w, h])
    cb = fig.colorbar(mappable, cax=cax, orientation=orientation, ticks=ticks)
    cb.ax.tick_params(labelsize=6.5, width=0.5, length=2)
    cb.outline.set_linewidth(0.5)
    if label: cb.set_label(label, fontsize=7, labelpad=2)
    return cb

def save(fig, name):
    os.makedirs(FIGDIR, exist_ok=True); os.makedirs('figcheck', exist_ok=True)
    fig.savefig(FIGDIR + name + '.pdf')
    fig.savefig('figcheck/' + name + '.png', dpi=200)
    plt.close(fig)
