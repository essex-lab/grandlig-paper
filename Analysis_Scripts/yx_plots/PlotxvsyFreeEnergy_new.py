# Modules to import
import argparse
import scipy.stats
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error
import numpy as np
from simtk.unit import *
from scipy.optimize import curve_fit, OptimizeWarning

# Set font sizes for the graph
small_font = 16
medium_font = 18
large_font = 20


plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=small_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=14)


# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-x", "--x", help="Input x data")
parser.add_argument("-y", "--y", help="Input y data")
parser.add_argument("-x_label", "--x_label", help="Input x label")
parser.add_argument("-y_label", "--y_label", help="Input y label")
parser.add_argument("-lims", "--lims", help="Input Insertion Work files", type=float, nargs='+', default=None)
parser.add_argument("-o", "--output", help="Input Output File Name")
parser.add_argument("--legend", help="Legend Placement", default=None, type=str)

args = parser.parse_args()

def func(x, m, c):
    y = (m*x) + c
    return y

def calc_stats(x, y):
    rmse = np.sqrt(mean_squared_error(x, y))
    mae = mean_absolute_error(x, y)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
        x, y
    )
    r_squared = r_value **2
    pearsons = scipy.stats.pearsonr(x, y)[0]
    kendall = scipy.stats.kendalltau(x, y)[0]
    spearman = scipy.stats.spearmanr(x, y)[0]

    return rmse, mae, slope, intercept, r_squared, pearsons, kendall, spearman


def boot_stats(x, y, x_err, y_err, ci=0.95, nboot=1000):
    s_rmse = np.zeros([nboot])
    s_mae = np.zeros([nboot])
    s_r_squared = np.zeros([nboot])
    s_pearsons = np.zeros([nboot])
    s_kendall = np.zeros([nboot])
    s_spearman = np.zeros([nboot])

    sample_size = len(y)
    for r in range(nboot):
        x_sample = np.zeros_like(x)
        y_sample = np.zeros_like(y)
        for i, j in enumerate(
            np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)
        ):
            stddev_true = np.fabs(x_err[j])
            stddev_pred = np.fabs(y_err[j])
            #stddev_true = 0
            #stddev_pred = 0

            x_sample[i] = np.random.normal(loc=x[j], scale=stddev_true, size=1)
            y_sample[i] = np.random.normal(loc=y[j], scale=stddev_pred, size=1)
        s_rmse[r], s_mae[r], _, _, s_r_squared[r], s_pearsons[r], s_kendall[r], s_spearman[r] = calc_stats(x_sample, y_sample)

    low_frac = (1.0 - ci) / 2.0
    high_frac = 1.0 - low_frac

    s_rmse = np.sort(s_rmse)
    s_mae = np.sort(s_mae)
    s_r_squared = np.sort(s_r_squared)
    s_pearsons = np.sort(s_pearsons)
    s_kendall = np.sort(s_kendall)
    s_spearman = np.sort(s_spearman)

    rmse_ci = [s_rmse[int(np.floor(nboot*low_frac))], s_rmse[int(np.ceil(nboot*high_frac))]]
    mae_ci = [s_mae[int(np.floor(nboot*low_frac))], s_mae[int(np.ceil(nboot*high_frac))]]
    r_squared_ci = [s_r_squared[int(np.floor(nboot*low_frac))], s_r_squared[int(np.ceil(nboot*high_frac))]]
    pearsons_ci = [s_pearsons[int(np.floor(nboot*low_frac))], s_pearsons[int(np.ceil(nboot*high_frac))]]
    kendall_ci = [s_kendall[int(np.floor(nboot*low_frac))], s_kendall[int(np.ceil(nboot*high_frac))]]
    spearman_ci = [s_spearman[int(np.floor(nboot*low_frac))], s_spearman[int(np.ceil(nboot*high_frac))]]
    return rmse_ci, mae_ci, r_squared_ci, pearsons_ci, kendall_ci, spearman_ci


frags = {}

with open(args.x, 'r') as fi:
    for line in fi.readlines():
        f = str(line.split()[2])
        dg = float(line.split()[0])
        dg_err = float(line.split()[1])
        frags[f] = [dg, dg_err]

with open(args.y, 'r') as fi:
    for line in fi.readlines():
        f = str(line.split()[2])
        dg = float(line.split()[0]) # - 0.68  # STD state coprrect for host guest testing
        dg_err = float(line.split()[1])
        frags[f].append(dg)
        frags[f].append(dg_err)

# frags[f] = [dg_x, err, dg_y, err]

data = np.asarray(list(frags.values()))
frags = list(frags.keys())

x_dg = data[:, 0]
x_error = data[:, 1]
y_dg = data[:, 2]
y_error = data[:, 3]


for i in range(len(x_dg)):
    print('x - y = ', x_dg[i] - y_dg[i], frags[i])

fig1 = plt.figure(figsize=(8.5, 8.5)) #, facecolor='#EDEDED')
ax2 = fig1.add_subplot(111)

# Getting all the stats
x_data_stats = x_dg
y_data_stats = y_dg

rmse, mae, slope, intercept, r_squared, pearsons, kendall, spearman = calc_stats(x_data_stats, y_data_stats)

rmse_ci, mae_ci, r_squared_ci, pearsons_ci, kendall_ci, spearman_ci = boot_stats(x_dg, y_dg, x_error, y_error)

# Annotate with statistics including ranges
stats_text = (
    f"{'R squared = ' : >12}{r_squared:.3f} (95%: {r_squared_ci[0]:.3f}, {r_squared_ci[1]:.3f})\n"
    f"{'RMSE = ' : >12}{rmse:.3f} (95%: {rmse_ci[0]:.3f}, {rmse_ci[1]:.3f})\n"
    f"{'MAE = ' : >12}{mae:.3f} (95%: {mae_ci[0]:.3f}, {mae_ci[1]:.3f})\n"
    f"{'Kendall Tau = ' : >12}{kendall:.3f} (95%: {kendall_ci[0]:.3f}, {kendall_ci[1]:.3f})\n"
    f"{'Pearson r = ' : >12}{pearsons:.3f} (95%: {pearsons_ci[0]:.3f}, {pearsons_ci[1]:.3f})\n"
    f"{'Spearman ρ = ' : >12}{spearman:.3f} (95%: {spearman_ci[0]:.3f}, {spearman_ci[1]:.3f})"
)


# Curve fit
ax_min = min((min(x_dg), min(y_dg))) - 1
ax_max = max((max(x_dg), max(y_dg))) + 1

if args.lims:
    lims = [args.lims[0], args.lims[1]]
else:
    lims = [ax_min, ax_max]

print(f"LIMITS OF THE PLOT: {lims[0]} {lims[1]}", lims)

x_fit = np.linspace(lims[0], lims[1], 5)
y_fit = (x_fit * slope) + intercept


ax2.scatter(x_dg, y_dg, marker='o', s=75, linestyle='None', c='r', zorder=2)
ax2.errorbar(x_dg, y_dg, xerr=x_error, yerr=y_error, fmt='none', c='black', zorder=1, capsize=2)
ax2.plot(x_fit, y_fit, c='black', lw=2, zorder=1)

ax2.set_xlim(lims)
ax2.set_ylim(lims)

# y=x line and shading

ax2.plot(lims, lims, linestyle='--', c='.3')

ax2.fill_between(
    lims,
    [lims[0] - 1, lims[1] - 1],
    [lims[0] + 1, lims[1] + 1],
    color="#92A9CA",
    alpha=0.2,
)

ax2.fill_between(
    lims,
    [lims[0] - 2, lims[1] - 2],
    [lims[0] + 2, lims[1] + 2],
    color="#92A9CA",
    alpha=0.2,
)

# Add the text at the bottom right
if args.legend == 'upper left':
    plt.text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
         verticalalignment='top', horizontalalignment='left',
         fontsize=14, bbox=dict(facecolor='white', alpha=0.8))
else:
    plt.text(0.95, 0.05, stats_text, transform=plt.gca().transAxes,
         verticalalignment='bottom', horizontalalignment='right',
         fontsize=14, bbox=dict(facecolor='white', alpha=0.8))

ax2.margins(0)
ax2.set_xlabel(args.x_label + ' / ' + r'kcal mol$^{-1}$')
ax2.set_ylabel(args.y_label + ' / ' + r'kcal mol$^{-1}$')
ax2.set_axisbelow(True)

plt.grid(color='#92A9CA', alpha=0.75, linestyle='--')

plt.tight_layout()
plt.savefig(f'{args.output}.png')
plt.savefig(f'{args.output}.pdf')
plt.show()
