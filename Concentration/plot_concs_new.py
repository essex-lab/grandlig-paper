
import argparse
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy.stats import gaussian_kde
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input data file')
parser.add_argument('-c', '--concentration', type=float, help='Target concentration (M)')
parser.add_argument('-e', '--equil', type=int, default=0, help='Target concentration (M)')
parser.add_argument('-s', '--skip', default=[1, 1], nargs='+', type=int, help='Sample every nth frame')
parser.add_argument('-x', '--xlim', type=float, default=None, help='x-axis upper limit limit')
parser.add_argument('-b', '--binwidth', type=float, default=None, help='Bin width')
parser.add_argument('--yzero', default=False, action='store_true', help='Set y-axis minimum to zero')
parser.add_argument('--title', default=None, help='Graph title')
parser.add_argument('-o', '--output', default=None, help='Output file')
args = parser.parse_args()


small_font = 20
medium_font = 22
large_font = 24

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=medium_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=16)

# Set up plot
# Plot data
fig = plt.figure(figsize=(12.8, 5)) #, facecolor='#EDEDED')


gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
ax1 = fig.add_subplot(gs[0])
if args.title is not None:
    plt.title(args.title)
ax2 = fig.add_subplot(gs[1])
axes = [ax1, ax2]

# Plot conc vs number of moves
axes[0].axhline(args.concentration, linestyle='--', color='black', label='Macroscopic concentration', lw=2, zorder=10)

if args.equil > 0:
    axes[0].axvline(args.equil, linestyle='--', color='green', label='Equilibration')

axes[1].axhline(args.concentration, linestyle='--', color='black', lw=2, zorder=10)

def get_data(file, equil):
    names = ['Moves', 'repeat_0', 'repeat_1', 'repeat_2', 'repeat_3']
    df = pd.read_csv(file, usecols=names)
    print(df.head(5))
    n_repeats = len(df.columns) - 1
    after_equil = df[df["Moves"] >= equil]
    #after_equil = after_equil.dropna()

    repeat_means = []
    repeat_medians = []
    repeat_upper = []
    repeat_lower = []
    for i in range(n_repeats):
        repeat_means.append(after_equil[f"repeat_{i}"].mean())
        repeat_medians.append(after_equil[f"repeat_{i}"].median())
        repeat_upper.append(after_equil[f"repeat_{i}"].quantile(0.75))
        repeat_lower.append(after_equil[f"repeat_{i}"].quantile(0.25))

    mean_means = np.mean(repeat_means)
    std_err_means = np.std(repeat_means) / np.sqrt(len(repeat_means))

    mean_quantiles = [np.mean(repeat_medians), np.mean(repeat_upper), np.mean(repeat_lower)]


    # Drop moves column from og df
    df.drop("Moves", axis=1, inplace=True)
    for i in range(n_repeats):  # Make the top row zeros
        df[f"repeat_{i}"][0] = df[f"repeat_{0}"][0]

    initial_conc = df[f"repeat_{0}"][0]

    #df = df.dropna()
    moves = range(len(df))
    move_means = df.mean(axis=1)
    move_errs = df.std(axis=1) / np.sqrt(n_repeats)

    after_equil.drop("Moves", axis=1, inplace=True)

    return moves, initial_conc, mean_means, std_err_means, mean_quantiles, move_means, move_errs, after_equil

def write_results(mode, initial_conc, mean_conc, err, quantiles):
    with open("Final_Results.txt", mode) as fo:
        fo.write("Desired Concentration = {:.2f} M\n".format(args.concentration))
        fo.write("Starting Concentration = {:.2f} M\n".format(initial_conc))

        fo.write('<c>  = {:.2f} ± {:.2f} M\n'.format(mean_conc, err))

        fo.write('Median = {:.2f} M, Interquartile range: [{:.2f}, {:.2f}] M\n'.format(quantiles[0],
                                                                                        quantiles[1], quantiles[2]))

def plot(file, equil, binwidth, color, mode, axes=axes):

    moves, initial_conc, mean_means, std_err_means, mean_quantiles, move_means, move_errs, df_ = get_data(file, equil)
    write_results(mode, initial_conc, mean_means, std_err_means, mean_quantiles)

    axes[0].plot(moves, move_means, '-', color=color, label=f"From {initial_conc:.2f} M")
    axes[0].fill_between(moves, move_means - move_errs, move_means + move_errs, color=color, alpha=0.3, lw=0)


    all_concs = np.hstack(df_.values)


    # ax2.axhline(np.quantile(all_concs, 0.5), linestyle='dotted', lw=3, color="k")
    if binwidth is not None:
        bins = np.arange(min(all_concs) - 0.5 * binwidth, max(all_concs) + 0.5 * binwidth, binwidth)
    else:
        bins = 'auto'


    axes[1].hist(move_means[5000:], bins=bins, color=color, density=True, orientation='horizontal', histtype='stepfilled',
             alpha=0.3, lw=0)
    axes[1].hist(move_means[5000:], bins=bins, color=color, density=True, orientation='horizontal', histtype='step', alpha=1.0,
             lw=1.5)




plot(args.input, args.equil, args.binwidth, "red", mode="w")

# Format graph
ax1.set_xlabel('Number of GCNCMC moves')
ax1.set_ylabel('Concentration / M')
if args.xlim is not None:
    ax1.set_xlim(0, float(args.xlim))
else:
    ax1.set_xlim(left=0)
if args.yzero:
    ax2.set_ylim(bottom=0, top=2.0)


ax1.set_ylim(ax2.get_ylim())
ax2.set_xticks([])
ax2.set_yticks([])
ax1.legend(loc='lower right', ncol=2, fancybox=True, shadow=True) #, facecolor='#EDEDED')
ax1.axis('on')
ax2.axis('off')
ax1.set_facecolor("white")
ax2.set_facecolor('white')
plt.tight_layout()
plt.subplots_adjust(wspace=0)
plt.tight_layout()
plt.subplots_adjust(wspace=0)



if args.output is not None:
    plt.savefig(args.output+".png")
    plt.savefig(args.output + ".pdf")

plt.show()


