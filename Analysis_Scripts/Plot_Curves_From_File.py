# Modules to import
import argparse
import os

import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from glob import glob
from scipy.optimize import curve_fit
import pickle


small_font = 16
medium_font = 18
large_font = 20

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=medium_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=17)

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-pkl", "--pkl", help="Input the pickle file", type=str)
parser.add_argument("-pro", '--pro', help='Input protein name', type=str)
parser.add_argument("-out", '--out', help='Input output name', type=str)

#parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')


args = parser.parse_args()

# Load in the data
with open(args.pkl, 'rb') as fi:
    plot_data = pickle.load(fi)

print(plot_data)

frags = list(plot_data.keys())
all_dgs = []
for h, frag in enumerate(frags):
    x, y, results = plot_data[frag]
    dG_mean, dG_std, tau = results
    all_dgs.append(dG_mean)

sorted_dGs, sorted_frags = zip(*sorted(zip(all_dgs, frags)))

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)
plt.xlabel('Ligand Concentration / M')
plt.ylabel('Site Occupancy')

cmap = matplotlib.cm.get_cmap("tab10")
colors = cmap(np.linspace(0, 1, len(all_dgs)))

all_dgs = []
linestyles = ["-", '--', ':']
ls_index = 0
linestyle = linestyles[ls_index]
c_index = 0
all_dg_stds = []
for h, frag in enumerate(sorted_frags):
    x, y, results = plot_data[frag]
    x = [10**i for i in x]
    dG_mean, dG_std, tau = results
    all_dgs.append(dG_mean)
    ax1.plot(x, y, ls=linestyle, c=colors[c_index], label=f'{frag.title()}: '+r'${:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}\ $'.format(dG_mean._value, dG_std._value)+r'($\tau={:.2f}$)'.format(tau))
    c_index += 1
  #  if c_index == 10:
  #      c_index =0
 #       ls_index += 1
#        linestyle = linestyles[ls_index]

plt.xlim(10**-10, 10**-3)
plt.xscale("log")
ax1.axhline(0.5)
handles, labels = ax1.get_legend_handles_labels()
#dGs = [float(label.split()[1]) for label in labels]
all_dgs, handles, labels = zip(*sorted(zip(all_dgs, handles, labels), key=lambda t: t[0]))
#dGs, labels = zip(*sorted(zip(dGs, labels)))

plt.legend(handles, labels, loc='lower right', ncol=2, fancybox=True, shadow=True)
plt.ylim(0.00, 1.01)

#plt.title(r'T4L99A Ligands')
#plt.title(r'$\beta$-Cyclodextrin Ligands')
plt.title(r'{} Titrations'.format(args.pro))
plt.savefig(f"{args.out}.pdf", format='pdf', bbox_inches='tight')
plt.savefig(f"{args.out}.png", bbox_inches='tight')

plt.show()





