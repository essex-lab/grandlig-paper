"""
This is a basic script to plot the various useful GCNCMC metrics/track variables
Data must be avaliable in the ncmc_data.pkl file (settings data=True in ncmc_mover.report)
"""

# Modules to be imported
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import seaborn as sns
import pickle as pkl

# Load pickle file
with open("ncmc_data.pkl", "rb") as f:
        ncmc_data = pkl.load(f)

def plot_moves(ax, ncmc_data=ncmc_data):
    categories = ['Insert', 'Delete']
    n_inserts = ncmc_data["n_inserts"]
    n_deletes = ncmc_data["n_deletes"]
    n_accepted_inserts = ncmc_data["n_accepted_inserts"]
    n_accepted_deletes = ncmc_data["n_accepted_deletes"]

    # Positions of the bars on the x-axis
    x = np.arange(len(categories))

    # Width of the bars
    width = 0.75


    # Plotting the stacked bars
    bars_inserts = ax.bar(x[0], n_inserts, width, alpha=0.7)
    bars_accepted_inserts = ax.bar(x[0], n_accepted_inserts, width, bottom=n_inserts, label='Accepted Insertions', alpha=0.7)

    bars_deletes = ax.bar(x[1], n_deletes, width, alpha=0.7)
    bars_accepted_deletes = ax.bar(x[1], n_accepted_deletes, width, bottom=n_deletes, label='Accepted Deletions', alpha=0.7)

    # Adding labels and title
    ax.set_xlabel('Move Type')
    ax.set_ylabel('Counts')
    #ax.set_title('Stacked Operations with Accepted Overlays')
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.legend(loc="lower right", ncols=2)

    
def plot_outcomes(ax, ncmc_data=ncmc_data):
    uni_outcomes = np.unique(ncmc_data["outcome"], return_counts=True)
    ax.bar(uni_outcomes[0], uni_outcomes[1])
    # Adding labels and title
    ax.set_ylabel('Counts')
    ax.set_title('Move outcomes - N.B some double counting')
    ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=-45, ha='left')

def plot_avgN(ax, ncmc_data=ncmc_data):

    ax.plot(range(ncmc_data["n_moves"]), ncmc_data["Ns"], ls='None', marker='o')
    avg_N = [np.mean(ncmc_data["Ns"][:i]) for i in range(ncmc_data['n_moves'])]

    ax.plot(range(ncmc_data["n_moves"]), avg_N, label="Avg. N in sphere")
    ax.legend()
    ax.set_xlabel("GCNCMC Moves")
    ax.set_ylabel("N")

def plot_work_dist(ax, accepted=False, ncmc_data=ncmc_data):
    from openmm.unit import Quantity
    if not accepted:
        insert_works = np.asarray([i._value if type(i) == Quantity else i for i in ncmc_data["insert_works"]])
        delete_works = np.asarray([i._value if type(i) == Quantity else i for i in ncmc_data["delete_works"]]) * -1
        ax.set_title('Work Distributions - All Moves')


    else:
        insert_works = np.asarray([i._value if type(i) == Quantity else i for i in ncmc_data["accepted_insert_works"]])
        delete_works = np.asarray([i._value if type(i) == Quantity else i for i in ncmc_data["accepted_delete_works"]]) * -1
        ax.set_title('Work Distributions - Accepted Moves')
        
    if not np.all(np.isnan(insert_works)):
        ax.hist(insert_works, bins='auto', density=True, alpha=0.5, histtype='stepfilled', color=sns.color_palette('muted')[0], label="Insertion") 
        ax.hist(insert_works, bins='auto', density=True, alpha=1.0, histtype='step', lw=1.5, color=sns.color_palette('muted')[0])
    
    if not np.all(np.isnan(delete_works)):
        ax.hist(delete_works, bins='auto', density=True, alpha=0.5, histtype='stepfilled', color=sns.color_palette('muted')[1], label="Deletion") 
        ax.hist(delete_works, bins='auto', density=True, alpha=1.0, histtype='step', lw=1.5, color=sns.color_palette('muted')[1])
    ax.set_xlabel('Work Done / kcal/mol')
    ax.set_ylabel('Counts')
    ax.legend(loc='upper right')
    
def plot_acc_probs(ax, type, ncmc_data=ncmc_data):
    if type == 'all':
        acc_probs = np.asarray(ncmc_data['acceptance_probabilities'])
        n = ncmc_data['n_moves']
        ax.set_xlabel('Num. Moves')
        ax.set_title("Acceptance Probabilities")
        #ax.set_ylabel('Acceptance Probability')
    elif type == 'insert':
        acc_probs = np.asarray(ncmc_data['insert_acceptance_probabilities'])
        n = ncmc_data['n_inserts']
        ax.set_xlabel('Num. Insertion Moves')
        #ax.set_ylabel('Insertion Acceptance Probability')
    elif type == 'delete':
        acc_probs = np.asarray(ncmc_data['delete_acceptance_probabilities'])
        n = ncmc_data['n_deletes']  
        ax.set_xlabel('Num. Deletion Moves')
        #ax.set_ylabel('Deletion Acceptance Probability')
   
    # Replace -1 with 0
    acc_probs[acc_probs == -1] = 0
    # Anything greater than 1 make equal to 1
    acc_probs[acc_probs > 1] = 1
    avg_acc_probs = [np.mean(acc_probs[:i]) for i in range(n)] 


    #plt.ylim(0, 1)
    ax.plot(range(n), avg_acc_probs)

def format_axes(fig):
    for i, ax in enumerate(fig.axes):
        ax.text(0.5, 0.5, "ax%d" % (i+1), va="center", ha="center")
        ax.tick_params(labelbottom=False, labelleft=False)



if __name__ == '__main__':

    
    sns.set_theme()


    fig = plt.figure(layout="constrained", figsize=(25, 15))

    gs = GridSpec(5, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])

    gs00 = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1:3, 0:])
    ax4 = fig.add_subplot(gs00[0:, 0])
    ax5 = fig.add_subplot(gs00[0, 1])

    gs01 = GridSpecFromSubplotSpec(3, 2, subplot_spec=gs[3:5, 0:])

    ax6 = fig.add_subplot(gs01[0, 0:])
    ax7 = fig.add_subplot(gs01[1, 0:])
    ax8 = fig.add_subplot(gs01[2, 0:])

    plot_moves(ax1)
    plot_outcomes(ax2)
    plot_avgN(ax3)
    plot_work_dist(ax4, accepted=False)
    plot_work_dist(ax5, accepted=True)

    plot_acc_probs(ax6, "all")
    plot_acc_probs(ax7, "insert")
    plot_acc_probs(ax8, "delete")


    fig.suptitle("GCNCMC Simulation Analysis")
    #format_axes(fig)
    plt.savefig("Test_analysis.pdf")
    #plt.show()
