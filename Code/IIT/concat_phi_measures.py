# -*- coding: utf-8 -*-
# /usr/bin/env python
# @author: link
# cd ~/Documents/thèse_onera/scripts_python/EXP_EEG/scripts

###############################################################################

import glob
import numpy as np
import pandas as pd
from os import chdir
from os import path as op
# import seaborn as sns
# import matplotlib.pyplot as plt

###############################################################################

sampfreq = 125
cluster = "temporo_sagittal"
sourcepath = "/home/link/Documents/thèse_onera"
eegpath = "experimentation/eeg/manipeeg/results/EEG/"
resultpath = "data_for_stats/hits_miss_around_detection/125Hz/"
path = op.join(sourcepath, eegpath)
chdir(path + "phi_%s/%dHz" % (cluster, sampfreq))

data_file_list_hits = [element for element in glob.glob("*.csv") if 'hits' in
                       element]
data_file_list_hits.sort()

data_file_list_miss = [element for element in glob.glob("*.csv") if 'miss' in
                       element]
data_file_list_miss.sort()

suj_hits = [int(''.join(filter(str.isdigit, data_file_list_hits[idx]))) for idx
            in range(len(data_file_list_hits))]

suj_miss = [int(''.join(filter(str.isdigit, data_file_list_miss[idx]))) for idx
            in range(len(data_file_list_miss))]

suj_both = [idx for idx in np.arange(1, 21, 1) if idx in suj_miss and idx in
            suj_hits]

###############################################################################

# suj = 2
detection = ["hits", "miss"]
# det = "hits"
phi_label = ['Phi_Geo', 'Phi_Star', 'Phi_H', 'Phi_MI']
tau = np.arange(0, 650, 1).tolist()
time = np.linspace(3, -2.16, len(tau)).round(2)
fenetre = np.repeat(np.arange(0, len(tau), int(len(tau)/65)).tolist(),
                    int(len(tau)/65)).tolist()
# x = np.repeat(np.arange(0, len(tau), int(len(tau)/65)).tolist(),
#               int(len(tau)/65)).tolist()
# y = np.repeat(np.arange(10,
#                         len(tau)+int(len(tau)/65),
#                         int(len(tau)/65)).tolist(),
#               int(len(tau)/65)).tolist()
# fenetre = [''.join(str(a) + '-' + str(b)) for a, b in zip(x, y)]

frames = []
for suj in suj_both:
    for det in detection:
        df = pd.read_csv(
            path + 'phi_%s/%dHz/phi_measures_mean_over_trials_%s_sujet_%d.csv'
            % (cluster, sampfreq, det, suj), names=tau)
        df = df.transpose()
        df.columns = phi_label
        df['Tau'] = tau
        df['Time'] = time
        df['Fenetre'] = fenetre
        df['Sujet'] = suj
        df['Detection'] = det
        cols = ['Sujet', 'Detection', 'Tau', 'Time', 'Fenetre',
                'Phi_Geo', 'Phi_Star', 'Phi_H', 'Phi_MI']
        df = df[cols]
        frames.append(df)
data = pd.concat(frames)

data.to_csv(path + resultpath + "phi_%s_measures.csv" % (cluster), index=False)

"""
data.boxplot(column='Phi_Geo', by=['Fenetre', 'Detection'], grid=False)
data.groupby(by=['Detection', 'Fenetre']).mean().boxplot(
    column='Phi_Geo', by=['Detection', 'Fenetre'])
plt.show()

g = sns.catplot(x="Time", y="Phi_Star", hue="Detection", col="Sujet",
                data=data, kind="point", palette="bright")
g.set_xticklabels(rotation=90)
# g.savefig(path + "figures/figures_epochs_%s/%s_moving_window_detection_clusters.png" %(analysis, var))
plt.show()


g = sns.FacetGrid(data=df, hue="Detection", col="Cluster", col_wrap=3, height=3, col_order=cluster_sets)
g.map(sns.pointplot, "Time", dep_var)
g.set_xticklabels(rotation=30)
plt.tight_layout()
plt.show()
g.savefig(path + "figures/figures_epochs_%s/%s_mowing_window_detection_clusters.png" %(analysis, dep_var), dpi=150)

g = sns.FacetGrid(data=data, hue="Detection", col="Sujet", col_wrap=4, height=3, col_order=None)
g.map(sns.pointplot, "Time", dep_var)
g.set_xticklabels(rotation=30)
plt.tight_layout()
plt.show()

sns.pointplot(data=data, x='Time', y=dep_var, hue='Detection', dodge=True, markers=['o', 's'], capsize=.1, errwidth=1, palette='colorblind')
plt.show()

df.groupby(['Time', 'Detection', 'Cluster'])[dep_var].agg(['mean', 'std']).round(3)

g = sns.FacetGrid(data, col="Sujet", col_wrap=4, height=2)

g = sns.FacetGrid(df, col="Cluster", hue="Detection", col_wrap=3, height=3)
g.map(sns.scatterplot, "Time", "Hurst", alpha=.7)
g.add_legend()
plt.show()

g = sns.FacetGrid(df, col="Cluster", hue="Detection", col_wrap=3, height=3, margin_titles=True)
g.map(sns.regplot, "Time", "Hurst", color=".3", fit_reg=False, x_jitter=.1)
g.add_legend()
plt.show()
"""

###############################################################################
