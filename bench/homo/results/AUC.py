"""Plots AUC of all methods on each benchmark.

Yuzhen Ye, Indiana University, Nov 2023
"""

from operator import itemgetter
import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt

benchlist = ["pfam_max50", "pfam_nomax50", "pfam_localpfam_nomax50", "supfam_nomax50", "gene3d_nomax50"]
benchlist_show = ["max50", "pfam_nomax50", "local", "supfam_nomax50", "gene3d_nomax50"]
approaches = ["DCTMCdomain", "DCTMCglobal", "csblast", "phmmer", "hhsearch", "blast", "fasta", "ublast", "usearch", "prost"]
approaches_show = ["DCTdomain", "DCTglobal", "CS-BLAST", "phmmer", "HHsearch", "BLAST", "FASTA", "UBLAST", "USEARCH", "PROST"]

for b in range(len(benchlist)):  #pylint: disable=consider-using-enumerate
    bench = benchlist[b]
    bench_show = benchlist_show[b]
    print(f"\nBench: {bench}")
    rocs = []
    for (ap, ap_show) in zip(approaches, approaches_show):
        file = f"bench/homo/results/{bench}/{ap}_results.csv"
        df = pd.read_csv(file)
        df.columns = ['label', 'pred']
        fpr, tpr, thresholds = metrics.roc_curve(df['label'], df['pred'], pos_label=1)
        print("fpr", fpr)
        print("tpr", tpr)
        print("thresholds", thresholds)
        auc = metrics.auc(fpr, tpr)
        print(f"{ap_show} {auc:.3f}")
        rocs.append([ap_show, fpr, tpr, auc])

    roc_sorted = sorted(rocs, key=itemgetter(3), reverse=True)
    plt.figure(b)
    for (ap, fpr, tpr, auc) in roc_sorted:
        plt.plot(fpr, tpr, label=f"{ap} (AUC={auc:.3f})")

    plt.plot([0,1], [0,1], "k--", label="Random (AUC = 0.500)")
    plt.axis("square")
    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.title(f"{bench_show}", fontsize=12)
    if "local" in bench:
        plt.legend(fontsize=10)
    else:
        plt.legend(fontsize=11)
    plt.savefig(f"bench/homo/results/{bench_show}-roc.pdf")
    plt.close()
