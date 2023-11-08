import pandas as pd
from sklearn import metrics
import matplotlib.pyplot as plt
from operator import itemgetter

benchlist = ["max50", "nomax50", "localpfam_nomax50"]
benchlist_show = ["max50", "nomax50", "local"]
approaches = ["DCTMCdomain", "DCTMCglobal", "csblast", "phmmer", "hhsearch", "blast", "fasta", "ublast", "usearch"]
approaches_show = ["DCTdomain", "DCTglobal", "CS-BLAST", "phmmer", "HHsearch", "BLAST", "FASTA", "UBLAST", "USEARCH"]

for b in range(len(benchlist)):  #pylint: disable=consider-using-enumerate
    bench = benchlist[b]
    bench_show = benchlist_show[b]
    print(f"\nBench: {bench}")
    rocs = []
    for (ap, ap_show) in zip(approaches, approaches_show):
        file = f"bench/homo/results/pfam_{bench}/{ap}_results.csv"
        df = pd.read_csv(file)
        df.columns = ['label', 'pred']
        fpr, tpr, thresholds = metrics.roc_curve(df['label'], df['pred'], pos_label=1)
        auc = metrics.auc(fpr, tpr)
        print(f"{ap_show} {auc:.3f}")
        rocs.append([ap_show, fpr, tpr, auc])

    roc_sorted = sorted(rocs, key=itemgetter(3), reverse=True)
    plt.figure(b)
    for (ap, fpr, tpr, auc) in roc_sorted:
        plt.plot(fpr, tpr, label=f"{ap} (AUC={auc:.3f})")

    plt.plot([0,1], [0,1], "k--", label="Random (AUC = 0.500)")
    plt.axis("square")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"{bench_show}")
    plt.legend()
    plt.savefig(f"{bench_show}-roc.pdf")
    plt.close()
