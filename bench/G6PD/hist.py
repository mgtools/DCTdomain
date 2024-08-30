import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

g6pd = pd.read_csv("G6PD-dctsim.txt", sep=" ")
label = pd.read_csv("G6PD.pair", sep=" ")
g6pd['label'] = label['label']

nomax50_dctsim = "./nonhomolog-dctsim.txt"
nomax50_labelfile = "./nonhomolog-pair.txt"
nomax50_data = pd.read_csv(nomax50_dctsim, sep=" ")
label = pd.read_csv(nomax50_labelfile, sep=" ")
nomax50_data['label'] = label['label']

data = g6pd
bench = "G6PD containing proteins"

simlist = ["domain", "global"]
for idx in range(2):
    sim = simlist[idx]
    tag = f'sim-{sim}'
    h = data[tag]
    hl = data.loc[data['label'] == 'local', tag]
    hg = data.loc[data['label'] != 'local', tag]
    nonhom = nomax50_data.loc[nomax50_data['label'] == 'nonhom', tag]

    num_bins = 10

    #local
    n, bins, patch = plt.hist(hl, alpha=0.5, label='local similarity', bins=num_bins, density=True, color="skyblue")
    mu, sigma = hl.mean(), hl.std()
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    plt.plot(bins, y, '--', color="skyblue")

    #global homolog
    n, bins, patch = plt.hist(hg, alpha=0.5, label='global homolog', bins=num_bins, density=True, color="orange")
    mu, sigma = hg.mean(), hg.std()
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    plt.plot(bins, y, '--', color="orange")

    #non homolog
    n, bins, patch = plt.hist(nonhom, alpha=0.5, label='non-homolog', bins=num_bins, density=True, color="gray")
    mu, sigma = nonhom.mean(), nonhom.std()
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    plt.plot(bins, y, '--', color="gray")

    plt.xlabel(f'DCT-{sim} similarity score', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    #plt.legend(title='Protein pairs', fontsize=12)
    plt.legend(fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=11)
    plt.savefig(f"G6PD-{tag}-hist.pdf")
    plt.close()

#fig.tight_layout()
#plt.title('DCT similarity score distribution (G6PD containing proteins)')
#plt.show
