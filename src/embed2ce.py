import torch
import esm
import numpy as np
import argparse
import sys
import os
#Yuzhen Ye, Indiana University, Nov 2023
def readfasta(filename):
    seqid, seqseq = [], []
    inf = open(filename, "r")
    for aline in inf:
        if aline[0] == '>':
            seqid.append(aline[1:-1])
            seqseq.append("")
        else:
            seqseq[-1] += aline.strip()
    return seqid, seqseq

#write in CE for domain segmentation 
#the top t * L contact pairs, L is the length (t is the alpha parameter in FUpred paper)
from operator import itemgetter
def writece_a(outfile, seqid, seqseq, a, t = 2.6):
    out_a = open(outfile, "w")
    a1 = np.array(a) 
    slen = len(seqseq)
    a2 = a1.reshape(slen, slen)
    out_a.write(f"INF   protein {slen}\n")
    out_a.write(f"SEQ   {seqseq}\n")
    out_a.write(f"SS    {'C' * slen}\n")
    ssep = 5
    data = []
    for i in range(slen - ssep):
        for j in range(i + ssep, slen):
            data.append([a2[i][j], i, j])
    ct_sorted = sorted(data, key=itemgetter(0), reverse=True)

    sout = ""
    #t: 0.5 to 5
    tot = int(t * slen)
    print(f"protein length {slen} contacts considered {tot} confidence {ct_sorted[tot - 1][0]}")
    for s in range(tot):
        i, j = ct_sorted[s][1], ct_sorted[s][2]
        if not sout:
            sout = f"CON   {i} {j} {a2[i][j]:.6f}"
        else:
            sout += f",{i} {j} {a2[i][j]:.6f}"
    out_a.write(sout + "\n")
    out_a.close()

#write in CE for domain segmentation (select contacts according to confidence)
from operator import itemgetter
def writece_t(outfile, seqid, seqseq, a, t = 0.1):
    out_a = open(outfile, "w")
    a1 = np.array(a) 
    slen = len(seqseq)
    a2 = a1.reshape(slen, slen)
    out_a.write(f"INF   protein {slen}\n")
    out_a.write(f"SEQ   {seqseq}\n")
    out_a.write(f"SS    {'C' * slen}\n")
    ssep = 5
    sout = ""
    for i in range(slen - ssep):
        for j in range(i + ssep, slen):
            if a2[i][j] >= t:
                if not sout:
                    sout = f"CON   {i} {j} {a2[i][j]:.6f}"
                else:
                    sout += f",{i} {j} {a2[i][j]:.6f}"
    out_a.write(sout + "\n")
    out_a.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="input file (fasta)", required=True)
    parser.add_argument("--npy", help="embedding files under the given folder (one for each protein)", required=True)
    parser.add_argument("--ce", help="directory for saving the contact map to files", required=True)
    parser.add_argument("--ce_t", help="CE threshold (default using ce_a approach; confidence defaults to 0.1)", required=False, type=float)
    parser.add_argument("--ce_a", help="CE threshold (default using ce_t approach; ce_a defaults to 2.6)", required=False, default=2.6, type=float)
    args = parser.parse_args()

    if not os.path.exists(args.ce):
        os.mkdir(args.ce) 
    if not os.path.exists(args.npy):
        sys.exit(args.npy, " doesn't exist") 

    seqid, seqseq = readfasta(args.input)
    for idx in range(len(seqid)):
        npyfile = args.npy + "/" + seqid[idx] + ".npz"
        #load data
        if not os.path.exists(npyfile):
            print(f"{npyfile} doesn't exist; please run embedding first")
            continue
        out_a = args.ce + "/" + seqid[idx] + ".ce"
        if os.path.exists(out_a):
            print(f"{out_a} exists; skip")
            continue
        data = np.load(npyfile)
        cta = data["ct"]
        out_a = args.ce + "/" + seqid[idx] + ".ce"
        if args.ce_t:
            writece_t(out_a, seqid[idx], seqseq[idx], cta, t=args.ce_t)
        else:
            writece_a(out_a, seqid[idx], seqseq[idx], cta, t=args.ce_a)

if __name__ == "__main__":
    main()
