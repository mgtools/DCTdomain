"""Writes most confident contacts to file after sequences have been embedded.

Yuzhen Ye, Indiana University, Nov 2023
"""

import argparse
import sys
import os
from operator import itemgetter
import numpy as np

def readfasta(filename: str) -> tuple:
    """Read fasta file and return a list of (seqid, seq) tuples.

    Args:
        filename (str): path to fasta file

    Returns:
        list of (seqid, seq) tuples
    """

    seqid, seqseq = [], []
    inf = open(filename, "r", encoding='utf-8')
    for aline in inf:
        if aline[0] == '>':
            seqid.append(aline[1:-1])
            seqseq.append("")
        else:
            seqseq[-1] += aline.strip()
    inf.close()
    return seqid, seqseq

def writece_a(outfile: str, seqid: str, seqseq: str, cta: np.ndarray, t: float = 2.6):
    """write in CE for domain segmentation
    the top t * L contact pairs, L is the length (t is the alpha parameter in FUpred paper)

    Args:
        outfile (str): path to output file
        seqid (str): sequence id
        seqseq (str): sequence
        cta (np.array): contact map
        t (float): threshold for contact map (0.5 to 5)
    """

    # Sort contacts by confidence
    slen = len(seqseq)
    cta1 = cta.reshape(slen, slen)
    ssep = 5
    data = []
    for i in range(slen - ssep):
        for j in range(i + ssep, slen):
            data.append([cta1[i][j], i, j])
    ct_sorted = sorted(data, key=itemgetter(0), reverse=True)

    # Get top t * L contacts
    sout = ""
    tot = int(t * slen)
    print(f"protein length {slen} contacts considered {tot} confidence {ct_sorted[tot - 1][0]}")
    for s in range(tot):
        i, j = ct_sorted[s][1], ct_sorted[s][2]
        if not sout:
            sout = f"CON   {i} {j} {cta1[i][j]:.6f}"
        else:
            sout += f",{i} {j} {cta1[i][j]:.6f}"

    # Write sequence info and top contacts to file
    with open(outfile, "w", encoding='utf8') as out_f:
        out_f.write(f"INF   {seqid} {slen}\n")
        out_f.write(f"SEQ   {seqseq}\n")
        out_f.write(f"SS    {'C' * slen}\n")
        out_f.write(sout + "\n")

def writece_t(outfile: str, seqid: str, seqseq: str, cta: np.ndarray, t: float = 0.1):
    """#write in CE for domain segmentation (select contacts according to confidence)
    """

    slen = len(seqseq)
    cta1 = cta.reshape(slen, slen)
    ssep = 5
    sout = ""
    for i in range(slen - ssep):
        for j in range(i + ssep, slen):
            if cta1[i][j] >= t:
                if not sout:
                    sout = f"CON   {i} {j} {cta1[i][j]:.6f}"
                else:
                    sout += f",{i} {j} {cta1[i][j]:.6f}"

    # Write sequence info and top contacts to file
    with open(outfile, "w", encoding='utf8') as out_f:
        out_f.write(f"INF   {seqid} {slen}\n")
        out_f.write(f"SEQ   {seqseq}\n")
        out_f.write(f"SS    {'C' * slen}\n")
        out_f.write(sout + "\n")

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
        sys.exit(f"{args.npy} doesn't exist")

    seqid, seqseq = readfasta(args.input)
    for idx in range(len(seqid)):  #pylint: disable=C0200

        # Load corresponding embedding file
        npyfile = args.npy + "/" + seqid[idx] + ".npz"
        if not os.path.exists(npyfile):
            print(f"{npyfile} doesn't exist; please run embedding first")
            continue
        out_a = args.ce + "/" + seqid[idx] + ".ce"
        if os.path.exists(out_a):
            print(f"{out_a} exists; skip")
            continue
        data = np.load(npyfile)

        # Get contact map and write most confident contacts to file
        cta = data["ct"]
        out_a = args.ce + "/" + seqid[idx] + ".ce"
        if args.ce_t:
            writece_t(out_a, seqid[idx], seqseq[idx], cta, t=args.ce_t)
        else:
            writece_a(out_a, seqid[idx], seqseq[idx], cta, t=args.ce_a)

if __name__ == "__main__":
    main()
