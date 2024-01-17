"""Embeds a FASTA file using multiple layers of ESM2_t33_650M_UR50D.

Usage:
    embed-mult.py --input <fasta> --npy <npy> [--gpu] [--batchsize <batchsize>] [--overwrite]

Yuzhen Ye, Indiana University, Nov 2023
"""

import argparse
import sys
import os
import torch
import esm
import numpy as np


def readfasta(filename: str, batchsize: int):
    """Read fasta file and return a list of (seqid, seq) tuples.

    Args:
        filename (str): path to fasta file
        batchsize (int): number of seqs to read at a time

    Returns:
        list of (seqid, seq) tuples
    """
    seqid, seqseq = [], []
    inf = open(filename, "r", encoding="utf-8")
    for aline in inf:
        if aline[0] == '>':
            if len(seqid) == batchsize:
                yield seqid, seqseq
                seqid, seqseq = [], []
            seqid.append(aline[1:-1])
            seqseq.append("")
        else:
            seqseq[-1] += aline.strip()
    yield seqid, seqseq
    inf.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="FASTA file", required=True)
    parser.add_argument("--npy", help="save embeddings and contact to files under the given folder (one for each protein)", required=True)
    parser.add_argument("--gpu", help="using gpu, otherwise using cpu", required=False)
    parser.add_argument("--batchsize", help="batchsize; default 1", default=1, required=False)
    parser.add_argument("--overwrite", help="redo if npy exists; otherwise skip if exists", type=bool, default=False, required=False)
    args = parser.parse_args()

    if not os.path.exists(args.npy):
        os.mkdir(args.npy)

    print("load model...")
    torch.cuda.empty_cache()
    device = torch.device('cuda:0' if torch.cuda.is_available() and args.gpu else 'cpu')

    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results

    if args.gpu:
        model.to(device)

    fasta = readfasta(args.input, args.batchsize)
    while True:
        try:
            seqid, seqseq = next(fasta)
        except StopIteration:
            break
        data = []

        # For each sequence in batch, prepare inputs for embedding
        for idx in range(len(seqid)):  #pylint: disable=C0200
            filen = args.npy + "/" + seqid[idx] + ".npz"
            if (not args.overwrite) and os.path.exists(filen):  # skip if exists
                print(f"seq {seqid[idx]} already exists")
                continue
            if len(seqseq[idx]) > 6000:  # skip if too long
                print(f"seq {seqid[idx]} too long {len(seqseq[idx])}")
                continue
            print(f"seq {seqid[idx]} {len(seqseq[idx])}")
            data.append((seqid[idx], seqseq[idx]))
            _, _, batch_tokens = batch_converter(data)  # batch_labels, batch_strs
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
        if len(data) == 0:  # skip if no inputs were processed
            continue

        # Extract per-residue representations (on CPU)
        #batch_tokens = batch_tokens.to('cuda')
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[13, 25, 33], return_contacts=True)
        r33a = results["representations"][33]
        r13a = results["representations"][13]
        r25a = results["representations"][25]
        cta = results["contacts"]

        # Save data into npz file
        np.set_printoptions(threshold=sys.maxsize)
        for (sid, seq), _, e13, e25, e33, ct in zip(data, batch_lens, r13a, r25a, r33a, cta):
            filen = args.npy + "/" + sid
            np.savez_compressed(filen, s=seq, e13=e13, e25=e25, e33=e33, ct=ct)

if __name__ == "__main__":
    main()
