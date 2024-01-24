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

def embed_batch(model: esm.pretrained, batch_tokens: torch.Tensor, data: list, npyfile: str):
    """Embeds a batch of protein sequences using the given model.

    Args:
        model (esm.pretrained): ESM model
        batch_tokens (torch.Tensor): batch of sequences to embed
        data (list): list of (seqid, seq) tuples
        npyfile (str): path to output file
    """

    # Extract per-residue representations (on CPU)
    try:
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[13, 25, 33], return_contacts=True)
    except RuntimeError:
        print("RuntimeError occured, try smaller --maxgpu")
        return
    r33a = results["representations"][33].cpu()
    r13a = results["representations"][13].cpu()
    r25a = results["representations"][25].cpu()
    cta = results["contacts"].cpu()

    # Save data into npz file
    np.set_printoptions(threshold=sys.maxsize)
    for (sid, seq), e13, e25, e33, ct in zip(data, r13a, r25a, r33a, cta):
        filen = npyfile + "/" + sid
        np.savez_compressed(filen, s=seq, e13=e13, e25=e25, e33=e33, ct=ct)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="FASTA file", required=True)
    parser.add_argument("--npy", help="save embeddings and contact to files under the given folder (one for each protein)", required=True)
    parser.add_argument("--gpu", help="using gpu, otherwise using cpu", default=False, required=False)
    parser.add_argument("--batchsize", help="batchsize; default 1", default=1, required=False)
    parser.add_argument("--overwrite", help="redo if npy exists; otherwise skip if exists", type=bool, default=False, required=False)
    parser.add_argument("--maxlen", help="max length of sequence to embed", type=int, default=1200, required=False)
    parser.add_argument("--maxgpu", help="max length of sequence to embed on GPU", type=int, default=600, required=False)
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

    # Get filename for .list file
    listfile = args.input.split('.')[0] + ".list"
    with open(listfile, "w", encoding='utf8') as listf:
        listf.write('#protein\n')

    fasta = readfasta(args.input, args.batchsize)
    while True:

        # Move back to GPU if on CPU
        if args.gpu and device == 'cpu':
            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            model.to(device)

        # Read batch of sequences
        try:
            seqid, seqseq = next(fasta)
        except StopIteration:
            break

        # For each sequence in batch, prepare inputs for embedding
        data = []
        for idx in range(len(seqid)):  #pylint: disable=C0200
            filen = args.npy + "/" + seqid[idx] + ".npz"

            # Skip sequence if embedding already exists or seq is too long
            if (not args.overwrite) and os.path.exists(filen):
                print(f"seq {seqid[idx]} already exists")
                continue
            if len(seqseq[idx]) > args.maxlen:
                print(f"seq {seqid[idx]} too long {len(seqseq[idx])}, skipping")
                continue
            if args.gpu and len(seqseq[idx]) > args.maxgpu:  # move to CPU if too long
                print(f"seq {seqid[idx]} too long {len(seqseq[idx])}, moving to CPU")
                device = 'cpu'
                model.to(device)

            # Add sequence to list file
            with open(listfile, "a", encoding='utf8') as listf:
                listf.write(f"{seqid[idx]}\n")

            # Add each sequence to list of sequences to embed
            print(f"seq {seqid[idx]} {len(seqseq[idx])}")
            data.append((seqid[idx], seqseq[idx]))
            _, _, batch_tokens = batch_converter(data)  # batch_labels, batch_strs
            batch_tokens = batch_tokens.to(device)
            #batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
        if len(data) == 0:  # skip if no inputs were processed
            continue

        # Embed batch of sequences
        embed_batch(model, batch_tokens, data, args.npy)

if __name__ == "__main__":
    main()
