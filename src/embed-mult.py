"""Embeds a FASTA file using multiple layers of ESM2.

Usage:
    embed-mult.py --input <fasta> --npy <npy> [--gpu] [--maxlen] [--overwrite] [--checkpoint]

Yuzhen Ye, Indiana University, Nov 2023
"""

import argparse
import sys
import os
import torch
import esm
import numpy as np

def readfasta(filename: str) -> iter:
    """Read fasta file and returns a generator of (seqid, sequence) tuples one at a time.

    Args:
        filename (str): path to fasta file

    Returns:
        generator: (seqid, sequence) tuples
    """
    seqid, seqseq = [], []
    inf = open(filename, "r", encoding="utf-8")
    for aline in inf:
        if aline[0] == '>':
            if len(seqid) == 1:
                yield seqid, seqseq
                seqid, seqseq = [], []
            seqid.append(aline[1:-1])
            seqseq.append("")
        else:
            seqseq[-1] += aline.strip()
    yield seqid, seqseq
    inf.close()

def split_seq(seqseq: str, seqid: str, maxlen: int, overlap: int) -> list:
    """Returns a list of (seqid, subseq) tuples, where subseq is a substring of seqseq[idx] of
    length maxlen.

    Args:
        seqseq (str): sequence
        seqid (str): sequence id
        maxlen (int): maximum length of subseq
        overlap (int): number of overlapping positions between adjacent subseqs

    Returns:
        list of (seqid, subseq) tuples
    """

    embed = []
    for i in range(0, len(seqseq), maxlen-overlap):
        subseq = seqseq[i:i+maxlen]
        if len(subseq) > overlap:
            embed.append((seqid, subseq))
    return embed

def combine_contacts(mat1: torch.Tensor, mat2: torch.Tensor, inc: int, times: int) -> torch.Tensor:
    """Returns a larger square matrix combining two smaller square matrices.

    Mat1 and mat2 are both square matrices. New matrix (mat3) is of size (n+m, n+m) where
    n and m are the dimensions of mat1 and mat2. Overlapping indices are averaged.

    Args:
        mat1 (torch.Tensor): Running matrix of contacts (n x n).
        mat2 (torch.Tensor): Matrix to be added to mat1 (m x m).
        inc (int): Number of indices to increase the size of the matrix by (inc x inc).
        times (int): Running total of times the function has been called.

    Returns:
        torch.Tensor: Matrix of combined contacts (n+inc x n+inc).
    """

    # Add space to mat1
    olp = inc * times
    mlen1, mlen2 = mat1.size(0), mat2.size(0)
    mlen3 = olp + mlen2
    new_mat = torch.zeros((mlen3, mlen3), device=mat1.device)
    new_mat[:mlen1, :mlen1] = mat1

    # Add mat2 to new_mat
    new_mat[olp:mlen3, olp:mlen3] = new_mat[olp:mlen3, olp:mlen3] + mat2
    new_mat[olp:mlen1, olp:mlen1] = new_mat[olp:mlen1, olp:mlen1] / 2

    return new_mat

def embed_seq(model: esm.pretrained, batch_converter: esm.Alphabet, seqid: str, seqseq: str, device: str) -> dict:
    """Returns a dictionary of embeddings and contacts for a single sequence.

    Args:
        model (esm.pretrained): ESM model
        batch_converter (esm.Alphabet): ESM tokenizer
        seqid (str): sequence id
        seqseq (str): sequence
        device (str): device to use for embedding

    Returns:
        dict: dictionary where key is output label and value is the embedding/contact map
    """

    # Prepare inputs for embedding
    data = [(seqid, seqseq)]
    _, _, batch_tokens = batch_converter(data)  # batch_labels, batch_strs
    batch_tokens = batch_tokens.to(device)

    # Embed sequence and throw error if sequence is too long
    try:
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[15, 21], return_contacts=True)
    except RuntimeError:
        print("RuntimeError occured, try smaller --maxlen")
        sys.exit()

    # indexing [0][1:-1] removes the start and end tokens
    r15a = results["representations"][15][0][1:-1]
    r21a = results["representations"][21][0][1:-1]
    cta = results["contacts"][0]
    return {'e15': r15a, 'e21': r21a, 'ct': cta}

def get_embeds(model: esm.pretrained, batch_converter: esm.Alphabet, seqid: str, seqseq: str, device: str, maxlen: int) -> dict:
    """Returns a dictionary of embeddings and contacts for a single sequence. If the sequence is
    too long, it is split into smaller, overlapping sequences. The overlapping regions are averaged
    and the non-overlapping regions are concatenated.

    Args:
        model (esm.pretrained): ESM model
        batch_converter (esm.Alphabet): ESM tokenizer
        seqid (str): sequence id
        seqseq (str): sequence
        device (str): device to use for embedding
        maxlen (int): maximum length of subseq

    Returns:
        dict: dictionary where key is output label and value is the embedding/contact map
    """

    # Split sequence if too long and embed individually
    overlap = 200
    if len(seqseq) > maxlen:
        print(f"seq {seqid} too long {len(seqseq)}, splitting")
        subseqs = split_seq(seqseq, seqid, maxlen, overlap)
    else:  # otherwise embed whole sequence
        print(f"seq {seqid} {len(seqseq)}")
        subseqs = [(seqid, seqseq)]

    # Embed each subsequence and combine them if necessary
    edata = {}
    for i, seq in enumerate(subseqs):
        embed = embed_seq(model, batch_converter, seq[0], seq[1], device)
        if not edata:
            edata = embed
            continue
        for key, value in embed.items():
            if key == 'ct':
                edata[key] = combine_contacts(edata[key], value, maxlen-overlap, i)
                continue
            edata[key][-overlap:] = (edata[key][-overlap:] + value[:overlap]) / 2
            edata[key] = torch.cat((edata[key], value[overlap:]), axis=0)

    # Convert to numpy arrays (on CPU)
    for key, value in edata.items():
        edata[key] = value.cpu().numpy()

    return edata

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="FASTA file", required=True)
    parser.add_argument("--npy", help="save embeddings and contact to files under the given folder (one for each protein)", required=True)
    parser.add_argument("--gpu", help="using gpu, otherwise using cpu", default=True, required=False)
    parser.add_argument("--maxlen", help="max length of sequence to embed", type=int, default=1000, required=False)
    parser.add_argument("--overwrite", help="redo if npy exists; otherwise skip if exists", type=bool, default=False, required=False)
    parser.add_argument("--checkpoint", help="ESM-2 model checkpoint", default="t30", required=False)
    args = parser.parse_args()

    if not os.path.exists(args.npy):
        os.mkdir(args.npy)

    #print("load model...")
    torch.cuda.empty_cache()
    device = torch.device('cuda:0' if torch.cuda.is_available() and args.gpu else 'cpu')

    # Load ESM-2 model
    if args.checkpoint == "t30":
        model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    elif args.checkpoint == "t33":
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    batch_converter = alphabet.get_batch_converter()
    model.eval()  # disables dropout for deterministic results
    model.to(device)

    # Get filename for .list file
    listfile = args.input.split('.')[0] + ".list"
    with open(listfile, "w", encoding='utf8') as listf:
        listf.write('#protein\n')

    fasta = readfasta(args.input)
    while True:

        # Read next sequence in file
        try:
            seqid, seqseq = next(fasta)
            seqid, seqseq = seqid[0], seqseq[0]
        except StopIteration:
            break

        # For each sequence in batch, prepare inputs for embedding
        filen = args.npy + "/" + seqid + ".npz"
        with open(listfile, "a", encoding='utf8') as listf:  # add to .list file
            listf.write(f"{seqid}\n")
        if (not args.overwrite) and os.path.exists(filen):  # skip if file exists
            print(f"seq {seqid} already exists")
            continue
        edata = get_embeds(model, batch_converter, seqid, seqseq, device, args.maxlen)
        np.savez_compressed(filen, s=seqseq, e15=edata['e15'], e21=edata['e21'], ct=edata['ct'])


if __name__ == "__main__":
    main()
