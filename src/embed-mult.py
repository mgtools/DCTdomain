"""Embeds a FASTA file using multiple layers of ESM2.

Usage:
    embed-mult.py --input <fasta> --npy <npy> [--gpu] [--overwrite]

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

def avg_emb(embeds: list, overlap: int) -> np.ndarray:
    """Returns a single embedding by averaging and concatenating a list of embeddings.

    Args:
        embeds (list): list of (seqid, np.ndarray) tuples
        overlap (int): number of overlapping positions between adjacent embeddings
    Returns:
        numpy array: single 1D embedding
    """

    full_emb = np.array([])
    for emb in embeds:
        if len(full_emb) == 0:
            full_emb = emb
            continue
        avg = (full_emb[-overlap:] + emb[:overlap]) / 2  # avg adjacent positions
        full_emb[-overlap:] = avg  # replace overlapping positions of previous emb with avg
        full_emb = np.concatenate((full_emb, emb[overlap:]), axis=0)  # concat rest of new emb
    return full_emb

def combine_contacts(mat1: np.ndarray, mat2: np.ndarray, inc: int, times: int) -> np.ndarray:
    """Returns a larger square matrix combining two smaller square matrices.

    mat1 has values starting from the top left corner, mat2 has values starting from the bottom
    right corner, and the overlapping indices are averaged. The number of overlapping indices is
    determine by the inc argument.

    Args:
        mat1 (numpy array): Running matrix of contacts (n x n).
        mat2 (numpy array): Matrix to be added to mat1 (m x m).
        inc (int): Number of indices to increase the size of the matrix by (inc x inc).
        times (int): Running total of times the function has been called.

    Returns:
        numpy array: Matrix of combined contacts (n+inc x n+inc).
    """

    # Create new matrices to store combined contacts
    mlen = len(mat1)
    size = mlen + inc
    zeros1 = np.zeros((size, size))
    zeros2 = np.zeros((size, size))

    # Add input matrices to new matrices
    olp = inc*times
    zeros1[:mlen, :mlen] = mat1
    zeros2[:len(mat2), :len(mat2)] = mat2
    zeros2 = np.roll(zeros2, olp, axis = 0)
    zeros2 = np.roll(zeros2, olp, axis = 1)

    # Average the overlapping indices
    zeros1[olp:mlen, olp:mlen] = (zeros1[olp:mlen, olp:mlen] + zeros2[olp:mlen, olp:mlen]) / 2

    # Add rest of zeros2 to zeros1
    zeros1[olp:size, mlen:size] = zeros2[olp:size, mlen:size]
    zeros1[mlen:size, olp:mlen] = zeros2[mlen:size, olp:mlen]

    # If any row or column is all 0's, remove it
    zeros1 = zeros1[~np.all(zeros1 == 0, axis=1)]
    zeros1 = zeros1[:, ~np.all(zeros1 == 0, axis=0)]

    return zeros1

def avg_ct(contacts: list, increase: int) -> np.ndarray:
    """Returns a single contact map by averaging and combining a list of contact maps.

    Args:
        contacts (list): list of (seqid, np.ndarray) tuples
        increase (int): number of indices to increase size of contact maps

    Returns:
        numpy array: single 2D contact map
    """

    full_contacts = np.array([])
    for i, cont in enumerate(contacts):
        if len(full_contacts) == 0:  # initialize full array
            full_contacts = cont
            continue
        full_contacts = combine_contacts(full_contacts, cont, inc = increase, times = i)
    return full_contacts

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
    r15a = results["representations"][15].cpu()[0][1:-1]
    r21a = results["representations"][21].cpu()[0][1:-1]
    cta = results["contacts"].cpu()[0]
    return {'e15': r15a, 'e21': r21a, 'ct': cta}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="FASTA file", required=True)
    parser.add_argument("--npy", help="save embeddings and contact to files under the given folder (one for each protein)", required=True)
    parser.add_argument("--gpu", help="using gpu, otherwise using cpu", default=True, required=False)
    parser.add_argument("--overwrite", help="redo if npy exists; otherwise skip if exists", type=bool, default=False, required=False)
    parser.add_argument("--maxlen", help="max length of sequence to embed", type=int, default=1000, required=False)
    args = parser.parse_args()

    if not os.path.exists(args.npy):
        os.mkdir(args.npy)

    #print("load model...")
    torch.cuda.empty_cache()
    device = torch.device('cuda:0' if torch.cuda.is_available() and args.gpu else 'cpu')

    # Load ESM-2 model
    model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
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

        # Split sequence if too long and embed individually
        if len(seqseq) > args.maxlen:
            print(f"seq {seqid} too long {len(seqseq)}, splitting")
            overlap = 200
            edata = {'e15': [], 'e21': [], 'ct': []}
            subseqs = split_seq(seqseq, seqid, args.maxlen, overlap)
            for seq in subseqs:
                embed = embed_seq(model, batch_converter, seq[0], seq[1], device)
                for key in edata:  #pylint: disable=C0206
                    edata[key].append(embed[key])
            edata['e15'] = avg_emb(edata['e15'], overlap)
            edata['e21'] = avg_emb(edata['e21'], overlap)
            edata['ct'] = avg_ct(edata['ct'], args.maxlen-overlap)
            np.savez_compressed(filen, s=seqseq, e15=edata['e15'], e21=edata['e21'], ct=edata['ct'])
            continue

        # Add each sequence to list of sequences to embed
        print(f"seq {seqid} {len(seqseq)}")
        edata = embed_seq(model, batch_converter, seqid, seqseq, device)
        np.savez_compressed(filen, s=seqseq, e15=edata['e15'], e21=edata['e21'], ct=edata['ct'])

if __name__ == "__main__":
    main()
