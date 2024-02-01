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

def readfasta(filename: str, batchsize: int) -> list:
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

def split_seq(seqseq: str, seqid: str, idx: int, maxlen: int) -> list:
    """Returns a list of (seqid, subseq) tuples, where subseq is a substring of seqseq[idx] of
    length maxlen.

    Args:
        seqseq (str): sequence
        seqid (str): sequence id
        idx (int): index of sequence in seqseq
        maxlen (int): maximum length of subseq

    Returns:
        list of (seqid, subseq) tuples
    """

    embed = []
    for i in range(0, len(seqseq[idx]), maxlen):
        if i == 0:
            subseq = seqseq[idx][i:i+maxlen]
        else:
            subseq = seqseq[idx][i-200:i+maxlen-200]
        if len(subseq) > 200:
            embed.append((seqid[idx], subseq))
    return embed

def split_emb(embed: list, batch_converter: esm.Alphabet, device: str, model: esm.pretrained):
    """Returns a dictionary of embeddings and contacts for each subsequence in a list of
    (seqid, subseq) tuples.

    Args:
        embed (list): list of (seqid, subseq) tuples
        batch_converter (esm.Alphabet): ESM alphabet
        device (str): device to use for embedding
        model (esm.pretrained): ESM model

    Returns:
        dict: dictionary where key is type of embedding and value is a list of overlapping
        embeddings or contacts
    """

    edata = {'e13': [], 'e25': [], 'ct': []}
    for emb in embed:
        _, _, batch_tokens = batch_converter([emb])  # batch_labels, batch_strs
        batch_tokens = batch_tokens.to(device)
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[13, 25], return_contacts=True)
        edata['e13'].append(results["representations"][13].cpu().numpy())
        edata['e25'].append(results["representations"][25].cpu().numpy())
        edata['ct'].append(results["contacts"].cpu().numpy())
    return edata

def avg_emb(embeds: list) -> np.ndarray:
    """Returns a single embedding by averaging and concatenating a list of embeddings.

    Args:
        embeds (list): list of (seqid, np.ndarray) tuples
    Returns:
        numpy array: single 1D embedding
    """

    full_emb = np.array([])
    for emb in embeds:
        emb = emb[0][1:-1]
        if len(full_emb) == 0:
            full_emb = emb
            continue
        avg = (full_emb[-200:] + emb[:200]) / 2  # avg adjacent positions
        full_emb[-200:] = avg  # replace last 200 positions of previous emb with avg
        full_emb = np.concatenate((full_emb, emb[200:]), axis=0)  # concat rest of new emb
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

def avg_ct(contacts: list) -> np.ndarray:
    """Returns a single contact map by averaging and combining a list of contact maps.

    Args:
        contacts (list): list of (seqid, np.ndarray) tuples

    Returns:
        numpy array: single 2D contact map
    """

    full_contacts = np.array([])
    for i, cont in enumerate(contacts):
        cont = cont[0]
        if len(full_contacts) == 0:  # initialize full array
            full_contacts = cont
            continue
        full_contacts = combine_contacts(full_contacts, cont, inc = 800, times = i)
    return full_contacts

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
            results = model(batch_tokens, repr_layers=[13, 25], return_contacts=True)
    except RuntimeError:
        print("RuntimeError occured, try smaller --maxgpu")
        return
    r13a = results["representations"][13].cpu()
    r25a = results["representations"][25].cpu()
    cta = results["contacts"].cpu()

    # Save data into npz file
    np.set_printoptions(threshold=sys.maxsize)
    for (sid, seq), e13, e25, ct in zip(data, r13a, r25a, cta):
        filen = npyfile + "/" + sid
        np.savez_compressed(filen, s=seq, e13=e13, e25=e25, ct=ct)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="FASTA file", required=True)
    parser.add_argument("--npy", help="save embeddings and contact to files under the given folder (one for each protein)", required=True)
    parser.add_argument("--gpu", help="using gpu, otherwise using cpu", default=False, required=False)
    parser.add_argument("--batchsize", help="batchsize; default 1", default=1, required=False)
    parser.add_argument("--overwrite", help="redo if npy exists; otherwise skip if exists", type=bool, default=False, required=False)
    parser.add_argument("--maxlen", help="max length of sequence to embed", type=int, default=1000, required=False)
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
    model.to(device)

    # Get filename for .list file
    listfile = args.input.split('.')[0] + ".list"
    with open(listfile, "w", encoding='utf8') as listf:
        listf.write('#protein\n')

    fasta = readfasta(args.input, args.batchsize)
    while True:

        # Read batch of sequences
        try:
            seqid, seqseq = next(fasta)
        except StopIteration:
            break

        # For each sequence in batch, prepare inputs for embedding
        data = []
        for idx in range(len(seqid)):  #pylint: disable=C0200
            filen = args.npy + "/" + seqid[idx] + ".npz"
            with open(listfile, "a", encoding='utf8') as listf:  # add to .list file
                listf.write(f"{seqid[idx]}\n")
            if (not args.overwrite) and os.path.exists(filen):  # skip if file exists
                print(f"seq {seqid[idx]} already exists")
                continue

            # Split sequence if too long and embed individually
            if len(seqseq[idx]) > args.maxlen:
                print(f"seq {seqid[idx]} too long {len(seqseq[idx])}, splitting")
                embed = split_seq(seqseq, seqid, idx, args.maxlen)
                edata = split_emb(embed, batch_converter, device, model)
                edata['e13'] = avg_emb(edata['e13'])
                edata['e25'] = avg_emb(edata['e25'])
                edata['ct'] = avg_ct(edata['ct'])
                np.savez_compressed(filen, s=seqseq[idx], e13=edata['e13'], e25=edata['e25'], ct=edata['ct'])
                continue

            # Add each sequence to list of sequences to embed
            print(f"seq {seqid[idx]} {len(seqseq[idx])}")
            data.append((seqid[idx], seqseq[idx]))
            _, _, batch_tokens = batch_converter(data)  # batch_labels, batch_strs
            batch_tokens = batch_tokens.to(device)
        if len(data) == 0:  # skip if no inputs were processed
            continue

        # Embed batch of sequences
        embed_batch(model, batch_tokens, data, args.npy)

if __name__ == "__main__":
    main()
