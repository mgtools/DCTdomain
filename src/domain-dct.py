"""Gets DCT fingerprints for a directory of protein embeddings and predicted domains.

Yuzhen Ye, Indiana University, Nov 2023
"""

import argparse
import os
import numpy as np
from scipy.fftpack import dct, idct

def iDCTquant(v,n):
    """From prosttools. Returns the iDCT quantization of a vector.
    """
    f = dct(v.T, type=2, norm='ortho')
    trans = idct(f[:,:n], type=2, norm='ortho')
    for i in range(len(trans)):
        trans[i] = scale(trans[i])
    return trans.T

def scale(v: np.ndarray):
    """From prosttools. Returns a scaled vector.

    Args:
        v (np.ndarray): vector
    """
    M = np.max(v)
    m = np.min(v)
    return (v - m) / float(M - m)

def quant2D(emb: np.ndarray,n=5,m=44) -> np.ndarray:
    """From prosttools. Returns the iDCT quantization of a 2D embedding.

    Args:
        emb (np.ndarray): 2D embedding
        n (int): first dimension size
        m (int): second dimension size

    Returns:
        1D iDCT quantization of a 2D embedding
    """

    dct = iDCTquant(emb,n)
    ddct = iDCTquant(dct.T,m).T
    ddct = ddct.reshape(n*m)
    return (ddct*127).astype('int8')

def get_domains(domstr: str, length: int) -> list:
    """Returns a list of domains from a string.

    Args:
        domstr (str): protein domains ex. "1-64;65-128"
        length (int): length of the protein

    Returns:
        list of domains
    """

    # Get each domain
    if not domstr:
        doms = []
    else:
        if domstr[-1] == ';':
            domstr = domstr[:-1]
        doms = domstr.split(";")
    if len(doms) != 1: #no domain given or multiple domain
        doms.append(f"1-{length}") #produce whole protein for multi domain proteins

    return doms

def domainEmb_help(a_list: list, s_list: list) -> np.ndarray:
    """Using multiple layers (one layer is just a special case).

    Args:
        a_list (list): Part of embedding selected for a domain
        s_list (list): Dimension sizes for DCT's

    Returns:
        DCT fingerprint for a domain
    """

    q = []
    for idx in range(len(a_list)):  #pylint: disable=C0200)
        a1 = a_list[idx]
        s1 = s_list[idx]
        q1 = quant2D(a1, s1[0], s1[1])
        if idx == 0:
            q = q1
        else:
            q = np.concatenate([q, q1])
    return q

def domainEmb(a_list: list, s_list: list, domstr: str) -> tuple:
    """Returns domains and DCTs in arrays.

    Args:
        a_list (list): list of embeddings
        s_list (list): list of dimension sizes for DCT's
        domstr (str): domain information
    
    Returns:
        domains and DCTs in arrays
    """

    num_layer = len(a_list)
    doms = get_domains(domstr, len(a_list[0]))

    # Get each domain's fingerprint
    fprints = []
    for one in doms:
        segs = one.split(",")
        embsel = []
        for seg in segs:
            subs = seg.split("-")
            beg, end = int(subs[0])-1, int(subs[1])  #python index starts from 0
            if not embsel:
                for l in range(num_layer):
                    embsel.append(a_list[l][beg:end, :])
            else:
                for l in range(num_layer):
                    #discontinous domains, embeddings of the segments will be merged
                    embsel[l] = np.concatenate((embsel[l], a_list[l][beg:end, :]), axis=0)

        fprints.append(domainEmb_help(embsel, s_list))

    return doms, fprints

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--npy", help="folders with embedding files", required=True)
    parser.add_argument("--domain", help="a list of proteins with/without domain information", required=True)
    parser.add_argument("--output", help="save domain based embedding to a file", required=True)
    args = parser.parse_args()

    inf = open(args.domain, "r")
    tosave_dct = []
    tosave_dom = []
    tosave_seqid =[]
    tosave_idx = []
    add = 0
    for aline in inf:
        if aline[0] == '#':
            continue
        subs = aline.strip().split()
        if len(subs) >= 3:
            (seqid, _, domstr) = subs[0:3]  # domsize is not used
        else:
            seqid, _, domstr = subs[0], 1, ""
        add += 1
        print(f"generate DCT fingerprint for protein {add} {seqid}")
        filen = args.npy + "/" + seqid + ".npz"
        if os.path.exists(filen):
            data = np.load(filen)
            (dom, dcta) = domainEmb([data['e15'], data['e21']], [[3, 80], [3, 80]], domstr)
            tosave_seqid.append(seqid)
            tosave_idx.append(len(tosave_dom))
            tosave_dom.extend(dom)
            tosave_dct.extend(dcta)
        else:
            print("npy not found; skip")
    tosave_idx.append(len(tosave_dom)) #add one more
    inf.close()
    np.savez(args.output, sid=tosave_seqid, idx=tosave_idx, dom=tosave_dom, dct=tosave_dct)

if __name__ == "__main__":
    main()
