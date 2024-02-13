"""Calculates similarity between proteins based on their DCT fingerprints.

Yuzhen Ye, Indiana University, Nov 2023
"""

import numpy as np
import argparse
from operator import itemgetter
import time


def prostSimilarity(emb1: np.ndarray, emb2: np.ndarray) -> float:
    """Returns similarity between two fingerprints using PROST distance.

    Args:
        emb1, emb2: two fingerprints, each is a 475-dimensional vector.

    Returns:
        A float between 0 and 1, where 0 means no similarity and 1 means
        identical fingerprints.
    """

    d = abs(emb1-emb2).sum()  #L1-distance
    d /= 17000 #normalize to 0-1 (works for 475 embedding)
    d=min(d, 1)
    return 1-d

def domain_sim(dct_i: np.ndarray, dct_j: np.ndarray) -> tuple:
    """Returns the max similarity between any two fingerprints in the two sets of fingerprints and
    the similarity between the global fingerprints.

    Args:
        dct_i, dct_j: two sets of fingerprints, each is a n x 475 matrix where n is the number of
        RecCut predicted domains in the protein (+1 for the whole protein).

    Returns:
        A tuple (maxs, s), where maxs is the max similarity between any two fingerprints
        (DCTdomain), and s is the similarity between the two global fingerprints (DCTglobal).
    """

    #print(f"dct_i {dct_i}, dct_j {dct_j}")
    di, s = dct_i.shape
    dj, s = dct_j.shape
    maxs = 0
    for pi in range(di):
        for pj in range(dj):
            s = prostSimilarity(dct_i[pi], dct_j[pj])
            if s > maxs:
                maxs = s
    return (maxs, s)

def load_dct(filename: str, asmap=True) -> tuple:
    """Returns a dictionary or a list of fingerprints from a npz file.

    Args:
        filename: the name of the npz file.
        asmap: if True, returns a dictionary mapping sequence ids to fingerprints;
        otherwise, returns a list of fingerprints.

    Returns:
        A tuple (dct, seqid), where dct is a dictionary or a list of fingerprints, and
        seqid is a list of sequence ids.
    """

    start = time.time()
    data_all0 = np.load(filename)
    seqid = data_all0['sid']
    domidx = data_all0['idx']
    dct_all = data_all0['dct']
    tot = len(seqid)
    if asmap:
        dct = {}
    else:
        dct = []
    for i in range(tot):
        s, e = domidx[i], domidx[i + 1]
        ai = dct_all[s:e, :]
        if asmap:
            dct[seqid[i]] = ai
        else:
            dct.append(ai)
    nowt = time.time()
    print(f"dct loaded for {len(seqid)} sequences, time used: {nowt - start:.1f}s")
    return (dct, seqid)

def pair_sim(npzfile: str, pairfile: str, pairfound: str, output: str):
    """Writes similarity between each protein pair to stdout or a file.

    Args:
        npzfile (str): file containing fingerprints.
        pairfile (str): file containing protein pairs.
        pairfound (str): file to save protein pairs with fingerprints.
        output (str): file to save similarity between each protein pair.
    """

    dct, _ = load_dct(npzfile, asmap=True)
    tot, totfound = 0, 0
    inf = open(pairfile, "r", encoding='utf8')
    if pairfound:
        out2 = open(pairfound, "w", encoding='utf8')
    for aline in inf:
        if aline[0] == '#':
            if pairfound:
                out2.write(aline)
            continue
        subs = aline.split()
        s1, s2 = subs[0], subs[1]
        tot += 1
        if (s1 in dct) and (s2 in dct):
            (maxs, s) = domain_sim(dct[s1], dct[s2])
            tmp = f"{s1} {s2} {maxs} {s}"
            if output:
                with open(output, "a", encoding='utf8') as out:
                    out.write(f'{tmp}\n')
            else:
                print(tmp)
            if pairfound:
                out2.write(aline)
            totfound += 1
    print(f"total pair {pairfile} found {totfound} (not found: {tot - totfound})")
    inf.close()
    if pairfound:
        print(f"pairs saved to file {pairfound}")
        out2.close()

def db_search(npzfile: str, dbfile: str, top: int, threshold: float, output: str):
    """Writes top hits for each protein in the query file to stdout or a file.

    Args:
        npzfile (str): file containing fingerprints.
        dbfile (str): file containing fingerprints of proteins in the database.
        top (int): report at most this many hits for each query protein.
        threshold (float): similarity threshold for reporting hits.
        output (str): file to save hits.
    """

    dct, seqid = load_dct(npzfile, asmap=False)
    totseq = len(seqid)
    db_dct, db_seqid = load_dct(dbfile, asmap=False)
    db_tot = len(db_seqid)
    for i in range(totseq):
        results = []
        for q in range(db_tot):
            (maxs, s) = domain_sim(dct[i], db_dct[q])
            results.append([db_seqid[q], maxs, s])
        results_sorted = sorted(results, key=itemgetter(2), reverse=True)
        for q in range(db_tot):
            if (q >= top) and (results_sorted[q][2] < threshold):
                break
            hit = results_sorted[q]
            tmp = f"{seqid[i]} {hit[0]} {hit[1]} {hit[2]}"
            if output:
                with open(output, "a", encoding='utf8') as out:
                    out.write(tmp + "\n")
            else:
                print(tmp)

def all_sim(npzfile: str, output: str):
    """Writes similarity between all protein pairs to stdout or a file.

    Args:
        npzfile (str): file containing fingerprints.
        output (str): file to save similarity between each protein pair.
    """

    dct, seqid = load_dct(npzfile, asmap=False)
    totseq = len(seqid)
    for i in range(totseq - 1):
        for j in range(i + 1, totseq):
            (maxs, s) = domain_sim(dct[i], dct[j])
            tmp = f"{seqid[i]} {seqid[j]} {maxs:.3f} {s:.3f}"
            if output:
                with open(output, "a", encoding='utf8') as out:
                    out.write(tmp + "\n")
            else:
                print(tmp)


def main():
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--dct", help="dct in a npz file", required=True)
    parser.add_argument("--output", help="save results to a file", required=False)
    parser.add_argument("--pair", help="calculate distance between the proteins in the given file", required=False)
    parser.add_argument("--pairfound", help="pairs of proteins with similarity computed", required=False)
    parser.add_argument("--db", help="search query dct against this db", required=False)
    parser.add_argument("--top", help="report at most this many hits for database search", default=5, type=int)
    parser.add_argument("--threshold", help="similarity threshold for reporting hits for database search", default=0.25, type=float)
    args = parser.parse_args()

    if args.output:
        out = open(args.output, "w", encoding='utf8')
        out.write("#prot1 prot2 sim-domain sim-global\n")
        out.close()
    else:
        print("#prot1 prot2 sim-domain sim-global")

    nowt = time.time()
    if args.pair: #given pair
        pair_sim(args.dct, args.pair, args.pairfound, args.output)
    elif args.db: #query against database
        db_search(args.dct, args.db, args.top, args.threshold, args.output)
    else: #all againt all
        all_sim(args.dct, args.output)
    if args.output:
        print("results saved to", args.output)
        out.close()

    end = time.time()
    print(f"total time used {end - start:.1f}s")
    print(f"distance calculation used {end - nowt:.1f}s")

if __name__ == '__main__':
    main()
