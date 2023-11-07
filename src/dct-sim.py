import numpy as np
import argparse
from operator import itemgetter
import time
#Yuzhen Ye, Indiana University, Nov 2023

def prostSimilarity(emb1,emb2):
    #return abs(emb1-emb2).sum()/2 (PROST distance)
    d = abs(emb1-emb2).sum()  #L1-distance
    d /= 17000 #normalize to 0-1 (works for 475 embedding)
    if d > 1:
        d = 1
    return (1-d)

def domain_sim(dct_i, dct_j):
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

def load_dct(filename, asmap=True):
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
    return (dct, seqid)

def main():
    start = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--dct", help="dct in a npz file", required=True)
    parser.add_argument("--output", help="save results to a file", required=False)
    parser.add_argument("--pair", help="calculate distance between the proteins in the given file", required=False)
    parser.add_argument("--pairfound", help="pairs of proteins with similarity computed", required=False)
    parser.add_argument("--db", help="search query dct against this db", required=False)
    parser.add_argument("--top", help="report at most this many hits for database search", default=5, type=int)
    parser.add_argument("--threshold", help="similarity threshold for reporting hits for database search", default=0.4, type=float)
    args = parser.parse_args()

    if args.output:
        out = open(args.output, "w")

    if args.pair:
        dct, seqid = load_dct(args.dct, asmap=True)
    else:
        dct, seqid = load_dct(args.dct, asmap=False)
    nowt = time.time()
    print(f"dct loaded for {len(seqid)} sequences, time used: {nowt - start:.1f}s")

    tot = len(seqid)
    if args.output:
        out.write(f"#prot1 prot2 sim-domain sim-global\n")
    else:
        print(f"#prot1 prot2 sim-domain sim-global")

    tot, totfound = 0, 0
    if args.pair: #given pair
        inf = open(args.pair, "r")
        if args.pairfound:
            out2 = open(args.pairfound, "w")
        for aline in inf:
            if aline[0] == '#':
                if args.pairfound:
                    out2.write(aline)
                continue
            subs = aline.split()
            s1, s2 = subs[0], subs[1]
            tot += 1
            if (s1 in dct) and (s2 in dct):
                (maxs, s) = domain_sim(dct[s1], dct[s2])
                tmp = f"{s1} {s2} {maxs} {s}"
                if args.output:
                    out.write(tmp + "\n")
                else:
                    print(tmp)
                if args.pairfound:
                    out2.write(aline)
                totfound += 1
        print(f"total pair {args.pair} found {totfound} (not found: {tot - totfound})")
        if args.pairfound:
            print(f"pairs saved to file {args.pairfound}")
            out2.close()
    elif args.db: #query against database
        db_dct, db_seqid = load_dct(args.db, asmap=False)
        db_tot = len(db_seqid)
        for i in range(tot):
            results = []
            for q in range(db_tot):
                (maxs, s) = domain_sim(dct[i], db_dct[q])
                results.append([db_seqid[q], s, maxs])
            results_sorted = sorted(results, key=itemgetter(2), reverse=True) 
            for q in range(db_tot):
                if (q >= args.top) and (results_sorted[q][2] < args.threshold):
                    break
                hit = results_sorted[q]
                tmp = f"{seqid[i]} {hit[0]} {hit[1]} {hit[2]}"
                if args.output:
                    out.write(tmp + "\n")
                else:
                    print(tmp)
    else: #all againt all
        for i in range(tot - 1):
            for j in range(i + 1, tot):
                (maxs, s) = domain_sim(dct[i], dct[j])
                #sim = similarity_level(ann, seqid[i], seqid[j])
                tmp = f"{seqid[i]} {seqid[j]} {maxs:.3f} {s:.3f}"
                if args.output:
                    out.write(tmp + "\n")
                else:
                    print(tmp)
    if args.output:
        print("results saved to", args.output)
        out.close()

    end = time.time()
    print(f"total time used {end - start:.1f}s")
    print(f"distance calculation used {end - nowt:.1f}s")

if __name__ == '__main__':
    main()
