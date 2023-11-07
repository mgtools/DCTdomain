import numpy as np
import argparse
import sys
from scipy.fftpack import dct, idct
import os
#Yuzhen Ye, Indiana University, Nov 2023

#ref protsttools for iDCTquant
def iDCTquant(v,n):
    f = dct(v.T, type=2, norm='ortho')
    trans = idct(f[:,:n], type=2, norm='ortho')
    for i in range(len(trans)):
        trans[i] = scale(trans[i])
    return trans.T

def scale(v):
    M = np.max(v)
    m = np.min(v)
    return (v - m) / float(M - m)

def quant2D(emb,n=5,m=44):
    dct = iDCTquant(emb[1:len(emb)-1],n)
    ddct = iDCTquant(dct.T,m).T
    ddct = ddct.reshape(n*m)
    return (ddct*127).astype('int8')

#using multiple layers (one layer is just a special case)
def domainEmb_help(a_list, s_list):
    q = []
    for idx in range(len(a_list)):
        a1 = a_list[idx]
        s1 = s_list[idx]
        #a1: excluding begining and ending positions
        q1 = quant2D(a1[:, 1:-1], s1[0], s1[1])  
        if idx == 0:
            q = q1
        else:
            q = np.concatenate([q, q1])
    return q

#return domains and DCTs in arrays
def domainEmb(a_list, s_list, domstr):
    layer = len(a_list)
    if not domstr:
        doms = []
    else:
        if domstr[-1] == ';':
            domstr = domstr[:-1]
        doms = domstr.split(";")
    if len(doms) != 1: #no domain given or multiple domain
        doms.append(f"1-{len(a_list[0]) - 2}") #produce whole protein for multi domain proteins
    q_list = []
    for one in doms:
        segs = one.split(",")
        embsel = []
        for seg in segs:
            subs = seg.split("-")
            #beg end start at index 1 (not 0); same as the embedding which has the beginning of the sequence as the first token
            beg, end = int(subs[0]), int(subs[1]) + 1
            if not embsel:
                for l in range(layer):
                    embsel.append(a_list[l][beg:end, :])
            else:
                for l in range(layer):
                    #discontinous domains, embeddings of the segments will be merged
                    embsel[l] = np.concatenate((embsel[l], a_list[l][beg:end, :]), axis=0)

        q_list.append(domainEmb_help(embsel, s_list))

    return doms, q_list

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
        if len(subs) == 3:
            (seqid, domsize, domstr) = subs[0:3] 
        else:
            seqid, domsize, domstr = subs[0], 1, "" 
        add += 1
        print(f"generate DCT fingerprint for protein {add} {seqid}")
        filen = args.npy + "/" + seqid + ".npz"
        if os.path.exists(filen):
            data = np.load(filen)
            #layer 14: 3x85, layer 26, 5x44
            #print(f"embedding shape {data['e25'].shape}")
            (dom, dcta) = domainEmb([data['e13'], data['e25']], [[3, 85], [5, 44]], domstr)
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
