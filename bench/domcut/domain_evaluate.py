import sys
import numpy as np
import argparse
import sys
import math
import statistics

parser = argparse.ArgumentParser()
parser.add_argument("--ref", help="ground truth", required=True)
parser.add_argument("--pred", help="prediction", required=True)
args = parser.parse_args()

seqid, dompred = [], []
in1 = open(args.pred, "r")
for aline in in1:
    if aline[0] == '#':
        continue
    subs = aline.split()
    sid, dom = subs[0], int(subs[1])
    seqid.append(sid)
    dompred.append(dom)
in1.close()

tot = len(seqid)

domref = [-1] * tot

in2 = open(args.ref, "r")
for aline in in2:
    if aline[0] == '#':
        continue
    subs = aline.split()
    sid, dom = subs[0], int(subs[1])
    if sid in seqid:
        domref[seqid.index(sid)] = dom
in2.close()

#report, acc, pre for single
pred_single, pred_multi = 0, 0
#ts: true single, tm: true multi, fs: false single, fm: false multi
ts, tm, fs, fm = 0, 0, 0, 0
for idx in range(tot):
    if dompred[idx] == 1:
        pred_single += 1 
        if domref[idx] == 1:
            ts += 1
    else:
        pred_multi += 1 
        if domref[idx] != 1:
            tm += 1
fs = pred_single - ts
fm = pred_multi - tm

print(f"Single Precision {ts/(ts + fs):0.3f} Recall {ts/(ts + fm):0.3f}")
print(f"Multi Precision {tm/(tm + fm):0.3f} Recall {tm/(tm + fs):0.3f}")
acc = (ts + tm) / (ts + fs + tm + fm)
mcc = (tm * ts - fm * fs) / math.sqrt((tm + fm)*(tm + fs)*(fm + ts)*(ts + fs))
print(f"All ACC {acc:0.3f} MCC {mcc:.3f}")
