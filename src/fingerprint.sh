#!/bin/bash
if test "$#" -ne 2; then
    echo "Usage: $0 fasta-file output-prefix"
    exit
fi
start=$(date +%s)
temp=$( realpath "$0" )
d="$(dirname "$0")"

#Calculate ESM-2 embedding
echo "step 1: compute embeddings.."
python $d/embed-mult.py --input $1 --npy NPY

#Extract contact maps
echo "step 2: extract contact maps.."
python $d/embed2ce.py --input $1 --npy NPY --ce CE

#Domain segmentation
echo "step 3: compute domains.."
listfile="${1%.*}.list"
bash $d/domaincut.sh < $listfile > $2.dom

#Prepare DCTdomain fingerprints (saved to *.npz)
echo "step 4: generate DCT fingerprints.."
python $d/domain-dct.py --npy NPY --domain $2.dom --output $2-dct.npz

end=$(date +%s)
echo "#Total time used: $(($end-$start)) seconds"
