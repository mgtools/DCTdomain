#python ../../src/embed-mult.py --input test.fasta --npy npy
#npy files are large and they can be removed after contact maps are extracted 

#python ../../src/embed2ce.py --input test.fasta --npy npy --ce CE

bash ../../src/domaincut.sh <test.info > test.dom

#test.info -- containing "ground-truth" domain segmentations (from SCOP)
#test.dom -- predicted domains from RecCut
python ./domain_evaluate.py --ref test.info --pred test.dom
