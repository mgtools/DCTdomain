#Inputs: G6PD.fasta, G6PD.pair
#Calculate ESM-2 embedding
python ../../src/embed-mult.py --input G6PD.fasta --npy NPY

#Extract contact maps
python ../../src/embed2ce.py --input G6PD.fasta --npy NPY --ce CE

#Domain segmentation
bash ../../src/domaincut.sh < G6PD.list > G6PD.dom

#Prepare DCTdomain fingerprints (saved to G6PD-dct.npz)
python ../../src/domain-dct.py --npy NPY --domain G6PD.dom --output G6PD-dct.npz

#Compute similarity scores using DCTdomain fingerprints (save the results to G6PD-dctsim.txt
python ../../src/dct-sim.py --pair G6PD.pair --dct G6PD-dct.npz --output G6PD-dctsim.txt

#to prepare plots showing the comparison of the score distribions (as those in the paper)
#python hist.py
