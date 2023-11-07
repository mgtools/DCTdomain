#Inputs: 
#  1. example.fasta -- sequence file
#  2. example.list -- protein list
#  3. example.pair -- pairs of proteins 

#generate DCTdomain fingerprints
bash ../src/fingerprint.sh example.fasta example.list example

#calculate pairwise similarity (for the pairs given by --pair) using DCTdomain fingerprints
python ../src/dct-sim.py --pair example.pair --dct example-dct.npz --output example-dctsim.txt
