# Run script to generate similarity values between all pairs in the corresponding benchmark datasets
# Make sure you are in the directory 'bench/homo' to run this script

#pfam-max50 bench
python ../../src/dct-sim.py --dct pfam_max50-dct.npz --pair pfam_max50.pair --output pfam_max50-dctsim.txt

#pfam-nomax50 bench
python ../../src/dct-sim.py --dct pfam_nomax50-dct.npz --pair pfam_nomax50.pair --output pfam_nomax50-dctsim.txt

#pfam-local bench
python ../../src/dct-sim.py --dct pfam_nomax50-dct.npz --pair pfam_localpfam_nomax50.pair --output pfam_localpfam_nomax50-dctsim.txt

#superfamily-nomax50 bench
python ../../src/dct-sim.py --dct supfam_nomax50-dct.npz --pair supfam_nomax50.pair --output supfam_nomax50-dctsim.txt

#gene3d-nomax50 bench
python ../../src/dct-sim.py --dct gene3d_nomax50-dct.npz --pair gene3d_nomax50.pair --output gene3d_nomax50-dctsim.txt
