#!/bin/bash

# parallel --keep-order --tag ./cstoy.py --overwrite --py-pthatmin 80 --nev 5000 --py-seed {} ::: 22 33 44 55 66 77
# hadd -f output_alpha_0_dRmax_0.0_SDzcut_0.1_seed_merged.root output_alpha_0_dRmax_0.0_SDzcut_0.1_seed_??.root

# parallel --keep-order --tag ./cstoy.py --overwrite --py-pthatmin 80 --nev 5000 --dRmax 0.25 --py-seed {} ::: 22 33 44 55 66 77
# hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_??.root

# parallel --keep-order --tag ./cstoy.py --overwrite --py-pthatmin 80 --nev 5000 --zcut 0.2 --dRmax 0.25 --py-seed {} ::: 22 33 44 55 66 77
# hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.2_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.2_seed_??.root

# embs="--embed /Users/ploskon/data/alice/LHC18qr/local_PbPb_file_list.txt"
# parallel --keep-order --tag ./cstoy.py --overwrite --py-pthatmin 80 --nev 10000 --zcut 0.2 --dRmax 0.25 --py-seed {} ${embs} ::: 22 33 44 55 66 77 88
# hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.2_emb_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.2_seed_??_emb.root

nev=2000
embs="--embed /Users/ploskon/data/alice/LHC18qr/local_PbPb_file_list.txt"
#parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.2 --dRmax 0.25 --py-seed {} ${embs} --efficiency ::: 22 33 44 55 66 77 88
#hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.2_emb_effi_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.2_seed_??_emb_effi.root

parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.2 --dRmax 0.25 --py-seed {} ${embs} --efficiency --SDsignal ::: 22 33 44 55 66 77 88
hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.2_emb_effi_SDsignal_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.2_seed_??_emb_effi_SDsignal.root

#parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.1 --dRmax 0.25 --py-seed {} ${embs} --efficiency ::: 22 33 44 55 66 77 88
#hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.1_emb_effi_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_??_emb_effi.root

parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.1 --dRmax 0.25 --py-seed {} ${embs} --efficiency --SDsignal ::: 22 33 44 55 66 77 88
hadd -f output_alpha_0_dRmax_0.25_SDzcut_0.1_emb_effi_SDsignal_seed_merged.root output_alpha_0_dRmax_0.25_SDzcut_0.1_seed_??_emb_effi_SDsignal.root

#parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.1 --py-seed {} ${embs} --efficiency --csjet ::: 22 33 44 55 66 77 88
#hadd -f output_alpha_0_dRmax_0.0_SDzcut_0.1_emb_effi_csjet_seed_merged.root output_alpha_0_dRmax_0.0_SDzcut_0.1_seed_??_emb_effi_csjet.root

parallel --joblog run_cstoy.log --keep-order --tag $PYJETTYDIR/pyjetty/cstoy/cstoy.py --py-pthatmin 80 --nev ${nev} --zcut 0.1 --py-seed {} ${embs} --efficiency --csjet --SDsignal ::: 22 33 44 55 66 77 88
hadd -f output_alpha_0_dRmax_0.0_SDzcut_0.1_emb_effi_SDsignal_csjet_merged.root output_alpha_0_dRmax_0.0_SDzcut_0.1_seed_??_emb_effi_SDsignal_csjet.root
