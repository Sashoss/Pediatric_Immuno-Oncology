

python ./src/arboreto_with_multiprocessing.py ./out/init_loom_harmony_SCPCP000001_50_2000_3000.loom \
                                               ./in/hs_hgnc_curated_tfs.txt --num_workers 8 \
                                               -o out/adjacency_matrix.csv \
                                               --method grnboost2
