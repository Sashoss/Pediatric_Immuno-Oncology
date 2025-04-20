pyscenic ctx ./out/adjacency_matrix.csv \
        ./in/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        ./in/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather \
        --annotations_fname ./in/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname ./out/init_loom_harmony_SCPCP000001_50_2000_3000.loom \
        --output ./out/regulons.csv \
        --mask_dropouts --num_workers 8

