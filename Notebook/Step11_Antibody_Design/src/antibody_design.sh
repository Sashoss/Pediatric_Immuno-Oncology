
mkdir -p ./out/seq_evolve_thread
antibody_designer.linuxgccrelease \
        -s ./in/SPP1_antobody_bound.pdb \
        -seq_design_cdrs H1 H2 H3 H4 L1 L2 L3 L4 \
        -nthreads 12 \
        -multithreading:total_threads 8 \
        -light_chain kappa \
        -nstruct 2500 \
        -mintype relax \
        -relax:default_repeats 1 \
        -seq_design_profile_samples 30 \
        -run:timer true \
        -ab_template_db_path ./db \
        -out:prefix spp1_rabd_ \
        -out:file:scorefile ./out/spp1_rabd_score.sc \
        -out:path:all ./out/seq_evolve_thread 

