# Pediatric_glioblastoma_singelcell


## Things to do
1. Run antibody modelling


# Antibody modelling
1. Dock SPP1 resi 120-14 to anti-osteopontin antibody 23C3 (PDB ID: 3CXD)
2. Identify 23C3 heavy and light chains using [PyIgClassify2](https://dunbrack2.fccc.edu/PyIgClassify2/Search/pdb.aspx). Make sure to apply for PyIgClassify2 database download request by filling their form. Each DB access is granted to a single user. Therefore, dont share it with others. Carefully read its usage policy. 
3. Run antibody designer 
antibody_designer.linuxgccrelease \
        -s ./in/SPP1_antobody_bound.pdb \
        -seq_design_cdrs H1 H2 H3 H4 L1 L2 L3 L4 \
        -light_chain kappa \
        -nstruct 2500 \
        -mintype relax \
        -relax:default_repeats 2 \
        -seq_design_profile_samples 30 \
        -run:timer true \
        -ab_template_db_path ./db \
        -out:prefix spp1_rabd_$i \
        -out:file:scorefile ./out/spp1_rabd_score_$i.sc \
        -out:path:all ./out/seq_evolve_$i &
 

# Results

## Fig 1. 
- a) Unsupervised clustering and UMAP projection identified eleven major cell populations, including malignant glial‐lineage states (MES-like, MES-AC-like, MES-AC-like Cycling, MES-APC-like, OPC-like), lymphoid cells (T cells, NK cells, Naïve B), and myeloid cells segregating into two macrophage subsets (MGD Macrophage, MGD TAM) plus a Microglia‐derived TAM cluster (MD Macrophage). A small “Undetermined” cluster likely represents doublets or rare stromal elements.

- b) Infercnv results
We next applied inferCNV to detect large-scale chromosomal copy‐number alterations at single‐cell resolution, using T cells, NK cells, and Naïve B cells as a diploid reference. The resulting CNV heatmap (Fig. 1b) revealed that all malignant glial clusters display characteristic alterations—such as chromosome 7 gain and chromosome 10 loss—while immune clusters show flat profiles consistent with a normal karyotype. Among tumor states, the MES-like and MES-AC-like Cycling populations exhibited the most pronounced copy gains (for example, focal amplification on 7p), underscoring their heightened genomic instability relative to OPC- and APC-like cells.

- c) Marker dotplot
Cell type marker expression dotplot. 

- d) GSEA normalized enrichment score dotplot

   - d-1) Immune cells - 

      - d-1.2) T cells show strong signatures for T cell receptor complex, MHC class II protein complex, and immunological synapse. NK cells are uniquely enriched in cytotoxic granule lumen, perforin complex, and natural killer cell–mediated cytotoxicity pathways, consistent with their tumor‐killing potential.  

      - d-1.4) Homeostatic/transitioning macrophages - Summary: These Macrophages are setting up their secretory functions to modulate the immune microenvironment—either in support or suppression of immune cells, which indicates that these cells are in highly plastic state, typical of macrophages in transition. 

         - d-1,4.1) Tertiary granule membrane, Azurophil granule membrane - Suggests mobilization and trafficking of immune granules, which may indicate that large number of macrophage cells have transitioned into M2 Macrophage. 
         - d-1.4.2) Golgi stack, golgi cisterna membrane, golgi apparatus subcompartment - Indicating high protein processing, which is essential during transition phase. 
         - d-1-4.3) Ficolin-1 Rich Granule Membrane - Ficolins are PRRs (pattern recognition receptors) involved in lectin pathway complement activation, which has a dual role in cancer. Its activation can lead to downstream pathway activation either resulting in the generation of C3a and C5a anaphylatoxins, which can suppress immune cells and promote tumor progression, or the formation of the membrane attack complex (MAC), which can directly lyse cancer cells.

      - d-1.5) Microglia-derived TAMs - Summary: Enriched pathway indicate following properties - Highly metabolically and translationally active, antigen presentation, actively remodeling their local environment via secretion and phagocytosis (key to cell-cell communication).

         - d-1.5.1) Ribosome subunit, cytosolic ribosome, large/small subunit pathways - indicating high protein synthesis activity, which is often seen in activated or metabolically active cells. This suggests that the MGD TAMs are actively producing proteins, possibly for cytokine secretion, antigen presentation, or immune modulation. They are definitly active and not in a quiescent state.
         - d.1.5.2) Primary lysosome, lytic vacuole pathways - Indicates active processes like autophagy, endocytosis, and phagocytosis, which are all important for maintaining cellular homeostasis and responding to various stimuli within the tumor microenvironment.
         - d.1.5.3) MHC class II protein complex - Classic marker of professional antigen-presenting cells (APCs). These TAMs may be presenting antigens to CD4⁺ T cells, which indicates presence of CD4⁺ T cells in the TME.
 
   - d-2) Tumor cells
      - d-2.1) MES-APC-like cells - Summary: Invasive, undergoing transition from APC to MES cell type. 
         - d-2.1.1) Combined enrichment of mitochondrial inner membrane, mitochondrial envelope, intermembrane space, membrane, NADH dehydrogenase, respiratory chain complex I, Proton transporting ATP synthase complex, Cytochrome complex, oxidoreductase complex indicates that these MES-like cells are undergoing high oxidative phosphorylation (OXPHOS) to meet elevated energy demands to meet likely high more invasive characteristics, which is common in glioblastoma. 
         - d-2.1.2) These pathway - Ribosome, mitochondrial ribosome, ribonucleoprotein complex, rough endoplasmic reticulum, peptidase complex, peribosome, Nucleolus and small subunit processome are collectively a translation powerhouse of cells, suggesting that these MES-like cells are producing large amounts of protein, to meet the resources required for their high invasive and proliferative behavior.
         - d-2.1.3) Activation of the above pathways is further supported by the activation of U2 type spliceosomal complex, U12 type spliceosomal complex, spliceosomal snRNP complex, pre-catalytic spliceosome, U4/U6.U5 tri-snRNP complex, small nuclear ribonucleoprotein complex, and Sm-like protein family complex which are associated with high level of RNA processing and alternative splicing. These MES-cells are likely fine-tuning gene expression by inducing alternate splicing events. This raises a question of how much does alternate splicing control cell type transition these glioma's primarily via alternate splicing mechinery, which could make it a strong target. 

      - d-2.2) MES-AC-like: This is one of the most intriguing cell types we have among all. This cluster does not express any classical neuron or neuron progenitor cell markers, and clearly expresses both MES-like and AC-like cell type markers. Although, it does express multiple canonical neuronal pathways. This indicates that these transitioning cells are likely partially acquiring neuronal mimicry.

      - d-2.3) MES-AC-like cyclcing: 
         - d-2.3.1) Condensed chromosome, spindle, spindle pole, mitotic spindle, spindle midzone, kinetochore, centromere, replication fork, midbody, cleavage furrow, centriole - all are hallmarks of active mitosis, indicating that these cells are actively dividing and we caught them mid-cycle — probably in S/G2/M phases.
         - d-2.3.2) Activation of kinesin complex, microtubule cytoskeleton, supramolecular fiber, microtubule organizing center are associated with chromosome segregation and cytokinesis during cell division. 
         - d-2.3.3) Activation of nuclear chromosome, nuclear replication fork, heterochromatin, DNA repair complex indicates S-phase DNA replication activity and DNA damage sensing/repair, common in rapidly proliferating cells under stress, such as tumor cells.


<img src="Manuscript/Fig1.png">


## Fig 2.
- a) Outgoing signal heatmap
- b) Incoming signal heatmap
- c) Outgoing signal dotplot
- d) Incoming signal dotplot
- e) Polarized microglia derived tumor associated macrophage outgoing signals
- f) SPP1 pathway signal communication heatmap

<img src="Manuscript/Fig2.png">

## Fig 3.
- a) TF Regulon activity
- b) TF expression
- c) MAFB feature plot 
- d) MAFB violin plot

<img src="Manuscript/Fig3.png">

## Fig 4.
- a) Root mean squared fluctuation of SPP1 protein residues under replica exchange molecular dynamics simulation at temperatures 283.15K, 303.15K, 333.15K, and 353.15K 
- b) Energy landscape based conformation sampling from replicate exchange runs
- c) Most stable conformation across all 4 replica exchange trajectories. CD44 binding interface (residue 121-140) is shown in orange, other known SPP1 antibody binding interface shown in blue, and cyan.
- d) SPP1 protein residue phosphorylation heatmap. Arrow highlights the top likely phosphorylation site 169 with score > 0.8.
- e) Sequence alignment of SPP1 across various mammals. Orange bar highlights SPP1-CD44 binding interface (residue 121-140), other known SPP1 antibody binding interface shown in blue, and cyan. Red arrow shows phosphorylation residue 169.
- f) Root mean squared fluctuation quadratic mean of normal SPP1 (black) and residue 169 phosphorylated SPP1 (red), showing increase in stability in SPP1-CD44 binding interface upon phosphorylation.
<img src="Manuscript/Fig4.png">

## Fig 5.
- a) UMAP showing 2500 new antibodies generated using Rosetta antibody design tool and the native 23C3 antibody. Color of the dots highlight free energy required to seperate SPP1 from the corresponding antibody variant. UMAP is calculated using cosine distance between sequence embedding of CDR sequence. 
- b) Hexplot showing 2500 new antibodies free energy required to seperate SPP1 from the corresponding antibody variant and distance between complete CDR sequence.
- c) Most stable antibody protein structure (grey) bonded to SPP1 protein structure (wheat). Blue spheres are the residue changes between the most stable Antibody and 23C3 in light chain and red spheres are the residue changes between the most stable Antibody and 23C3 in heavy chain.
<img src="Manuscript/Fig5.png">




### Cuurated molecules that can cross blood brain barrier
- Manuscript - [https://www.nature.com/articles/s41597-021-01069-5](https://www.nature.com/articles/s41597-021-01069-5)
- Github - [https://github.com/theochem/B3DB/tree/main](https://github.com/theochem/B3DB/tree/main)