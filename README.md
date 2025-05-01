# Project long-term goal:
1. Collect and process publicly available single cell, bulk transcriptomic, and structural bioinformatics data to explore, organize, and annotate key pathways and cytokines associated with immune suppression in Pediatric High Grade Glioma (pHGG). 
2. Create a database of small molecules that may cross blood-brain barrier and help target certain cytokine-receptor pairs.
3. Apply directed evolution of antibodies to effectively target cytokines of interest and minimize off-target hits. 

Given the high mortality rate of pHGG, multiple labs around the world are making their transcriptomic data available for everyone to use and explore. Us bioinformatician's can organize and annotate these datasets in a way that can help wet lab and clinical researchers to easily explore their hypothesis without worrying too much about programming or tool usage limitations. Feel free to reach out to us is you find this project interesting, and if you would want to contribute. This is a complete open source project with all code, results, and datasets freely available to public.   

## Project 1: Targeting the SPP1-CD44 Axis in Pediatric High-Grade Glioma through Integrated Single-Cell and Structural Bioinformatics 

### <b>Manuscript under review</b> 
- [Manuscript draft link: For Nationwide Children's Hospital employees](https://nationwidechildrens-my.sharepoint.com/:w:/r/personal/ambuj_kumar_nationwidechildrens_org/Documents/SPP1_Project/Manuscript.docx?d=w6709113cea954b4794edb8017e88f115&csf=1&web=1&e=Fe9XyH) 

### Authors
- Shiwani Limbu, University of California Merced 
- Ambuj Kumar, Nationwide Children's Hospital


## Scripts and workflows
1. <b>Collect and clean pHGG single cell data</b> - Here we collect single cell data from public database and preprocess it for our specific study.   
   - 1.1. Single cell data collection - [Notebook](/Notebook/Step1_Dataset/pull_dataset.ipynb)
   - 1.2. Create seurat object - [Notebook](/Notebook/Step2_Generate_Seurat_Object/create_seurat.ipynb)
   - 1.3. Doublet finder - [Notebook](/Notebook/Step3_Preprocessing/doublet_finder.ipynb) 
   - 1.4. Preprocessing - [Notebook](/Notebook/Step3_Preprocessing/filter_cells.ipynb)

2. <b>Identify cell types present in pHGG single cell data</b> - Here, multiple clean preprocessed sample data is then integrated, clustered, annotated to identify various tumor and immune cell types present within it. 
   - 2.1. Harmony Integration - [Notebook](/Notebook/Step4_Integration/harmony_integration.ipynb)
   - 2.2. Clustering and cell annotation - [Notebook](/Notebook/Step5_Clustering/clustering.ipynb)
<img src="Notebook/Step5_Clustering/out/umap_base.png" alt="UMAP plot" width="500" height="400">

   - 2.3. InferCNV - [Notebook](/Notebook/Step6_Infercnv/infercnv.ipynb) 
<img src="Notebook/Step6_Infercnv/out/infercnv.png" alt="UMAP plot" width="600" height="400">


3. <b>Identify SPP1 cytokine activity and its regulatory features in pHGG.</b> This step enables us to identify active immunosupressive cytokine-receptor signals in pHGG, and transcription factors that may likely regulate it. Identifying these TFs can help us design specific antibody or small molecules to inhibit it.
   - 3.1. Cellchat cell-cell communication - [Notebook](Notebook/Step7_Cellchat/run_cellchat.ipynb)
   - 3.2. Pyscenic transcription factor prediction - [Notebook](Notebook/Step8_Pyscenic/pyscenic_plots.ipynb)

4. <b>Build antibody to target SPP1</b> - Here we aim to build an antibody that can specifically target CD44 binding motif of SPP1 cytokine. 

   - 4.1. <i><b>Molecular dynamics simulation</b></i> - Build and optimize SPP1 protein structure.
      - [Notebook](Notebook/Step9_SPP1_Molecular_Dynamics/simulate.ipynb)
      - [Slurm automation script](Notebook/Step9_SPP1_Molecular_Dynamics/src/simulation_sbatch_script.sh) 

   - 4.2. <i><b>Rosetta antibody design (build 23C3-v1)</b></i> - Build an antibody to target CD44 binding region of SPP1 protein.
      - [Rosetta script](Notebook/Step11_Antibody_Design/src/antibody_design.sh)
      - [esm2 embedding based AB variant selection](Notebook/Step11_Antibody_Design/esm_embedding_workflow.ipynb)

   - 4.3. <i><b>Humanization of 23C3-v1 (Hu23C3-v1)</b></i> - [Notebook](Notebook/Step11_Antibody_Design/humanize_epitope.ipynb) - Modify 23C3-v1 class I and class II epitope region amino acid sequences in 23C3-v1 to minimize immune response against it. This step requires careful selection of point mutations to eliminate epitope hits, specifically on the surface exposed regions on the protein, without loosing its SPP1 binding affinity. Below are the steps outline implemented in this study.  
      - 4.3.1. Pull class I and class II epitopes in 23C3-v1
      - 4.3.2. Select epitope regions that are exposed on antibody surface with <b>SASA >= 50%</b>
      - 4.3.3. Compare epitope regions of 23C3-v1 with other closely related humanized antibodies
      - 4.3.4. Incorporate relevant point mutations in epitopes found in CDR regions (if any) to neutralize it without loosing SPP1 binding.
      - 4.3.5. Incorporate replace 23C3-v1 residues with corresponding alignment position residue of humanized AB in non-CDR epitope regions.
      - 4.3.6. Test the modified 23C3-v1 binding affinity changes with SPP1 using molecular dynamics. 

   - 4.4. <i><b>Hu23C3-v1-SPP1 binding affinity test</b></i> - Check if we significantly loose SPP1 binding affinity upon humanization of 23C3-v1 
      - [MD simulation of Hu23C3-v1 docked with SPP1 and 23C3-v1 docked with SPP1](Notebook/Step12_Antibody_Molecular_Dynamics/simulate.ipynb)
      - [Binding affinity test](Notebook/Step13_Binding_Affinity_MMPBSA/mmpbsa.sh)

5. <b>Utility scripts<b> - [dir](Notebook/Utility)

<hr>

## Results


<b>Insights from GSEA results</b>

- d-1) Immune cells - 

   - T cells show strong signatures for T cell receptor complex, MHC class II protein complex, and immunological synapse. NK cells are uniquely enriched in cytotoxic granule lumen, perforin complex, and natural killer cell–mediated cytotoxicity pathways, consistent with their tumor‐killing potential.  

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


<hr>

# Other resources
## 1. Molecules capable of passing through blood brain barrier
- Manuscript - [https://www.nature.com/articles/s41597-021-01069-5](https://www.nature.com/articles/s41597-021-01069-5)
- Github - [https://github.com/theochem/B3DB/tree/main](https://github.com/theochem/B3DB/tree/main)