# Exploring role of SPP1 in regulating T cell exhaustion in pediatric high grade glioma (pHGG)
We had a recent dealth of our family member due to pediatric high grade glioma, and we aim to slowly clean and analyze all available pediatric high grade glioma transcriptomic, structural bioinformatics, and proteomic data, and release it to the community, to contribute to pHGG research as much as we can from our end.  

## Intro
T cell exhaustion is a leading cause contributing to failure of immunotherapy in pHGG. SPP1 released from Tumor associated macrophage (TAM) inducing T cell exhaustion. Here we first aim to,
- 1. Validate SPP1 mediated communication between Tumor associated macrophage (TAM) and T cells, leading to T cell exhaustion in pHGG, using pHGG single cell data.
- 2. Explore methods to block SPP1-CD44 interation
   - 2.1. via inhibiting TFs regulating SPP1 expression in TAM. This can be obtained by 
      - 2.1.1. Gene expression based prediction of likely SPP1 TFs in TAM using pyscenic.
      - 2.1.2. Gene sequence based prediction of likely SPP1 TFs in TAM using MEME+Tomtom workflow.
      - 2.1.3. Validate predicted TFs using gene expression spearman correlation score.
      - 2.1.4. Identify gene therapy vectors efficient in crossing blood brain barrier and targetting SPP1 tracriptional mechanism in TAMs in pHGG.
   - 2.2. via evolving existing SPP1 antibody molecule to specifically target SPP1 protein motif associated with SPP1-CD44 interation.
      - 2.2.1. We have antibodies targetting residues 43-48, as well as residues 162-168. Light and Heavy chains of these antibodies can be evolved further to direct its binding towards its CD44 bidning motif using Rosetta antibody designer and directed evolution methods.  
      - 2.2.2. Designing bivalent peptides to target SPP1 and its career protein.
      - 2.2.3. Target SPP1 CD44 binding motif using small molecules that can cross blood brain barrier. 
         - 2.2.3.1. List of these molecules can be obtained from [https://www.nature.com/articles/s41597-021-01069-5](https://www.nature.com/articles/s41597-021-01069-5), [https://github.com/theochem/B3DB/tree/main](https://github.com/theochem/B3DB/tree/main). 
         - 2.2.3.2. Subset of these molecules can be selected using molecular weight, hydrogen donor and acceptor counts, and number of rotation bonds filter.
         - 2.2.3.3. Virtual screening of the compunds can help identify high rank compounds that may have higher binding affinity to SPP1 CD44 binding motif.
         - 2.2.3.4. Functional groups of high rank compounds can be further modulated to reduce off-targets.   

## SPP1 workflow 
- 1. Data download - [Notebook](Notebook/Step1_dataset/pull_dataset.ipynb)
- 2. Generate_Seurat_Object - [Notebook](Notebook/Step2_Seurat_Object/create_seurat.ipynb)
- 3. Data preprocessing - 
   - 3.1. Doubletfinder - [Notebook](Notebook/Step3_Preprocessing/doublet_finder.ipynb)
   - 3.2. Base filters applied using cell and gene counts - [Notebook](Notebook/Step3_Processing/filter_cells.ipynb)
- 4. Sample integration - <b>Harmony</b> - [Notebool](Notebook/Step4_Integration/harmony_integration.ipynb)
- 5. Cell clustering and annotation - [Notebook](Notebook/Step5_Clustering/clustering.ipynb)
- 6. Copy number variant prediction - <b>Infercnv</b> - [Notebook](Notebook/Step6_Infercnv/infercnv.ipynb)
- 7. Cell-cell communication - <b>Cellchat</b> - [Notebook](Notebook/Step7_Cellchat/run_cellchat.ipynb)
- 8. Transcription factor prediction - <b>Pyscenic</b> 
   - 8.1. [Pyscenic TF prediction scripts](Notebook/Step8_Pyscenic_plots/src) 
   - 8.2. Plots - [Notebook](Notebook/Step8_Pyscenic/pyscenic_plots)
- 9. Model SPP1 protein structure using Alphafold, as well as ab-initio modelling.
- 10. Improve structure model using
   - 10.1. Replica exhange to identify regions sensitive to subtle change in tumor microenvironment.
   - 10.2. Long term production MD to stabilize highly dissordered regions before and after phosphorylation.
   - 10.3. Pull stable conformations using free energy landscape.  
- 11. Antibody directed evolution 
   - 11.1. Dock known SPP1 antibody with its CD44 targetting motif (residue 121-140) with SPP1 antibody its light and heavy chains.
   - 11.2. Apply Rosetta antobody design workflow to model ~2500 AB variants.
   - 11.3. Calculate light chain and heavy chain diversity amoung 2500 AB variants to create AB clusters. 
   - 11.4. Apply Protein binding energy estimator PBEE to calculate binding affinity for all 2500 AB variants.
   - 11.5. Pick top high rank candidate using high Rosetta score, high binding affinity, and low sequence diversity.
   - 11.6. Apply Molecular mechanics-poisson-boltzman surface area based binding affinity calculation and comaprison between top AB candidates. 
   - 11.7. Scan AB surface residues for likely off-target binding within pHGG tumor microenvironment. 
   - 11.8. If potential off-targets predicted - 
      - 11.8.1. Pull known sequences of corresponding AB from Therian species, to be conservative when picking  mutations for evolving regions associated with off-target binding to create mutation library 1 (lib 1) for directed evolution. 
      - 11.8.2. Apply Rosetta mutation scan to identify AB off-targetting region surface conservative mutations to pick the ones that may not lead to obvious deleterious impact on the structure.
      - 11.8.3. Apply ESM inveserfold based mutation scan to identify conservative mutations on AB off-targetting region surface that may not lead to obvious deleterious impact on the structure.
      - 11.8.4. Apply elastic network model to identify surface residues that may induce allosteric deleterious response to heavy and light chain of AB upon mutations, leading to dead antobody variants in the sample. 
      - 11.8.5 Remove Cystine mutations to avoid oxidative damage.
      - 11.8.6. Pick remaining mutations to explore directed evolution of AB to reduce off-target binding.

## Results
- 1. Figure 1.

<img src="Manuscript/Fig1.png">

- 2. Figure 2. 

<img src="Manuscript/Fig2.png">

-3. Figure 3.

<img src="Manuscript/Fig3.png">
