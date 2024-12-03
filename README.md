# Analysis Code for "Immune Checkpoint Blockade Induces Distinct Alterations in the Microenvironment of Primary and Metastatic Brain Tumors" (Sun et al., Journal of Clinical Investigation, 2023)

# scRNAseq Analysis Steps and Codes
1.	scRNAseq.Seurat.CombinedAnalysis.rmd: Integration, normalization, and cell type annotation of scRNAseq data using both our own and public datasets.
2.	scRNAseq.Seurat.MyeloidAnalysis.rmd: Sub-clustering and subtype annotation of myeloid cells.
3.	scRNAseq.Seurat.LymphoidAnalysis.rmd: Sub-clustering and subtype annotation of lymphoid cells.
4.	scRNAseq.DiffusionMap.LymphoidAnalysis.rmd: DiffusionMap analysis of lymphoid cell subsets.
5.	PlotDiffusionMap.R: Plotting 2D and 3D diffusion maps.
6.	scVelo.LymphoidAnalysis.py: RNA velocity analysis of lymphoid cell subsets.
7.	scRNAseq.Cellchat.Step1.CreateMetaData.rmd: Creating metadata for cell-cell interaction analysis using the Cellchat package.
8.	scRNAseq.Cellchat.Step2.CreatingSingleCellChatObject.rmd: Creating Cellchat objects for each treatment condition.
9.	scRNAseq.Cellchat.Step3.ComparativeAnalysis.rmd: Comparing cell-cell interactions across different treatment groups.
    
# Visium Spatial Transcriptomics (ST) Analysis Steps and Codes
1.	Visium.Seurat.CombinedAnalysis.rmd: Integration and normalization of ST data.
2.	Visium.RCTD.CombinedAnalysis.rmd: Cell type decomposition of the ST data.
3.	GenerateGMTFromScRNAseq.rmd: Generating gene signature GMT files for each cell type based on scRNAseq analysis. These GMT files are used in the subsequent spot-based gene set enrichment analysis of the Visium ST data.
4.	Visium.AUCell.NeighborhoodAnalysis.rmd: Performing gene set enrichment analysis to map each cell subtype on Visium ST, followed by cell neighborhood analysis.

# Notes
1. Obtain and unzip the [BM_Input.zip](https://drive.google.com/file/d/17vyuxs7EcM6WQYvOo93aTxmh2neRiLaW/view?usp=drive_link) file to access the input files.
2. The code available here represents the functional code utilized throughout the project. We are in the process of enhancing the ST analysis code and transforming it into a more user-friendly and easily executable version.
3. Requirements: R (tested with version 4.2.0).
