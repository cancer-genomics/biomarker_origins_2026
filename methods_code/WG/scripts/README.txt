WGS-CRE
Jason Zavras
Cancer Genomics Lab

Overview: The function of this pipeline is to extract footprint features from low coverage cfDNA WGS data using regulatory element coordinates.
Key Data:
    - CRE (TFBS, cell type specific ATAC peaks) coordinates
    - cfDNA (bed, starch, rds) patient files

Key Methodology:
    - Apply strict QC of both CRE and cfDNA data
    - Apply GC correction (weighted fragments from Cristiano et al 2019) to generate accurate coverage profiles
    - Extract many different feature famiies
        - Rel_cov, coverage-based nucleosome depletion
        - FLD, fragment-length-distribution profiling at center positions compared to flanks
        - SWFSE, fragment length bin entropy (sliding window, 100-120bp, 120-140bp, etc.)