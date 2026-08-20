# De novo mutational signatures in lung adenocarcinoma (TCGA-LUAD)
**About:** This project extracts mutational signatures from lung cancer somatic mutations by NMF and matches them to known COSMIC signatures.

### Question
Which mutation processes have operated in LUAD, and do de novo signatures extracted by NMF match known COSMIC signatures?

### Data
| Input                | Source                                      | Build         | Notes                        | Reference                                                                                                                                                                 |
|----------------------|---------------------------------------------|---------------|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Somatic mutations    | TCGA-LUAD Masked Somatic Mutation MAF (GDC) | GRCh38 / hg38 | one MAF per case, aggregated | GDC: [Heath et al. (2021)](https://doi.org/10.1038/s41588-021-00791-5), TCGA-LUAD: [The Cancer Genome Atlas Research Network (2014)](https://doi.org/10.1038/nature13385) |
| Reference genome     | UCSC hg38 (`hg38.fa`)                       | hg38          | for trinucleotide context    | [Casper et al. (2025)](https://doi.org/10.1093/nar/gkaf1250)                                                                                                              |
| Reference signatures | COSMIC v3.4 SBS (GRCh38)                    | —             | for cosine matching          | [Tate et al. (2018)](https://doi.org/10.1093/nar/gky1015), Website: https://cancer.sanger.ac.uk/cosmic/                                                                   |

**Note:** Raw data is not committed (see .gitignore); the table above is sufficient to obtain every input. The processed SBS-96 matrix is small and is included.

### Software
| Tool                                  | Version | Purpose                         | Reference                                                             |
|---------------------------------------|---------|---------------------------------|-----------------------------------------------------------------------|
| scikit-learn                          | 1.6.1   | NMF signature extraction        | [Pedregosa et al. (2011)](https://scikit-learn.org/stable/about.html) |
| SigProfilerExtractor (Alexandrov lab) | 1.4.1   | Validation (de novo extraction) | [Islam et al. (2022)](https://doi.org/10.1016/j.xgen.2022.100179)     |



### Method
**Step 1 - Building the SBS-96 Catalogue:** TCGA-LUAD somatic-mutation MAFs were retrieved from the GDC API and filtered to single-base substitutions on the standard chromosomes. Each substitution was classified by trinucleotide context (e.g. A[C>A]G) under the pyrimidine convention, producing a 96-channel profile per tumor. Contexts were read from the hg38 genome reference, with the middle base validated against the recorded reference allele. The aggregate profile showed a strong C>A component, the expected the tobacco-smoking signature.

**Step 2 - NMF:** Non-negative matrix factorization from scikit-learn decomposes the 96 x tumor catalogue (V) into signatures (H) and per-tumor exposures (W), such that V ≈ W x H. Non-negativity is required because mutation counts cannot be negative. For each candidate rank k, we analyze reconstruction error, which measures fit, and stability, which measures reproducibility. The NMF function uses Kullback-Leibler divergence.

**Step 3 - Rank Selection:** Across k = 2-20 (25 restarts each), reconstruction error and restart stability indicated the range of best k-values to be between 2 and 6. Error's steep decline flattened in that range and stability remained above 0.90. Upon plotting the extracted signatures per k, k = 2 and k = 3 show distinct, interpretable signatures; beyond k = 3 the additional signatures increasingly resemble split or noisy versions of the same processes rather than new biology. Thus, k = 3 is chosen as the appropriate rank.

**Step 4 - COSMIC Matching:** Cosine similarity between the extracted signatures and COSMIC v3.4 associated Signature 3 with SBS4 (tobacco smoking), Signature 2 with SBS2/13 (APOBEC), and Signature 1 with SBS5 ('clock').

**Step 5 - Exposures:** Per-tumor signature fractions (W) showed smoking and SBS5/SBS6 broadly distributed, while APOBEC activity was concentrated in a subset of tumors.

**Step 6 - Validation:** Signatures were cross-checked against SigProfilerExtractor, which independently recovered three signatures. Signature 1 matched to SBS5, Signature 2 matched to SBS2 followed by SBS13, and Signature 3 to SBS4. The fact that an established, more heavily optimized tool converges on the same three processes reinforces the validity of the NMF pipeline.  

### Findings
The NMF model extracted three signatures. We selected the rank k=3 from two diagnostics, reconstruction error and stability, computed for k values from 2 to 20. Error declined sharply and slowed near k=3, while stability remained high (above 0.90) through this value. Inspecting the signature profile plots at nearby ranks, k=3 yielded distinct, interpretable singatures, while higher ranks began splitting single processes into near-duplicate signatures. Signature 1 most closely matched SBS5 (with a cosine similarity of 0.768), followed by SBS6 (0.763). SBS5 is a clock-like mutational signature whose accumulation correlates with age; SBS6 is linked to defective DNA mismatch repair. Signature 2 mapped to SBS2 (with a cosine similarity of 0.762), with SBS13 the next-beset match (0.635); both are linked to APOBEC, a family of proteins known to cause hypermutations in DNA. Finally, Signature 3 matched SBS4 (with a cosine similarity of 0.979), which is strongly associated with tobacco smoking. This is a confident assignment, considering both the cosine value and the characteristic C>A profile. Per-tumor exposures showed smoking and background signatures broadly distributed, while APOBEC activity was concentrated in a subset of tumors. To validate these signatures, we ran SigProfilerExtractor — an established, heavily-optimized tool for de novo signature extraction — on the same catalogue. It retrieved three signatures, matching our rank, whose COSMIC identities agreed with ours: SBS5 (0.767, with SBS6 next at 0.618) for Signature 1, SBS2 (0.762) then SBS13 (0.635) for Signature 2, and SBS4 (0.979) for Signature 3. The cross-tool agreement supports the validity of our pipeline. We conclude that three processes underlie the TCGA-LUAD mutational landscape: tobacco smoking (SBS4), APOBEC activity (SBS2/13), and an SBS5/6-like background. 

### Limitations
- Signatures were derived from one cancer cohort (TCGA-LUAD) of bulk tumor samples, so they may not generalize to other cohorts, and rare processes might be missed.
- scikit-learn's general-purpose NMF was used rather than a tool built specifically for signature extraction. However, results were cross-validated against SigProfilerExtractor; bootstrap resampling would provide further validation. 
- Only dominant processes were recovered; finer sub-signatures were not differentiated.
- The rank of k=3 was chosen from reconstruction error, stability, and analysis of signature profile graphs. It remains a bit of a subjective choice; a different criterion could support another k value. 

### Repository Structure
    .
    ├── README.md
    ├── requirements.txt
    ├── luad_signatures.ipynb  # the full analysis, narrated
    ├── results/
    │   ├── figures/                    # profiles, rank selection, exposures
    │   └── sigprofiler/                # SigProfiler output
    └── data/
        ├── luad_mafs/                  
        ├── raw/                        
        └── processed/
            ├── sbs96_matrix.csv
            └── sbs96_matrix.txt        

### Reproduce
    # 1. Environment
    python -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt

    # 2. Obtain data (see Data section) into data/raw/:
    #    - TCGA-LUAD MAF from the GDC portal
    #    - hg38.fa from UCSC
    #    - COSMIC v3.4 SBS reference

    # 3. Run the notebook top-to-bottom
    jupyter lab notebooks/luad_signatures.ipynb

### Citations + Attributions
This project uses publicly available datasets and open-source software. Please refer to the linked publications and resources for the appropriate citations and attribution.