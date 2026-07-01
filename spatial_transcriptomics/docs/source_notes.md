# Source Notes

This scaffold was designed from three sources of evidence.

## HattieChungLab/ovarian_aging

Repository inspected on 2026-06-29 from `https://github.com/HattieChungLab/ovarian_aging`, main branch `82406337d5c3fd8785372e13b1d102db9e13649b`.

Reusable design points absorbed into this pipeline:

- Keep raw spot-level AnnData objects and write processed `.h5ad` checkpoints.
- Add coordinates from platform-specific barcode/coordinate files into `obsm["spatial"]`.
- Run sample-aware clustering after normalization, HVG selection, PCA and optional Harmony integration.
- Use spatial plotting early to catch smear/registration/cropping problems.
- Build object/segment-level summaries from spot-level labels, then assign object labels from cell-type composition.
- Preserve both spot-level and segment-level outputs.
- Compute object composition, object progression/pseudotime, and scaled radial distance within tissue objects.
- Keep cell-cell communication as a downstream interpretation layer, not as a prerequisite for core QC.
- Keep cNMF/NMF-style gene programs as a reusable downstream module.

## AnySearch results

Searches were run with the AnySearch skill on 2026-06-29.

Key sources returned:

- Orchestrating Spatial Transcriptomics Analysis with Bioconductor, intermediate processing chapter: normalization, HVG/SVG, PCA/UMAP, BayesSpace and marker workflows.
- Single-cell best practices, spatial deconvolution chapter: Cell2location, RCTD, SPOTlight, Stereoscope and deconvolution practical guidance.
- Squidpy Visium tutorial: image features, spatial neighbors, neighborhood enrichment, co-occurrence, ligand-receptor testing and Moran's I.
- Scanpy/Seurat spatial tutorials: Visium object loading, plotting and spatial workflow conventions.
- Chung lab publications and Nature Aging entry for the ovarian aging spatial study.

Additional supplementation was added on 2026-06-29 after the user requested a maximal method shelf.

Method and workflow categories recorded:

- Containers and workflow ecosystems: SpatialData paper and user guide (`https://doi.org/10.1038/s41592-024-02212-x`, `https://spatialdata.scverse.org/en/latest/user_guide.html`), Sopa documentation and paper (`https://prism-oncology.github.io/sopa/`, `https://doi.org/10.1038/s41467-024-48981-z`), Giotto Suite (`https://doi.org/10.1038/s41592-025-02817-w`), OSTA/Bioconductor infrastructure, Voyager, Squidpy and Seurat spatial v5/Visium HD vignettes.
- Spatial QC and normalization: SpotSweeper (`https://www.nature.com/articles/s41592-025-02713-3`), SpaNorm (`https://doi.org/10.1186/s13059-025-03565-y`), OSTA normalization and QC recommendations.
- Spatially variable genes: nnSVG (`https://doi.org/10.1038/s41467-023-39748-z`), SVG method categorization (`https://doi.org/10.1038/s41467-025-56080-w`), SVG benchmarking (`https://doi.org/10.1093/bioinformatics/btaf131`), SpatialDE/SPARK-X/trendsceek/SOMDE review material.
- Spatial domain and niche methods: BANKSY (`https://doi.org/10.1038/s41588-024-01664-3`), GraphST (`https://doi.org/10.1038/s41467-023-36796-3`), BayesSpace, SpaGCN, STAGATE, PRECAST, SpaceFlow, domain-detection benchmarks (`https://doi.org/10.1002/imt2.70084`, `https://doi.org/10.3390/cells14141060`), MISTy (`https://doi.org/10.1186/s13059-022-02663-5`), NicheCompass (`https://doi.org/10.1038/s41588-025-02120-6`) and scNiche.
- Deconvolution and reference mapping: Cell2location, RCTD/spacexr, SPOTlight, Stereoscope, Tangram, SpatialScope (`https://doi.org/10.1038/s41467-023-43629-w`), SpatialcoGCN (`https://doi.org/10.1093/bib/bbae130`), STdGCN (`https://doi.org/10.1186/s13059-024-03353-0`), SpaJoint (`https://doi.org/10.1093/bib/bbag158`) and practical benchmark (`https://doi.org/10.1038/s41467-023-37168-7`).
- High-resolution and image-aware analysis: Seurat Visium HD vignette, Bin2cell (`https://doi.org/10.1093/bioinformatics/btae546`), Sopa, SPArrOW (`https://doi.org/10.1101/2024.07.04.601829`), SpatialData, STPath (`https://doi.org/10.1038/s41746-025-02020-3`), ST-Align (`http://arxiv.org/abs/2411.16793`) and H&E/spatial transcriptomics survey material.
- Multi-sample alignment and spatial multiomics: STAligner (`https://www.nature.com/articles/s43588-023-00543-x`, `https://github.com/zhoux85/STAligner`), PASTE/PASTE2, moscot (`https://doi.org/10.1101/2023.05.11.540374`), MEFISTO/MOFA2 (`https://biofam.github.io/MOFA2/MEFISTO.html`) and Panpipes (`https://doi.org/10.1186/s13059-024-03322-7`).
- Cell-cell communication: COMMOT (`https://doi.org/10.1038/s41592-022-01728-4`), CellChat (`https://doi.org/10.1038/s41467-021-21246-9`), spatial CellChat preprint, LIANA/MISTy tutorials, MultiNicheNet and reviews/comparisons of ligand-receptor methods (`https://doi.org/10.1038/s41467-022-30755-0`, `https://doi.org/10.1038/s41576-023-00685-8`).
- Trajectory, object and geometry analysis: Spateo (`https://doi.org/10.1101/2022.12.07.519417`), StPedf (`https://doi.org/10.1371/journal.pcbi.1014346`), moscot, radial/object summaries inspired by `ovarian_aging`, topology and spatial geometry review material.

## Pipeline choices

- Default execution backend is Python/scverse because it handles `h5ad`, Visium, Slide-seq-like coordinate data and image-aware Squidpy analysis.
- R/Rmd is used for review plots and reports so outputs follow the local `wgbs_sperm` result-reading style.
- Heavy or context-specific tools are optional: Cell2location, RCTD, SPOTlight, BayesSpace, BANKSY, SpaGCN, MultiNicheNet and COMMOT can be enabled once their dependencies and references are available.
- The default path remains tissue-agnostic: input validation, object construction, QC, normalization, integration, spatial graph/statistics, spatial domains, marker annotation, marker-score deconvolution fallback, generic segmentation, gene programs, communication and summary.
- The expanded method shelf is intentionally over-complete. It is meant for selection and sensitivity analysis, not for enabling every branch in one run.
