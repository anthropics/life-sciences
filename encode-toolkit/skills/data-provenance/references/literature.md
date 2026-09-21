# Data Provenance — Literature References

**Last updated:** 2026-03-07
**Purpose:** Reference catalog for the data-provenance skill — papers on data management principles, reproducibility standards, workflow engines, and containerization that underpin provenance tracking in genomics.

---

### Wilkinson et al. 2016 — FAIR Guiding Principles for scientific data

- **Citation:** Wilkinson MD, Dumontier M, Aalbersberg IJ, Appleton G, Axton M, Baak A, Blomberg N, Boiten JW, da Silva Santos LB, Bourne PE, Bouwman J, Brookes AJ, Clark T, Crosas M, Dillo I, Dumon O, Edmunds S, Evelo CT, Finkers R, Gonzalez-Beltran A, Gray AJG, Groth P, Goble C, Grethe JS, Heringa J, 't Hoen PAC, Hooft R, Kuhn T, Kok R, Kok J, Lusher SJ, Martone ME, Mons A, Packer AL, Persson B, Rocca-Serra P, Roos M, van Schaik R, Sansone SA, Schultes E, Sengstag T, Slater T, Strawn G, Swertz MA, Thompson M, van der Lei J, van Mulligen E, Velterop J, Waagmeester A, Wittenburg P, Wolstencroft K, Zhao J, Mons B. The FAIR Guiding Principles for scientific data management and stewardship. *Scientific Data*, 3, 160018, 2016.
- **DOI:** [10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18)
- **PMID:** 26978244 | **PMC:** PMC4792175
- **Citations:** ~8,000
- **Key findings:** Established the FAIR principles (Findable, Accessible, Interoperable, Reusable) as the standard for scientific data management. The principles emphasize machine-actionability and persistent identifiers, providing the conceptual foundation for data provenance tracking. ENCODE data follows FAIR principles through DOI assignment, standardized metadata schemas, and open API access, making this paper the cornerstone reference for any provenance workflow.

---

### Baker 2016 — Reproducibility crisis in science

- **Citation:** Baker M. 1,500 scientists lift the lid on reproducibility. *Nature*, 533(7604), 452-454, 2016.
- **DOI:** [10.1038/533452a](https://doi.org/10.1038/533452a)
- **PMID:** 27225100
- **Citations:** ~3,200
- **Key findings:** Reported results from a Nature survey of 1,576 researchers, finding that over 70% had tried and failed to reproduce another scientist's experiments. The survey highlighted that selective reporting, low statistical power, and inadequate methodological documentation are major contributors. This paper provides the empirical motivation for rigorous provenance tracking: without detailed records of tools, versions, parameters, and inputs, computational results cannot be reliably reproduced.

---

### Stodden et al. 2016 — Empirical study of computational reproducibility

- **Citation:** Stodden V, McNutt M, Bailey DH, Deelman E, Gil Y, Hanson B, Heroux MA, Ioannidis JPA, Taufer M. Enhancing reproducibility for computational methods. *Science*, 354(6317), 1240-1241, 2016.
- **DOI:** [10.1126/science.aah6168](https://doi.org/10.1126/science.aah6168)
- **PMID:** 27940837
- **Citations:** ~1,800
- **Key findings:** Outlined concrete recommendations for improving computational reproducibility, including sharing code and data, using version control, recording computational environments, and adopting workflow systems. The authors argued that journals, funders, and institutions must all contribute to a culture of reproducibility. These recommendations directly inform the provenance standard used in the ENCODE toolkit, where every operation logs tool versions, commands, input accessions, and output descriptions.

---

### Grüning et al. 2018 — Bioconda for reproducible bioinformatics

- **Citation:** Grüning B, Dale R, Sjodin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Koester J, Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. *Nature Methods*, 15(7), 475-476, 2018.
- **DOI:** [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7)
- **PMID:** 29967506
- **Citations:** ~1,200
- **Key findings:** Introduced Bioconda, a distribution of bioinformatics software through the Conda package manager, with over 3,000 software packages and automatic container generation. Bioconda ensures that specific tool versions can be reliably installed across different computing environments. For provenance tracking, Bioconda-managed environments provide exact version pinning, enabling bit-for-bit reproducibility of ENCODE analysis pipelines.

---

### Di Tommaso et al. 2017 — Nextflow workflow engine

- **Citation:** Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35(4), 316-319, 2017.
- **DOI:** [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
- **PMID:** 28398311
- **Citations:** ~2,000
- **Key findings:** Presented Nextflow, a reactive workflow framework that enables scalable and reproducible scientific workflows using software containers (Docker, Singularity). Nextflow decouples the pipeline logic from execution environment, supporting local, HPC, and cloud execution. The ENCODE toolkit pipeline skills use Nextflow as the standard execution engine, and its built-in provenance features (trace files, execution reports, DAG visualization) are integral to the data-provenance workflow.

---

### Amstutz et al. 2016 — Common Workflow Language specification

- **Citation:** Amstutz P, Crusoe MR, Tijanic N, Chapman B, Chilton J, Heuer M, Kartashov A, Leehr D, Menager H, Nedeljkovich M, Scales M, Soiland-Reyes S, Stojanovic L. Common Workflow Language, v1.0. *Specification*, Common Workflow Language working group, 2016.
- **DOI:** [10.6084/m9.figshare.3115156.v2](https://doi.org/10.6084/m9.figshare.3115156.v2)
- **Citations:** ~600
- **Key findings:** Defined the Common Workflow Language (CWL), an open standard for describing analysis workflows in a portable, interoperable manner. CWL specifications are vendor-neutral and can be executed across different workflow engines. ENCODE pipelines were originally specified in CWL/WDL formats, and the standard provides a formal framework for capturing the complete computational provenance of a workflow, from input data to final outputs.

---

### Goecks et al. 2010 — Galaxy for accessible, reproducible genomics

- **Citation:** Goecks J, Nekrutenko A, Taylor J, Galaxy Team. Galaxy: a comprehensive approach for supporting accessible, reproducible, and transparent computational research in the life sciences. *Genome Biology*, 11(8), R86, 2010.
- **DOI:** [10.1186/gb-2010-11-8-r86](https://doi.org/10.1186/gb-2010-11-8-r86)
- **PMID:** 20738864 | **PMC:** PMC2945788
- **Citations:** ~4,000
- **Key findings:** Described Galaxy as a web-based platform that automatically tracks provenance by recording every computational step, parameter, and input/output relationship. Galaxy's history system provides a complete audit trail for all analyses, and workflows can be shared and re-executed. Galaxy's approach to built-in provenance tracking serves as a model for the ENCODE toolkit's own provenance standard, where every operation is logged with sufficient detail to regenerate results.

---
