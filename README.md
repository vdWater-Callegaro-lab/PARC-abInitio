# PARC-abInitio

Shared repository for the **PARC ab initio case study**, containing scripts and analyses contributed by consortium partners.

This repository supports the generation, quality control, reporting, and manuscript-level visualization of Benchmark Dose (BMD) analyses performed as part of the PARC ab initio case study.

## Repository structure

```text
PARC-abInitio/
├── 1_Analyses/
├── 2_BMD_Reporting_Template/
└── 3_Manuscript_Figures/
```

## Folder contents

### `1_Analyses/`

This folder contains scripts provided by consortium partners to generate input files for the Benchmark Dose (BMD) analysis and to support quality control.

- scripts used for differential expression analysis;
- partner-specific analysis workflows.

For analyses using **DRomics**, this folder may also include:

- scripts used to run the BMD analysis

### `2_BMD_Reporting_Template/`

This folder contains the template used to evaluate BMD analysis results.

Running the script generates an HTML report that can be used to visually inspect the data and assess the BMD results.

Details on the required input files are provided at the beginning of the script.

### `3_Manuscript_Figures/`

This folder contains scripts used to generate figures included in the manuscript.

## Intended use

This repository is intended to provide a transparent and reproducible record of the analysis workflows used in the PARC ab initio case study.

Users can use the repository to:

- review partner-specific analysis scripts;
- inspect the generation of BMD input files;
- generate HTML reports for BMD result evaluation;
- reproduce manuscript figures.


## Citation and archived version

The version of this repository associated with the manuscript submission is archived on Zenodo.


- GitHub release: `v1.0.0-submission`
- Zenodo version DOI: https://doi.org/10.5281/zenodo.20413136 


Please cite the Zenodo version DOI when referring to the repository version used for the submitted manuscript.

The GitHub repository may be updated over time. Archival versions associated with the manuscript are preserved on Zenodo. Each GitHub release archived on Zenodo receives a version-specific DOI.


## License

This repository is licensed under the MIT License. See the [`LICENSE`](LICENSE) file for details.
