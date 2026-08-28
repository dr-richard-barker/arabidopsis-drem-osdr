# Licensing

**Code** (`scripts/`, build files): MIT, see [LICENSE](LICENSE).

**Data, figures and prose** (`results/`, `figures/`, `manuscript/`, `docs/`): CC-BY-4.0.

## Third-party components

| Component | Licence | How it is used |
|---|---|---|
| DREM 2.0.7 (`drem.jar`) | GPL-3.0, Ernst laboratory | **Downloaded at run time by `scripts/08_run_drem.py`, not redistributed here.** The pipeline invokes it as an external program; it is not linked into or derived from this code. |
| NASA OSDR study data | NASA Open Data | Re-fetched by `scripts/01_fetch_osdr.py`, not redistributed. Cite the original investigators listed in `MANIFEST.tsv`. |
| DAP-seq peaks (GEO GSE60143) | NCBI GEO terms | Downloaded, not redistributed. Cite O'Malley et al. 2016. |
| SOG1 ChIP-seq (GEO GSE112529) | NCBI GEO terms | Downloaded, not redistributed. Cite Bourbousse et al. 2018. |
| Single-nucleus atlas (GEO GSE226097) | NCBI GEO terms | Used via derived signature matrices from the sibling Tropism project. |
| Ensembl Plants release-54 annotation | EMBL-EBI, no restriction | Downloaded, not redistributed. |
| PlantTFDB TF list | PlantTFDB terms | Downloaded, not redistributed. |
| TAIR gene aliases | TAIR terms | Downloaded, not redistributed. |

Every file the repository *does* ship is listed in `MANIFEST.tsv` with
`redistributed_here = yes`, its source, its licence and its SHA-256.
