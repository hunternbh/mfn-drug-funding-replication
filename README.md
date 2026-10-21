# U.S. most-favored-nation pricing and foreign drug funding

**Replication package for the paper on SSRN.**  
Hunter Ng · Shou-Ming Chang · Yung-Cheng Lin · October 21, 2026

[Read the paper on SSRN](https://ssrn.com/abstract=7652061) · [Paper DOI](https://doi.org/10.2139/ssrn.7652061) · [Reproduction instructions](documentation/REPRODUCE.md) · [Data dictionary](documentation/DATA_DICTIONARY.md)

[![Replication checks](https://github.com/hunternbh/mfn-drug-funding-replication/actions/workflows/reproduce.yml/badge.svg)](https://github.com/hunternbh/mfn-drug-funding-replication/actions/workflows/reproduce.yml)

## What this repository provides

This package reproduces the paper's numerical calculations and figures **from final aggregate analysis tables**. The study compares documented public-funding entry in Australia, Switzerland and Britain with Brazil, Hong Kong and Singapore around the May 2025 U.S. MFN announcement, following the same 361 medicines.

The public files contain country-period totals, country-month figure values, aggregate regression cross-products and reference results. **Raw scraped archives, document text, medicine-level records, linkage masters and individual medicine weights are not distributed.** Collection, medicine-level matching and weight construction precede the starting point of this public package. See [Data access and replication scope](documentation/DATA_ACCESS.md).

## Reproduce the results

Requires Python 3.12. From this folder:

```bash
python -m pip install -r requirements.txt
python reproduce.py
```

The analysis runs offline: no accounts, API keys or source downloads are needed. Outputs appear in `build/tables/`, `build/figures/` and `build/replication_report.json`. Run the same commands on Windows, macOS or Linux.

The checks recompute:

- All 2,688 alternative specifications and their 2,100 distinct designs, retaining every result.
- 107 main and historical regressions, including country omissions, alternative weights and historical comparisons, from country-level regression cross-products.
- 20 calendar/date checks, including both admissible British date conventions.
- The five-year, three-country DiD and all four controlled projection specifications, earlier prediction checks and annual sensitivities.
- The five manuscript figures from their final plotted tables.

The [results map](documentation/RESULTS_MAP.md) connects files to manuscript exhibits. The [methods note](documentation/METHODS.md) explains the estimators and the distinction between observed data and prediction draws.

## Folder guide

| Folder | Contents |
|---|---|
| `data/` | Final aggregate tables needed to rerun calculations and draw figures |
| `code/` and `reproduce.py` | Offline numerical and plotting code |
| `expected/` | Reference numerical results used for verification |
| `figures/` | Ready-to-view manuscript figures |
| `documentation/` | Instructions, dictionary, methods and data-access information |

## Cite the paper

Ng, Hunter, Shou-Ming Chang, and Yung-Cheng Lin. 2026. “U.S. most-favored-nation pricing and foreign drug funding.” SSRN, October 21. https://doi.org/10.2139/ssrn.7652061.

```bibtex
@techreport{ng2026mfn,
  author = {Ng, Hunter and Chang, Shou-Ming and Lin, Yung-Cheng},
  title = {U.S. most-favored-nation pricing and foreign drug funding},
  year = {2026},
  institution = {SSRN},
  number = {7652061},
  doi = {10.2139/ssrn.7652061},
  url = {https://ssrn.com/abstract=7652061}
}
```

GitHub's **Cite this repository** button also supplies the paper citation through `CITATION.cff`. This repository is the companion replication package, not a separate paper or journal publication.
