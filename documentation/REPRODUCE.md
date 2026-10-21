# Reproduction instructions

1. Download this repository as a ZIP, or clone it.
2. Open a terminal in the extracted repository folder.
3. Create an optional virtual environment: `python -m venv .venv`. Activate it using the command appropriate to your operating system.
4. Install the pinned packages: `python -m pip install -r requirements.txt`.
5. Run `python reproduce.py`.

A successful run ends with `"status": "passed"`. It verifies every final-input checksum, recomputes coefficients and uncertainty, compares them with the numerical reference files, and regenerates figures. All generated files go into `build/`. The supplied final data tables are not modified. Small floating-point differences across numerical libraries are allowed by explicit tolerances.

## What is estimated again

The alternative-specification calculations use country-specific pre/post means. The main and older regression checks use aggregate cluster cross-products, which are sufficient for OLS/WLS coefficients and the retained clustered covariance formula. The timing checks use final country-month outcomes. The five-year projection models are fitted again using only their observed pre-period aggregate table. Bootstrap draws quantify uncertainty; they are not additional observed years.

## What is read as a final table

The graph data contain the final exposure-stratum rates, FDA-cohort rates, country means and historical comparison series. The documentation-share curve starts from final numerators and denominators. These inputs permit exact arithmetic and graphical replication without disclosing the underlying medicine-level matching archive. Recreating extraction, linkage and the individual medicine weights requires materials outside the public package.
