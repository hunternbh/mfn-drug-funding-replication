# Methods and interpretation of the final tables

## Main comparison

The fixed cohort contains 361 medicines linked to 2023 US Medicare Part D spending. Three OECD countries are compared with three non-OECD markets. The primary before period is January 2023–April 2025; the after period is June 2025–May 2026. May 2025 is omitted. The outcome is documented entry in the captured funding programmes, not an unfunded risk-set hazard. The higher-minus-lower spending comparison uses the fixed FDA-age weights already incorporated into the final aggregates.

## Alternative specifications

The full list of 2,688 labelled specifications and all results are retained. They represent 2,100 distinct effective designs. The 96% negative-sign figure is 288 of 300 distinct designs in the original FDA-age family; it is not a confidence level or 300 independent replications. The complete set is reproduced from the country-specific pre/post table and the specification definitions. Country-level CR1 and CR2/Satterthwaite calculations are rerun.

## Regression cross-products

For each model and country, the package supplies X'X and X'y. The metadata provide sample size, design rank and the coefficient column of interest. Summing the country totals identifies the OLS/WLS coefficient vector. The country score is X'y minus X'X times the estimated coefficient. Its outer products reproduce the cluster covariance calculation, including the retained finite-sample correction. No medicine identities or source passages enter this step. Some original algebraic checks use phase means; the cross-products are computational sufficient statistics, not newly observed monthly histories.

## Five-year comparison and projections

Australia and Switzerland are compared with Brazil using May 2020–April 2025 and June 2025–May 2026. Brazil is one comparison country; single-donor synthetic DiD reduces to ordinary DiD here. Four controlled gap projections are retained. Their selection uses earlier prediction errors, not the post-announcement sign. Predictive uncertainty uses circular 12-month residual blocks with model refitting. All outcomes remain conditional on the retrospective fixed medicine cohort. Five pre-period years do not extend the observed post-period beyond one year.

## Figures and dates

Main figures use country-period means and approval-cohort aggregate results. Historical smoothed series restart after the omitted announcement month. Date conventions preserve the range of admissible British records. Unmatched medicines are not automatically classified as unfunded; the FDA-age documentation graph gives documented positive shares, not a complete national waiting-time distribution.
