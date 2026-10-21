# Final-table data dictionary

All tables are aggregate or model-level. No row is an individual medicine record.

## `age_balanced_country_period_outcomes.csv`

Rows: 6. Columns: `country`, `oecd`, `pre`, `post`, `change`.

## `all_six_country_results.csv`

Rows: 20. Columns: `test`, `estimate`, `cluster_SE`, `lower95`, `upper95`, `country_clusters`, `approx_MDE80`, `units`.

## `alternative_country_period_outcomes.csv`

Rows: 17,920. Columns: `spec_id`, `country`, `oecd`, `pre`, `post`, `change`.

## `alternative_specification_definitions.csv`

Rows: 2,688. Columns: `spec_id`, `country_set`, `outcome`, `window`, `measurement_rule`, `UK_endpoint`, `pre_start`, `pre_end`, `post_start`, `post_end`, `pre_months`, `post_months`, `medicines`, `units`, `effective_design_id`.

## `Approval_cohort_country_outcomes.csv`

Rows: 18. Columns: `cohort`, `country`, `pre`, `post`, `change`, `cluster_influence`.

## `Approval_cohort_results.csv`

Rows: 3. Columns: `cohort`, `High`, `Low`, `pooled_cohort_share`, `estimate`, `cluster_SE`, `lower95`, `upper95`, `contribution_to_main`, `OECD_change`, `comparison_change`, `country_clusters`.

## `calendar_country_month_panel.csv`

Rows: 492. Columns: `country`, `month`, `aggregate`, `High`, `Low`, `exposure`, `endpoint`, `OECD`, `announcement_month_omitted`.

## `calendar_event_bins.csv`

Rows: 28. Columns: `endpoint`, `outcome`, `bin`, `label`, `start`, `end`, `reference`, `estimate`, `lower95`, `upper95`, `df`, `reference_bin`, `inference`.

## `Comparison_country_changes.csv`

Rows: 3. Columns: `country`, `change`.

## `controlled_projection_paths.csv`

Rows: 242. Columns: `outcome`, `method`, `month`, `observed_gap`, `predicted_gap`, `observed_OECD`, `realised_Brazil`, `projected_OECD_conditional_on_Brazil`, `omitted_month`.

## `Country_contributions_main_estimate.csv`

Rows: 3. Columns: `country`, `change`, `comparison_mean_change`, `relative_change`, `contribution_to_pooled_estimate`.

## `Exposure_split_country_group_period_rates.csv`

Rows: 24. Columns: `country`, `spending_group`, `phase`, `value`, `group`, `monthly_rate`.

## `Exposure_split_group_period_rates.csv`

Rows: 8. Columns: `group`, `spending_group`, `phase`, `rate`.

## `FDA_age_documented_funding_curves.csv`

Rows: 1,128. Columns: `cohort`, `country`, `months_after_FDA`, `denominator`, `documented_by_horizon`, `documented_percent`, `unresolved_by_horizon`, `upper_bound_percent`.

## `long_history_country_month.csv`

Rows: 219. Columns: `country`, `month`, `aggregate`, `exposure`.

## `long_history_source_coverage.csv`

Rows: 3. Columns: `country`, `required`, `observed`, `missing`.

## `main_country_counts.csv`

Rows: 6. Columns: `country`, `pre_entries`, `post_entries`, `pre_monthly`, `post_monthly`, `change`.

## `regression_cluster_moments.csv.gz`

Rows: 71,666. Columns: `model_id`, `country`, `term_i`, `term_j`, `crossproduct`.

## `regression_model_definitions.csv`

Rows: 107. Columns: `model_id`, `family`, `label`, `N`, `parameters`, `rank`, `target_column`, `country_clusters`.

## `Synthetic_control_country_weights.csv`

Rows: 5. Columns: `outcome`, `country_method`, `country`, `weight`.

## `Synthetic_control_monthly_series.csv`

Rows: 41. Columns: `date`, `Australia`, `Brazil`, `Hong Kong SAR`, `Singapore`, `Switzerland`, `OECD_age_balanced`, `synthetic_age_balanced`, `gap`.

## `weighting_comparison_results.csv`

Rows: 16. Columns: `outcome`, `country_method`, `unit`, `estimate`, `cluster_SE`, `lower95`, `upper95`, `countries_with_positive_weight`, `status`.

## Key meanings

- `pre` and `post`: country-specific mean monthly outcome in the stated period; `change` is post minus pre.
- `oecd`: treatment-group indicator for that specification; alternative country assignments are explicit.
- `effective_design_id`: identifies specifications generating the same effective contrast.
- `term_i`, `term_j`, `crossproduct`: upper-triangle entries of each country's X-transpose-X table; `term_j = -1` denotes X-transpose-y. Omitted entries are exactly zero.
- `N`, `rank`, `parameters`, `target_column`: design sample size, rank, column count and zero-based coefficient index.
- `aggregate`: recorded funding additions; `exposure`: the preweighted higher-minus-lower spending outcome.
- `documented_by_horizon` / `denominator`: final numerator and denominator for the documentation curve; the unresolved group is not labelled unfunded.
- Source-coverage rows record availability of monthly editions for the captured programmes, not a medicine-level unfunded status.
