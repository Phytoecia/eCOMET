# Perform ANOVA and Tukey's HSD test on the `mmo`

This function performs ANOVA and Tukey's HSD test on the feature data of
the `mmo`, Returns a list of ANOVA results, Tukey's HSD results, Tukey's
significance letters, and Dunnett's test results.

## Usage

``` r
anova_tukey_dunnett(df, formula)
```

## Arguments

  - df:
    
    The data frame containing the feature data and metadata

  - formula:
    
    The formula for the ANOVA test, e.g., "feature \~ group"

## Value

A list containing the ANOVA results, Tukey's HSD results, Tukey's
significance letters, and Dunnett's test results

## Examples

``` r
if (FALSE) {
anova_results <- anova_tukey_dunnett(df = merged_data, formula = "feature_value ~ group")
}
```
