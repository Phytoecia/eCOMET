# Write results of anova\_tukey\_dunnett to a CSV file

This function writes the results of ANOVA and Tukey's HSD test to a CSV
file.

## Usage

``` r
write_anova(anova_data, outdir, way = "oneway")
```

## Arguments

  - anova\_data:
    
    A list containing the results of ANOVA and Tukey's HSD test

  - outdir:
    
    The output directory where the results will be saved

  - way:
    
    The type of ANOVA test to perform. Options are 'oneway' or 'twoway'

## Examples

``` r
if (FALSE) {
write_anova(anova_data = anova_results, outdir = "anova_tukey_results.csv", way = 'oneway')
}
```
