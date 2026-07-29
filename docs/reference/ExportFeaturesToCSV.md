# ExportFeaturesToCSV

This function exports selected features, their annotations, and pairwise
comparisons to a CSV file.

## Usage

``` r
ExportFeaturesToCSV(mmo, id_list, normalization = "None", output_dir)
```

## Arguments

  - mmo:
    
    The `mmo` containing feature data, annotations, and pairwise
    comparisons

  - id\_list:
    
    A list of feature names to filter and export

  - normalization:
    
    The normalization method to use for feature data. Options are
    'None', 'Log', 'Meancentered', or 'Z' (default: 'None')

  - output\_dir:
    
    The output directory to save the CSV file

## Examples

``` r
if (FALSE) {
ExportFeaturesToCSV(mmo,
                    id_list = Glucosinolates,
                     normalization = 'Z',
                     output_dir = 'output.csv')
ExportFeaturesToCSV(mmo, id_list = DAMs_up$control_vs_treatment1.up,
                         normalization = 'None',
                          output_dir = 'output.csv')
}
```
