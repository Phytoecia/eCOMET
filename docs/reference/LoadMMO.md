# Load a `mmo` previously saved with SaveMMO

This function returns the loaded `mmo` (visible return). By default it
prints basic information about the R version and recorded packages that
were present when the object was saved.

## Usage

``` r
LoadMMO(file, check_session = TRUE, verbose = TRUE)
```

## Arguments

  - file:
    
    Path to an RDS file created with SaveMMO

  - check\_session:
    
    Logical; if TRUE and save-time session info is present, print a
    short summary (default: TRUE)

  - verbose:
    
    Logical; print messages about saved session info when available
    (default: TRUE)

## Value

The loaded `mmo`
