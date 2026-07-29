# Save entire mmo object to a file (RDS)

Save entire mmo object to a file (RDS)

## Usage

``` r
SaveMMO(mmo, file = "mmo.rds", compress = "xz", include_session = TRUE)
```

## Arguments

  - mmo:
    
    The `mmo` to save

  - file:
    
    File path to write (default: "mmo.rds")

  - compress:
    
    Compression type passed to saveRDS ("gzip", "bzip2", "xz", or
    logical) (default: "xz")

  - include\_session:
    
    Logical; if TRUE attach sessionInfo() as an attribute to the saved
    object (default: TRUE)

## Value

Invisibly returns the file path
