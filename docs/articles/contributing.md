# How to Contribute

## Contributing to eCOMET

eCOMET is a community-driven tool built by ecologists for ecologists.
The package grows most usefully when the people using it in real
workflows contribute ideas, report problems, and share new methods. This
page explains how to do that.

### Design philosophy

Before contributing, it helps to understand what eCOMET is trying to be:

  - **The `mmo`, the S3 Class Object, is the center of everything in the
    eCOMET package** All contributed functions should work with an `mmo`
    as its primary input and either return a modified `mmo` or a
    plot/data.frame or list. This keeps workflows composable and
    readable.
  - **Reproducibility over interactivity.** Functions should work in
    scripts, vignettes, and automated pipelines — not just interactively
    in RStudio.
  - **ggplot2 for visualization.** We designed eCOMET around ggplot so
    that users could use ggplot conventions to keep modifying plots
    after they are produced by eCOMET functions. New plot functions
    should return a `ggplot`/`ggtree` object so users can layer
    modifications with `+` without rewriting the function.

-----

### Conventions for new functions

#### Function structure

All new analysis or visualization functions should follow this pattern:

``` r
MyNewFunction <- function(mmo, ..., group_col = NULL) {
  # validate inputs
  # do work
  # return modified mmo OR a ggplot object OR a data.frame
}
```

  - Functions that add new results to the dataset should store them
    inside the `mmo` object and return the updated `mmo`.
  - Functions that produce a plot should return a `ggplot` or `ggtree`
    object — never call `print()` or open a graphics device inside the
    function. You can however save the plot as an output.
  - Functions that produce a summary table should return a plain
    `data.frame`.

#### Plot functions

  - Use `ggplot2` (or `ggtree` for tree layouts) as the base.
  - Accept a `palette` argument using `colorspace::qualitative_hcl()`
    for color scales.
  - Return the plot object so users can modify it with `+`.
  - Default to showing group-level information derived from
    `mmo$metadata`.

<!-- end list -->

``` r
# Good — user can keep layering
p <- MyPlotFunction(mmo, group_col = "treatment")
p + ggplot2::theme_minimal() + ggplot2::labs(title = "My title")
```

#### Documentation

Every exported function needs roxygen2 documentation with at minimum:

``` r
#' Short one-line title
#'
#' One paragraph description of what the function does and when to use it.
#'
#' @param mmo An mmo object created by \code{GetMZmineFeature()}.
#' @param ... Other parameters
#' @return Description of what is returned.
#' @export
#' @examplesIf FALSE
#' result <- MyNewFunction(mmo)
```

Use `@examplesIf FALSE` rather than `@examples` unless the example can
run in under a few seconds using only built-in demo data.

#### Dependencies

  - Before adding a new package dependency, check whether the
    functionality already exists in a package eCOMET already imports
    (see `DESCRIPTION`).
  - Heavy or specialized packages (e.g. those requiring system libraries
    or Bioconductor) should go in `Suggests`, not `Imports`, and be
    loaded with `.require_pkg("pkgname")` inside the function so they
    are only required when actually called.

-----

### Sharing a complete new workflow

If you have developed a full analysis workflow using eCOMET — a new type
of study design, a combination of functions for a specific ecological
question, or an integration with another tool — consider contributing it
as a tutorial vignette.

Vignettes live in `vignettes/` and are built into the website. A good
tutorial vignette:

  - Loads data using `system.file("extdata/...", package = "ecomet")` so
    it is portable.
  - Walks through a complete, realistic analysis from raw data to
    interpretation.
  - Explains the ecological motivation for each step, not just the code.
  - Uses demo data small enough to be bundled in the package (target \<
    5 MB total per tutorial).

-----

### Questions

If you are unsure whether an idea fits, or want to discuss an approach
before writing code, open a
[Discussion](https://github.com/Phytoecia/eCOMET/discussions) on GitHub.
This is the best place for open-ended questions about directions for the
package.
