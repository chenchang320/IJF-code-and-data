# Replication package for "High-dimensional portfolio optimization using GARCH-EVT, R-vine Copula, and SCAD penalty with weight constrains"

Lv jian, Chenxu Wang and Wenyong Yuan

## Overview & contents

The replication code provided herein is designed to generate all 56 figures and 9 tables presented in the paper “High-dimensional portfolio optimization using GARCH-EVT, R-vine Copula, and SCAD penalty with weight constraints”. The empirical strategy remains consistent throughout, with the sole variation being the selection of truncation points. To accommodate this, all graphical outputs are produced within a unified script, where the key distinctions are governed by adjustable parameters—specifically, the input objects and the length of the estimation window—within the respective functions. Regarding the tabular results, while they can be directly outputted from the code, we opt for manual transcription to ensure the utmost accuracy and to allow for flexible post-hoc formatting adjustments.===

The following list specifies the correspondence between each figure/table and the code file from which it is generated：===

- `plots/`: folder of generated plots as PDF files
- `tables/`: folder of generated tables as txt files
- `data-raw/`: folder of raw data files and the functions for processing them
- `data/`: folder of processed data files
- `Figure_[xx]_*.R`: R scripts to create the respective figures
- `Table_[xx]_*.R`: R scripts to create the respective tables


## Instructions & computational requirements.

For consistency, all file paths in the code are specified as absolute paths referencing our local working environment. Users are therefore advised to download the data files and reconfigure these paths to match their own directory structures.===

Additionally, as the routine generates intermediate datasets during execution, we strongly recommend running the scripts sequentially in the order they are provided to prevent potential runtime errors.===

These analyses were run on R 4.3.1, and we explicitly use the following packages in the analysis files: `triptych` (0.1.2), `ggplot2` (3.4.3), `patchwork` (1.1.3), `dplyr` (1.1.3), `tidyr` (1.3.0), `purrr` (1.0.2), `grid` (base R), `lubridate` (1.9.2).

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Data availability and provenance

All data necessary for replicating the empirical results are provided in tabular format within this replication package, with no access restrictions or usage limitations. These data are ready for direct use in the accompanying code.===


