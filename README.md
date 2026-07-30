# Replication package for "High-dimensional portfolio optimization using GARCH-EVT, R-vine Copula, and SCAD penalty with weight constrains"

Lv jian, Chenxu Wang and Wenyong Yuan

## Overview & contents

The code in this replication material generates the 56 figures and 9 tables for
the paper "High-dimensional portfolio optimization using GARCH-EVT, R-vine Copula, and SCAD penalty with weight constrains".
However, since the main research strategy itself has no other differences except for the cut-off point, the image generation is placed in the same code file. The specific differences can be simply adjusted by modifying the input objects and time lengths of the functions within the file.The data in the table can be presented through output display. However, to ensure accuracy, we use manual input for this part. The processing form for this section can be freely selected.

In terms of content, the specific positions of the pictures and tables are matched as follows：

- `plots/`: folder of generated plots as PDF files
- `tables/`: folder of generated tables as txt files
- `data-raw/`: folder of raw data files and the functions for processing them
- `data/`: folder of processed data files
- `Figure_[xx]_*.R`: R scripts to create the respective figures
- `Table_[xx]_*.R`: R scripts to create the respective tables


## Instructions & computational requirements.

For the sake of uniformity, all data file paths have been set as the absolute paths of the work computer. Please download the data files and then set your absolute path accordingly.===

Since intermediate data will be generated during the operation, to avoid any errors during the process, it is more recommended to conduct the tests in the order of the uploaded code.===

These analyses were run on R 4.3.1, and we explicitly use the following packages in the analysis files: `triptych` (0.1.2), `ggplot2` (3.4.3), `patchwork` (1.1.3), `dplyr` (1.1.3), `tidyr` (1.3.0), `purrr` (1.0.2), `grid` (base R), `lubridate` (1.9.2).

A comprehensive list of dependencies can be found in the `renv.lock` file. For a convenient setup in a (local) R session, we recommend using the `renv` package. The following steps are required once:
```
# install.packages("renv")
renv::activate()
renv::restore() # install dependencies
renv::status() # check environment
```

## Data availability and provenance

All the data are presented in tabular form, without any restrictions or limitations. Just use them directly.===


