# vCcCR: Vegetation/Canopy Cover Calculator in R

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A powerful R package for calculating vegetation and/or canopy cover ratios from raster data for polygon features. `vCcCR` is designed to handle large polygons efficiently through intelligent splitting, preventing memory issues and speeding up analysis.

---

## Installation

You can install the development version of `vCcCR` directly from GitHub using `devtools`:

```R
# Install devtools if you haven't already
# install.packages("devtools")

# Install vCcCR
devtools::install_github("DijoG/vCcCR")
```

## Usage example

Import dependencies:

```R
require(terra)        # For raster data handling
require(sf)           # For spatial vector data handling
require(exactextractr) # For fast raster exact extraction
require(dplyr)        # For data manipulation (e.g., mutate, filter, group_by, summarise)
require(progress)     # For progress bars
require(cli)          # For formatted command-line output
require(rlang)        # For quosures and unquoting (e.g., !!sym)
require(purrr)        # For functional programming (e.g., map, map_dfr)
require(tictoc)       # For timing code execution
require(tibble)       # For add_column() 
```

Compute vegetation/canopy cover ratio (%) using 10m resolution annual composites (with singe-feature input vector):

```R
# Input raster is an annual composite of mounthly mosaics (value 1 for vegetation/canopy, 0 for anything else) 
# Output is the updated inputSHAPE (here a single-feature vector file) with the computed VCr_date attributes)
vCcCR::get_VCr(inputRAST = ".../VC_Annual_2024_thr_0_15.tif",
               inputSHAPE = ".../02032025_Riyadh_METROPOLITAN.geojson", 
               outputSHAPE = ".../test/22052025_Riyadh_METROPOLITAN.geojson")
```

Compute vegetation/canopy cover ratio (%) using 10m resolution annual composites (with multi-feature input vector):

```R
vCcCR::get_VCr(inputRAST = ".../VC_2024/VC_Annual_2024_thr_0_15.tif",
        inputSHAPE = ".../0_2_Green Riyadh Project Boundaries/05112024_GRP_ARABIC — 20241105_GRP_ARABIC_DISSsel02.geojson", 
        outputSHAPE = ".../test/05112024_GRP_ARABIC — 20241105_GRP_ARABIC_DISSsel02.geojson",
        id_field = "NAME_ENGLI")
```

Compute vegetation/canopy cover ratio (%) using 0.35/0.3 resolution binarized raster file (with multi-feature complex vector):

```R
# A multi-featured mixed (large and small) vector file whose features are POLYGON
POLY <- sf::st_read(".../TestPoly.geojson") %>%
  sf::st_make_valid() %>%      
  sf::st_cast("POLYGON") 

# Binarized raster with pixel values 1 and NA/0, 1 indicating vagetetion or canopy
VCR <- terra::rast(".../VC_EPSG32638.tif")


# Run the vegetation/canopy ratio (%) computation:
tic("Total Vegetation Processing Time")
vCcCR::get_VEGETATION(
  polygons = POLY, 
  veg_raster = VCR,
  output_path = ".../CC_resultest.gpkg",
  split_threshold = 1850000,   # Polygons > 1.85 km^2 will be split in n_areas
  n_areas = 4,                 # Large polygons into n_areas sub-areas
  id_field = "NAME_ENGLI"      # Column name for polygon IDs
)
toc()
# Total Vegetation Processing Time: 1688.43 sec elapsed
```

  