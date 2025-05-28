# vCcCR: Vegetation Canopy Cover Calculator in R

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

Import data, POLY with multiple features, and VCR the binarized raster file:

```R
# A multi-feautured vector file whose featires are POLYGON
POLY <- sf::st_read(".../TestPoly.geojson") %>%
  sf::st_make_valid() %>%      
  sf::st_cast("POLYGON") 

# Vegetation Cover (raster with pixel values 1 and NA/0, 1 indicating vagetetion or canopy)
VCR <- terra::rast(".../VC_EPSG32638.tif")
```

Run the vegetation processing for mixed (large and small) or only large or small polygons:

```R
# Example: Run the vegetation processing for large polys
tic("Total Vegetation Processing Time")
get_VEGETATION(
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

  