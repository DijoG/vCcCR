#' Vegetation Canopy Cover Calculator in R (vCcCR)
#'
#' A package for calculating vegetation and/or canopy  cover ratios from raster data for polygon features.
#' Supports processing of large polygons through intelligent splitting.
#'
#' @keywords internal
"_PACKAGE"

#' Calculates Vegetation Cover Ratio (VCr) for each month in a raster stack (annual composite of monthly mosaics)
#'
#' @param inputRAST Path to the input raster stack (12 monthly vegetation cover layers)
#' @param inputSHAPE Path to the input shapefile (polygons for analysis)
#' @param outputSHAPE Path and name (without file extension!) where to save the output spatial files in gpkg and shp formats (default = NULL ~ same as inputSHAPE with '_VCr.gpkg' extension)
#' @param id_field Attribute name or id to keep, all other fields are removed (default = NULL ~ keeping first attribute)
#' @return An sf object with added VCr columns (invisibly), writes out a .geojson and a .shp 
#' @export 
get_VCratio <- function(inputRAST, 
                        inputSHAPE, 
                        outputSHAPE = NULL, 
                        id_field = NULL) {
  
  # Load required packages
  require(terra)
  require(sf)
  require(exactextractr)
  require(dplyr)
  require(purrr)
  require(tictoc)
  
  # Load data
  message("Loading raster data...")
  VCRAST = rast(inputRAST)
  
  message("Loading shapefile...")
  SHP = st_read(inputSHAPE, quiet = TRUE)
  if (!is.null(id_field)) {
    SHP = SHP %>%
      dplyr::select(all_of(id_field))
  } else {
    SHP = SHP %>%
      dplyr::select(1)
  }
  
  # Check CRS match
  if(crs(VCRAST) != st_crs(SHP)$wkt) {
    message("Reprojecting shapefile to match raster CRS...")
    SHP = st_transform(SHP, crs(VCRAST))
  }
  
  # Obtain pixel size
  pixel_size = res(VCRAST)[1]
  
  # Create output object
  result_shp = SHP
  
  # Process each monthly layer
  message("\nProcessing monthly vegetation cover areas:")
  for(i in 1:nlyr(VCRAST)) {
    month_name = names(VCRAST)[i]
    vcr_colname = paste0("VCr", month_name)
    
    tictoc::tic(paste("Month", month_name))
    
    # Extract vegetation pixels (value = 1)
    veg_counts = exact_extract(
      VCRAST[[i]], 
      SHP,
      fun = function(values, coverage_fractions) {
        sum(values == 1, na.rm = TRUE)
      }
    )
    
    # Calculate VCr and add to result
    result_shp[[vcr_colname]] = round(
      (veg_counts * (pixel_size^2)) / as.numeric(st_area(SHP)) * 100, 2
    )
    
    # Add processing time message
    tictoc::toc()
  }
  
  # Save results
  if (is.null(outputSHAPE)) {
    outputSHAPEgpkg = paste0(sub("\\.(gpkg|geojson|shp)$", "", inputSHAPE), "_VCr.gpkg")
    outputSHAPEshp = paste0(sub("\\.(gpkg|geojson|shp)$", "", inputSHAPE), "_VCr.shp")
  } else {
    outputSHAPEgpkg = paste0(sub("\\.(gpkg|geojson|shp)$", "", outputSHAPE), "_VCr.gpkg")
    outputSHAPEshp = paste0(sub("\\.(gpkg|geojson|shp)$", "", outputSHAPE), "_VCr.shp")
  }
  message("\nSaving output files...")
  st_write(result_shp, outputSHAPEgpkg, delete_layer = TRUE, quiet = TRUE)
  st_write(result_shp, outputSHAPEshp, delete_layer = TRUE, quiet = TRUE)
  message("Successfully saved to:\n", normalizePath(outputSHAPEgpkg), "\n", normalizePath(outputSHAPEshp))
  
  # Return the sf object invisibly
  invisible(result_shp)
}

#' Calculates Vegetation Cover Area (VCa) for each month in a raster stack (annual composite of monthly mosaics)
#'
#' @param inputRAST Path to the input raster stack (12 monthly vegetation cover layers)
#' @param inputSHAPE Path to the input shapefile (polygons for analysis)
#' @param outputSHAPE Path and name (without file extension!) where to save the output spatial files in gpkg and shp formats (default = NULL ~ same as inputSHAPE with '_VCr.gpkg' extension)
#' @param id_field Attribute name or id to keep, all other fields are removed (default = NULL ~ keeping first attribute)
#' @return An sf object with added VCa columns (invisibly), writes out a .geojson and a .shp 
#' @export 
get_VCarea <- function(inputRAST, 
                       inputSHAPE, 
                       outputSHAPE = NULL, 
                       id_field = NULL) {
  
  # Load required packages
  require(terra)
  require(sf)
  require(exactextractr)
  require(dplyr)
  require(purrr)
  require(tictoc)
  
  # Load data
  message("Loading raster data...")
  VCRAST = rast(inputRAST)
  
  message("Loading shapefile...")
  SHP = st_read(inputSHAPE, quiet = TRUE)
  if (!is.null(id_field)) {
    SHP = SHP %>%
      dplyr::select(all_of(id_field))
  } else {
    SHP = SHP %>%
      dplyr::select(1)
  }
  
  # Check CRS match
  if(crs(VCRAST) != st_crs(SHP)$wkt) {
    message("Reprojecting shapefile to match raster CRS...")
    SHP = st_transform(SHP, crs(VCRAST))
  }
  
  # Obtain pixel size
  pixel_size = res(VCRAST)[1]
  
  # Create output object
  result_shp = SHP
  
  # Process each monthly layer
  message("\nProcessing monthly vegetation cover ratios:")
  for(i in 1:nlyr(VCRAST)) {
    month_name = names(VCRAST)[i]
    vcr_colname = paste0("VCa", month_name)
    
    tictoc::tic(paste("Month", month_name))
    
    # Extract vegetation pixels (value = 1)
    veg_counts = exact_extract(
      VCRAST[[i]], 
      SHP,
      fun = function(values, coverage_fractions) {
        sum(values == 1, na.rm = TRUE)
      }
    )
    
    # Calculate VCr and add to result
    result_shp[[vcr_colname]] = round(
      veg_counts * (pixel_size^2), 2)
    
    # Add processing time message
    tictoc::toc()
  }
  
  # Save results
  if (is.null(outputSHAPE)) {
    outputSHAPEgpkg = paste0(sub("\\.(gpkg|geojson|shp)$", "", inputSHAPE), "_VCa.gpkg")
    outputSHAPEshp = paste0(sub("\\.(gpkg|geojson|shp)$", "", inputSHAPE), "_VCa.shp")
  } else {
    outputSHAPEgpkg = paste0(sub("\\.(gpkg|geojson|shp)$", "", outputSHAPE), "_VCa.gpkg")
    outputSHAPEshp = paste0(sub("\\.(gpkg|geojson|shp)$", "", outputSHAPE), "_VCa.shp")
  }
  message("\nSaving output files...")
  st_write(result_shp, outputSHAPEgpkg, delete_layer = TRUE, quiet = TRUE)
  st_write(result_shp, outputSHAPEshp, delete_layer = TRUE, quiet = TRUE)
  message("Successfully saved to:\n", normalizePath(outputSHAPEgpkg), "\n", normalizePath(outputSHAPEshp))
  
  # Return the sf object invisibly
  invisible(result_shp)
}

#' Validate and transform CRS between polygons and raster
#'
#' Ensures that the Coordinate Reference Systems (CRSs) of the input polygons
#' and raster match. If they don't, polygons are transformed to the raster's CRS.
#'
#' @param poly An sf object (polygons).
#' @param raster A SpatRaster object.
#' @return The sf object with its CRS validated and potentially transformed.
#' @noRd  
.validate_crs <- function(poly, raster) {
  poly_crs = sf::st_crs(poly)
  raster_crs = sf:: st_crs(terra::crs(raster))
  
  if (is.na(poly_crs) || is.na(raster_crs)) {
    stop("Either polygons or raster has no CRS defined. Please ensure both have valid CRS information.")
  }
  
  if (!identical(poly_crs, raster_crs)) {
    cli::cli_alert_warning("CRS mismatch detected: polygons ({poly_crs$input}) vs. raster ({raster_crs$input}).")
    cli::cli_alert_info("Transforming polygons to match raster CRS...")
    poly = sf::st_transform(poly, raster_crs)
    cli::cli_alert_success("Polygons successfully transformed.")
  } else {
    cli::cli_alert_info("CRS of polygons and raster match ({poly_crs$input}). No transformation needed.")
  }
  
  return(poly)
}

#' Splits a large polygon into smaller sub-polygons using ADAPTIVE grid sizing method.
#'
#' This is useful for processing very large polygons that might cause memory
#' issues or be too slow for direct raster extraction.
#'
#' @param poly An sf object containing a single polygon.
#' @param n_areas An integer, the desired number of sub-polygons to create, 
#' @param id_field A character string, the name of the ID field in `poly`.
#' @return An sf object containing the smaller split polygons, inheriting attributes from the original polygon.
#' @noRd  
.split_POLY <- function(poly, n_areas, id_field) {
  if (nrow(poly) != 1) stop("split_POLY expects a single polygon")
  
  orig_id = poly[[id_field]]
  
  # Adaptive grid sizing based on polygon complexity
  bbox = sf::st_bbox(poly)
  width = bbox$xmax - bbox$xmin
  height = bbox$ymax - bbox$ymin
  
  # Aspect ratio aware grid - better for long/narrow polygons
  aspect_ratio = width / height
  if (aspect_ratio > 3) {
    # Long and narrow - more cells along the long axis
    grid_dims = c(ceiling(sqrt(n_areas * aspect_ratio)), ceiling(sqrt(n_areas / aspect_ratio)))
  } else if (aspect_ratio < 0.33) {
    # Tall and narrow
    grid_dims = c(ceiling(sqrt(n_areas * aspect_ratio)), ceiling(sqrt(n_areas / aspect_ratio)))
  } else {
    # Roughly square
    grid_dims = c(ceiling(sqrt(n_areas)), ceiling(sqrt(n_areas)))
  }
  
  # Ensure minimum dimensions
  grid_dims = pmax(grid_dims, 2)
  
  cli::cli_alert_info("Using {grid_dims[1]}x{grid_dims[2]} grid for complex polygon")
  
  grid = sf::st_make_grid(poly, n = grid_dims)
  grid_sf = sf::st_sf(geometry = grid)
  sf::st_crs(grid_sf) = sf::st_crs(poly)
  
  # Robust intersection with complex polygons
  splitted_poly = tryCatch({
    result = sf::st_intersection(grid_sf, poly) %>%
      filter(sf::st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>%
      mutate(
        PolyArea = as.numeric(sf::st_area(.)),
        !!id_field := orig_id
      ) %>%
      select(all_of(id_field), PolyArea) %>%
      sf::st_make_valid() %>%
      sf::st_collection_extract("POLYGON")
    
    # Filter out tiny slivers that can occur with complex shapes
    area_threshold = max(result$PolyArea) * 0.01  # 1% of largest sub-polygon
    result = result %>% filter(PolyArea >= area_threshold)
    
    result
  }, error = function(e) {
    cli::cli_alert_danger("Grid intersection failed: {e$message}")
    cli::cli_alert_info("Falling back to bounding box grid")
    
    # Ultimate fallback - use bounding box
    bbox_poly = sf::st_as_sfc(sf::st_bbox(poly))
    grid = sf::st_make_grid(bbox_poly, n = c(3, 3))
    grid_sf = sf::st_sf(geometry = grid)
    sf::st_crs(grid_sf) = sf::st_crs(poly)
    
    sf::st_intersection(grid_sf, poly) %>%
      filter(sf::st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>%
      mutate(
        PolyArea = round(as.numeric(sf::st_area(.)), 2),
        !!id_field := orig_id
      ) %>%
      select(all_of(id_field), PolyArea) %>%
      sf::st_make_valid() %>%
      sf::st_collection_extract("POLYGON")
  })
  
  return(splitted_poly)
}

#' Processes vegetation/canopy cover for a set of polygons using a vegetation/canopy raster
#' derived from very high resolution multispectral instrument.
#' 
#' Handles large polygons by splitting them into smaller chunks. Runs in sequential mode only.
#'
#' @param polygons An sf object containing polygons.
#' @param veg_raster A SpatRaster object representing vegetation/canopy cover with pixel values 1.
#' @param output_path A character string, the path to save the results.
#' @param split_threshold A numeric value (m^2), polygons larger than this will be split into n_areas.
#' @param n_areas An integer, the number of sub-polygons to split large polygons into.
#' @param id_field A character string, the name of the unique ID field in the polygons that has to be kept (default = NULL, taking the first field)
#' @param by_ROW Logical, TRUE triggers feature by feature processing to avoid RAM overhead (default = FALSE).
#' @param return Logical, default = FALSE.
#' @return An sf object with id_field plus calculated VegArea and VegRatio attributes added.
#' @export 
get_VEGETATION <- function(polygons,
                           veg_raster,
                           output_path,
                           split_threshold = 1850000,
                           n_areas = 8,
                           id_field = NULL,
                           by_ROW = FALSE,
                           return = FALSE) {
  
  # Start processing message
  cli::cli_h1("Vegetation Cover Analysis")
  cli::cli_alert_info("Starting processing at {Sys.time()}")
  
  # Input validation and data loading
  if (is.character(polygons)) {
    polygons = sf::st_read(polygons, quiet = TRUE)
  }
  if (!inherits(polygons, "sf")) {
    stop("The 'polygons' object must be an sf object or a path to a spatial file.")
  }
  
  if (is.character(veg_raster)) {
    veg_raster = terra::rast(veg_raster)
  }
  if (!inherits(veg_raster, "SpatRaster")) {
    stop("The 'veg_raster' object must be a SpatRaster or a path to a raster file.")
  }
  
  # Validate CRS
  polygons = .validate_crs(polygons, veg_raster)
  
  # Check if id_field exists
  if (!is.null(id_field)) {
    if (!id_field %in% names(polygons)) {
    stop("The specified id_field '", id_field, "' does not exist in the polygons.")
    }
  } else {
    id_field = names(polygons)[1]
  }
  
  # Create output directory if needed
  if (!dir.exists(dirname(output_path))) {
    dir.create(dirname(output_path), recursive = TRUE)
    cli::cli_alert_info("Created output directory: {dirname(output_path)}")
  }
  
  # Calculate polygon areas if not already present
  if (!"PolyArea" %in% names(polygons)) {
    polygons$PolyArea = round(as.numeric(sf::st_area(polygons)), 2)
  }
  
  # Obtain pixel size
  pix_size = terra::res(veg_raster)[1]
  
  # Processing logic
  if (by_ROW) {
    cli::cli_alert_info("Processing features one by one with individual splitting checks")
    pb = progress::progress_bar$new(
      format = "Processing: [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = nrow(polygons),
      clear = FALSE,
      width = 60
    )
    
    # Pre-allocate result columns for efficiency
    polygons$VegArea = NA_real_
    polygons$VegRatio = NA_real_
    
    for (row in seq_len(nrow(polygons))) {
      pb$tick()
      current_poly = polygons[row, ]
      parea = current_poly$PolyArea
      
      # Check if splitting is needed for this feature
      if ((parea > split_threshold) | (current_poly[[id_field]] == "Sports Boulevard")) {
        cli::cli_alert_info("Splitting polygon: {current_poly[[id_field]]}")
        poly_parts = .split_POLY(current_poly, n_areas, id_field)
      } else {
        poly_parts = current_poly
      }
      
      # Extract vegetation data
      veg_pixels = exactextractr::exact_extract(veg_raster, poly_parts) %>%
        purrr::map_int(~ sum(.x$value, na.rm = TRUE)) %>%
        sum(.)
      
      # Calculate and store results
      polygons$VegArea[row] = round(veg_pixels * (pix_size * pix_size), 2)
      polygons$VegRatio[row] = round(polygons$VegArea[row] / parea * 100, 2)
      
      cli::cli_alert_success("Processed {current_poly[[id_field]]}: Area = {round(parea)}, VegArea = {polygons$VegArea[row]}, VegRatio = {polygons$VegRatio[row]}")
    }
    
  } else {
    cli::cli_alert_info("Processing all features at once with uniform splitting criteria")
    
    # Separate large and small polygons
    large_polys = polygons %>% dplyr::filter(PolyArea > split_threshold)
    small_polys = polygons %>% dplyr::filter(PolyArea <= split_threshold)
    
    cli::cli_alert_info("Polygons exceeding split threshold ({split_threshold} m²): {nrow(large_polys)}")
    cli::cli_alert_info("Polygons below threshold: {nrow(small_polys)}")
    
    all_results = list() # To store results from both large and small polygons
    
    # Process large polygons with splitting
    if (nrow(large_polys) > 0) {
      cli::cli_h2("Processing Large Polygons")
      pb_large = progress::progress_bar$new(
        format = "Large Polygons: [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
        total = nrow(large_polys),
        clear = FALSE,
        width = 60
      )
      
      large_computed_results = purrr::map_dfr(
        1:nrow(large_polys),
        function(row) {
          pb_large$tick()
          current_poly = large_polys[row, ]
          poly_parts = .split_POLY(current_poly, n_areas, id_field)
          
          veg_pixels = exactextractr::exact_extract(veg_raster, poly_parts) %>%
            purrr::map_int(~ sum(.x$value, na.rm = TRUE)) %>%
            sum(.)
          
          tibble::tibble(
            !!id_field := current_poly[[id_field]],
            PolyArea_calc = current_poly$PolyArea, # Store original PolyArea for later ratio calculation if needed
            VegArea = round(veg_pixels * (pix_size * pix_size), 2)
          )
        }
      )
      all_results[[length(all_results) + 1]] = large_computed_results
    }
    
    # Process small polygons without splitting
    if (nrow(small_polys) > 0) {
      cli::cli_h2("Processing Small Polygons")
      pb_small = progress::progress_bar$new(
        format = "Small Polygons: [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
        total = nrow(small_polys),
        clear = FALSE,
        width = 60
      )
      
      small_computed_results = purrr::map_dfr(
        1:nrow(small_polys),
        function(row) {
          pb_small$tick()
          current_poly = small_polys[row, ]
          
          veg_pixels = exactextractr::exact_extract(veg_raster, current_poly) %>%
            purrr::map_int(~ sum(.x$value, na.rm = TRUE)) %>%
            sum(.)
          
          tibble::tibble(
            !!id_field := current_poly[[id_field]],
            PolyArea_calc = current_poly$PolyArea, # Store original PolyArea
            VegArea = round(veg_pixels * (pix_size * pix_size), 2)
          )
        }
      )
      all_results[[length(all_results) + 1]] = small_computed_results
    }
    
    # Combine all results and calculate VegRatio
    if (length(all_results) > 0) {
      combined_results = dplyr::bind_rows(all_results) %>%
        dplyr::mutate(VegRatio = round(VegArea / PolyArea_calc * 100, 2)) %>%
        dplyr::select(-PolyArea_calc) # Remove temporary PolyArea_calc column
      
      # Join the results back to the original polygons object
      polygons = polygons %>%
        dplyr::left_join(combined_results, by = id_field)
    }
  }
  
  # Save output
  ext = tools::file_ext(output_path)
  if (tolower(ext) == "gpkg") {
    sf::st_write(polygons, output_path, delete_layer = TRUE, quiet = TRUE) 
  } else if (tolower(ext) %in% c("geojson", "json")) {
    sf::st_write(polygons, output_path, delete_dsn = TRUE, quiet = TRUE) 
  } else {
    stop("Unsupported output format. Please use '.gpkg' or '.geojson'")
  }
  
  cli::cli_alert_success("Results saved to: {output_path}")
  
  # Final summary
  total_veg_area = sum(polygons$VegArea, na.rm = TRUE)
  cli::cli_alert_info("Total vegetation area calculated: {round(total_veg_area / 1000000, 2)} km²")
  
  if (return) return(polygons)
}
