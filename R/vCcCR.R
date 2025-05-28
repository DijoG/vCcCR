#' Vegetation Canopy Cover Calculator in R (vCcCR)
#'
#' A package for calculating vegetation and/or canopy  cover ratios from raster data for polygon features.
#' Supports processing of large polygons through intelligent splitting.
#'
#' @docType package
#' @name vCcCR
NULL

#' Calculate Vegetation Cover Ratio (VCr) for each month in a raster stack (annual composite of monthly mosaics)
#'
#' @param inputRAST Path to the input raster stack (12 monthly vegetation cover layers)
#' @param inputSHAPE Path to the input shapefile (polygons for analysis)
#' @param outputSHAPE Path where to save the output shapefile with VCr attributes (default = NULL ~ same as inputSHAPE with '_VCr.geojson' extension)
#' @param id_field Attribute name or id to keep, all other fields are removed (default = NULL ~ keeping first attribute)
#' @return An sf object with added VCr columns (invisibly), writes out a .geojson and a .shp 
#' @export 
get_VCr <- function(inputRAST, 
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
    vcr_colname = paste0("VCr_", month_name)
    
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
    outputSHAPEgeoj = paste0(sub("\\.(geojson|shp)$", "", inputSHAPE), "_VCr.geojson")
    outputSHAPEshp = paste0(sub("\\.(geojson|shp)$", "", inputSHAPE), "_VCr.shp")
  } else {
    outputSHAPEgeoj = paste0(sub("\\.(geojson|shp)$", "", outputSHAPE), "_VCr.geojson")
    outputSHAPEshp = paste0(sub("\\.(geojson|shp)$", "", outputSHAPE), "_VCr.shp")
  }
  message("\nSaving output shapefile...")
  st_write(result_shp, outputSHAPEgeoj, delete_layer = TRUE, quiet = TRUE, append = FALSE)
  st_write(result_shp, outputSHAPEshp, delete_layer = TRUE, quiet = TRUE, append = FALSE)
  message("Successfully saved to:\n", normalizePath(outputSHAPEgeoj), "\n", normalizePath(outputSHAPEshp))
  
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

#' Splits a large polygon into smaller sub-polygons using k-means and Voronoi tessellation.
#'
#' This is useful for processing very large polygons that might cause memory
#' issues or be too slow for direct raster extraction.
#'
#' @param poly An sf object containing a single polygon.
#' @param n_areas An integer, the desired number of sub-polygons to create.
#' @param id_field A character string, the name of the ID field in `poly`.
#' @return An sf object containing the smaller split polygons, inheriting attributes from the original polygon.
#' @noRd  
.split_POLY <- function(poly, n_areas, id_field) {
  # Ensure we have exactly one polygon
  if (nrow(poly) != 1) {
    stop("split_POLY expects a single polygon")
  }
  
  # Store original ID
  orig_id = poly[[id_field]]
  
  # Sample points and create voronoi polygons
  points = sf::st_sample(poly, size = min(200, n_areas * 50))
  points_df = as.data.frame(sf::st_coordinates(points))
  
  k_means = kmeans(points_df, centers = n_areas)
  voronoi_polys = dismo::voronoi(k_means$centers, ext = poly)
  
  # Convert to sf and set CRS
  voronoi_sf = sf::st_as_sf(voronoi_polys)
  sf::st_crs(voronoi_sf) = sf::st_crs(poly)
  
  # Intersect with original polygon
  splitted_poly = sf::st_intersection(voronoi_sf, poly) %>%
    mutate(
      PolyArea = as.numeric(sf::st_area(.)),
      !!id_field := orig_id  # Preserve original ID
    ) %>%
    select(all_of(id_field), PolyArea)
  
  # Ensure valid geometries and return as single polygons
  splitted_poly = splitted_poly %>%
    sf::st_make_valid() %>%
    sf::st_collection_extract("POLYGON")
  
  return(splitted_poly)
}

#' Calculates vegetation area and ratio for input polygons.
#'
#' Uses exact_extract to count raster pixels within polygons and calculates
#' vegetation area and ratio based on pixel size. Includes robust error handling.
#'
#' @param poly An sf object (single or multi-feature polygon).
#' @param rast A SpatRaster object.
#' @param id_field A character string, the name of the ID field in `poly`.
#' @return A tibble with the ID field, calculated VegArea, and VegRatio.
#' @noRd  
.calculate_vegetation <- function(poly, rast, id_field, split = FALSE, n_areas = 10) {
  # Progress bar tick handled by calling function
  
  # Validate inputs
  if (!inherits(poly, "sf")) stop("`poly` must be an sf object.")
  if (!inherits(rast, "SpatRaster")) rast = terra::rast(rast) # Ensure it's a SpatRaster
  
  # Ensure id_field exists in poly
  if (!(id_field %in% names(poly))) {
    warning(paste0("ID field '", id_field, "' not found in polygon attributes for some features. Returning NA results."))
    # Return a tibble with NA values but correct structure and number of rows
    return(tibble(
      !!rlang::sym(id_field) := rep(NA_character_, nrow(poly)), 
      VegArea = rep(NA_real_, nrow(poly)),
      VegRatio = rep(NA_real_, nrow(poly))
    ))
  }
  # Perform exact extraction
  if (split) {
    counts = exactextractr::exact_extract(
      rast,
      .split_POLY(poly, n_areas, id_field) 
    ) %>%
      map_int(~ sum(.x$value == 1, na.rm = TRUE)) %>%
      sum()
  } else {
    counts = exactextractr::exact_extract(
      rast,
      poly
    ) %>%
      map_int(~ sum(.x$value == 1, na.rm = TRUE))
  }
  
  # Initialize veg_area and veg_ratio with zeros, matching the number of polygons
  veg_area = round(counts*(.35*.35), 2)
  poly_area = round(as.numeric(sf::st_area(poly)), 2)
  veg_ratio = round(veg_area/poly_area * 100, 2)
  
  polyout = 
    poly %>%
    dplyr::mutate(PolyArea = poly_area,
                  VegArea = veg_area,
                  VegRatio = veg_ratio)
  
  return(polyout)
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
#' @param id_field A character string, the name of the unique ID field in the polygons that has to be kept.
#' @param return Logical, default = FALSE
#' @return An sf object with id_field plus calculated VegArea and VegRatio attributes added.
#' @export 
get_VEGETATION <- function(polygons, 
                           veg_raster, 
                           output_path,
                           split_threshold = 1850000, 
                           n_areas = 8,
                           id_field = "NAME_ENGLI",
                           return = FALSE) {
  
  cli::cli_h1("Vegetation/Canopy Cover Analysis")
  cli::cli_alert_info("Starting processing at {Sys.time()}")
  
  # Load data if paths are provided
  if (is.character(polygons)) {
    polygons = sf::st_read(polygons, quiet = TRUE)
  }
  if (is.character(veg_raster)) {
    veg_raster = terra::rast(veg_raster)
  }
  
  cli::cli_alert_info("Total input polygons: {nrow(polygons)} features")
  
  # Validate CRS and ensure input is sf object
  polygons = .validate_crs(polygons, veg_raster)
  if (!inherits(polygons, "sf")) {
    polygons = sf::st_as_sf(polygons)
  }
  
  # Ensure the id_field exists
  if (!(id_field %in% names(polygons))) {
    stop(paste0("Error: The specified `id_field` '", id_field, "' is not found in the input polygons."))
  }
  
  # Pre-process polygons
  initial_poly_count = nrow(polygons)
  polygons = polygons %>%
    dplyr::mutate(PolyArea = as.numeric(st_area(geometry))) %>%
    dplyr::filter(PolyArea > 0)
  
  if (nrow(polygons) < initial_poly_count) {
    cli::cli_alert_warning("Filtered out {initial_poly_count - nrow(polygons)} polygon(s) with zero or invalid area.")
  }
  
  # Separate large and small polygons
  large_polys = polygons %>% dplyr::filter(PolyArea > split_threshold)
  small_polys = polygons %>% dplyr::filter(PolyArea <= split_threshold)
  
  cli::cli_alert_info("Polygons exceeding split threshold ({split_threshold} m²): {nrow(large_polys)}")
  cli::cli_alert_info("Polygons below or equal to split threshold: {nrow(small_polys)}")
  
  # Process polygons with error handling
  results = list()
  
  if (nrow(large_polys) > 0) {
    cli::cli_h2("Processing Large Polygons")
    pb_large = progress::progress_bar$new(
      format = "Large Polygons: [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = nrow(large_polys),
      clear = FALSE,
      width = 60
    )
    
    # Process each large polygon individually with proper aggregation
    large_results = purrr::map_dfr(
      1:nrow(large_polys),
      function(i) {
        pb_large$tick()
        poly_single = large_polys[i, , drop = FALSE]
        
        # Calculate stats for this polygon (with splitting)
        result = .calculate_vegetation(
          poly_single, 
          veg_raster, 
          id_field, 
          split = TRUE, 
          n_areas = n_areas
        )
        
        # Ensure we return exactly one row per original polygon
        result %>%
          dplyr::group_by(!!rlang::sym(id_field)) %>%
          dplyr::summarise(
            VegArea = sum(VegArea, na.rm = TRUE),
            VegRatio = mean(VegRatio, na.rm = TRUE),
            .groups = "drop"
          )
      }
    )
    
    results$large = large_results
    cli::cli_alert_success("Completed large polygon processing.")
  }
  
  if (nrow(small_polys) > 0) {
    cli::cli_h2("Processing Small Polygons")
    pb_small = progress::progress_bar$new(
      format = "Small Polygons: [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
      total = nrow(small_polys),
      clear = FALSE,
      width = 60
    )
    
    small_results = purrr::map_dfr(
      1:nrow(small_polys),
      function(i) {
        pb_small$tick()
        poly_single = small_polys[i, , drop = FALSE]
        
        # Calculate stats for this polygon (no splitting)
        .calculate_vegetation(
          poly_single, 
          veg_raster, 
          id_field, 
          split = FALSE
        )
      }
    )
    
    results$small = small_results
    cli::cli_alert_success("Completed small polygon processing.")
  }
  
  # Combine results and ensure proper format
  all_results = bind_rows(results) %>%
    sf::st_drop_geometry() %>%  
    dplyr::distinct(!!rlang::sym(id_field), .keep_all = TRUE) 
  
  # Join results back to original polygons
  final_polygons = polygons %>%
    dplyr::left_join(all_results, by = id_field) %>%
    dplyr::mutate(
      PolyArea = round(as.numeric(st_area(geometry)), 2),
      VegArea = ifelse(is.na(VegArea), 0, VegArea),
      VegRatio = ifelse(is.na(VegRatio), 0, VegRatio)
    ) %>%
    dplyr::select(!!rlang::sym(id_field), PolyArea, VegArea, VegRatio)
  
  # Save results
  cli::cli_alert_info("Saving results to {output_path}")
  sf::st_write(final_polygons, output_path, quiet = TRUE, delete_layer = TRUE)
  cli::cli_alert_success("Analysis completed successfully!")
  
  # Final summary
  total_veg_area = sum(final_polygons$VegArea, na.rm = TRUE)
  cli::cli_alert_info("Total vegetation area calculated: {round(total_veg_area / 1000000, 2)} km²")
  
  if (return) {
    return(final_polygons)
  }
}




