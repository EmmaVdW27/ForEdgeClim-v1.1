#' Generate a voxelized point density grid from a LAS file from TLS
#'
#' This function reads a LAS file and creates a 3D voxel grid where each voxel represents
#' a 1x1x1m^3 volume (or a user-defined size). The density of each voxel is computed
#' based on the number of points in the point cloud. The output is a data.table containing
#' rescaled X, Y, and Z coordinates and a normalized density value.
#'
#' Additionally, a 3D plot of the voxel grid and a digital terrain model (DTM) are generated and saved as PNG files in
#' the "Output" directory, i.e. 'output_path' as defined by the user (see script 'testing_foredgeclim.R')
#'
#' @param las_file Path to the LAS file
#' @param voxel_size Numeric value defining the voxel size (default = 1m).
#' @return List containing DTM and voxel grid. Voxel grid is a data.table with columns: X, Y, Z, and density.
#' @importFrom lidR readLAS classify_ground rasterize_terrain filter_poi voxel_metrics csf knnidw normalize_height
#' @importFrom data.table as.data.table setnames CJ :=
#' @export
generate_DTM_grid_TLS <- function(las_file, voxel_size = 1) {

  start = Sys.time()

  # Read the LAS file
  # Voxel density (voxel point cloud value) was determined via octree filtering in RiSCAN:
  # Laser beam divergence of RIEGL VZ-400i = 0.02°, at 38m distance (~height canopy), this translates
  # into a distance error of 0.0256m: tan(0.02°) = (x/2)/38m => x = 0.0256m
  # Because of this, a 5cm downsampling should do to correct for distance to scan positions.
  # Remark that this does not correct for occlusion (nor for the very dense ground layer).
  las <- readLAS(las_file)

  #######
  # DTM #
  #######

  # Use a filter to retain only the ground points (create DTM)
  dtm <- classify_ground(las, csf())
  # Create a raster model of the ground (DTM)
  dtm_raster <- rasterize_terrain(dtm, res = 1, algorithm = knnidw(k=10L))
  # Convert the raster to a dataframe
  dtm_df <- as.data.frame(dtm_raster, xy = TRUE)
  colnames(dtm_df) <- c("X", "Y", "Z")

  # Normalize height of the LAS file (everywhere set ground level to Z = 0)
  las <- normalize_height(las, dtm_raster)
  # Remove the points below ground level (for points that were e.g. wrongly not classified as ground points)
  las <- filter_poi(las, Z >= 0)

  ###################
  # STRUCTURAL GRID #
  ###################

  # Compute the 3D voxel grid with point count per voxel
  voxels <- voxel_metrics(las, ~length(Z), res = voxel_size)

  # Convert to a data.table and rename columns
  voxel_grid <- as.data.table(voxels)
  setnames(voxel_grid, c("X", "Y", "Z", "n"))

  # Remove voxels on min and max X & Y boundaries
  voxel_grid <- voxel_grid[voxel_grid$X > min(voxel_grid$X) & voxel_grid$X < max(voxel_grid$X), ]
  voxel_grid <- voxel_grid[voxel_grid$Y > min(voxel_grid$Y) & voxel_grid$Y < max(voxel_grid$Y), ]

  # Rescale coordinates so that the minimum value for each axis is 1
  voxel_grid$X <- voxel_grid$X - min(voxel_grid$X, na.rm = TRUE) + 1
  voxel_grid$Y <- voxel_grid$Y - min(voxel_grid$Y, na.rm = TRUE) + 1
  voxel_grid$Z <- voxel_grid$Z - min(voxel_grid$Z, na.rm = TRUE) + 1

  # Fill missing voxels with density = 0
  full_grid <- CJ(
    X = seq(min(voxel_grid$X), max(voxel_grid$X), by = voxel_size),
    Y = seq(min(voxel_grid$Y), max(voxel_grid$Y), by = voxel_size),
    Z = seq(min(voxel_grid$Z), max(voxel_grid$Z), by = voxel_size)
  )
  full_grid <- merge(full_grid, voxel_grid, by = c("X", "Y", "Z"), all.x = TRUE)
  full_grid[is.na(n), n := 0]  # Replace NA with 0 for voxels without points

  # Normalize the density values
  # Log transformation: Reduces dominance of high-density areas (mainly ground level) (necessary even after downsampling*) -> make distribution less extreme
  # *downsampling is done with the RiSCAN octree filter previous to the import of the las file
  # Exponential scaling: Boosts contribution of higher voxels even with lower density -> sort of correction for occlusion
  # The factor 0.5 determines how strong higher layers get more weight. 0.5 seems reasonable for a temperate forest.
  # For this application (a structural proxy in a microclimate model), this transformation is likely sufficient,
  # especially since the goal is not to quantify each voxel precisely, but rather to capture a spatial pattern.
  scaling <- log1p(full_grid$n) * exp(0.5*full_grid$Z/max(full_grid$Z))
  full_grid$density <- scaling/max(scaling)

  finish = Sys.time()
  print(paste0('Time to create DTM and structural grid from TLS data = ', round(as.numeric(finish - start, units = "secs"), 2), ' s'))

  # Return the voxel grid
  return(list(dtm = dtm_df,
              grid = as.data.frame(full_grid)))
}


#' Generate a voxelized point density grid from a LAS file from ALS
#'
#' @param las_file Path to the LAS file
#' @param voxel_size Numeric value defining the voxel size (default = 1m).
#' @return List containing DTM and voxel grid. Voxel grid is a data.table with columns: X, Y, Z, and density.
#' @importFrom lidR readLAS classify_ground rasterize_terrain filter_poi voxel_metrics csf knnidw normalize_height is.empty
#' @importFrom data.table as.data.table setnames CJ :=
#' @export
generate_DTM_grid_ALS <- function(las_file, voxel_size = 1) {

  start <- Sys.time()

  # Read ALS point cloud
  las <- readLAS(las_file)

  if (is.empty(las)) {
    stop("The LAS file contains no points.")
  }

  #######
  # DTM #
  #######

  # Classify ground points
  dtm <- classify_ground(las, csf())

  # Generate digital terrain model
  dtm_raster <- rasterize_terrain(
    dtm, res = 1,
    algorithm = knnidw(k = 10L)
  )

  dtm_df <- as.data.frame(dtm_raster, xy = TRUE)
  colnames(dtm_df) <- c("X", "Y", "Z")

  # Normalize heights relative to terrain
  las <- normalize_height(las, dtm_raster)
  las <- filter_poi(las, Z >= 0)

  ###################
  # STRUCTURAL GRID #
  ###################

  # Count points in each voxel
  voxels <- voxel_metrics(
    las, ~length(Z), res = voxel_size
  )

  voxel_grid <- as.data.table(voxels)
  setnames(
    voxel_grid, c("X", "Y", "Z", "n")
  )

  # Remove outer X and Y boundaries
  voxel_grid <- voxel_grid[
    X > min(X) & X < max(X) &
      Y > min(Y) & Y < max(Y)
  ]

  # Keep physical voxel heights for weighting
  voxel_grid[, height := Z]

  # Rescale coordinates to start at 1
  voxel_grid[, X := X - min(X) + 1]
  voxel_grid[, Y := Y - min(Y) + 1]
  voxel_grid[, Z := Z - min(Z) + 1]

  # Fill empty voxels
  full_grid <- CJ(
    X = seq(min(voxel_grid$X), max(voxel_grid$X),
            by = voxel_size),
    Y = seq(min(voxel_grid$Y), max(voxel_grid$Y),
            by = voxel_size),
    Z = seq(min(voxel_grid$Z), max(voxel_grid$Z),
            by = voxel_size)
  )

  full_grid <- merge(
    full_grid,
    voxel_grid[, .(X, Y, Z, n)],
    by = c("X", "Y", "Z"),
    all.x = TRUE
  )

  full_grid[is.na(n), n := 0]

  ###################
  # ALS NORMALIZATION
  ###################

  # Relative height (0 = lowest, 1 = highest)
  z_rel <- (full_grid$Z - min(full_grid$Z)) /
    (max(full_grid$Z) - min(full_grid$Z))

  # Higher weighting towards the ground
  height_weight <- exp(0.5 * (1 - z_rel))

  # Log-transformed point counts
  # multiplied by height-dependent weight
  scaling <- log1p(full_grid$n) * height_weight

  # Normalize to [0, 1]
  full_grid[, density := scaling / max(scaling)]

  finish <- Sys.time()

  message(
    "Time to create DTM and structural grid from ALS data = ",
    round(as.numeric(finish - start, units = "secs"), 2),
    " s"
  )

  return(list(
    dtm = dtm_df,
    grid = as.data.frame(full_grid)
  ))
}



#' Generate a voxelized density grid from TLS and ALS, equally weighted
#'
#' Each source is voxelised on the SAME grid. Its voxel counts are separately
#' log-transformed and normalised to [0,1], then averaged 50/50.
#' The combined point cloud is still used to derive the DTM.
#'
#' @param tls_file Path to TLS LAS/LAZ file.
#' @param als_file Path to ALS LAS/LAZ file.
#' @param voxel_size Voxel size in metres (default 1).
#' @return List with dtm and grid (X,Y,Z,density, plus diagnostics).
#' @importFrom lidR readLAS classify_ground rasterize_terrain filter_poi voxel_metrics csf knnidw normalize_height is.empty
#' @importFrom data.table as.data.table copy CJ :=
#' @export
generate_DTM_grid_TLS_ALS <- function(tls_file, als_file, voxel_size = 1) {
  start <- Sys.time()
  stopifnot(is.numeric(voxel_size), length(voxel_size) == 1L,
            is.finite(voxel_size), voxel_size > 0)

  tls <- readLAS(tls_file)
  als <- readLAS(als_file)
  if (is.empty(tls) || is.empty(als))
    stop("One or both LAS files contain no points.")

  # Align ALS coordinates with TLS (same offsets as the original function).
  # data.table::copy prevents .internal.selfref warnings when adding columns.
  tls@data <- data.table::copy(tls@data)
  als@data <- data.table::copy(als@data)
  als@data[, X := X - 620000]
  als@data[, Y := Y - 1000000]

  n_tls_original <- nrow(tls@data)
  n_als_original <- nrow(als@data)

  # LAS files may contain different extra attributes. Retain the intersection
  # before merging, and add the SAME source flag to both point clouds.
  common_fields <- intersect(names(tls@data), names(als@data))
  if (!all(c("X", "Y", "Z") %in% common_fields))
    stop("Both LAS files must contain X, Y and Z coordinates.")
  message("TLS-only attributes omitted: ",
          paste(setdiff(names(tls@data), common_fields), collapse = ", "))
  message("ALS-only attributes omitted: ",
          paste(setdiff(names(als@data), common_fields), collapse = ", "))

  # [ keeps the LAS point table as a data.table; use data.table::copy to avoid
  # reference-related warnings. Same attribute names AND order in both files.
  tls@data <- data.table::copy(tls@data[, ..common_fields])
  als@data <- data.table::copy(als@data[, ..common_fields])
  tls@data[, source_TLS := 1L]
  als@data[, source_TLS := 0L]

  # lidR may warn about differing LAS scales/offsets. Check that the merged
  # point cloud remains accurately aligned in the desired coordinate system.
  las <- rbind(tls, als)

  # Same DTM approach as in the original function.
  ground <- classify_ground(las, csf())
  dtm_raster <- rasterize_terrain(
    ground, res = 1, algorithm = knnidw(k = 10L)
  )
  dtm_df <- as.data.frame(dtm_raster, xy = TRUE)
  colnames(dtm_df) <- c("X", "Y", "Z")

  las <- normalize_height(las, dtm_raster)
  las <- filter_poi(las, Z >= 0)
  if (is.empty(las)) stop("No points remain after height normalisation.")

  # Voxelise both sources TOGETHER, using the same voxel origin and boundaries.
  # This avoids subtle shifts that can occur with separate voxel_metrics() calls.
  voxels <- voxel_metrics(
    las,
    ~list(n_TLS = sum(source_TLS == 1L),
          n_ALS = sum(source_TLS == 0L)),
    res = voxel_size
  )
  voxel_grid <- as.data.table(voxels)
  if (!all(c("X", "Y", "Z", "n_TLS", "n_ALS") %in% names(voxel_grid)))
    stop("Unexpected voxel_metrics output; check lidR version and columns.")

  # Keep original removal of outer XY boundaries.
  voxel_grid <- voxel_grid[
    X > min(X) & X < max(X) & Y > min(Y) & Y < max(Y)
  ]
  if (!nrow(voxel_grid)) stop("No voxels remain after trimming XY boundaries.")

  # Keep original 1-based coordinate convention.
  voxel_grid[, X := X - min(X) + 1]
  voxel_grid[, Y := Y - min(Y) + 1]
  voxel_grid[, Z := Z - min(Z) + 1]

  full_grid <- CJ(
    X = seq(min(voxel_grid$X), max(voxel_grid$X), by = voxel_size),
    Y = seq(min(voxel_grid$Y), max(voxel_grid$Y), by = voxel_size),
    Z = seq(min(voxel_grid$Z), max(voxel_grid$Z), by = voxel_size)
  )
  full_grid <- merge(full_grid, voxel_grid,
                     by = c("X", "Y", "Z"), all.x = TRUE)
  full_grid[is.na(n_TLS), n_TLS := 0]
  full_grid[is.na(n_ALS), n_ALS := 0]

  # EQUAL weighting: normalise each source independently before averaging.
  # Counts from one source cannot swamp the other source's normalisation.
  scale_tls <- log1p(full_grid$n_TLS)
  scale_als <- log1p(full_grid$n_ALS)
  max_tls <- max(scale_tls)
  max_als <- max(scale_als)
  if (max_tls <= 0 || max_als <= 0)
    stop("One source has no points in the retained voxel grid.")

  full_grid[, density_TLS := scale_tls / max_tls]
  full_grid[, density_ALS := scale_als / max_als]
  full_grid[, density := 0.5 * density_TLS + 0.5 * density_ALS]
  full_grid[, n := n_TLS + n_ALS]

  message(sprintf("Original point counts: TLS = %s; ALS = %s",
                  format(n_tls_original, big.mark = ","),
                  format(n_als_original, big.mark = ",")))
  message(sprintf("Original point shares: TLS = %.1f%%; ALS = %.1f%%",
                  100 * n_tls_original / (n_tls_original + n_als_original),
                  100 * n_als_original / (n_tls_original + n_als_original)))
  message("Density weighting after separate normalisation: TLS = 50%; ALS = 50%")
  message(sprintf("Time to create equal-weight TLS-ALS grid = %.2f s",
                  as.numeric(difftime(Sys.time(), start, units = "secs"))))

  return(list(dtm = dtm_df, grid = as.data.frame(full_grid)))
}
