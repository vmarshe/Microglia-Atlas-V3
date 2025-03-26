#-------------------------------------------------------------------------------
#
# ReadVizgen
# EDITS FROM https://github.com/satijalab/seurat/pull/7190
#
#-------------------------------------------------------------------------------
ReadVizgen <- function(
    data.dir,
    transcripts = NULL,
    spatial = NULL,
    molecules = NULL,
    type = c('centroids', 'segmentations'),
    mol.type = 'microns',
    metadata = NULL,
    filter = NA_character_,	
    z = 3L,
    use.BiocParallel = TRUE,
    workers.MulticoreParam = 12,
    DTthreads.pct = NULL,
    min.area = 5,
    verbose = FALSE
) {
  # TODO: handle multiple segmentations per z-plane
  # NOTE: this is only needed when segmentations differ between z-planes
  
  # packages that need to be installed a priori
  pkgs <- c("data.table", "arrow", "sfarrow",
            "tidyverse", "sf", "BiocParallel", "Matrix")
  lapply(pkgs %>% length %>% seq, function(i) 
  { !requireNamespace(pkgs[i], quietly = TRUE) } ) %>% 
    unlist %>% 
    { if (c(which(.) > 0) %>% any()) 
    { c("Please install ->", "\n",
        paste0("'", pkgs[which(.)], "'", collapse = ", "), " for this function") %>% 
        stop(., call. = FALSE) } }
  
  # setting workers to use for parallelization - `BiocParallel` ----
  if (use.BiocParallel) {
    if (verbose) { message("Using parallelization with: `BiocParallel`") }
    if (is.null(workers.MulticoreParam)) { 
      workers.MulticoreParam <- quantile(BiocParallel::multicoreWorkers() %>% seq)[4] %>% ceiling
      if (verbose) { message(workers.MulticoreParam, " of total workers available will be used") }
    } else { if (verbose) { message("Setting total workers to: ", workers.MulticoreParam) } }   
  } else { if (verbose) { message("Using parallelization with: `future`", "\n",
                                  "NOTE: set workers for parallelization, eg: `future::plan('multisession', workers = 10)`") }
  }
  # support parallelization on unix and windows
  # credit to https://github.com/Bioconductor/BiocParallel/issues/98
  BPParam <-
    if (.Platform$OS.type == "windows") {
      if (verbose) { message("Using parallelization for Windows with: `BiocParallel::SerialParam`") }
      BiocParallel::SerialParam(force.GC = FALSE, progressbar = TRUE)
    } else {
      BiocParallel::MulticoreParam(workers.MulticoreParam, tasks = 50L, 
                                   force.GC = FALSE, progressbar = TRUE)
    }
  
  # (optional) 
  # setting additional cores to use for parallelization in `data.table`
  # might allow to read large files faster
  if (!is.null(DTthreads.pct)) {
    if (verbose) { message("Using parallelization with: `data.table`", "\n",
                           "..for `data.table::fread`") }
    data.table::setDTthreads(threads = 0) # all cores
    DTthreads <- data.table::getDTthreads() # max cores
    DTthreads <- c((DTthreads * DTthreads.pct) / 100) %>% ceiling # percentage from total threads
    if (verbose) { message("Setting DTthreads to: ", DTthreads, " (", paste0(DTthreads.pct, "%"), ")") }
    data.table::setDTthreads(threads = DTthreads) # set
  }
  
  # hdf5r is only used for loading polygon boundaries
  # Not needed for all Vizgen input
  hdf5 <- requireNamespace("hdf5r", quietly = TRUE)
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('segmentations', 'centroids', 'boxes'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels', 'microns'),
    several.ok = TRUE
  )
  if (!is.null(x = metadata)) {
    metadata <- match.arg(
      arg = metadata,
      choices = c('volume', 'fov'),
      several.ok = TRUE
    )
  }
  if (!z %in% seq.int(from = 0L, to = 6L)) {
    stop("The z-index must be in the range [0, 6]")
  }
  use.dir <- all(vapply(
    X = c(transcripts, spatial, molecules),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Vizgen directory ", data.dir)
  } else { 
    if (verbose) { message("Reading data from:" , "\n", data.dir) }
  }
  
  # Use segmentation output from ".parquet" file ----
  # check if file exists
  parq <-
    # look in the current directory
    list.files(data.dir, 
               pattern = ".parquet$",
               full.names = TRUE) %>%
    { if (length(.) == 0) { 
      # look in the sub directory (if nothing is found)
      list.files(data.dir,
                 pattern = ".parquet$",
                 full.names = TRUE,
                 recursive = TRUE)
    } else { (.) }}
  # set to use .parquet" file if present
  use.parquet <- 
    parq %>% length() %>% any
  
  if (use.parquet) {
    if (verbose) { message("Cell segmentations are found in `.parquet` file", "\n", 
                           ">>> using ", gsub("s", "", mol.type), " space coordinates") }}
  
  # Identify input files..
  # if no files are found in the current directory..
  #..look for them in the sub directory
  files <- c(transcripts = NULL %||% 
               list.files(data.dir, 
                          pattern = "cell_by_gene",
                          full.names = TRUE) %>%
               { if (length(.) == 0) { 
                 list.files(data.dir,
                            pattern = "cell_by_gene",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }} %>% 
               { if (length(.) > 1) { 
                 stop("There are > 1 `cell_by_gene` files",
                      "\n", "make sure only 1 file is read")
               } else if (!length(.)) {
                 stop("No `cell_by_gene` file is available")
               } else { (.) }},
             
             spatial = NULL %||% 
               list.files(data.dir, 
                          pattern = "cell_metadata",
                          full.names = TRUE) %>%
               { if (length(.) == 0) { 
                 list.files(data.dir,
                            pattern = "cell_metadata",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }} %>% 
               { if (length(.) > 1) { 
                 stop("There are > 1 `cell_metadata` files",
                      "\n", "make sure only 1 file is read")
               } else if (!length(.)) {
                 stop("No `cell_metadata` file is available")
               } else { (.) }},
             
             molecules = NULL %||% 
               list.files(data.dir, 
                          pattern = "detected_transcripts",
                          full.names = TRUE) %>%
               { if (length(.) == 0) { 
                 list.files(data.dir,
                            pattern = "detected_transcripts",
                            full.names = TRUE,
                            recursive = TRUE)
               } else { (.) }} %>%
               { if(length(.) > 1) { 
                 stop("There are > 1 `detected_transcripts` files", 
                      "\n", "make sure only 1 file is read")
               } else if (!length(.)) {
                 stop("No `detected_transcripts` file is available")
               } else { (.) }}
  )
  
  files[is.na(x = files)] <- NA_character_
  h5dir <- file.path(
    ifelse(
      test = dirname(path = files['spatial']) == '.',
      yes = data.dir,
      no = dirname(path = files['spatial'])
    ),
    'cell_boundaries'
  )
  zidx <- paste0('zIndex_', z)
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  if (all(is.na(x = files))) {
    stop("Cannot find Vizgen input files in ", data.dir)
  }
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['spatial']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    sp <- data.table::fread(
      file = files[['spatial']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE
      # showProgress = progressr:::progressr_in_globalenv(action = 'query')
      # showProgress = verbose
    )
    pprecoord(type = 'finish')
    rownames(x = sp) <- as.character(x = sp[, 1])
    #sp <- sp[, -1, drop = FALSE]
    if ((names(sp) == "transcript_count") %>% any) {
      if (verbose) { message(">>> filtering `cell_metadata` - keep cells with `transcript_count` > 0") }
      sp %<>% select(-1) %>% filter(transcript_count > 0)
    } else { sp %<>% select(-1) }
    
    # Check to see if we should load segmentations
    if ('segmentations' %in% type) {
      poly <- if (isFALSE(x = hdf5) && !use.parquet) {
        warning(
          "Cannot find hdf5r; unable to load segmentation vertices",
          immediate. = TRUE
        )
        FALSE
      } else if (!dir.exists(paths = h5dir) && !use.parquet) {
        warning("Cannot find cell boundary H5 or `.parquet` file(s)", immediate. = TRUE)
        FALSE
      } else if (use.parquet) { # for non .hdf5 files
        if (length(parq)) {
        }
        TRUE
      }
      else {
        TRUE
      }
      if (isFALSE(x = poly) && isFALSE(use.parquet)) {
        type <- setdiff(x = type, y = 'segmentations')
      }
    }
    spatials <- rep_len(x = files[['spatial']], 
                        length.out = length(x = type))
    names(x = spatials) <- type
    files <- c(files, spatials)
    files <- files[setdiff(x = names(x = files), y = 'spatial')]
  } else if (!is.null(x = metadata)) {
    warning(
      "metadata can only be loaded when spatial coordinates are loaded",
      immediate. = TRUE
    )
    metadata <- NULL
  }
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    
    # optionally - load molecules
    mx <- data.table::fread(
      file = files[['molecules']],
      sep = ',',
      verbose = FALSE
      # showProgress = verbose
    )
    mx <- mx[mx$global_z == z, , drop = FALSE]
    if (!is.na(x = filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = filter, x = mx$gene), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules')]
  }
  files <- files[!is.na(x = files)]
  
  # Read input data ----
  outs <- vector(mode = 'list', length = length(x = files))
  names(x = outs) <- names(x = files)
  if (!is.null(metadata)) {
    outs <- c(outs, list(metadata = NULL))
  }
  for (otype in names(x = outs)) {
    outs[[otype]] <- 
      switch(EXPR = otype, 
             transcripts = {
               ptx <- progressor()
               ptx(message = "Reading counts matrix", class = "sticky",
                   amount = 0)
               tx <- data.table::fread(file = files[[otype]], sep = ",",
                                       data.table = FALSE, verbose = FALSE)
               rownames(x = tx) <- as.character(x = tx[, 1])
               # keep cells with `transcript_count` > 0
               if (verbose) { message(">>> filtering `cell_by_gene` - keep cells with counts > 0") }
               tx %<>% select(-1) %>%
                 filter_all(any_vars(. > 0)) %>% {
                   if ((names(sp) == "transcript_count") %>% any) {
                     # match cells to filtered data from spatial (cell_metadata)
                     filter(., rownames(.) %in% rownames(sp)) 
                   } else { (.) } } %>% # return filtered count data
                 as.sparse() %>% # convert to sparse matrix
                 Matrix::t() # transpose sparse matrix
               
               # filter cell metadata df
               # match filtered cell IDs from count matrix to cell metadata df
               if (!(names(sp) == "transcript_count") %>% any) {
                 sp %<>% filter(rownames(.) %in% colnames(tx))
               }
               
               ratio <- getOption(x = "Seurat.input.sparse_ratio",
                                  default = 0.4)
               if ((sum(tx == 0)/length(x = tx)) > ratio) {
                 ptx(message = "Counts are converted to sparse matrix",
                     class = "sticky", amount = 0)
                 #tx <- as.sparse(x = tx)
               }
               if (!is.na(x = filter)) {
                 ptx(message = paste("Filtering genes with pattern",
                                     filter), class = "sticky", amount = 0)
                 tx <- tx[!grepl(pattern = filter, x = rownames(x = tx)),
                          , drop = FALSE]
               }
               ptx(type = "finish")
               tx
               
             }, 
             centroids = {
               pcents <- progressor()
               pcents(message = "Creating centroid coordinates",
                      class = "sticky", amount = 0)
               pcents(type = "finish")
               data.frame(x = sp$center_x, y = sp$center_y, cell = rownames(x = sp),
                          stringsAsFactors = FALSE)
             }, 
             segmentations = {
               # use segmentations from ".parquet"
               if (use.parquet) {
                 if (length(parq) > 1) {
                   # eg, if two files are present:
                   # `cellpose_micron_space.parquet`
                   # `cellpose_mosaic_space.parquet`
                   parq %<>% { 
                     if (mol.type == "pixels") {
                       # for `cellpose_mosaic_space.parquet`
                       grep("mosaic", ., value = TRUE)
                     } else if (mol.type == "microns") {
                       # for `cellpose_micron_space.parquet`
                       grep(gsub("s", "", mol.type), ., value = TRUE)
                     }
                   } %>%
                     { if (length(.) == 0) {
                       # only if single ".parquet" file present
                       # eg, `cell_boundaries.parquet`
                       parq %>%
                         grep(".parquet$", ., value = TRUE)
                     } else { (.) } }
                 }
                 
                 # Read .parquet file ----
                 parq %<>% sfarrow::st_read_parquet(.) %>%
                   # keep only selected z-plane
                   filter(., ZIndex == z)
                 
                 # Sanity checks on segmentation polygons
                 # Using `sf` for filtering polygons ----
                 # keep multiple polygons pieces per single cell id
                 parq %<>%
                   .filter_polygons(., 
                                    min.area = min.area, 
                                    BPPARAM = 
                                      if (BPParam$workers > 2) {
                                        # use 2 workers less than total workers
                                        BiocParallel::MulticoreParam(14 - 2, tasks = 50L, 
                                                                     force.GC = FALSE, progressbar = TRUE)
                                      } else if (!BPParam$workers > 2)
                                      { BPParam } else { SerialParam() }, 
                                    verbose = verbose)
                 
                 # Get filtered segmentation geometries/polygons
                 segs <- 
                   parq %>%
                   st_geometry() %>%
                   # add cell ID
                   set_names(pull(parq, EntityID) %>% as.character)
                 
                 # Extract cell boundaries per cell id ----
                 # TODO: (optionally) resample & make cell boundaries equidistant? 
                 if (use.BiocParallel) {
                   gc() %>% invisible() # free up memory
                   if (verbose) { message("Extracting cell segmentations - using `BiocParallel`") } 
                   segs_list <-
                     BiocParallel::bplapply(segs %>% seq,
                                            function(i) {
                                              # keep multiple polygons per cell
                                              segs[[i]] %>% 
                                                unique() %>% # remove any duplicates
                                                purrr::lmap(.f = data.table::as.data.table) %>% 
                                                data.table::rbindlist() %>%
                                                mutate(cell = names(segs)[i]) },
                                            BPPARAM = BPParam)
                 } else {
                   if (verbose) { message("Extracting cell segmentations - using `future`") }    
                   segs_list <-
                     future.apply::future_lapply(segs %>% seq,
                                                 function(i) {
                                                   segs[[i]] %>% 
                                                     unique() %>% # remove any duplicates
                                                     purrr::lmap(.f = data.table::as.data.table) %>% 
                                                     data.table::rbindlist() %>%
                                                     mutate(cell = names(segs)[i])
                                                 }
                     )
                 }
                 
                 # df of all cell segmentations
                 segs <- 
                   data.table::rbindlist(segs_list) %>% 
                   data.table::setnames(c("x", "y", "cell"))
                 if (verbose) { message("All cell segmentations are loaded..") }
                 npg <- length(x = unique(x = segs$cell))
                 if (npg < nrow(x = sp)) {
                   warning(nrow(x = sp) - npg, " cells missing polygon information",
                           immediate. = TRUE) }
                 segs  # final segmentaions out..
               } else {
                 # else use ".hdf5" files from ./cell_boundaries (older version)
                 ppoly <- progressor(steps = length(x = unique(x = sp$fov)))
                 ppoly(message = "Creating polygon coordinates", class = "sticky",
                       amount = 0)
                 # use `BiocParallel` or `future`
                 if (use.BiocParallel) {
                   if (verbose) { message("Reading '.hdf5' files..") }
                   pg <- BiocParallel::bplapply(X = unique(x = sp$fov), FUN = function(f, ...) {
                     fname <- file.path(h5dir, paste0("feature_data_", f, ".hdf5"))
                     if (!file.exists(fname)) {
                       warning("Cannot find HDF5 file for field of view ",
                               f, immediate. = TRUE)
                       return(NULL)
                     }
                     # reading hdf5 files
                     hfile <- hdf5r::H5File$new(filename = fname, 
                                                mode = "r")
                     on.exit(expr = hfile$close_all())
                     cells <- rownames(x = subset(x = sp, subset = fov == f))     
                     # creating df for cell boundaries     
                     df <- lapply(X = cells, FUN = function(x) {
                       return(tryCatch(expr = {
                         cc <- hfile[["featuredata"]][[x]][[zidx]][["p_0"]][["coordinates"]]$read()
                         cc <- as.data.frame(x = t(x = cc))
                         colnames(x = cc) <- c("x", "y")
                         cc$cell <- x
                         cc
                       }, error = function(...) {
                         return(NULL)
                       }))
                     })   
                     ppoly()
                     #return(do.call(what = "rbind", args = df))
                     return(data.table::rbindlist(df))
                   }, BPPARAM = BPParam)
                 } else { 
                   pg <- 
                     future.apply::future_lapply(X = unique(x = sp$fov), 
                                                 FUN = function(f, ...) {
                                                   fname <- file.path(h5dir, paste0("feature_data_",
                                                                                    f, ".hdf5"))
                                                   if (!file.exists(fname)) {
                                                     warning("Cannot find HDF5 file for field of view ",
                                                             f, immediate. = TRUE)
                                                     return(NULL)
                                                   }
                                                   hfile <- hdf5r::H5File$new(filename = fname,
                                                                              mode = "r")
                                                   on.exit(expr = hfile$close_all())
                                                   cells <- rownames(x = subset(x = sp, subset = fov == f))
                                                   df <- lapply(X = cells, FUN = function(x) {
                                                     return(tryCatch(expr = {
                                                       cc <- hfile[["featuredata"]][[x]][[zidx]][["p_0"]][["coordinates"]]$read()
                                                       cc <- as.data.frame(x = t(x = cc))
                                                       colnames(x = cc) <- c("x", "y")
                                                       cc$cell <- x
                                                       cc
                                                     }, error = function(...) {
                                                       return(NULL)
                                                     }))
                                                   })
                                                   ppoly()
                                                   return(do.call(what = "rbind", args = df))
                                                 }
                     )
                 }
                 ppoly(type = "finish")
                 # cell polygons
                 pg <- do.call(what = "rbind", args = pg)
                 npg <- length(x = unique(x = pg$cell))
                 if (npg < nrow(x = sp)) {
                   warning(nrow(x = sp) - npg, " cells missing polygon information",
                           immediate. = TRUE)
                 }
                 pg # final segmentaions out..
               }
             }, 
             boxes = {
               pbox <- progressor(steps = nrow(x = sp))
               pbox(message = "Creating box coordinates", class = "sticky",
                    amount = 0)
               # use parallel or future
               if (use.BiocParallel) {
                 if (verbose) { message("Creating box coordinates..") }
                 bx <- BiocParallel::bplapply(X = rownames(x = sp), FUN = function(cell) {
                   row <- sp[cell, ]
                   # faster version for grid construction
                   df <- data.table::CJ(x = c(row$min_x, row$max_x),
                                        y = c(row$min_y, row$max_y), 
                                        cell = cell) %>% 
                     slice(c(1, 3, 4, 2))
                   #df <- expand.grid(x = c(row$min_x, row$max_x),
                   #                 y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                   #                stringsAsFactors = FALSE)
                   #df <- df[c(1, 3, 4, 2), , drop = FALSE]
                   pbox()
                   return(df)
                 }, BPPARAM = BPParam)
               } else {
                 bx <- future.apply::future_lapply(X = rownames(x = sp), FUN = function(cell) {
                   row <- sp[cell, ]
                   df <- data.table::CJ(x = c(row$min_x, row$max_x),
                                        y = c(row$min_y, row$max_y), cell = cell) %>%
                     slice(c(1, 3, 4, 2))
                   #df <- expand.grid(x = c(row$min_x, row$max_x),
                   #y = c(row$min_y, row$max_y), cell = cell, KEEP.OUT.ATTRS = FALSE,
                   #stringsAsFactors = FALSE)
                   #df <- df[c(1, 3, 4, 2), , drop = FALSE]
                   pbox()
                   return(df)
                 })
               }
               pbox(type = "finish")
               do.call(what = "rbind", args = bx)
             }, 
             metadata = {
               pmeta <- progressor()
               pmeta(message = "Loading metadata", class = "sticky",
                     amount = 0)
               pmeta(type = "finish")
               sp[, metadata, drop = FALSE]
             }, 
             pixels = {
               ppixels <- progressor()
               ppixels(message = "Creating pixel-level molecule coordinates",
                       class = "sticky", amount = 0)
               df <- data.frame(x = mx$x, y = mx$y, gene = mx$gene,
                                stringsAsFactors = FALSE)
               ppixels(type = "finish")
               df
             }, 
             microns = {
               pmicrons <- progressor()
               pmicrons(message = "Creating micron-level molecule coordinates",
                        class = "sticky", amount = 0)
               df <- data.frame(x = mx$global_x, y = mx$global_y,
                                gene = mx$gene, stringsAsFactors = FALSE)
               pmicrons(type = "finish")
               df
             }, stop("Unknown MERFISH input type: ", type))
  }
  
  # add z-slice index for cells ----
  outs$zIndex <- 
    data.frame(z = rep_len(z, length.out = outs$centroids %>% pull(cell) %>% length), 
               cell = outs$centroids %>% pull(cell))
  
  return(outs)
}
#-------------------------------------------------------------------------------
#
# This is the .filter_polygons function from Seurat.
#
#-------------------------------------------------------------------------------
#' @import dplyr
#' @import sf
#' @import magrittr
#' Helper function to filter geometryy polygons using \code{sf} package
#' modified function from \code{SpatialFeatureExperiment:::.filter_polygons}
.filter_polygons <-
  function(polys, # object of class \code{sf} and \code{data.frame}
           min.area, # minimal polygon area to use as a threshold
           BPPARAM = BiocParallel::SerialParam(),
           verbose) {
    # Sanity check on nested polygon lists
    test.segs <- lapply(polys %>% st_geometry() %>% seq, function(i) {
      polys %>%
        st_geometry() %>% .[[i]] %>% length() }) %>% unlist()
    if (any(test.segs > 1)) {
      segs.art.index <- which(test.segs > 1)
      if (verbose) { log_info("Sanity checks on cell segmentation polygons:", "\n",
                              ">>> ..found ", segs.art.index %>% length,
                              " cells with (nested) polygon lists", "\n",
                              ">>> ..applying filtering") }
    }
    # remove empty elements
    polys.orig <- polys
    polys %<>% filter(!st_is_empty(.))
    empty.inds <- which(!polys.orig$ID %in% polys$ID)
    if (verbose && length(empty.inds)) { log_info(">>> ..removing ",
                                                  length(empty.inds), " empty polygons") }
    
    if (st_geometry_type(polys, by_geometry = FALSE) == "MULTIPOLYGON") {
      polys_sep <- lapply(st_geometry(polys), function(x) {
        st_cast(st_sfc(x), "POLYGON")
      })
      areas <- lapply(polys_sep, st_area)
      
      if (!is.null(min.area)) {
        which_keep <- lapply(areas, function(x) which(x > min.area))
        multi_inds <- which(lengths(which_keep) > 1L)
        if (length(multi_inds)) {
          if (verbose) { log_info("There are ", length(multi_inds), " cells with multiple",
                                  " pieces in cell segmentation larger than min.area,",
                                  " whose first 10 indices are: ",
                                  paste(multi_inds %>% head(10), # print only first 10 indices
                                        collapse = ", "),
                                  ". The largest piece is kept.") }
          which_keep[multi_inds] <- lapply(areas[multi_inds], which.max)
        }
        inds <- lengths(which_keep) > 0L
        polys <- polys[inds,]
        # using parallelization, else can take a while when `which_keep` length is towards 100K
        which_keep <- unlist(which_keep[inds])
        geo <- st_geometry(polys)
        new_geo <-
          BiocParallel::bplapply(seq_along(which_keep), function(i) {
            geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]] %>%
              unique() %>% # remove any duplicates
              st_polygon()
          }, BPPARAM = BPPARAM) |> st_sfc()
        st_geometry(polys) <- new_geo
      } else if (is.null(min.area)) {
        # use only maximal area, ie the largest polygon
        if (verbose) { log_info(">>> ..keeping polygons with the largest area only") }
        which_keep <-
          lapply(areas, function(x) which.max(x)) %>% unlist()
        geo <- st_geometry(polys)
        new_geo <-
          BiocParallel::bplapply(seq_along(which_keep), function(i) {
            geo[[i]] <- st_cast(geo[i], "POLYGON")[[which_keep[i]]] %>%
              unique() %>% # remove any duplicates
              st_polygon()
          }, BPPARAM = BPPARAM) |> st_sfc()
        st_geometry(polys) <- new_geo
      }
    } else {
      inds <- st_area(st_geometry(polys)) > min.area
      polys <- polys[inds,]
    }
    polys
  }
#-------------------------------------------------------------------------------
#
# LoadVizgen
# EDITS FROM https://github.com/satijalab/seurat/pull/7190
#
#-------------------------------------------------------------------------------
LoadVizgen <- function(
    data.dir, 
    fov = 'vz', 
    assay = 'Vizgen',
    mol.type = 'microns',
    filter = '^Blank-',
    z = 3L,
    add.zIndex = TRUE, 
    update.object = TRUE,
    add.molecules = TRUE,
    min.area = 5,
    verbose,
    ...)
{
  # reading data..
  data <- ReadVizgen(data.dir = data.dir,
                     mol.type = mol.type,
                     filter = filter,
                     z = z,
                     min.area = min.area,
                     verbose = verbose,
                     ...)
  
  if (verbose) { message("Creating Seurat object..") }  
  obj <- CreateSeuratObject(counts = data[["transcripts"]], assay = assay)
  
  # in case no segmentation is present, use boxes
  if (!"segmentations" %in% names(data)) {
    if ("boxes" %in% names(data)) {
      bound.boxes <- CreateSegmentation(data[["boxes"]])
      cents <- CreateCentroids(data[["centroids"]])
      bound.boxes.data <- list(centroids = cents, 
                               boxes = bound.boxes)
      if (verbose) { 
        message("Creating FOVs..", "\n", 
                if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
                "\n",
                ">>> using box coordinates instead of segmentations") 
      }
      coords <- 
        CreateFOV(coords = bound.boxes.data,
                  type = c("boxes", "centroids"),
                  molecules = 
                    if (add.molecules) {
                      data[[mol.type]] } else { NULL }, 
                  assay = assay) %>%
        subset(x = .,
               cells = intersect(x = Cells(x = .[["boxes"]]),
                                 y = Cells(x = obj)))
    } else { 
      # in case no segmentation & no boxes are present, use centroids only
      cents <- CreateCentroids(data[["centroids"]])
      if (verbose) { 
        message("Creating FOVs..", "\n", 
                if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
                "\n", 
                ">>> using only centroids") 
      }
      coords <-
        CreateFOV(coords = list(centroids = cents),
                  type = c("centroids"),
                  molecules = 
                    if (add.molecules) {
                      data[[mol.type]] } else { NULL }, 
                  assay = assay) %>%
        subset(x = ., 
               cells = intersect(x = Cells(x = .[["centroids"]]),
                                 y = Cells(x = obj)))
    }
  } else if ("segmentations" %in% names(data)) {
    segs <- CreateSegmentation(data[["segmentations"]])
    cents <- CreateCentroids(data[["centroids"]])
    segmentations.data <- list(centroids = cents, segmentation = segs)
    if (verbose) { 
      message("Creating FOVs..", "\n", 
              if (!add.molecules) { ">>> `molecules` coordidates will be skipped" }, 
              "\n", 
              ">>> using segmentations") 
    }
    coords <-
      CreateFOV(coords = segmentations.data, 
                type = c("segmentation", "centroids"), 
                molecules = 
                  if (add.molecules) {
                    data[[mol.type]] } else { NULL }, 
                assay = assay) %>%
      # only consider the cells we have counts and a segmentation.
      # Cells which don't have a segmentation are probably found in other z slices.
      subset(x = .,
             cells = intersect(x = Cells(x = .[["segmentation"]]),
                               y = Cells(x = obj)))
  }
  
  # add z-stack index for cells
  if (add.zIndex) { obj$z <- data$zIndex %>% pull(z) }
  
  # add metadata vars
  if (verbose) { message(">>> adding metadata infos") }
  if (c("metadata" %in% names(data))) {
    metadata <- match.arg(arg = "metadata", choices = names(data), several.ok = TRUE)
    meta.vars <- names(data[[metadata]])
    for (i in meta.vars %>% seq) {
      obj %<>% AddMetaData(metadata = data[[metadata]][[meta.vars[i]]], 
                           col.name = meta.vars[i])
    }
  }
  
  # sanity on fov name
  fov %<>% gsub("_|-", ".", .)
  
  if (verbose) { message(">>> adding FOV") }
  obj[[fov]] <- coords
  
  ## filter - keep cells with counts > 0
  # helper function to return metadata
  callmeta <- function (object = NULL) { return(object@meta.data) }
  nCount <- grep("nCount", callmeta(obj) %>% names, value = TRUE)
  if (any(obj[[nCount]] == 0)) {
    if (verbose) { message(">>> filtering object - keeping cells with counts > 0") }
    obj %<>% subset(subset = !!base::as.symbol(nCount) > 0)
  } else { if (verbose) { message(">>> all counts are > 0") } }
  
  if (update.object) { 
    if (verbose) { message("Updating object:") 
      obj %<>% UpdateSeuratObject()
    } else { 
      obj %<>% 
        UpdateSeuratObject() %>% 
        suppressMessages() } }
  
  if (verbose) { message("Object is ready!") } 
  return(obj)
}

#-------------------------------------------------------------------------------
#
# load_merscope
#
#-------------------------------------------------------------------------------
load_merscope <- function(folder_path = NULL,
                         fov_name = "fov1",
                         filter = "^Blank-"){
  
  require(progressr)
  require(tidyverse)
  require(magrittr)
  require(arrow)
  require(sfarrow)
  require(sf)
  require(BiocParallel)
  require(data.table)
  require(logger)
  require(SeuratObject)
  require(Seurat)
  require(Matrix)
  
  obj <- LoadVizgen(data.dir = folder_path,
                    filter = filter,
                    verbose = T,
                    fov = fov_name,
                    min.area = 1)
  
  obj@meta.data <- obj@meta.data %>%
    mutate(cell_names = as.character(rownames(obj@meta.data)))
  
  keep_cells <- obj@meta.data$cell_names
  
  # ADD ADDITIONAL METADATA
  #-----------------------------------------------------------------------------
  if(file.exists(file.path(folder_path, "cell_metadata.csv"))){
    metadata <- fread(file.path(folder_path, "cell_metadata.csv")) %>%
      rename_at(vars(any_of(c("...1", "EntityID", "V1"))), ~"cell_names") %>%
      mutate(cell_names = as.character(cell_names)) %>%
      filter(cell_names %in% keep_cells) %>%
      as_tibble()
    
    for(vars in c("fov", "solidity", "perimeter_area_ratio", "volume")){
      tmp_var <- metadata %>% pull(vars)
      names(tmp_var) <- metadata$cell_names
      obj <- AddMetaData(obj, tmp_var, vars)
    }
  }
  
  # ADD INFO ABOUT BLANKS
  #-----------------------------------------------------------------------------
  cell_by_gene <- data.table::fread(file.path(folder_path, "cell_by_gene.csv"))
  names(cell_by_gene)[1] <- "cell"
  cell_by_gene <- cell_by_gene %>% as.data.frame() %>%
    mutate(cell = as.character(cell)) %>%
    filter(cell %in% keep_cells) %>%
    dplyr::select(cell_names = cell, all_of(starts_with("Blank-")))
  
  blanks <- rowSums(cell_by_gene[,-1])
  names(blanks) <- cell_by_gene$cell_names
  
  obj <- AddMetaData(obj, blanks, col.name = "nCount_Vizgen_blanks")
  
  # SUBSET
  #-----------------------------------------------------------------------------
  invisible(lapply(capture.output(obj), log_info))
  
  cell_w_centroid <- obj@images[[fov_name]]$centroids@cells
  cell_w_segmentation <- names(obj@images[[fov_name]]$segmentation@polygons)
  keep_cells <- intersect(cell_w_centroid, cell_w_segmentation)
  
  log_info(length(keep_cells), " cells have centroid & segmentation out of ", ncol(obj))
  
  obj <- subset(obj, subset = cell_names %in% keep_cells)
  log_info("=== Final object: ")
  invisible(lapply(capture.output(obj), log_info))
  
  return(obj)
  
}


#===============================================================================
#
# GET CELL GEOMETRY (FORM FACTOR)
#
#===============================================================================
get_cell_geometry <- function(data = NULL, fov_names = NULL){
  
  require(sf)
  require(doFuture, include.only = c("registerDoFuture", "%dofuture%"))
  require(foreach)
  require(tidyverse)
  require(future)
  
  registerDoFuture()
  plan(multisession)
  options(future.globals.maxSize = 8000 * 1024^2)
  
  res <- data.frame()
  for(fov in fov_names){
    
    cell_names <- data@images[[fov]]@boundaries$centroids@cells
    polygons <- lapply(cell_names, function(x) st_polygon(list(data@images[[fov]]@boundaries$segmentation@polygons[[x]]@Polygons[[1]]@coords)))
    
    shape_features <- foreach(i = 1:length(cell_names), .combine = rbind, .options.future = list(packages = c("sf", "lwgeom"))) %dofuture% {
      
      cell_sfc <- st_sfc(polygons[i])
      
      # For definition of form factor, see CellProfiler
      # https://cellprofiler-manual.s3.amazonaws.com/CellProfiler-3.0.0/modules/measurement.html
      # https://github.com/CellProfiler/CellProfiler/blob/67a758ef2c74b07e3b8e97d9bd08ebce9c0f9087/cellprofiler/modules/measureobjectsizeshape.py
      
      form_factor <- 4.0 * pi * st_area(cell_sfc) / (lwgeom::st_perimeter_2d(cell_sfc)^2)
      form_factor <- ifelse(length(form_factor) == 0, NA, form_factor)
      
      cbind.data.frame(cell_names = cell_names[i], fov, form_factor)
    }
    res <-rbind(res, shape_features)
  }
  
  plan(sequential)
  
  return(res)
  
}

#-------------------------------------------------------------------------------
#
# qc_spatial
#
#-------------------------------------------------------------------------------
qc_spatial <- function(data = NULL,
                       min_umi = 0,
                       max_umi = 2500,
                       min_genes = 0,
                       min_volume = 0,
                       max_volume = 2500,
                       fov_name = "fov1",
                       res_path = NULL,
                       dpi = 300,
                       assay = "Vizgen",
                       size_feature = "volume"){
  
  require(Seurat)
  require(logger)
  require(patchwork)
  require(Matrix)
  
  lapply(capture.output(data), log_trace)
  
  log_info("UMI threshold: [", min_umi, ", ", max_umi, "]")
  log_info("Gene count threshold: ", min_genes)
  log_info("Volume/area threshold: [", min_volume, ", ", max_volume, "]")
  
  # UMI
  #-----------------------------------------------------------------------------
  log_info("# UMI distribution ")
  umi_sum <- summary(data@meta.data[,paste0("nCount_", assay)])
  log_info(paste0(names(umi_sum), ": ", sapply(umi_sum, function(x) round(x, 1)), collapse = ", "))
  
  pdf(file.path(res_path, "qc_hist_UMI.pdf"), height = 3, width = 4)
  p <- data@meta.data %>%
    rename_at(all_of(paste0("nCount_", assay)), ~"feature") %>%
    ggplot(aes(x = feature)) +
    geom_histogram(fill = proj_cols("light blue"), color = "black", alpha = .8, bins = 50) +
    geom_vline(xintercept = min_umi, col = proj_cols("pink"))+
    geom_vline(xintercept = max_umi, col = proj_cols("pink"))+
    geom_vline(xintercept = umi_sum["Mean"], col = proj_cols("yellow"))+
    labs(x = "UMI Count", y = expression(paste(italic('N'),' Cells')))+
    theme_proj()+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "in"))
  print(p)
  dev.off()
  
  pdf(file.path(res_path, "qc_fov_UMI.pdf"), height = 5.5, width = 5)
  p <- qc_fov(data, 
              min_feature = min_umi, 
              max_feature = max_umi, 
              assay = assay, 
              fov = fov_name, 
              legend_title = "UMI Threshold", 
              feature = paste0("nCount_", assay), 
              feature_label = " ")
  print(p)
  dev.off()
  
  keep_umi <- colnames(data)[which(colSums(data[[assay]]) >= min_umi &  colSums(data[[assay]]) <= max_umi)]
  
  # GENE COUNT
  #-----------------------------------------------------------------------------
  log_info("# Gene count distribution ")
  gene_sum <- summary(data@meta.data[,paste0("nFeature_", assay)], na.rm = T)
  log_info(paste0(names(gene_sum), ": ", sapply(gene_sum, function(x) round(x, 1)), collapse = ", "))
  
  n_cells <- length(data@meta.data %>% filter(!!sym(paste0("nFeature_", assay)) < min_genes))
  log_warn("There are ", n_cells, " (", round(n_cells/ncol(data)*100, 1), "%) cells with a gene count < ", min_genes)
  keep_genes <- data@meta.data %>% filter(!!sym(paste0("nFeature_", assay)) >= min_genes) %>% pull(cell_names)
  
  # VOLUME
  #-----------------------------------------------------------------------------
  n_cells_min_vol <- nrow(data@meta.data %>% filter(!!sym(size_feature) < min_volume))
  n_cells_max_vol <- nrow(data@meta.data %>% filter(!!sym(size_feature) > max_volume))
  log_warn("There are ", n_cells_min_vol, " (", round(n_cells_min_vol/ncol(data)*100, 1), "%) cells with a volume/area < ",
           min_volume, " and ", n_cells_max_vol, " (", round(n_cells_max_vol/ncol(data)*100, 1), "%) with a volume/area > ",
           max_volume)
  
  keep_volume <- data@meta.data %>% 
    filter(!!sym(size_feature) >= min_volume & !!sym(size_feature) <= max_volume) %>% 
    pull(cell_names)
  
  log_info("# Volume/Area distribution ")
  size_sum <- summary(data@meta.data[,size_feature], na.rm = T)
  log_info(paste0(names(size_sum), ": ", sapply(size_sum, function(x) round(x, 1)), collapse = ", "))
  
  pdf(file.path(res_path, "qc_hist_volume.pdf"), height = 3, width = 4)
  p <- data@meta.data %>%
    rename_at(all_of(size_feature), ~"feature") %>%
    ggplot(aes(x = .data$feature)) +
    geom_histogram(fill = proj_cols("light blue"), color = "black", alpha = .8, bins = 50) +
    geom_vline(xintercept = min_volume, col = proj_cols("pink"))+
    geom_vline(xintercept = max_volume, col = proj_cols("pink"))+
    geom_vline(xintercept = size_sum["Mean"], col = proj_cols("yellow"))+
    labs(x = "Volume/Area", y = expression(paste(italic('N'),' Cells')))+
    theme_proj()+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "in"))
  print(p)
  dev.off()
  
  pdf(file.path(res_path, "qc_fov_volume.pdf"), height = 5.5, width = 5)
  p <- qc_fov(data, 
              min_feature = min_volume, 
              max_feature = max_volume, 
              assay = "Vizgen", 
              fov = fov_name, 
              legend_title = "Volume/Area", 
              feature = size_feature, 
              feature_label = " ")
  print(p)
  dev.off()
  
  # SUBSET CELLS THAT PASS QC
  #-----------------------------------------------------------------------------
  keep_cells <- names(which(table(c(keep_umi, keep_genes, keep_volume)) == 3))
  
  n_cells_remove <- ncol(data) - length(keep_cells)
  pct_cells_remove <- n_cells_remove/ncol(data)
  log_warn(n_cells_remove, " (", round(pct_cells_remove*100, 2), "%) cells will be removed.")
  
  p1 <- data@meta.data %>%
    ggplot(aes(x = !!sym(paste0("nCount_", assay)), y = !!sym(size_feature), color = status)) +
    geom_point(color = proj_cols("light blue"), alpha = 0.7) +
    geom_vline(xintercept = min_umi, col = proj_cols("pink"))+
    geom_vline(xintercept = max_umi, col = proj_cols("pink")) +
    geom_hline(yintercept = min_volume, col = proj_cols("blue"))+
    geom_hline(yintercept = max_volume, col = proj_cols("blue"))+
    labs(x = "UMI Count", y = "Cell Size")+
    theme_proj()+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "in"))
  
  p2 <- data@meta.data %>%
    ggplot(aes(x = !!sym(paste0("nFeature_", assay)), y = !!sym(size_feature), color = status)) +
    geom_point(color = proj_cols("light blue"), alpha = 0.7) +
    geom_vline(xintercept = min_genes, col = proj_cols("pink"))+
    geom_hline(yintercept = min_volume, col = proj_cols("blue"))+
    geom_hline(yintercept = max_volume, col = proj_cols("blue"))+
    labs(x = "Gene Count", y = "Cell Size")+
    theme_proj()+
    theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "in"))
  
  pdf(file.path(res_path, "qc_combined.pdf"), height = 3, width = 7)
  p <- wrap_plots(p1, p2)
  print(p)
  dev.off()
  
  ## BEFORE QC
  p5 <- qc_fov_feature(data, feature = paste0("nCount_", assay), title = "UMI (Before QC)", fov = fov_name)
  p6 <- qc_fov_feature(data, feature = paste0("nFeature_", assay), title = "Genes (Before QC)", fov = fov_name)
  p7 <- qc_fov_feature(data, feature = size_feature, title = "Volume/Area", fov = fov_name)
  
  ## QC
  data <- subset(x = data, subset = cell_names %in% keep_cells)
  lapply(capture.output(data), log_trace)
  
  ## AFTER QC
  p8 <- qc_fov_feature(data, feature = paste0("nCount_", assay), title = "UMI (Before QC)", fov = fov_name)
  p9 <- qc_fov_feature(data, feature = paste0("nFeature_", assay), title = "Genes (Before QC)", fov = fov_name)
  p10 <- qc_fov_feature(data, feature = size_feature, title = "Volume/Area", fov = fov_name)
  
  pdf(file.path(res_path, "qc_final_plots.pdf"), height = 9, width = 7)
  p <- wrap_plots(p5, p8, p6, p9, p7, p10) + plot_layout(nrow = 3)
  print(p)
  dev.off()
  
  return(data)
}

qc_fov_feature <- function(data, feature = NULL, title = NULL, fov = NULL){
  
  plot_data <- data.frame(cell_names = data@images[[fov]]@boundaries$centroids@cells,
                          data@images[[fov]]@boundaries$centroids@coords) %>%
    left_join(data@meta.data) %>%
    rename_at(all_of(feature), ~"feature")
  
  p <- ggplot()+
    geom_point(data = plot_data, aes(x = x, y = y, color = feature), size = 0.1, stroke = 0)+
    labs(title = title)+
    coord_flip() +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_umap()+
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          legend.margin = margin(l = -0.15, unit = "in"),
          legend.key.width = unit(0.2, "in"),
          legend.key.height = unit(0.25, "in"),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1.5),
          plot.title = element_text(face = "bold", size = 14, margin = margin(b = 0.05, unit = "in")),
          axis.line.x = element_line(linewidth = 0.2, color = "black"),
          axis.line.y = element_line(linewidth = 0.2, color = "black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank()) 
  
  return(p)
  
}

qc_fov <- function(data, 
                   min_feature, 
                   max_feature, 
                   assay = "Vizgen", 
                   fov = NULL, 
                   legend_title = NULL, 
                   feature = NULL, 
                   feature_label = NULL){
  
  colors <- unname(c(proj_cols_grey("light grey"), proj_cols("blue", "pink")))
  names(colors) <- c("Retained", paste0(feature_label, " < ", min_feature), paste0(feature_label, " > ", max_feature))
  
  data@meta.data <- data@meta.data %>%
    mutate(threshold = ifelse(!!sym(feature) < min_feature, paste0(feature_label, " < ", min_feature), "Retained"),
           threshold = ifelse(!!sym(feature) > max_feature, paste0(feature_label, " > ", max_feature), threshold)) %>%
    mutate(threshold = factor(threshold, levels = names(colors)))
  
  plot_data <- data.frame(cell_names = data@images[[fov]]@boundaries$centroids@cells,
                          data@images[[fov]]@boundaries$centroids@coords) %>%
    left_join(data@meta.data)
  
  p <- ggplot(data = plot_data, aes(x = x, y = y, color = threshold))+
    geom_point(data = plot_data %>% filter(threshold=="Retained"), size = 0.1, stroke =0)+
    geom_point(data = plot_data %>% filter(threshold!="Retained"), size = 0.1, stroke = 0)+
    guides(color = guide_legend(title = legend_title, override.aes = aes(size = 3), title.position="top", title.hjust = 0.5))+
    coord_flip() +
    scale_color_manual(values = colors) +
    theme_umap()+
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key = element_rect(color = NA),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10, margin = margin(l = 0.05, unit = "in")),
          panel.border = element_rect(fill = NA, color = "black", linewidth = 1),
          axis.line.x = element_line(color = "black", linewidth = 0.25),
          axis.line.y = element_line(color = "black", linewidth = 0.25),
          axis.ticks.x =  element_line(color = "black", linewidth = 0.5),
          axis.ticks.y =  element_line(color = "black", linewidth = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10, family = "Helvetica"),
          axis.text.y = element_text(size = 10, family = "Helvetica"))
  
  return(p)
}
