split_tx_files <- function(output_dir, max_tx_per_grid) {
    ## Check if dir exists
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    ## Check if dir is writeable 
    if (file.access(output_dir, mode = 2) == -1L) {
        stop('Specified output_dir does not have write permissions.')
    }

    grid <- split_tx(
        tx, 
        max_tx = max_tx_per_grid, 
        max_voxels = round(nrow(tx) / max_tx_per_grid) * 5 ## to prevent worst case behavior of too much splitting for 
    ) 


    err_status <- grid %>%
        purrr::map('tx') %>% 
        purrr::imap(function(tx_grid, grid_id) {
            dirname <- file.path(output_dir, paste0('g', grid_id))        
            if (!dir.exists(dirname)) dir.create(dirname)
            fname <- as.character(glue::glue('{dirname}/tx_baysor.csv'))
            data.table::fwrite(tx_grid, fname, sep = ',') 
        }) 
    return(length(grid))
}

split_grid <- function(grid_point) {
    bbox <- grid_point$bbox
    xmid <- mean(bbox[1:2])
    ymid <- mean(bbox[3:4])
    ## NOTE: the 1e-10 term avoids duplicates 
    bboxes_new <- list(
        c(bbox[1], xmid - 1e-10, bbox[3], ymid - 1e-10),
        c(xmid, bbox[2], bbox[3], ymid - 1e-10),
        c(xmid, bbox[2], ymid, bbox[4]),
        c(bbox[1], xmid - 1e-10, ymid, bbox[4])
    )
    
    ## Make four new grid points 
    res <- purrr::map(bboxes_new, function(bbox_test) {
        res <- list(
            tx = dplyr::filter(grid_point$tx, between(x, bbox_test[1], bbox_test[2]) & between(y, bbox_test[3], bbox_test[4])), 
            # tx = grid_point$tx[between(x, bbox_test[1], bbox_test[2]) & between(y, bbox_test[3], bbox_test[4])], 
            bbox = bbox_test,
            bbox_geom = st_rectangle(bbox_test[1], bbox_test[2], bbox_test[3], bbox_test[4])
        )
        res$n <- nrow(res$tx)
        return(res)
    })
}


split_tx <- function(tx, max_tx, max_voxels) {
    ## Initialize grid with all transcripts
    grid <- list(
        list(
            tx = tx, 
            bbox = c(min(tx$x), max(tx$x), min(tx$y), max(tx$y))
        )
    )
    grid[[1]]$n <- nrow(grid[[1]]$tx)
    grid[[1]]$bbox_geom <- st_rectangle(grid[[1]]$bbox[1], grid[[1]]$bbox[2], grid[[1]]$bbox[3], grid[[1]]$bbox[4])

    ## Keep splitting grid points until each has at most max_tx transcripts
    .i <- 0
    while (TRUE) {
        .i <- .i + 1
        grids_split <- which(purrr::map_int(grid, 'n') > max_tx)
        if (length(grid) >= max_voxels) {
            break
        } else if (length(grids_split) > 0) {
            i <- grids_split[1]
            grid <- append(grid, split_grid(grid[[i]]))
            grid[[i]] <- NULL
            grids_split <- which(purrr::map_int(grid, 'n') > max_tx)
        } else {
            break
        }
    }
    ## For QC purposes, compute the transcript density of each region 
    for (i in seq_len(length(grid))) {
        grid[[i]]$density <- grid[[i]]$n / sf::st_area(grid[[i]]$'bbox_geom')
    }    
    return(grid)
}


