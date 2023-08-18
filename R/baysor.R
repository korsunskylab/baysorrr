try_cmd = function(cmd, attempts_left, verbose=FALSE) {
    if (attempts_left==0) {
        msg = glue::glue('Failed to run {cmd}.')
        if (verbose) message(msg)
        # return(msg)
        return(1L)
    }
    tryCatch({
        msg = system(cmd, intern = TRUE)
        if (verbose) message(msg)
        # return(msg)
        return(0L)
    }, error = function(msg) {
        if (verbose) message(glue::glue('Try to run again: {cmd}'))
        system('sleep 10')
        try_cmd(cmd, attempts_left-1)
    }, warning = function(msg) {
        if (verbose) message(glue::glue('Try to run again: {cmd}'))
        system('sleep 10')
        try_cmd(cmd, attempts_left-1)
    })    
}



baysor.collect_tx = function(dir) {
    tx <- list.files(dir, pattern = "segmentation.csv$", recursive = TRUE, full.names = TRUE) %>% 
        purrr::map(data.table::fread) %>% data.table::rbindlist(idcol = "tile")
    if (is(tx$cell, "character")) {
        ## Assumes cell is of format "{name}-{number}"
        cell_id = stringr::str_split(tx$cell, '-') %>% purrr::map(2)
        cell_id[purrr::map_lgl(cell_id, is.null)] = 0
        tx$cell = as.integer(cell_id)
        
        ## OLD VERSION: assign random integers
        # tx[, cell_int := as.integer(factor(cell)) - 1, by = tile] ## renaming directly keeps cell as chr - why? 
        # tx[, cell := cell_int][, cell_int := NULL]
        # tx$cell = as.integer(tx$cell)
    }
    if (length(unique(tx$tile)) > 1) {
        tx <- unique(tx[cell != 0, .(tile, cell)])[
            , .N, by = tile
        ][
            , 
            `:=`(offset, cumsum(N))
        ][
            , `:=`(offset, c(0L, offset[1:(.N - 1)]))
        ][
            tx, on = "tile"
        ][
            cell != 0, `:=`(cell, cell + offset)
        ][, `:=`(N = NULL, offset = NULL)][]
    }
    # tx[, `:=`(N, sum(!grepl("Blank", gene))), by = cell][, `:=`(cell, case_when(N >= mintx ~ cell, TRUE ~ 0L))] ## OLD OPTION: MIN TRANSCRIPTS
    return(tx)
}
    
baysor.read_shapes = function(dir) {
    fnames = list.files(dir, pattern='segmentation_polygons.json', recursive=TRUE, full.names=TRUE)
    json_list = fnames %>% purrr::map(jsonlite::read_json) %>% purrr::map('geometries')
    shapes = purrr::map(json_list, function(json) {
        cell_shapes = json %>% 
            purrr::map('coordinates') %>%purrr::map(1) %>% 
            purrr::map(function(y) cbind(purrr::map_dbl(y, 1), purrr::map_dbl(y, 2))) %>% 
            purrr::map(function(cell) {
                if (nrow(cell) < 3) {
                    ## some Baysor cells get returned as lines? 
                    sf::st_polygon(list())
                } else {
                    sf::st_polygon(list(rbind(cell, tail(cell, 1))))  
                }
            })
        cell_ids = purrr::map_int(json, 'cell') 
        sf::st_sfc(cell_shapes[order(cell_ids)])
    }) %>% 
        purrr::reduce(c)    
}

mean_hex = function(str) {
    paste0(
        '#',
        sprintf('%02X', as.integer(mean(as.hexmode(substr(str, 2, 3))))),
        sprintf('%02X', as.integer(mean(as.hexmode(substr(str, 4, 5))))),
        sprintf('%02X', as.integer(mean(as.hexmode(substr(str, 6, 7)))))
    )    
}

baysor.collect_cells = function(dir, no_ncv_estimation) {
    cells = sf::st_sf(shape=baysor.read_shapes(dir))
    cells = cbind(cells, sf::st_coordinates(sf::st_centroid(cells$shape))) %>% dplyr::rename(x = X, y = Y)
    cells$area <- sf::st_area(cells$shape)

    if (no_ncv_estimation) {
        cell_summary = tx[
            cell != 0, 
            .(
                segmentation_tile = unique(.SD$tile), 
                n_transcripts = .N, 
                avg_confidence = mean(.SD$assignment_confidence), 
                # ncv_color = mean_hex(ncv_color), 
                cluster = .SD[, .N, by = cluster][order(-N)][1, cluster]
            ), 
            by = cell
        ][order(cell)]
    } else {
        cell_summary = tx[
            cell != 0, 
            .(
                segmentation_tile = unique(.SD$tile), 
                n_transcripts = .N, 
                avg_confidence = mean(.SD$assignment_confidence), 
                ncv_color = mean_hex(.SD$ncv_color), 
                cluster = .SD[, .N, by = cluster][order(-N)][1, cluster]
            ), 
            by = cell
        ][order(cell)]
    }
        
    cells = cbind(cells, cell_summary) 
    return(cells)    
}

baysor.finish = function(remove_temp_files) {
    if (remove_temp_files) {
        for (fname in list.files(output_dir, full.names = TRUE)) {
            unlink(fname, recursive = TRUE)
        }
    }
    data.table::fwrite(tx, file.path(output_dir, "transcripts.csv"))
    writeLines(rownames(counts), file.path(output_dir, "genes.tsv"))
    spatula::writeMM(counts, file.path(output_dir, "counts.mtx"))
    suppressWarnings({
        sfarrow::st_write_parquet(dplyr::select(cells, shape), file.path(output_dir, "shapes.parquet"))
    })
    data.table::fwrite(sf::st_drop_geometry(cells), file.path(output_dir, "cells.csv"), sep = ",")
}

#' @export 
baysor.run = function(
    tx, 
    baysor_binpath, 
    output_dir, 
    n_clusters = 10,
    prior_segmentation_confidence = 0.7, 
    max_tx_per_tile = 5e6, 
    scale = 5, ## need to specify in v0.6.0, even with cellpose priors 
    no_ncv_estimation = FALSE, 
    min_molecules_per_cell = 10, ## doesn't seem to work? 
    remove_temp_files = TRUE, ## only set to FALSE for debugging purposes 
    max_attempts = 2 ## retry Baysor up to max_attempts-1 times. max_attempts=1 means no retries.
) {
    message('----- This Function is Designed to Work with Baysor v0.6.0 -----')
    stopifnot(all(colnames(tx) == c('x', 'y', 'gene', 'cell')))
    stopifnot(file.exists(baysor_binpath))
    
    ## split transcripts into tiles
    environment(split_tx_files) <- environment()
    ntiles <- split_tx_files(output_dir, max_tx_per_tile) 

    ncv_str = ''
    if (no_ncv_estimation) ncv_str = '--no-ncv-estimation'
    
    ## run baysor in each tile
    cmds = purrr::map_chr(glue::glue('{output_dir}/g{1:ntiles}/'), function(outdir) {
        as.character(glue::glue('{baysor_binpath} run --x-column x --y-column y --gene-column gene --scale={scale} --save-polygons=GeoJSON --min-molecules-per-cell={min_molecules_per_cell} {ncv_str} --prior-segmentation-confidence={prior_segmentation_confidence} --n-clusters={n_clusters}  -o {outdir} {outdir}/tx_baysor.csv :cell'))
    }) 
    baysor_err = furrr::future_map_int(cmds, try_cmd, attempts_left=max_attempts)
    if (any(baysor_err == 1L)) {
        failed_tiles = paste(which(baysor_err == 1L), collapse = ' ')
        warning(glue::glue('baysor failed to complete on tiles {failed_tiles}'))
    }
    
    ## stitch and summarize results 
    tx <- baysor.collect_tx(output_dir)
    counts <- tx_to_counts(tx$gene, tx$cell, remove_bg = TRUE)
    environment(baysor.collect_cells) <- environment()
    cells = baysor.collect_cells(output_dir, no_ncv_estimation)  
    
    ## QC: remove low count cells 
    ##     why doesn't Baysor do this internally? 
    i_keep = which(cells$n_transcripts >= min_molecules_per_cell)
    tx$cell[!tx$cell %in% i_keep] = 0
    cells = cells[i_keep, ]
    counts = counts[, i_keep]    

    ## cache results
    environment(baysor.finish) <- environment()
    baysor.finish(remove_temp_files) 
}

