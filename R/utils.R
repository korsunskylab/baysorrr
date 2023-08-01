#' @export 
st_assign_pts_to_polygon = function(pts, polygons) {
    ores = st_intersects(pts, polygons) 
    ores = map(ores, head, 1)
    ores = as.integer(ores)
    ores[which(is.na(ores))] = 0
    return(ores)
}


writeMM <- function(X, fname) {
    stopifnot(is(X, 'dgCMatrix'))
    nelem <- length(X@i)
    nrow <- X@Dim[1]
    ncol <- X@Dim[2]    
    writeLines(paste0(c('%%MatrixMarket matrix coordinate real general\n', nrow, ' ', ncol, ' ', nelem), collapse = ''), fname)    
    fwrite(data.table(X@i+1, rep(seq_len(ncol), times = diff(X@p)), X@x), fname, append = TRUE, sep = ' ')        
}

readMM <- function(fname, max_header_size = 100, nthreads = NULL) {
    ## First, figure out how many lines to skip
    ## Assumes that MM format has comments that start with %
    nlines_skip <- 0
    con <- file(fname, open = 'r')
    for (i in seq_len(max_header_size)) {
        line <- readLines(con, 1)
        nlines_skip <- nlines_skip + 1
        if (!grepl('^\\W*\\%', line)) {
            ## This is the line with dimension information 
            ## We need dimension information to handle empty rows and columns 
            nrow <- as.integer(strsplit(line, ' ')[[1]][1])
            ncol <- as.integer(strsplit(line, ' ')[[1]][2])
            break
        }
    }
    close(con)

    if (is.null(nthreads)) {
        nthreads <- data.table::getDTthreads()
    }
    ## Then, read the file and make a matrix 
    with(
        fread(fname, skip = nlines_skip, nThread = nthreads),
        Matrix::sparseMatrix(i = V1, j = V2, x = V3, dims = c(nrow, ncol))
    )
}

st_rectangle <- function(xmin, xmax, ymin, ymax) {
    res <- st_polygon(list(
        rbind(
            c(xmin, ymin),
            c(xmax, ymin), 
            c(xmax, ymax), 
            c(xmin, ymax), 
            c(xmin, ymin))
    ))
    return(res)    
}

tx_to_counts <- function(genes, cells, remove_bg = TRUE) {
    if (remove_bg) {
        idx <- which(cells != 0)
        cells <- cells[idx]
        genes <- genes[idx]
    }
    genes <- factor(genes)
    cells <- factor(cells)
    counts <- Matrix::sparseMatrix(
        i = as.integer(genes), 
        j = as.integer(cells), 
        x = rep(1, length(genes)),
        dims = c(length(levels(genes)), length(levels(cells)))
    )
    rownames(counts) <- levels(genes)
    colnames(counts) <- levels(cells)
    return(counts)
}