

### Write function for getting WPS score for a region
get.wps <- function(fragments, windows=NULL, size=120) {

    gr <- range(windows)
    fragments2 <- get.fragments2(fragments, gr)
    wdth <- width(fragments2)
    fragments2 <- fragments2[wdth >= 120 & wdth < 200]
#     windows <- get.windows(gr, size=size)
    wps <- 2*countOverlaps(windows, fragments2, type="within") - countOverlaps(windows, fragments)
#     wps <- 2*countOverlaps(windows, fragments2, type="within") - countOverlaps(windows, fragments)
    wps
    ### median center?
}

### Write function for getting WPS score for a region
get.wps2 <- function(fragments, windows=NULL, size=120) {

    gr <- range(windows)
    fragments2 <- get.fragments2(fragments, gr)
    wdth <- width(fragments2)
    fragments2 <- fragments2[wdth >= 120 & wdth < 200]
#     windows <- get.windows(gr, size=size)
    fragments.spanning <- countOverlaps(windows, fragments2, type="within")
    fragments.any <- countOverlaps(windows, fragments2, type="any")
    fragments.notspanning <- fragments.any - fragments.spanning
    wps <- fragments.spanning/(fragments.spanning + fragments.notspanning)
#     wps <- 2*countOverlaps(windows, fragments2, type="within") - countOverlaps(windows, fragments)
    wps
    ### median center?
}

get.starts <- function(fragments, windows=NULL, size=120) {
    gr <- range(windows)
    fragments2 <- get.fragments2(fragments, gr)
    start(fragments2)

}


get.reads <- function(file, gr, bamparams) {
    # Read bam
    bam <- scanBam(file, param=bamparams)[[1]]
    ## find fragments from these
    chr <- as.character(seqnames(gr))
    if(length(bam[["pos"]]) == 0) return(GRanges())
    else
        reads <- GRanges(seqnames=chr, IRanges(bam[["pos"]], width=bam[["qwidth"]]), strand=bam[["strand"]])
    reads <- reads[!duplicated(reads)]
    reads
}

get.disjoin <- function(fragments, gr) {
    rg <- range(fragments)[1]

    chr <- as.character(seqnames(gr))
    ## ignore strand or not?
    dj <- disjoin(fragments, ignore.strand=TRUE)
    ## for now add gaps back in
    dj <- c(dj, gaps(dj)[-1])
    dj$cov <- countOverlaps(dj, fragments)
    dj <- sort(dj)

    dj2 <- dj[queryHits(findOverlaps( dj, gr, type="within"))]

    dj2
}

get.fragments <- function(file, gr, bamparams) {
    chr <- as.character(seqnames(gr))
#     if(length(reads)==0) return(reads)
#     reads2 <- reads[queryHits(findOverlaps( reads, gr, type="within"))]
    ##galignmentpairs

#     fbamparams <- ScanBamParam(flag = bamflags, which=gr, what=scanBamWhat())
#     fragments <- readGAlignmentPairs(file, param=bamparams)
    fragments <- readGAlignmentPairs(BamFile(file, asMates=TRUE), param=bamparams)
    if(length(fragments) == 0) return(GRanges())
#     starts <- rowMins(cbind(start(first(fragments)), start(last(fragments))))
#     ends <- rowMaxs(cbind(end(first(fragments)), end(last(fragments))))
#     fragments <- GRanges(seqnames=chr, IRanges(starts, ends))

    ##### REMOVE DUPLICATES
    fragments <- GRanges(fragments)
    fragments <- unique(fragments)

#     fragments <- fragments[width(fragments) < 200]
    fragments
}

get.fragments2 <- function(fragments, gr) {
    fragments2 <- fragments[queryHits(findOverlaps( fragments, gr, type="within"))]
    fragments2
}

get.windows <- function(gr, size=120) {
    chr <- as.character(seqnames(gr))
    window.start <- seq(start(gr), end(gr)-size)
    window.end <- window.start + size
    windows <- GRanges(Rle(chr, length(window.start)), IRanges(window.start, window.end))

}

##### 
integrate.summaries <- function(windows, cov) {
    mid <- (start(windows)  + end(windows))/2
    windows2 <- GRanges(seqnames(windows), IRanges(mid, width=1))
    mcols(windows2) <- mcols(windows)
    windows <- windows2
    hits <- findOverlaps(cov, windows)
#     windows2 <- windows[subjectHits(hits)]

    gr3 <- c(granges(windows), granges(cov))
    gr3 <- disjoin(gr3)
    hits1 <- findOverlaps(gr3, windows)
    hits2 <- findOverlaps(gr3, cov)

    gr3$cov <- NA
    i <- queryHits(hits2)
    j <- subjectHits(hits2)
    gr3$cov[i] <- cov$cov[j]

    gr3$wps <- NA
    gr3$wps.kz <- NA
#     gr3$wps.kz2 <- NA
    gr3$gc <- NA

    k <- queryHits(hits1)
    m <- subjectHits(hits1)
    gr3$wps[k] <- windows$wps[m]
    gr3$wps.kz[k] <- windows$wps.kz[m]
#     gr3$wps.kz2[k] <- windows$wps.kz2[m]
    gr3$gc[k] <- windows$gc[m]

    gr3
}

read.piles <- function(reads) {
    height <- 1
    sep <- 0.5
    bins <- disjointBins(IRanges(start(reads), end(reads) +1 ))
    ybottom <- bins*(sep + height) - height
    ytop <- ybottom + height

    df <- data.frame(xleft=start(reads) - 0.5, xright=end(reads) + 0.5,
                     ybottom = ybottom, ytop = ytop)
    df
}

get.nsome <- function(gr.summary, gr) {
    gr2 <- gr.summary
    rl <- rle(gr2$wps.kz > 0)

    rl$values[((rl$lengths >= 50 & rl$lengths <= 450) & rl$value == TRUE)] 

    chr <- as.character(seqnames(gr))
    ind <- which(rl$values == TRUE)
    if(ind[1] == 1) ind <- ind[-1]
    st <- rep(NA, length(ind))
    for(s in seq_along(st)) st[s] <- sum(rl$lengths[1:(ind[s]-1)])+1
    st <- start(gr2)[st]
    en <- st + rl$lengths[ind] 
    nsomes <- GRanges(seqnames=chr, IRanges(st, en))
    
}


nsome.peak <- function(gr.summary, gr) {
    nas <- is.na(gr.summary$wps.kz)
    g <- gr.summary[gr.summary$wps.kz > 0 & !nas]
    g.red <- reduce(g)
    g.red <- g.red[width(g.red) >= 50 & width(g.red) <= 450]

    hits <- findOverlaps(g.red, g)
    g.list <- split(g[subjectHits(hits)], queryHits(hits))
    
    glist.2 <- sapply(g.list, function(x) {
                      med <- median(x$wps.kz)
                      x2 <- x[x$wps.kz > med]
                      x.red <- reduce(x2)
                      x.red <- x.red[which.max(width(x.red))]
                      ## Sometimes this is empty, throwing error when trying
                      ## to save metadata
                      x3 <- x2[queryHits(findOverlaps(x2, x.red))]
                      if(length(x3) == 0 | all(is.na(x3$cov))) GRanges()
                      else{
                          x.red$medcov <- median(x3$cov)
                          x.red$maxcov <- max(x3$cov)
                          x.red$maxwps <- round(max(x3$wps.kz), 2)
                          x.red$whichmaxcov <- start(x3[which.max(x3$cov)])
                          x.red$whichmaxwps <- start(x3[which.max(x3$wps)])
                          x.red$meangc <- as.integer(round(mean(x3$gc)))
                          x.red
                      }
                             })

#     nsome.peaks <- unlist(GRangesList(sapply(g.list, function(x) x[which.max(x$wps.kz)])))
    unlist(GRangesList(glist.2)) 
}

## todo: handle when windows is NULL
## currently ignores fragment input if both fragments and file provided
nsome.track <- function(file=NULL, windows=NULL, fragments=NULL){
    if(is.null(fragments) & is.null(file)) stop("No bamfile or fragments GRanges provided")
    gr <- range(windows)
    if(is.null(windows$gc)) {
        gc <- as.numeric(GCcontent(Hsapiens, windows))
        windows$gc <- as.integer(round(gc * 100))
    }
    if(!is.null(file)) {
        bamflags <- scanBamFlag(isProperPair = TRUE, isDuplicate = FALSE, isUnmappedQuery=FALSE)
        bamparams <- ScanBamParam(flag = bamflags, which=gr, what=scanBamWhat())

#         reads <- get.reads(file=file, gr=gr, bamparams=bamparams)
        fragments <- get.fragments(file, gr, bamparams)

        ## If there are no reads/fragments in a window then return empty GRanges object
        if(length(fragments)==0) return(fragments)
        #     fragments2 <- get.fragments2(fragments, gr)
    }
    else {
        fragments <- fragments[subjectHits(findOverlaps(fragments, windows))]
        fragments <- fragments
    }
    wps <- get.wps(fragments, windows)
#     wps2 <- get.wps2(fragments, windows)
    wps.kz <- kz(wps - runmed(wps, 201, "constant"), w=21)
#     wps.kz2 <- kz(wps2 - runmed(wps2, 201, "constant"), w=21)

    windows$wps <- wps
    windows$wps.kz <- wps.kz
#     windows$wps.kz2 <- wps.kz2

    cov <- get.disjoin(fragments, gr)
    gr2 <- integrate.summaries(windows, cov)

    peaks2 <- nsome.peak(gr2, gr)

    peaks2
}
