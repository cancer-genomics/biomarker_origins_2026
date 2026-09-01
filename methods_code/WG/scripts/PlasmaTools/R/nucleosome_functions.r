bedSegment <- function(fp, chr, start, end) {
    syscmd <- paste0("awk '($1 == \"", chr,"\") && ($2 >= ",start,
                     ") && ($3 <= ",end,")' ", fp)
    dt <- fread(cmd=syscmd)
    setnames(dt, c("chr", "start", "end", "mapq"))
    dt[,start:=start+1]
    gr <- makeGRangesFromDataFrame(dt, keep.extra.columns=TRUE)
    gr
}

getCoverage <- function(fragments) {
    if(class(fragments) != "GRanges") {
        fragments <- makeGRangesFromDataFrame(fragments)
    }
    dj <- disjoin(fragments, ignore.strand=TRUE)
    ## for now add gaps back in
    dj <- c(dj, gaps(dj)[-1])
    dj$cov <- countOverlaps(dj, fragments)
    dj <- sort(dj)
    dj
}

getWPS <- function(fragments, windows=NULL, size=120) {
    if(class(fragments) != "GRanges") {
        fragments <- makeGRangesFromDataFrame(fragments)
    }
    if(is.null(windows)) {
        gr <- range(fragments)
        windows <- get.windows(gr, size=size)
    }
    wdth <- width(fragments)
    fragments2 <- fragments[wdth >= 120 & wdth < 200]
    wps <- 2*countOverlaps(windows, fragments2, type="within") -
        countOverlaps(windows, fragments)
    wps
}

getWindows <- function(gr, size=120) {
    library(BSgenome.Hsapiens.UCSC.hg19)
    chr <- as.character(seqnames(gr))
    window.start <- seq(start(gr), end(gr)-size)
    window.end <- window.start + size
    windows <- GRanges(Rle(chr, length(window.start)),
                       IRanges(window.start, window.end))
    windows$gc <- GCcontent(Hsapiens, windows)
    windows
}

integrateSummaries <- function(windows, cov) {
    mid <- (start(windows)  + end(windows))/2
    windows2 <- GRanges(seqnames(windows), IRanges(mid, width=1))
    mcols(windows2) <- mcols(windows)
    windows <- windows2
    hits <- findOverlaps(cov, windows)

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
    gr3$gc <- NA

    k <- queryHits(hits1)
    m <- subjectHits(hits1)
    gr3[k]$wps <- windows[m]$wps
    gr3[k]$wps.kz <- windows[m]$wps.kz
    gr3[k]$gc <- windows[m]$gc

    gr3
}

pileFragments <- function(frags) {
    height <- 1
    sep <- 0.5
    bins <- disjointBins(IRanges(start(frags), end(frags) +1 ))
    ybottom <- bins*(sep + height) - height
    ytop <- ybottom + height

    ## Add motif color
    frags2 <- GRanges(seqlevels(frags), IRanges(start=start(frags)-1,
                                                end=end(frags)+1))

    seq <- getSeq(Hsapiens, frags2)
    forward <- substring(seq, 1, 3)
    reverse <- substring(reverseComplement(seq), 1, 3)
    tcc <- ifelse(forward == "TCC" ,#| reverse == "TCC",
                  "TCC", "Not TCC")
    ###

    dt <- data.table(xleft=start(frags) - 0.5, xright=end(frags) + 0.5,
                     ybottom = ybottom, ytop = ytop, tcc=tcc)
    dt
}

ggNucleosomes <- function(win, fragment.dt) {
    gr <- range(win)
    dt = as.data.table(win)
    dt <- dt[!is.na(cov), ]
    dt[,cov2:=as.numeric(filterv(cov, filter=rep(1/20, 20)))]

    mytheme <- theme_classic() %+replace%
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face="bold",angle=90),
          plot.margin=unit(c(-0.1, 0.1, 0, 0.1), "cm"),
          panel.grid.major.x = element_line(colour = "lightgray"),
          legend.position = "none")

#     rect <- data.frame(xmin=start(peaks2)-0.5, xmax=end(peaks2)+0.5, ymin=0.6, ymax=1.0)
#     rect2 <- data.frame(xmin=start(peaks)-0.5, xmax=end(peaks)+0.5, ymin=0.1, ymax=0.5)
    m <- nrow(dt)
    xlim <- gr
    xlim <- c(start(ranges(reduce(xlim, ignore.strand = TRUE))),
              end(ranges(reduce(xlim, ignore.strand = TRUE))))
    fig1 <- ggplot(dt, aes(start, cov2)) + geom_ribbon(aes(ymin=rep(0, m), ymax=cov2, fill="green")) + scale_fill_manual(values=c("green" = rgb(49, 163, 84, maxColorValue=255)), guide=FALSE) +  mytheme
    fig1  <- fig1 + xlim(xlim) + ylab("Coverage") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    fig1 <- fig1 + geom_line(color=rgb(28, 125, 57, maxColorValue=255), lwd=0.5)
    fig1 <- fig1 + scale_y_continuous(expand=c(0,0))
#     fig1 <- fig1 + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
#                              color="gray20",
#                              alpha=0.5,
#                              inherit.aes=FALSE)
#     fig1 <- fig1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
#                              color="gray20",
#                              alpha=0.5,
#                              inherit.aes=FALSE)
 
    wps.fig <- ggplot(dt, aes(start, wps.kz)) + geom_line(color="tomato2", size=1.1) + geom_hline(yintercept=0, color="gray40", linetype="dashed") + mytheme
    wps.fig <- wps.fig + xlim(xlim) + ylab("WPS")  + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

    gc.fig <- ggplot(dt, aes(start, gc)) + geom_line(color="gray", size=1.1) + mytheme
    gc.fig <- gc.fig + xlim(xlim) + ylab("GC") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())


    frags.fig <- ggplot(fragment.dt, aes(xmin=xleft, xmax=xright, ymin=ybottom, ymax=ytop, fill=tcc)) + geom_rect(color="gray40") + mytheme
    frags.fig <- frags.fig + xlim(xlim) + labs(x=seqlevels(gr), y="Fragments")
    frags.fig <- frags.fig + theme( plot.margin=unit(c(-0.1, 0.1, 0.5, .1), "cm"), axis.title.x = element_text(face="bold"))
    frags.fig <- frags.fig + scale_y_continuous(expand=c(0,0))

    plot_grid(fig1, wps.fig, gc.fig, frags.fig, ncol=1, align="v", rel_heights=c(2.2,1,1,2))
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
