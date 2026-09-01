# Plot coverage using base R.
# Pass in granges object of reads
plotCoverage <- function(gr, col="blue", xlab="Index", ylab="Coverage", main="") {
    rg <- range(gr)
    start <- start(rg)
    end <- end(rg)
#     xWindow <- as.vector(window(cv, start, end))
#      smcov <- round(runmean(gr$cov, 21, endrule="constant"))
#      smcov <- round(runmed(gr$cov, 21, endrule="constant"))
    xWindow <- rep(gr$cov, width(gr))
    ## how to filter? kz?
    xWindow <- filterv(xWindow, filter=rep(1/20, 20))
    x <- start:end
    x <- x[!is.na(xWindow)]
    xWindow <- xWindow[!is.na(xWindow)]
#     xWindow[!is.na(xWindow)] <- 0
    xlim <- c(start, end)
    ylim <- c(0, max(xWindow))
#     ylim <- c(0, 30)
    plot(x = start, y = 0, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, type = "n", xaxt="n", frame.plot=FALSE)
    polygon(c(start, x, end), c(0, xWindow, 0), col = col, border=col) 
}

plotCoverage2 <- function(gr, col="blue", xlab="Index", ylab="Coverage", main="") {
    rg <- range(gr)
    start <- start(rg)
    end <- end(rg)
#     xWindow <- as.vector(window(cv, start, end))
#      smcov <- round(runmean(gr$cov, 21, endrule="constant"))
#      smcov <- round(runmed(gr$cov, 21, endrule="constant"))
    xWindow <- rep(gr$cov, width(gr))
    ## how to filter? kz?
    xWindow <- filterv(xWindow, filter=rep(1/20, 20))
    x <- start:end
    x <- x[!is.na(xWindow)]
    xWindow <- xWindow[!is.na(xWindow)]
#     xWindow[!is.na(xWindow)] <- 0
    xlim <- c(start, end)
    ylim <- c(0, max(xWindow))
#     ylim <- c(0, 30)
    plot(x = start, y = 0, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, type = "n", xaxt="n", frame.plot=FALSE)
    polygon(c(start, x, end), c(0, xWindow, 0), col = col, border=col) 
}
# Plot reads or fragments.
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, axis=TRUE, ...) {
    height <- 1
    if (is(xlim, "GRanges"))
        xlim <- c(min(start(xlim)), max(end(xlim)))
    bins <- disjointBins(IRanges(start(x), end(x) + 1))
    plot.new()
    plot.window(xlim, c(0, max(bins)*(height + sep)))
    ybottom <- bins * (sep + height) - height
    rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
    title(main)
    if(axis) axis(1) 
}

plot.nsome.track <- function(df, fragment.df, gr, peaks, peaks2) {
    mytheme <- theme_classic() %+replace%
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(face="bold",angle=90),
          plot.margin=unit(c(-0.1, 0.1, 0, 0.1), "cm"),
          panel.grid.major.x = element_line(colour = "lightgray"))

# rect(start(peaks)-0.5, 0.1, end(peaks)+0.5, 0.2, col = "darkgray")
#     rect <- data.frame(xmin=start(peaks2)-0.5, xmax=end(peaks2)+0.5, ymin=0.6, ymax=1.0)
    rect2 <- data.frame(xmin=start(peaks)-0.5, xmax=end(peaks)+0.5, ymin=0.1, ymax=0.5)
    m <- nrow(df)
    xlim <- gr
    xlim <- c(start(ranges(reduce(xlim, ignore.strand = TRUE))),
              end(ranges(reduce(xlim, ignore.strand = TRUE))))
    fig1 <- ggplot(df, aes(start, cov2)) + geom_ribbon(aes(ymin=rep(0, m), ymax=cov2, fill="green")) + scale_fill_manual(values=c("green" = rgb(49, 163, 84, maxColorValue=255)), guide=FALSE) +  mytheme
    fig1  <- fig1 + xlim(xlim) + ylab("Coverage") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    fig1 <- fig1 + geom_line(color=rgb(28, 125, 57, maxColorValue=255), lwd=0.5)
    fig1 <- fig1 + geom_rect(data=rect2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                             color="gray20",
                             alpha=0.5,
                             inherit.aes=FALSE)
#     fig1 <- fig1 + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
#                              color="gray20",
#                              alpha=0.5,
#                              inherit.aes=FALSE)
 
    wps.fig <- ggplot(df, aes(start, wps)) + geom_line(color="tomato2", size=1.1) + geom_hline(yintercept=0, color="gray40", linetype="dashed") + mytheme
    wps.fig <- wps.fig + xlim(xlim) + ylab("WPS")  + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

#     gc.fig <- ggplot(df, aes(start, gc)) + geom_line(color="gray", size=1.1) + mytheme
#     gc.fig <- gc.fig + xlim(xlim) + ylab("GC") + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
#     wps.fig2 <- ggplot(df, aes(start, wps2)) + geom_line(color="tomato2", size=1.1) + geom_hline(yintercept=0, color="gray40", linetype="dashed") + mytheme
#     wps.fig2 <- wps.fig2 + xlim(xlim) + ylab("WPS2")  + theme(axis.ticks = element_blank(), axis.text.x = element_blank())


    frags.fig <- ggplot(fragment.df, aes(xmin=xleft, xmax=xright, ymin=ybottom, ymax=ytop)) + geom_rect(color="black", fill="skyblue2") + mytheme
    frags.fig <- frags.fig + xlim(xlim) + labs(x="Position", y="Fragments")
    frags.fig <- frags.fig + theme( plot.margin=unit(c(-0.1, 0.1, 0.5, .1), "cm"), axis.title.x = element_text(face="bold"))

    plot_grid(fig1, wps.fig, gc.fig, frags.fig, ncol=1, align="v", rel_heights=c(2.2,1,1,2))
#     plot_grid(fig1, wps.fig, wps.fig2, frags.fig, ncol=1, align="v", rel_heights=c(2.2,1,1,2))
}



plotNucleosomes <- function(file, gr, windows=NULL, nucleosomes=NULL, refpeaks=NULL) {
    if(is.null(windows))
        windows <- get.windows(gr, size=120)
    if(is.null(windows$gc)) {
        gc <- GCcontent(Hsapiens, windows)
        windows$gc <- gc
    }
    bamflags <- scanBamFlag(isProperPair = TRUE, isDuplicate = FALSE, isUnmappedQuery=FALSE)
    bamparams <- ScanBamParam(flag = bamflags, which=gr, what=scanBamWhat())

#     reads <- get.reads(file=file, gr=gr, bamparams=bamparams)
    fragments <- get.fragments(file, gr, bamparams)
    fragments2 <- get.fragments2(fragments, gr)
#     wps <- get.wps(fragments, gr)
     wps <- get.wps(fragments, windows)
#      wps2 <- get.wps2(fragments, windows)
#     wps.kz <- kz(wps - median(wps), w=21)
    wps.kz <- kz(wps - runmed(wps, 201, "constant"), w=21)
#     wps.kz2 <- kz(wps2 - runmed(wps2, 201, "constant"), w=21)

    windows$wps <- wps
    windows$wps.kz <- wps.kz
#     windows$wps.kz2 <- wps.kz2

#     cov <- get.disjoin(reads, gr)
    cov <- get.disjoin(fragments, gr)
    gr2 <- integrate.summaries(windows, cov)
#     df = data.frame(wps=gr2$wps.kz, wps2=gr2$wps.kz2, cov=gr2$cov, gc=gr2$gc,
#                     start=start(gr2), end=end(gr2))
    df = data.frame(wps=gr2$wps.kz, cov=gr2$cov, gc=gr2$gc,
                    start=start(gr2), end=end(gr2))
    df <- df[!is.na(df$cov), ]
    ## filter will feel if remove too many NAs?
    df$cov2 <- as.numeric(filterv(df$cov, filter=rep(1/20, 20)))

    frags.df <- read.piles(fragments2)

#     if(is.null(peaks)) peaks <- get.nsome(gr2, gr)
#    peaks2 <- get.nsome(gr2, gr)
    if(is.null(nucleosomes)) nucleosomes <- nsome.peak(gr2, gr)
    plot.nsome.track(df, frags.df, gr, nucleosomes, refpeaks)
}

