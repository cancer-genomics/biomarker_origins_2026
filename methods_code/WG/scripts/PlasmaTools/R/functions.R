### TODO:
#

#' @export
frag.density <- function(fragments, x="w", ...,  groups=NULL) {
    d.list <- function(x, ...) {
        d <- density(x, ...)
        list(x=d$x, y=d$y)
    }

    ## Add grouping arguments (gc etc) to pass into keyby
    dens <- fragments[,c(d.list(get(x), ...)), keyby=groups]
    dens[,ynorm:=y/max(y)][]
}

#' @export
frag.stats <- function(fragments, cutoff=150, groups=NULL) {
    Mode <- function(x){
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    fragments[,.(mode=Mode(w), median=as.integer(median(w)),
                 mean=mean(w), size_iqr=iqr(w),
                 nfrags=.N, mononucs=sum(w>=100 & w<=250),
                 multinucs = sum(w>=250), ultrashort = sum(w<100),
                 slratio = sum(w<=cutoff)/sum(w > cutoff & w <= 250),
                 meangc=mean(gc), mediangc=median(gc),
                 gc_iqr=iqr(gc), chrMrep=-log(sum(chr=="chrM")/.N)),
            keyby=groups][]
}

#' @export
gcCorrectTarget <- function(fragments, ref, bychr=TRUE, rescale=TRUE){
    fragments[, gc := round(gc, 2)]
    if(bychr) {
        DT.gc <- fragments[,.(n=.N), by=.(gc, chr)]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc, chr)]
    } else {
        DT.gc <- fragments[,.(n=.N), by=gc]
        DT.gc <- DT.gc[gc >= .20 & gc <= .80]
        DT.gc <- DT.gc[order(gc)]
    }
# setkey(mediandt, gc, seqnames)

    if(bychr) {
        setkey(DT.gc, gc, chr)
        setkey(ref, gc, chr)
    } else {
        setkey(DT.gc, gc)
        setkey(ref, gc)
    }
#     DT.gc <- DT.gc[ref][order(chr, gc)]
    DT.gc <- DT.gc[ref]
    DT.gc[,w:=target/n]
    if(bychr) {
        fragments[DT.gc, on= .(chr, gc), weight := i.w]
    }
    else fragments[DT.gc, on= .(gc), weight := i.w]
    fragments <- fragments[!is.na(weight)]
    if(rescale) fragments[,weight := weight * .N/sum(weight)]
    fragments[]
}

#' @export
gcCorrectTargetSize <- function(fragments, ref){
    fragments[, gc := round(gc, 2)]
    DT.gc <- fragments[,.(n=.N), by=.(size, gc)]
    DT.gc <- DT.gc[gc >= .20 & gc <= .80]
    DT.gc <- DT.gc[order(size, gc)]

    setkey(DT.gc, size, gc)
    setkey(ref, size, gc)
    DT.gc <- DT.gc[ref]
    DT.gc[,wt:=target/n]
    fragments[DT.gc, on= .(size, gc), weight := wt]
    fragments <- fragments[!is.na(weight)]
    fragments[]
}

getChrM <- function(bed) {
    dt <- fread(cmd=paste("awk '$1 == \"chrM\" { print $0 }'", bed))
    setnames(dt, paste0("V", 1:4), c("chr", "start", "end", "mapq"))
    dt
}

## Get ChrM Representation

#' @export
binFrags <- function(fragments, bins, cutoff=150,
                     chromosomes=paste0("chr",c(1:22, "X"))) {
    setkey(bins, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)
    bins2 <- fragbins[,.(arm=unique(arm), gc=gc[1], map=map[1],
                         short = sum(w >= 100 & w <= cutoff ),
                         long = sum(w > cutoff & w <= 250),
                         short.cor = sum(weight[w >= 100 & w <= cutoff]),
                         long.cor = sum(weight[w > cutoff & w <= 250]),
                         ultrashort = sum(w < 100),
                         ultrashort.cor = sum(weight[w < 100]),
                         multinucs = sum(w > 250),
                         multinucs.cor = sum(weight[w > 250]),
                         mediansize = as.integer(median(w)),
                         frag.gc = mean(fraggc)),
            by=.(chr, start, end)]

    setkey(bins2, chr, start, end)
    bins2 <- bins2[bins]
    bins2 <- bins2[is.na(i.gc), which(grepl("short|long|multi", colnames(bins2))):=0]
    bins2[,`:=`(gc=i.gc, map=i.map, arm=i.arm)]
    bins2[,which(grepl("^i.", colnames(bins2))):=NULL]
    bins2[, bin:=1:.N]
    setcolorder(bins2, c("chr", "start", "end", "bin"))
    bins2[]
}


#' @export
binKLDiv <- function(fragments, bins, chromosomes=paste0("chr",c(1:22, "X"))) {
    setkey(bins, chr, start, end)
    fragbins <- foverlaps(fragments[chr %in% chromosomes],
                          bins, type="within", nomatch=NULL)

    fragbins <- fragbins[w>=100 & w <= 250][,.(n=.N), keyby=.(chr, bin, w)]
    fragbins[,p:=n/sum(n), keyby=.(bin)]

    v <- dcast(fragbins, chr + bin ~ w, value.var=c( "p"))
    for (j in names(v))
        set(v,which(is.na(v[[j]])),j,0)

    v[,chr:=factor(chr, chromosomes)]
    vv <- v[, KLmat(.SD), .SDcols=!"bin", keyby=chr][,bin:=1:.N]
#     vv <- v[, ks.matrix(.SD), .SDcols=!"bin", keyby=chr][,bin:=1:.N]
#     vv <- v[,eigen:=.standardize(eigen), by=.(chr)]
    vv[]
#     v[]
}

## Function will append corrected measure to bins object
#' @export
removePCs <- function(bins, measure, pcs=1) {
    Y <- dcast(bins, id ~ bin , value.var=measure)
    svd <- svd(Y[,-1])
    subspace <- svd$u[,c(1,pcs),drop=FALSE] %*%
        diag(svd$d[c(1,pcs)]) %*%
        t(svd$v[,c(1,pcs),drop=FALSE])
    Ystar <- Y[,-1] - subspace
    Ystar[,id:=Y$id]
    Ystar.long <- melt(Ystar, id.vars="id",
                       value.name="corrected",
                       variable.name="bin")
    Ystar.long[,bin:=as.integer(bin)]
    setkey(Ystar.long, id, bin)
    setkey(bins, id, bin)
    bins <- bins[Ystar.long]
    bins[]
}

### Heatmap functions for CNN

#' @export
KLDmat <- function(v) {
    .kl.div <- function(x,y) {
        p <- (x+1)/sum(x+1)
        q <- (y+1)/sum(y+1)
        ## set 0 * log(0/0) = 0

#         ts <- sum(p * log(p/q), na.rm=TRUE)
        ## Symmetrized KL
        ts <- sum(p * log(p/q), na.rm=TRUE) + sum(q * log(q/p), na.rm=TRUE)
        ts
    }

    mat<-apply(v, 1, function(x) {
                    apply(v, 1, function(y) .kl.div(y, x))
           })

}

## KS
# TODO: need better way to handle 0 counts

#' @export
KSmat <- function(v) {
    .ks.distance <- function(x, y) {
        nx <- sum(x)
        ny <- sum(y)
        x1 <- cumsum(x)/nx
        y1 <- cumsum(y)/ny
        d <- max(abs(x1 - y1))
        ts <- d/sqrt((nx+ny)/(nx*ny))
        ts
    }

    mat<-apply(v, 1, function(x) {
               apply(v, 1, function(y) .ks.distance(y, x))
           })
}

## Anderson Darling
#' @export
ADmat <- function(v) {
    .ad <- function(x,y) {
        nx <- sum(x+1)
        ny <- sum(y+1)
        N <- nx+ny
        x1 <- (cumsum(x)+1)/nx
        y1 <- (cumsum(y)+1)/ny
        h <- (cumsum(x+y)+2)/(nx+ny)
        d <- (nx*ny/N) * sum((x1 - y1)^2/(h * (1-h)))
        d
    }

    mat<-apply(v, 1, function(x) {
                    apply(v, 1, function(y) .ad(y, x))
           })

}

## 
