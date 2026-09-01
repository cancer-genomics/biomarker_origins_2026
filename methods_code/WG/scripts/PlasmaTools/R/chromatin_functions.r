## eigen decompose matrix of measurement values across each grouping variable.
## Needs work to be made more generalized (ie grouping variables beyond just
## source and chr).
compartments <- function(tib, val, ids, ...) {
    eigenlist <- tib %>% group_by(.dots=lazyeval::lazy_dots(...)) %>%
    select(ids, val) %>% ungroup() %>%
    split(list(.$source, .$chr)) %>%
    map(function(x) {
        x %>% pivot_wider(id, names_from=bin, values_from=val) %>%
        select(-id) %>%
        .svd.qr(ss=FALSE)  %>%
        .standardize()})
    
    names=map_dfr(strsplit(names(eigenlist), "\\."),
                  function(x) tibble(source=x[1], chr=x[2]))
    times=sapply(eigenlist, length)

    tibble(source=rep(names$source, times),
           chr=rep(names$chr, times),
           eigen=unlist(eigenlist))
}

compartments.ss <- function(tib, val, ids, ...) {
    eigenlist <- tib %>% group_by(.dots=lazyeval::lazy_dots(...)) %>%
    select(ids, val) %>% ungroup() %>%
    split(list(.$id, .$chr)) %>%
    map(function(x) {
        x %>% pivot_wider(id, names_from=bin, values_from=val) %>%
        select(-id) %>%
        .svd.qr(TRUE)  %>%
        .standardize()})
    
    names=map_dfr(strsplit(names(eigenlist), "\\."),
                  function(x) tibble(id=x[1], chr=x[2]))
    times=sapply(eigenlist, length)

    tibble(id=rep(names$id, times),
           chr=rep(names$chr, times),
           eigen=unlist(eigenlist))
}

### Fast PCA using QR decomposition (inspired by minfi)
.svd.qr = function(A, l=2, ss=FALSE) {
    if(ss) { 
#         A <- as.matrix(A)
#         v <- (A - mean(A))/sd(A)
        mat <- crossprod(t(A), t(1/A)) #         mat <- crossprod(t(A))
        mat[mat > 1] <- 1/mat[mat>1]
        corr <- cor(mat, method="pearson")
#         corr <- mat
    }
    else corr <- cor(A, method="pearson")
#     corr <- cor(A, method="spearman")

    center <- rowMeans2(corr, na.rm = TRUE)
    corr <- sweep(corr, 1L, center, check.margin = FALSE) 
#     scale <- rowMeans2(corr, na.rm = TRUE)
#     corr <- sweep(corr, 1L, scale, FUN="/", check.margin = FALSE) 
#     l <- 5
    n <- ncol(corr)
    m <- nrow(corr)

    ## Fast implementation (Projection onto lower dimensional subspace)
    if(TRUE) {
    G <- matrix(rnorm(n * l, 0, 1), nrow = n, ncol = l)
    h.list <- vector("list", 2L)
#     h.list <- vector("list", i + 1L)
    h.list[[1]] <- corr %*% G
    for (j in seq(2, 1 + 1L)) {
        h.list[[j]] <- corr %*% (crossprod(corr, h.list[[j - 1L]]))
    }
#     h.list[[2]] <- corr %*% (crossprod(corr, h.list[[1L]]))
    H <- do.call("cbind", h.list) # n x [(1+1)l] matrix
    # QR algorithm
    Q <- qr.Q(qr(H, 0))
    TT <- crossprod(Q, corr)
    svd <- svd(TT)
    u <- Q %*% svd$u[,1, drop = TRUE]
    }
    ## Slow
#     svd <- svd(corr)
#     u <- svd$u[, 1, drop = TRUE]
    u
}

.standardize <- function(u) {
    u <- u - median(u)
    u <- u/sqrt(sum(u^2))
    u
}

##  For smoothing out vectors for nicer plotting of AB compartments.
##  Credit goes to JP Fortin.
.meanSmoother <- function(x, k=1, iter=2, na.rm=TRUE){
    meanSmoother.internal <- function(x, k=1, na.rm=TRUE){
        n <- length(x)
        y <- rep(NA,n)
        
        window.mean <- function(x, j, k, na.rm=na.rm){
            if (k>=1){
                return(mean(x[(j-(k+1)):(j+k)], na.rm=na.rm))
            } else {
                return(x[j])
            }    
        }
        
        for (i in (k+1):(n-k)){
            y[i] <- window.mean(x,i,k, na.rm)
        }
        for (i in 1:k){
            y[i] <- window.mean(x,i,i-1, na.rm)
        }
        for (i in (n-k+1):n){
            y[i] <- window.mean(x,i,n-i,na.rm)
        }
        y
    }
    
    for (i in 1:iter){
        x <- meanSmoother.internal(x, k=k, na.rm=na.rm)
    }
    x
}

### Plot AB Compartments (similar to JP and KH paper, but using ggplot)
### geom_ribbon instead of bar? large bins & small chromosomes look bad
geom_ab <- function(tib, chromosome) { 
    tib  <- tib %>% filter(chr==chromosome) #%>%
#         mutate(eigen = .meanSmoother(eigen),
#                color=ifelse(eigen < 0, "gray50", "brickred"))
p <- ggplot(tib, aes(x=bin, y=eigen2, color=color, fill=color)) +
    geom_bar(stat="Identity", width=1)
p <- p + facet_grid(`source2`~.,  scale="free")
p <- p +  scale_x_continuous(expand = c(0, 0)) + xlab(chromosome)
p
}

#create heatmap and cluster eigenvectors of samples from one or many chr
heatmap_ab <- function(tib, chr=FALSE){
    
    if (chr[1]){
        tib <- tib %>% filter(chr==chromosome)
    }
    
    df <- tib %>% select(bin, eigen2, source2) %>% spread(bin, eigen2)
    df <- tibble::column_to_rownames(df, source2)
        
    clust <- as.dendrogram(hclust(dist(df, method = "euclidean"), method = "ward.D2" ))
    
    heatmap(data.matrix(clust), scale="none", Rowv=h_clust, Colv = NA, col=brewer.pal(11,"RdBu"), cexRow=0.75, labCol=FALSE, margins=c(1,7), main=chromosome)
    
    legend("topleft",legend=c("min", "0", "max"),fill=brewer.pal(3,"RdBu"),cex=0.6)
    
}


.fragsize.ks <- function(x, y) {
    xx <- rep(100:220, x)
    yy <- rep(100:220, y)
    suppressWarnings(ks.test(xx, yy)$p.value)
    }

ks.matrix <- function(m) {

    d.matrix<-apply(m, 1, function(x) {
                    apply(m, 1, function(y) .fragsize.ks(y, x))
           })
    u <- c(.svd.qr(d.matrix))
    list(u=u)
}
