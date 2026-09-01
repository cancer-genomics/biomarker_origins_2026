#TODO: 1) get rid of ichor dependency, 2) introduce grouping options in plot to allow for many density profiles

#have something that makes wig files
#then tell user to go to dump to ichor
#writing our own segmentation program is not high ROI

#turn coverage across bins into log r ratio across bins
#direction = 1 if each sample is a column, or 2 if each sample is a row in df
cov2logr <- function(df, direction=2){
    
    #healthy median coverage per 500 kb bin taken from 140 healthy plasma samples
    #from DELFI paper
    #needs to be able to generalize to any size bin
    healthy.bin.medians <- readRDS("../healthymedians.rda")
    
    df <- as_tibble(sweep(df, direction, healthy.bin.medians, "/"))
    df <- log2(df)
    
    return(df)
    
    
}


#the dataframes for the following functions are the direct result of collecting ichor outputs
#if ichor dependency still exist, we can incorporate a preprocessing of ichor results step
#columns are chr, start, end, pos (start + end /2), sample.id, logr, corr.copy.number)

#call outliers on log r bin level data
#column names are sample and logr and pos (noting a bin number of location in genome)
callOutliers <- function(df, source="all", long=TRUE){
    
    #if you have many types of samples in df and want to pick a certain type to call outliers on
    if (source != "all"){
        df <- df %>% filter(source==source)
    }
    
    if (long){
        call.outliers <- spread(df, sample, logr)
    }
    
    call.outliers <- call.outliers %>% mutate_at(vars(-pos), funs( (. - median(.)) ))
    call.outliers <- call.outliers %>% mutate_at(vars(-pos),funs( (.-mean(.))/sd(.)) )

    call.outliers <- abs(call.outliers)
    
    position <- which(call.outliers=="pos")
    call.outliers[-position] <- +(call.outliers[-positino] > 3)
    call.outliers$counts <- rowSums(call.outliers[-position])
    excludeBins <- call.outliers$pos[call.outliers$counts >= 2]
    
    #TRUE means bin is to be excluded
    df <- df %>% mutate(logr = ifelse(pos %in% excludeBins, NA , logr ), excludeBins=ifelse(pos %in% excludeBins, TRUE, FALSE))
    
    return(df)
}




#pick a sample
plotPloidy <- function(df.ploidy, df.dens, outliers=TRUE, density=TRUE, sample){
    
    if (outliers){
        df.ploidy <- callOutliers(df.ploidy)
    }
    if (density){
        df.dens <- frag.density(df.dens)
    }
    
    df.ploidy <- df.ploidy %>% filter(sample.id=sample)
    
    mytheme_ploidy <- theme(axis.title.x=element_blank(), panel.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    mytheme_dens <- theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.background = element_rect(fill = "transparent"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x=element_blank())
    
     p <- ggplot(df.ploidy, aes(x=pos,y=logr,color=color)) + geom_point(size=1) + ylim(-2,2) + geom_vline(xintercept=as.numeric(df.ploidy$pos[cumsum(rle(df.ploidy$chr)$lengths)]), linetype="dotted",size=0.25) + ylab("Copy Number (log r ratio)") + geom_segment(aes(x = df.ploidy$pos[1], y = df.ploidy$draw.line[1], xend = tail(df.ploidy$pos,n=1), yend = df.ploidy$draw.line[1]),color="grey") + mytheme_ploidy
     
     k <- ggplot(l, aes(x=x, y=y)) + geom_area() +  xlim(min(df.dens$x), 325) + mytheme_dens

     grid.arrange(p,k,nrow=1, top=paste(unique(t$annot), collapse="\n"), widths=2:1)
     
}
