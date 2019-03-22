library(grid)
library(ggplot2)

cb_palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gwas.manhattan.chr <- function(d, genomewideline=-log10(1e-4)) {
  d <- na.omit(d)
  d <- d[d$P>0&d$P<=1,]
  d$logp <- -log10(d$P)

  uchr <- unique(d$CHR)
  stopifnot(length(uchr)==1)
  
  d$pos <- d$BP

  p <- ggplot(d,aes(pos,logp)) +
       geom_point(size=rel(0.7)) +
       geom_hline(yintercept=genomewideline,linetype=2,size=rel(0.25)) +
	   scale_x_continuous(expand=c(0.01,0.01)) + 
       scale_y_continuous(expand=c(0.01,0.01)) +
       scale_colour_manual(values=chr.color) + 
       xlab(paste("Chromosome",uchr)) + 
       ylab(expression(-log[10]*italic(P))) + 
       theme_bw() + 
       theme(legend.position="none",panel.grid=element_blank())

  #topsnps <- d[d$logp >= genomewideline,]
  #if (nrow(topsnps) > 0) {
  #  p <- p + geom_point(data=topsnps,size=rel(1)) # + geom_text(data=topsnps,aes(label=SNP),vjust=1,angle=0,size=rel(3))
  #}
  
  p
}
##############################################################################################################################
gwas.manhattan <- function(d, genomewideline=-log10(1e-4)) {
  d <- na.omit(d)
  d <- d[d$P>0&d$P<=1,]

  uchr <- unique(d$CHR)
  nchr <- length(uchr)

  if (nchr == 1) return(gwas.manhattan.chr(d,genomewideline,genomewideline2))

  chr.color <- rep(c("gray10","gray50"),nchr)
  chr.color <- rep(cb_palette[c(4,6)],nchr)
  
  d$logp <- -log10(d$P)

  d$index <- NA
  ind <- 0
  for (i in uchr) {
    ind <- ind + 1
    d[d$CHR==i,]$index <- ind
  }

  d$pos <- NA
  ticks <- rep(NA,nchr+1)
  ticks[1] <- 0
  for (i in 1:max(d$index)) {
    d[d$index==i,]$pos <- (d[d$index==i,]$BP-d[d$index==i,]$BP[1]) + ticks[i]
    ticks[i+1] <- max(d[d$index==i,]$pos) + 1
  }
  ticks <- ticks[-1] - diff(ticks)/2
  tick.labs <- uchr

  p <- ggplot(d,aes(pos,logp,color=factor(CHR))) + 
       geom_point(size=rel(0.05)) + 
       geom_hline(yintercept=genomewideline,linetype=2,size=rel(0.3),color="black") + 
	   geom_hline(yintercept=-log10(1e-9),linetype=2,size=rel(0.25),color="gray50") +
	   scale_x_continuous(breaks=ticks,labels=tick.labs,expand=c(0.01,0)) + 
	   scale_y_continuous(expand=c(0.01,0)) + 
       scale_colour_manual(values=chr.color) + 
       xlab("Chromosome") + 
       ylab(expression(-log[10]*italic(P))) + 
       theme_bw() + 
       theme(legend.position="none", panel.grid=element_blank())

  #topsnps <- d[d$logp >= genomewideline,]
  #if (nrow(topsnps) > 0) {
  #  p <- p + geom_point(data=topsnps,size=rel(1)) # + geom_text(data=topsnps,aes(label=SNP),vjust=1,angle=0,size=rel(3))
  #}
  
  p
}

gwas.qq <- function(ps) {
  ps <- ps[!is.na(ps)]
  ps <- ps[ps>0 & ps<1]
  ps <- ps[order(ps,decreasing=F)]

  d <- data.frame(o=-log10(ps), e=-log10(ppoints(length(ps))))

  p <- ggplot(d,aes(e,o)) +
       geom_point(size=rel(0.7),color=cb_palette[6]) + 
       geom_abline(intercept=0,slope=1,color=cb_palette[1]) + 
       scale_x_continuous(expand=c(0,0)) + 
       scale_y_continuous(expand=c(0,0)) + 
       xlab(expression(Expected~~-log[10]*italic(P))) + 
       ylab(expression(Observed~~-log[10]*italic(P))) + 
       theme_bw() + 
       theme(panel.grid=element_blank())
  p
}

gwas.manqq <- function(d, f="manqq.tif", w=24, h=6) {
  man <- gwas.manhattan(d)
  qq <- gwas.qq(d$P)

  tiff(f,width=w,height=h,units="cm",compression="lzw",res=600)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,4)))

  print(man,vp=viewport(layout.pos.row=1,layout.pos.col=1:3))
  print(qq,vp=viewport(layout.pos.row=1,layout.pos.col=4))

  dev.off()

  return()
}

gwas.plot.all <- function() {
  files <- dir("data")
  for (f in files) {
    d <- read.table(paste("data/",f,sep=""), header=T)
    gwas.manqq(d,paste(f,"_manqq.tif",sep=""))
  }
}
