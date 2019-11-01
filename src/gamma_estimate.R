library(scales)
library(flexmix)

dat <- my.data$segments_raw
dat$seg.mean <- rowSums(dat[,c('nAraw', 'nBraw')])
dat$width <- with(dat, endpos-startpos)
dat$min.width <- floor(dat$width / quantile(dat$width, 0.05))

##################
#### Plotting ####
#### nABraw ####
## Adjusted
set.seed(1234)
par(mfrow=c(2,1))
AB.comp <- plotHistComponents(rep(dat$seg.mean, dat$min.width), col='black', lwd=2, 
                             add=F, reduce.dat=TRUE, reduce.to=100000, max.comp=6)
ab <- getCIdist(list(parameters(AB.comp)), 
                list(table(clusters(AB.comp))), multiplier=1)
.addLine(ab$adj, model.col='blue', lwd=3)


ABraw.comp <- plotHistComponents(dat$seg.mean, col='black', lwd=2, add=F, 
                                 max.comp=6, minprior=0.001, niter=1000)
ab.raw <- getCIdist(list(parameters(ABraw.comp)), 
                list(table(clusters(ABraw.comp))), multiplier=1)
.addLine(ab.raw$adj, model.col='blue', lwd=3)


#### nAraw and nBraw ####
set.seed(1234)
par(mfrow=c(2,1))
A.comp <- plotHistComponents(rep(dat$nAraw, dat$min.width), col='red', lwd=2, 
                             add=F, reduce.dat=TRUE, reduce.to=10000, max.comp=6)
B.comp <- plotHistComponents(rep(dat$nBraw, dat$min.width), col='black', lwd=2, 
                             add=T, reduce.dat=TRUE, reduce.to=10000, max.comp=6)
ab <- getCIdist(lapply(list(A.comp, B.comp), parameters), 
                lapply(list(A.comp, B.comp), function(i) table(clusters(i))), multiplier=1)
.addLine(ab$adj, model.col='blue', lwd=3)
# mean(ab$adj['mag.cv',]) / mean(ab$raw['mag.cv',])

## Raw
set.seed(1234)
par(mfrow=c(2,1))
Amin.comp <- plotHistComponents(dat$nAraw, col='red', lwd=2, add=F, max.comp=6)
Bmin.comp <- plotHistComponents(dat$nBraw, col='black', lwd=2, add=T, max.comp=6)
ab <- getCIdist(parameters(Amin.comp), parameters(Bmin.comp), 
                table(clusters(Amin.comp)), table(clusters(Bmin.comp)), multiplier=2)
.addLine(ab$adj, model.col='green', lwd=3)


#### Fit Theoretical Data ####
psi <- 2.1
ab <- getCIdist(list(parameters(AB.comp)), 
                list(table(clusters(AB.comp))), multiplier=1)
comp.param <- ab$adj 

A=comp.param['magnitude', order(comp.param['mean',]),drop=F]
A = round(A/max(A),2)
psi=2.1 * ncol(A)
n.states <- ncol(A)-1


  
A <- matrix(c(0.02,0.34,0.44,1,0.32), ncol=5)
showEqn(A, b)


###################
#### Functions ####

getCIdist <- function(ab.l, rep.l, ...){
  .eucldist <- function(x,y){
    as.numeric(unlist(sapply(x, function(x0) x0-y)))
  }
  
  ab <- do.call(cbind, ab.l)
  magnitude <- unlist(lapply(rep.l, function(r0) r0 / sum(r0)))
  ab <- rbind(ab, 'magnitude'=magnitude)
  ab <- rbind(ab, "coef.var"=round(ab[2,]/ab[1,],3))
  ab <- rbind(ab, "mag.cv"=round(ab[3,] / ab[4,],3))
  
  ## Calculate the 95 CI distance between gaussian distributions
  y <- apply(ab, 2, function(i, multiplier=1){
    j.all <- apply(ab, 2, function(j){
      d <- .eucldist(rnorm(n = 1000, mean = i[1], sd=i[2]),
                     rnorm(n = 1000, mean = j[1], sd=j[2]))
      
      d.ci <- c((mean(d) - (sd(d)*multiplier)), (mean(d) + (sd(d)*multiplier)))
      d.ci
    })
    return(data.frame(j.all))
  })
  
  ## Remove based on lowest magnitude, highest CV
  rm.idx <- removeOverlappingComponents(ab, y)
  list("adj"=ab[,-rm.idx], "raw"=ab)
}

## Check if the CI spans 0
genOverlapMat <- function(y, rm.idx){
  .spans0 <- function(x){
    max(x) < 0  | min(x) > 0
  }
  
  ov.mat <- sapply(y, function(i) apply(i[-rm.idx], 2, .spans0))
  diag(ov.mat) <- T
  ov.mat
}

removeOverlappingComponents <- function(ab, y){
  overlaps.exist<-TRUE
  rm.idx <- 100
  mag.cv <- ab['magnitude',] / ab['coef.var',]
  mag.cv.ord <- order(mag.cv)
  idx <- 1
  
  ## NEEEEEED to add checks for while loops!
  while(overlaps.exist){
    # check if any distributions overlap
    ov <- genOverlapMat(y[-rm.idx], rm.idx)
    if(all(ov)){
      # If all distributions are unique, exit
      overlaps.exist <- FALSE
    } else {
      # Cycle through distributions based on [magnitude/CV] order
      ov.same <- TRUE
      while(ov.same & (idx <= length(mag.cv))){
        min.idx <- mag.cv.ord[idx]
        new.rm.idx <- c(rm.idx, min.idx)
        # Check if there's a difference with the potential component removal
        ov.same <- (sum(!ov) == sum(!genOverlapMat(y[-new.rm.idx], new.rm.idx)))
        idx <- idx + 1
      }
      
      # IF there's a difference, cement in the new removed component
      print(mag.cv[min.idx])
      mag.cv[min.idx] <- NA
      rm.idx <- c(rm.idx, min.idx)
    }
  }
  rm.idx[-1]
}


plotHistComponents <- function(dat, col, breaks=seq(0,10, by=0.01),
                               lwd=2, add=F,...){
  hist(dat, breaks=breaks, xlim=c(0,4), col=alpha(col, 0.3), border = NA, 
       freq = FALSE, las=1, xlab="Copy-number", add=add)
  fit.bic.dat <- fitMixedModel(model.dat=dat, ...)
  addMixComp(fit.bic.dat, model.col=col, lwd=lwd)
  fit.bic.dat
}

addMixComp <- function(fit.model, ...){
  fit.param <- parameters(fit.model)
  lam <- table(clusters(fit.model))
  fit.param <- rbind(fit.param, lam)
  
  .addLine(fit.param, ...)
}

.addLine <- function(fit.param, ...){
  apply(fit.param, 2, function(comp, model.col='black', ...){
    lines(x=breaks, y=(comp[3]/sum(fit.param[3,]) * dnorm(breaks, mean=comp[1], sd=comp[2])), 
          col=model.col, ...)
  }, ...)
}


fitMixedModel <- function(model.dat, min.comp=1, max.comp=10, 
                          model.selection='BIC', reduce.dat=FALSE, 
                          reduce.to=1000, minprior=0.01, niter=200){
  if(reduce.dat){
    model.dat <- sample(model.dat, size = reduce.to, replace=F)
  }
  
  control<-new("FLXcontrol")
  control@minprior <- minprior
  control@iter.max <- niter
  fit<-flexmix::stepFlexmix(model.dat ~ 1, model=flexmix::FLXMCnorm1(),
                            k=min.comp:max.comp, control=control)
  fit.bic<-flexmix::getModel(fit,which=model.selection)
  fit.bic
}
