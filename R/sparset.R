sparset <- function(gridEnvs, gridIndsPerEnv, lb=0, ub=Inf,
                    nCrosses=30, nProgeny=10, 
                    G=NULL, h2=0.5, p=0.1, 
                    nItersMacs=10, nItersTrt=5, trts=NULL,
                    verbose=TRUE){
  if(is.null(trts)){
    if(missing(gridEnvs)){stop("Please provide a grid of number of environments to test", call. = FALSE)}
    if(missing(gridEnvs)){stop("Please provide a grid of number of individuals per environment to test", call. = FALSE)}
  }
  
  result <- list() # to store results
  counter <- 1
  
  for( iMacs in 1:nItersMacs){
    
    founderPop = runMacs2(
      nInd=50,  nChr = 10,  Ne = 30,
      segSites = rep(1000,10), # from runMacs how many segSites we leave (this affects the speed performance due to memory issues). Normally nSNPs+nQTLs is a good value to use.
      genLen = rep(1, 10), # 1 Morgan contains all recombinations (normally 1 Morgan is equivalent to 1 million bp in humans, species-dependent). A specie with 1.4 Morgans is equal to a length of 143 cM
      # bp = c(1e+08,1e+06), #  1e+06 = 1 x 10^6 = 10^6 is a million bp, as many values as chromosomes
      # recombination rate = genLen/bp = 1 Morgan / 1e+08 = 1e-08 that is 100 and 140 cM along 100 million and 1 million chromosomes
      mutRate = 2.5e-08, # 2.5 point mutations per 100 million bp. In a 1 billion bp would be 25 mutations. Natural units are: per base pair per cell division, per gene per generation, or per genome per generation
      histNe = c(500, 1500, 6000, 12000, 1e+05), # effective population size in previous generations
      histGen = c(100, 1000, 10000, 1e+05, 1e+06), # number of generations ago for effective population sizes given in histNe
      inbred = TRUE, # outbred or inbred
      split = NULL,  ploidy = 2L,  returnCommand = FALSE,  nThreads = NULL
    )
    if(is.null(G)){
      nEnvs=500#round(min(c(300,max(gridEnvs)*10))) # number of farms to comprise the TPE
      G = simGECorMat0(nEnv =nEnvs,nMegaEnv = 1,mu=0.5, v=0.2); # we assume 1 cluster, mu.rg=0.3, sd.rg=0.54,  
      # hist(G[lower.tri(G)])
      # lattice::levelplot(G)
    }else{
      nEnvs <- ncol(G)
    }
    
    SP = SimParam$new(founderPop)
    
    suppressWarnings(
      SP$addTraitAD(
        nQtlPerChr=40, # number of QTLs per chromosome. Can be a single value or nChr values.
        mean = abs(rnorm(nEnvs, 11, 4)), # a vector of desired mean genetic values for one or more traits
        # mean = rep(0,nEnvs), # a vector of desired mean genetic values for one or more traits
        var = rep(1,nEnvs), # a vector of desired genetic variances for one or more traits
        meanDD = rep(0,nEnvs), # mean dominance degree
        corA = G, # a matrix of correlations between additive effects
        useVarA = TRUE, gamma = FALSE, # should a gamma distribution be used instead of normal
        shape = 1, # the shape parameter for the gamma distribution
        force = FALSE
      ), classes = "warning"
    )
    
    SP$addSnpChip(
      nSnpPerChr=150, # number of SNPs per chromosome. Can be a single value or nChr values.
      minSnpFreq = NULL, # minimum allowable frequency for SNP loci. If NULL, no minimum frequency is used.
      refPop = NULL # reference population for calculating SNP frequency. If NULL, the founder population is used.
    )
    
    SP$setTrackPed(
      isTrackPed=TRUE, # should pedigree tracking be on.
      force = FALSE
    )
    
    pop = newPop(founderPop, simParam=SP) # 20 lines to start
    pop = randCross2(pop,pop,nCrosses,1,simParam=SP) # pick 30 crosses and obtain F1 hybrids
    pop <- self(pop,nProgeny,simParam=SP) # self each F1 and produce 10 F2 per cross
    pop = makeDH(pop=pop,simParam=SP) # 300 DH lines
    ###########$$$$$$$$$$$$$$$$$$$$
    ###########$$$$$$$$$$$$$$$$$$$$
    pop = setPheno(pop,H2=rep(h2,nEnvs),reps=1,simParam=SP) # produce phenotypes at low h2 independently of the trial size
    ###########$$$$$$$$$$$$$$$$$$$$
    ###########$$$$$$$$$$$$$$$$$$$$
    
    # calculate true GV
    trueGV = apply(pop@gv,1,mean)
    names(trueGV) <- pop@id
    # define the treatments
    if(is.null(trts)){
      trts <- expand.grid(gridEnvs, gridIndsPerEnv) # no.farms by no.inds 
      colnames(trts) <- c("nFarms","nIndsPerFarm")
      trts$total <- trts$nFarms*trts$nIndsPerFarm # total number of plots needed
      trts <- trts[ order(trts[,1], trts[,2]), ]
      trts$availableInds <- nInd(pop)
      # add balanced designs
      trts2 <- trts
      trts2$availableInds <- trts2$nIndsPerFarm
      trts <- rbind(trts,trts2)
    }else{
      if( ! all(colnames(trts) %in% c("nFarms","nIndsPerFarm","total", "availableInds")) ){
        stop("Please make sure that your 'trts' argument has column names 'nFarms', 'nIndsPerFarm', 'total', 'availableInds'   ", call. = FALSE)
      }
      gridEnvs <- trts[,1]
      gridIndsPerEnv <- trts[,2]
    }
    
    # warnings or controls
    if(max(gridEnvs) > nEnvs){stop("The gridEnvs cannot include more environments than defined in G. Please correct", call. = FALSE)}
    if(max(gridIndsPerEnv) > (nProgeny*nCrosses) ){stop("The gridIndsPerEnv cannot include more individuals than nCrosses x nProgeny. Please correct", call. = FALSE)}
    
    
    trts <- trts[which(trts$total <= ub),]
    trts <- trts[which(trts$total >= lb),]
    
    ##############################
    ##############################
    ## bigger fields less reps
    ##############################
    ##############################
    
    # nItersTrt = 5 # only 5 reps since this takes a long time
    sampledEnvsList <- apply(trts,1,function(x){sample(1:ncol(pop@gv), x[1])})
    for(iTrt in 1:nrow(trts)){ # iTrt = 1
      
      for(iRep in 1:nItersTrt){ # iRep=1
        #################################################################
        #################################################################
        # EXPERIMENTAL DESIGN
        # 1) extract the entry list and add the number of seed packets availableInds
        
        sampledEnvs <- sampledEnvsList[[iTrt]]
        X <- Matrix(0, nrow=nrow(pop@gv), ncol=ncol(pop@gv)); rownames(X) <- pop@id
        availableInds <- sample(1:nrow(X), trts[iTrt,"availableInds"])
        for(iEnv in sampledEnvs){
          
          picked <- sample(availableInds,trts[iTrt,"nIndsPerFarm"] )
          X[picked,iEnv] <- 1
        } # Matrix::image(X)
        nIndsAcrossFarms <- length(which(apply(X,1,function(x){sum(x)}) > 0))
        
        overlap <- (t(X[,sampledEnvs])%*%X[,sampledEnvs])/ trts[iTrt,"availableInds"] # nrow(X)
        if(nrow(overlap) > 1){
          overlapMu <- mean(overlap[lower.tri(overlap)])
        }else{overlapMu <- 1}
        sparsityMu <- 1 - overlapMu
        
        #################################################################
        #################################################################
        #################################################################
        ## MODEL FITTING
        # 1) Testing (erase phenotypes from this specific design)
        Y = pop@pheno # store in new object
        Y = Y * (X/X) # delete phenos
        # 2) Harvestiing (move to long format to use sommer)
        dt <- data.frame(id=rep(pop@id,ncol(X)),pheno=as.vector(Y), env=rep(1:ncol(X), nInd(pop)))
        dt[which(is.nan(dt$pheno)),"pheno"]=NA
        dt$envF <- as.factor(dt$env)
        # 3) Modeling. we fit a model with the GRM
        M <- pullQtlGeno(pop, simParam=SP) # pull true QTL genotypes
        M2 <- M[dt[which(!is.na(dt$pheno)),"id"],]
        fit <- mrr(Y=as.matrix(dt[which(!is.na(dt$pheno)),"pheno"]),X=M2)
        
        if(!inherits(fit, "try-error")){
          # calculate accuracy
          sparseGV <- M%*%fit$b
          common <- intersect(names(trueGV) , rownames(sparseGV)) # currently sommer doesn't complete the matrix, I'll work on it
          rt <- cor(trueGV[common],sparseGV[common,])
          # calculate realized gain
          sparseGV <- as.data.frame(sparseGV)
          sparseGV <- sparseGV[with(sparseGV, order(-V1)),,drop=FALSE ]
          
          # sparseGV <- sparseGV[ order(sparseGV[,1] ), ]
          
          expGain <- mean( trueGV[ rownames(sparseGV)[1:(round(nInd(pop)*p))] ] ) - mean(trueGV)
          # popSel <- selectInd(pop, nInd = round(nInd(pop)*p), trait= function(Y,b){return(b)},  simParam = SP, b=sparseGV )
          # popS1 <- randCross(popSel, nCrosses = nCrosses, nProgeny = nProgeny)
          # realizedGain <- mean(apply(popS1@gv,1,mean)) - mean(trueGV)
          
        }else{rt <- expGain <- NA}
        # save results in new data.frame
        result[[counter]] <- data.frame(accuracy=c(rt), expGain=expGain, nFarms=c(trts[iTrt,"nFarms"]),
                                        nIndsPerFarm=c(trts[iTrt,"nIndsPerFarm"]), nIndsAcrossFarms=nIndsAcrossFarms, 
                                        nIndsAvail=trts[iTrt,"availableInds"], overlapMu=overlapMu,
                                        sparsityMu=sparsityMu,
                                        repTrt=iRep, repMacs=iMacs )
        counter <- counter+1
        
      }
      
      if(verbose){
        message(paste("Macs rep",iMacs,"with treatment",trts[iTrt,"nFarms"], "envs and",trts[iTrt,"nIndsPerFarm"] ,"inds per env completed."))
      }
    }
    
  }
  
  final <- do.call(rbind,result)
  final$nPlots <- final$nFarms * final$nIndsPerFarm
  
  final$propMaxPlot <- round( (final$nPlots / max(final$nPlots))*100, 1)
  final$propMaxInds <- round( (final$nIndsAcrossFarms / max(final$nIndsAcrossFarms))*100, 1)
  
  final$trt <- paste( final$propMaxPlot, "% (", final$nPlots,";",final$nIndsAvail,")", sep="")
  
  final$overlapInds <- ( final$nIndsPerFarm / max(final$nIndsPerFarm) ) * 100
  final$sparsityInds <- 100 - final$overlapInds
  
  final$relativeAccuracy <- final$accuracy/mean(final$accuracy)  # * 100
  final$relativeGain <- final$expGain/mean(final$expGain)  # * 100
  
  
  class(final) <- c(class(final),"sparsetMod")
  return(final)
  
}
