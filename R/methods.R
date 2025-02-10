##' @importFrom stats update
##' @S3method update evolaMod
plot.sparsetMod <- function(object, y="expGain", lb=0, ub=Inf, box=TRUE, color="nIndsPerFarm", ...) {
  # if (is.null(call <- getCall(object))){stop("object should contain a 'call' component")}
  # call <- getCall(object)
  if(!inherits(object,"sparsetMod")){stop("not an object of class sparsetMod")}
  
  # restrain plot
  object <- object[which(object$nPlots <= ub),]
  object <- object[which(object$nPlots >= lb),]
  
  object$propMaxPlot <- round( (object$nPlots / max(object$nPlots))*100, 0)
  object$propMaxInds <- round( (object$nIndsAvail / max(object$nIndsAvail))*100, 0)
  object$trt <- paste0( addZeros(object$nIndsAvail),"ia_",  addZeros(object$nIndsPerFarm),"ipf_" , addZeros(object$nPlots ) , "pt ")#,  object$propMaxInds,"% ;", object$propMaxPlot ,"%)", sep="")
  
  object$colorBy <- object[,color]
  # make sure we color by averages
  for(iTrt in unique(object$trt)){ # iTrt = unique(object$trt)[1]
    object[which(object$trt == iTrt),"colorBy"] <- round(mean(object[which(object$trt == iTrt),"colorBy"]))
  }
  
  object$yVar <- object[,y]
  myYint <- mean(object[,y])
  
  p <- ggplot(object, aes(x=factor(trt), y=yVar, fill=colorBy) )
  if(box){
    p <- p + geom_boxplot( alpha = 0.7)
  }else{
    p <- p + geom_violin( alpha = 0.7) 
  }
  
  if(y=="accuracy"){myYlab <- "Accuracy [cor(true,estimated)]"}else if(y=="expGain"){myYlab <- "Expected genetic gain (accuracy and intensity relationship)"}else{myYlab <- y}
  if(color=="overlapMu"){colorTitle<-"entry \noverlap"}else if(color=="sparsityMu"){colorTitle <- "entry \nsparsity"}else if(color=="nIndsAcrossFarms"){colorTitle <- "#entries \nacross \nfarms"}else{colorTitle <- color}
  
  p <- p +  geom_point(position = position_jitter(seed = 1, width = 0.2), alpha=0.1) +
    scale_x_discrete(breaks=factor(object$trt), labels=object$trt) +
    facet_wrap(~nFarms, scales = "free_x") + #ylim(c(0,1)) +
    guides(x= guide_axis(angle=45), fill = guide_colourbar(barwidth = NULL, 
                                                           title=colorTitle,
                                                           nbins=length(unique(object$colorBy)),
                                                           barheight = 25, title.position = "bottom", 
                                                           title.hjust = 0.5)) +
    geom_hline(yintercept=myYint, linetype='dashed', color=c('red')) +
    xlab("#Inds available _ #Inds per farm _ #Plots total") + 
    ylab(myYlab) +
    labs(fill = color ) +
    scale_fill_gradient(low="forestgreen", high="red") 
  
  
  return(p)
}