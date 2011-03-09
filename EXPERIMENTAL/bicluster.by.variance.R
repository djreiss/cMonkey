require(foreach)
require(ref)
require(multicore)
require(doMC)


#' Sort a ratios matrix based on the times given in the colMaps
#' 
#' @param ratios  The ratios matrix
#' @param colMaps  The colMap object used by inferelator
#' @export
#' @usage ratios <- sortRatiosByTime(ratios, colMap)
sortRatiosByTime <- function(ratios, colMap) {
	expList<-colMap$time
	names(expList)<-rownames(colMap)
	expList<-sort(expList)
	sortIdx<-match(names(expList),colnames(ratios))
	sortIdx<-sortIdx[!is.na(sortIdx)]
	return(ratios[,sortIdx])
}


#' OBSOLETE Given the genes in a BiCluster and the ratios matrix, calculate the conditions that belong
#' 
#' @param geneNames  The names of the relevant genes
#' @param ratios  The ratios matrix
#' @param numSamples  The number of samples (DEFAULT: 10,000)
#' @param pVal  The pVal cutoff (DEFAULT: 0.05)
#' @export
#' @usage varCutoff <- getCondMembers(geneNames, ratios, numSamples = 10000, pVal = 0.05)
getCondMembers <- function(geneNames, ratios, numSamples = 10000, pVal = 0.05) {
	cutOffs <- getVarianceCutoff(ratios, length(geneNames),numSamples = 10000, pVal = 0.05)
	condVars <- apply(ratios[rownames(ratios) %in% geneNames,],2,var,na.rm=T)
	
	exp2include <- names(condVars)[condVars < cutOffs]
	return(exp2include)
}


#' OBSOLETE Use the distribution of variances to find the variance below a pValue cutoff
#' 
#' @param ratios  The ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: 10,000)
#' @param pVal  The pVal cutoff (DEFAULT: 0.05)
#' @export
#' @usage varCutoff <- getVarianceCutoff(ratios, n,numSamples = 10000, pVal = 0.05)
getVarianceCutoff <- function(ratios, n, numSamples = 10000, pVal = 0.05) {
	varDist <- getVarianceDist(ratios, n,numSamples)
	cutOffs <- apply(varDist,2, function(x) qnorm(pVal,mean(x),sd(x)) )
	return(cutOffs)
}

#' Given a ratios matrix and a number of genes, figure out the expected distribution of variances
#'   Old version that will calculate number of samples using 31 samples
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1, i.e. auto calculate 99% confidence within 1/30s of maximum sampled standard deviation)
#' @export
#' @usage varDist <- getVarianceDist(ratios, n,numSamples = 10000)
getVarianceDist.bk <- function(ratios, n, numSamples = -1) {
	ratios <- derefdata(ratios)

	#Optionally, determine the number of samples you should take with boot strapping
	if (numSamples <= 0) {
		sampleN <- 31 #Based on Statistics by Trioloa, this should be good assuming a normal distribution
		vars <- matrix(0,nrow=sampleN,ncol=ncol(ratios))
		means <- matrix(0,nrow=sampleN,ncol=ncol(ratios))
		for (i in 1:sampleN) {
			vars[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,var,na.rm=T)
			#means[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,mean,na.rm=T)
		}

		alpha <- 0.01 #99% confidence interval
		
		#n = [(Za/2 * sd) / E]
		#E = margin of error, set to 1/30 * variance
		#  Why 30?  Because in test case that would set the error range to .1 (i.e. sd ~= 3)
		#  Why not just use E = .1?  Data can be scaled differently so .1 may not always be best
		#n = [(Za/2 * sd) / (sd^2 / 30)]
		#n = [(30*Za/2)/sd]
		numSamples <- round(( (30*qnorm(alpha/2))/sqrt(max(vars)) )^2)
	}

	vars <- matrix(0,nrow=numSamples,ncol=ncol(ratios))
	for (i in 1:numSamples) {
		vars[i,] <- apply(ratios[sample.int(nrow(ratios),n),],2,var,na.rm=T)
	}
	colnames(vars)<-colnames(ratios)
	return(vars)
}

#' Given a ratios matrix and a number of genes, figure out the expected distribution of variances
#'   Will sample background until the mean and sd converge or the operation times out
#'   Non-parallel due to speed
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param tolerance  The fraction tolance to use as a stopping condition (DEFAULT: 0.01)
#' @param maxTime  The approximate maximum time to run in seconds (DEFAULT: 600)
#' @param chunkSize  The number of samples to add between test (DEFAULT: 200)
#' @param verbose  Set to false to supress output (DEFAULT: F)
#' @param numSamples  Does nothing.  Temp included for backwards compatability (DEFAULT: NULL)
#' @export
#' @usage varDist <- getVarianceMeanSD(ratios, n, tolerance = 0.05 ,maxTime=600, chunkSize=200, verbose=F, numSamples=NULL)
getVarianceMeanSD <- function(ratios, n, tolerance = 0.01 ,maxTime=600, chunkSize=200, verbose=F, numSamples=NULL) {
	ratios <- derefdata(ratios)

	keepRunning <- rep(T,ncol(ratios))
	prevMeans <- rep(0,ncol(ratios))
	prevSds <- rep(0.001,ncol(ratios))
	curMeans <- rep(Inf,ncol(ratios))
	curSds <- rep(Inf,ncol(ratios))
	keepRunning <- rep(T,ncol(ratios))
	names(keepRunning) <- names(curMeans) <- names(curSds) <- names(prevMeans) <- names(prevSds) <- colnames(ratios)
	
	rawScores <- list()
	
	startTime <- proc.time()[3]
	curRep<-1
	
	if (verbose) { cat('\nSample', n, "genes from", ncol(ratios), "conditions") }
	while(any(keepRunning) & ((proc.time()[3]-startTime) < maxTime)) {
		if (verbose) { cat('\nRep = ', curRep," ( ",round(proc.time()[3]-startTime),"s ) ",sep="") }
		prevMeans <- curMeans
		prevSds <- curSds
		newScores <- foreach (i = which(keepRunning)) %do% {
			#if (verbose) { cat(i,",",sep="") }
			curReg <- names(keepRunning)[i]
			curScores <- rep(0,chunkSize)
			for (j in 1:chunkSize) {
				curScores[j] <-var(ratios[sample.int(nrow(ratios),n),curReg],na.rm=T)
			}
			curScores
		}

		#browser()

		for (ctr in 1:sum(keepRunning)) {
			i<-which(keepRunning)[ctr]
			curReg <- names(keepRunning)[i]
			rawScores[[curReg]] <- c(rawScores[[curReg]],newScores[[ctr]])
			curMeans[i] <- mean(rawScores[[curReg]])
			curSds[i] <- sd(rawScores[[curReg]])
		}

		#browser()
		keepRunning <- (abs(curMeans-prevMeans) >= tolerance*abs(prevMeans)) | (abs(curSds-prevSds) >= tolerance*prevSds)
		curRep <- curRep+1
	}
	if(verbose) {cat('\n')}
	
	return(list(means=curMeans,sds=curSds))
}


#' Return the means and SDs for the variances for each experimental condition
#'   Uses: getVarianceDist.bk
#' 
#' @param ratios  A refdata pointing to the ratios matrix
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1 i.e. autodetect)
#' @export
#' @usage varCutoff <- getVarianceMeanSD(ratios, n,numSamples = 10000)
getVarianceMeanSD.bk <- function(ratios, n, numSamples = -1) {
	varDist <- getVarianceDist(ratios, n,numSamples)
	means <- apply(varDist,2, mean )
	sds <- apply(varDist,2, sd )
	return(list(means=means,sds=sds))
}

#' Return the mean and SD for the variances of all experimental conditions
#' 
#' @param ratios  A a list of ratios matrices
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1 i.e. autodetect)
#' @export
#' @usage varCutoff <- getVarianceMeanSD.all(ratios, n,numSamples = 10000)
getVarianceMeanSD.all <- function(ratios, n, numSamples = -1) {
	varDist <- NULL
	for (i in 1:length(ratios)) {
		varDist<-c(varDist,as.numeric(getVarianceDist(refdata(e$ratios[[i]]), n,numSamples)))
	}
	return(data.frame(means=mean(varDist),sds=sd(varDist)))
}


#' Return a dataframe containing the means and standard deviations for all ratios matrixes in list ratios
#' 
#' @param ratios  A list of ratios matrixes
#' @param n  The number of genes
#' @param numSamples  The number of samples (DEFAULT: -1, i.e. autodetect)
#' @param all  Set to false to get a different means and sd for each experimental condition (DEFAULT: T)
#' @export
#' @usage varCutoff <- getVarianceMeanSD(ratios, n,numSamples = 10000, all = T)
getVarianceMeanSD.DF <- function(ratios, n, numSamples = -1, all = T) {
    expNames <- as.character(unlist(lapply(ratios, colnames)))
    if (n == 1) {
        sds <- means <- rep(0,length(expNames))
    } else {
    	if (all == T) {
    		means.sds<-getVarianceMeanSD.all(ratios, n, numSamples)
    		means<-rep(means.sds$means, sum(sapply(ratios,ncol)))
        	sds<-rep(means.sds$sds, sum(sapply(ratios,ncol)))
        } else {
        	means.sds<-sapply(ratios,function(x) {getVarianceMeanSD(refdata(x), n, numSamples = numSamples)})
        	means<-unlist(means.sds['means',],use.names=F)
        	sds<-unlist(means.sds['sds',],use.names=F)
        }
    }
    names(means) <- names(sds) <- expNames
    return(data.frame(means=means,sds=sds))
}

#' UNTESTED Return a pVal or vector of pVals for cluster K.
#' Intended to replace cluster.resid
#' Assumes "rows","rats","cols" `as part of "..."
#' 
#' @param k  The cluster number to calculate the pValue for
#' @param rats.inds  If not set to combined, will return one pVal for each member of "ratios" (DEFAULT: COMBINED)
#' @param means.sds  a list containing one data.frame for each number of genes.  Each data frame has means and sds. (DEFAULT: list())
#' @export
#' @usage resids <- cluster.ratPval( k, rats.inds="COMBINED", means.sds=means.sds ) 
cluster.ratPval <- function( k, rats.inds="COMBINED", means.sds=list(), clusterStack = get("clusterStack"), ... ) {
	clust <- clusterStack[[k]]
	numGenesInClusters<- clust$nrows
	
	#Add to means.sds as necessary
	numGeneList<-numGenesInClusters[! numGenesInClusters %in% names(means.sds)]
	lGeneList<-length(numGeneList)
	if (lGeneList > 0 ) {
		numGeneLists<-split(1:lGeneList,cut(1:lGeneList,floor(lGeneList/e$parallel.cores)))
		cat("Calculating variance means and sd's for",lGeneList,"gene counts\n")
		for (geneListNums in numGeneLists) { #Embed in for loop to have a status monitor
			curNumGenes<-numGeneList[geneListNums]
			cat(curNumGenes,'',sep=",")
			flush.console()
			new.means.sds<- foreach (n=curNumGenes, .inorder=TRUE) %do% {
				getVarianceMeanSD.DF(rats, n) 
			} #for (n in numGeneList)
			names(new.means.sds)<-as.character(curNumGenes)
			means.sds<-c(means.sds,new.means.sds)
		}
		cat('\n')
	}
	
	#Get the average pValue
	curVarList <- means.sds[[as.character(numGenesInClusters)]]
	
	#Get the pValues for each experiment given the number of genes
	rw<-get( "row.weights" )
	pVals <- rep(NA,length(rw))
	names(pVals)<-names(rw)
	for (i in 1:length(pVals) ) {
		rats <- get ("ratios")
		curCols <- clust$cols[clust$cols %in% colnames(rats[[i]])]
		relRats<- rats[[i]][rownames(rats[[i]]) %in% clust$rows,curCols]
		vars <- apply(relRats,2,var,na.rm=T)
		curPs<-NULL
		for (x.idx in 1:length(vars)) {
			x<-vars[x.idx]
			curPs<-c(curPs, pnorm(x,curVarList[names(x),'means'],curVarList[names(x),'sds']) )
		}
		pVals[i]<-mean(curPs)
	}
	
	if ( rats.inds[ 1 ] == "COMBINED" ) pVals <- weighted.mean( pVals, row.weights[ inds ], na.rm=T )

	return(pVals)
}

#' Return a list of pVals for a list of genes under all conditions
#' 
#' @param e  The cMonkey environment.  Used for the ratios matrix
#' @param geneList  ne data.frame for each number of genes.
#' @param mean.sd  A single elements of means.sds.  A DF with means and sds, one for each experimental condition.  
#' @export
#' @usage pValList <- getRatPvals(e, geneList, mean.sd=means.sds[["6"]])
getRatPvals <- function(e, geneList, mean.sd) {
	col.pVals<-NULL
	for (ratIdx in 1:length(e$ratios)) {
		curExps <- colnames(e$ratios[[ratIdx]])
		relRats<- e$ratios[[ratIdx]] [rownames(e$ratios[[ratIdx]]) %in% geneList,]
		vars <- apply(relRats,2,var,na.rm=T)

		pVals<-vars
		for (x.idx in 1:length(vars)) {
			x<-vars[x.idx]
			pVals[x.idx]<-pnorm(x,mean.sd[names(x),'means'],mean.sd[names(x),'sds']) 
		}
		
		col.pVals <- c(pVals,col.pVals)
	}
	col.pVals
}

#' Given a cMonkey environment, build a new clusterStack with different cols & pValues based on variance
#' Returns "newClusterStack" and "means.sds".  "means.sds" may be used in subsequent calls to avoid recomputation
#' 
#' @param e  The cMonkey environment
#' @param means.sds  a list containing one data.frame for each number of genes.  Each data frame has means and sds. (DEFAULT: generate all)
#' @param numSamples  The number of times to sample the background distribution to determine the pValues (DEFAULT: -1, i.e. autodetect)
#' @param pValCut  The pValue at which to cut-off bicluster inclusion.  (DEFAULT: 0.1)
#' @param bonferroni  Set to TRUE to turn bonferroni correction on.  Prelimiary tests show this to be too strict. (DEFAULT: FALSE)
#' @param aveNumCond  The average number of conditions to include in each cluster.  If set, will ignor pValCut   (DEFAULT: NULL)
#' @export
#' @usage clusterStack <- resplitClusters(e, means.sds=list(), numSamples = 10000, bonferroni = F, aveNumCond=NULL))
resplitClusters <- function(e, means.sds=list(), numSamples = -1, pValCut = 0.1, bonferroni = F, aveNumCond=NULL, all = T) {

	#Calculate the background distribution
	numGenesInClusters<- unique(sort(sapply(e$clusterStack,function(x) {x$nrows} )))

	#Load means.sds as necessary
	numGeneList<-numGenesInClusters[! numGenesInClusters %in% names(means.sds)]
	lGeneList<-length(numGeneList)
	if (lGeneList > 0 ) {
		numGeneLists<-split(1:lGeneList,cut(1:lGeneList,floor(lGeneList/e$parallel.cores)))
		cat("Calculating variance means and sd's for",lGeneList,"gene counts\n")
		for (geneListNums in numGeneLists) { #Embed in for loop to have a status monitor
			curNumGenes<-numGeneList[geneListNums]
			cat(curNumGenes,'',sep=",")
			flush.console()
			new.means.sds<- foreach (n=curNumGenes, .inorder=TRUE) %dopar% {
				getVarianceMeanSD.DF(e$ratios, n, numSamples = numSamples, all = all)
			} #for (n in numGeneList)
			names(new.means.sds)<-as.character(curNumGenes)
			means.sds<-c(means.sds,new.means.sds)
		}
		cat('\n')
	}
	
	#Dynamically calculate the pValue cutoffs for each numGenesInCluster
	pValCuts <- rep (pValCut,length(numGenesInClusters))
	names(pValCuts) <- as.character(numGenesInClusters)
	
	if (! is.null(aveNumCond)) {
	
		for ( i in 1:length(numGenesInClusters) ) {
			numGenes <- numGenesInClusters[i]
			curIdxs <- sapply(e$clusterStack,function(x) {x$nrows} ) == numGenes
			if (any(curIdxs) & (numGenes > 1)){
				mean.sd <- means.sds[[as.character(numGenes)]]
				varLists<-lapply(which(curIdxs), function(x) {getRatPvals(e, e$clusterStack[[x]]$rows, mean.sd)})
				pValCuts[as.character(numGenes)]<-as.numeric(sort(c(varLists,recursive=T))[aveNumCond*sum(curIdxs)])
			} else {
				pValCuts[as.character(numGenes)] <- 1
			}
		}
	} #if (pValCut <= 0)
	
	#Build the new clusterStack
	newClustStack <- foreach (clust=e$clusterStack, .inorder=T) %do% {
		if (clust$nrows > 1) {
			pValCut <- pValCuts[as.character(clust$nrows)]
			
			#Get the pValues for each experiment given the number of genes
			col.pVals <- getRatPvals(e, clust$rows, means.sds[[as.character(clust$nrows)]])

			#Bonferroni correction
			if (bonferroni == T) {pValCut <- pValCut / sum(sapply(e$ratios,ncol))}

			newClust<-clust
			newClust$cols <- names(col.pVals)[col.pVals < pValCut]
			if (length(newClust$cols) <= 2) {  newClust$cols <- names(col.pVals)[order(col.pVals)[1:2]] }
			newClust$ncols <- length(newClust$cols)

			#Update the residuals to be the mean pValue of the included clusters
			for (idx in 1:length(newClust$resid)) {
				relCols <- newClust$cols[newClust$cols %in% colnames(e$ratios[[names(clust$resid[idx])]])]
				newClust$resid[idx] <- mean(col.pVals[relCols])
			}
		} else {
			newClust<-clust
		}
		newClust
	}
	
	return(list(newClusterStack=newClustStack, means.sds=means.sds))
}



#These two functions update cMonkey to pull the rows & cols directly from the clusterStack
runnit <- function( e ) {
	##infFile <- 'DR.Inferelator-run-sce-2011-01-05.RData'
	##load(infFile)
	##colMap<-DrData$colMap
	##cMonkeyFile <- 'cmonkey-run-sce-2011-03-04.RData'
	##cMonkeyFile <- 'e.resplit.RData'
	##if (exists('env')) { e <- env }
	
	##load(cMonkeyFile)
	##for (i in 1:length(e$ratios)){
	##	e$ratios[[i]]<-sortRatiosByTime(e$ratios[[i]], colMap)
	##}
	##attr(e$ratios,"cnames")<-unlist(sapply(e$ratios,colnames))
	##environment(e$ratios) <- e

	##e$get.rows <- function (k, cs = get("clusterStack")) 
	##{
	##    return(cs[[k]]$rows)
	##}
	##environment(e$get.rows) <- e

	##e$get.cols <- function (k, cs = get("clusterStack")) 
	##{
	##    return(cs[[k]]$cols)
	##}
	##environment(e$get.cols) <- e

	##e$parallel.cores<-multicore:::detectCores()
	##environment(e$parallel.cores) <- e

	##e$foreach.register.backend(multicore:::detectCores())
	
	newClusterStack <- resplitClusters(e, pValCut = 0.1, aveNumCond=round(sum(sapply(ratios,ncol))/2) )
	
	##e$clusterStack <- newClusterStack$newClusterStack
	##environment(e$clusterStack) <- e
	
	e$means.sds <- newClusterStack$means.sds
	##environment(e$means.sds) <- e
	
	row.col.membership.from.clusterStack <- function( cs ) {
		row.memb <- row.membership * 0
		col.memb <- col.membership * 0
		for ( k in 1:length( cs ) ) {
	  		if ( k > ncol( row.memb ) ) row.memb <- cbind( row.memb, rep( 0, nrow( row.memb ) ) )
	  		rows <- cs[[ k ]]$rows; rows <- rows[ ! is.na( rows ) ]
	  		row.memb[ rows, k ] <- k
	  		if ( k > ncol( col.memb ) ) col.memb <- cbind( col.memb, rep( 0, nrow( col.memb ) ) )
	  		cols <- cs[[ k ]]$cols; cols <- cols[ ! is.na( cols ) ]
	  		col.memb[ cols, k ] <- k
	  	}
	  	row.memb <- t( apply( row.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
	  	row.memb <- row.memb[ ,apply( row.memb, 2, sum ) != 0, drop=F ]
	  	colnames( row.memb ) <- NULL
	  	col.memb <- t( apply( col.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
	  	col.memb <- col.memb[ ,apply( col.memb, 2, sum ) != 0, drop=F ]
	  	colnames( col.memb ) <- NULL
	  	list( r=row.memb, c=col.memb )
	}
	
	##environment( row.col.membership.from.clusterStack ) <- e
	tmp <- row.col.membership.from.clusterStack( newClusterStack$newClusterStack ) ##e$clusterStack )
	e$row.membership <- tmp$r
	e$col.membership <- tmp$c
	e$row.memb <- t( apply( e$row.membership, 1, function( i ) 1:k.clust %in% i ) )
	e$col.memb <- t( apply( e$col.membership, 1, function( i ) 1:k.clust %in% i ) )
	e$clusterStack <- get.clusterStack( force=T )
	e$stats <- rbind( stats, get.stats() )
	stop()
	  
	save(e,file="e.resplit.RData")
	
	#browser()	
	source('write.project.R')
	e$write.project()
}

if ( TRUE ) {
  runnit( e )
}