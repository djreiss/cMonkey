library(foreach)
if (require(doMC, quietly = T) & !getDoParRegistered() & require(multicore) ) { registerDoMC(cores = multicore:::detectCores()) }

#' Run tomtom on a cMonkey environment
#'   tomtom requires meme4.6.1 because that's the version that the downloadable databases work on.
#'   That can be a real pain in the ass.  This will require a lot of programming to integrate into main cMonkey
#'   NOTE: This function doesn't seem to work unless you copy the contents and paste them.
#'     I guess it's not easy to pass an environment as a parameter
#' 
#' @param env
#' @return a list of TFs and scores
#' @export
#' @usage regList <- tomtom.env(env)
#tomtom.env <- function( env )
#{
	environment(tomtom.clusters) <- env
	environment(tomtom.one.cluster) <- env
	regList <- tomtom.clusters(verbose=F)
#}

#' Run tomtom on a cMonkey cluster to get putative regulators
#'   tomtom requires meme4.6.1 because that's the version that the downloadable databases work on.
#'   That can be a real pain in the ass.  This will require a lot of programming to integrate into main cMonkey
#' 
#' @param ks  The clusters to run tomtom on (DEFAULT: 1:k.clust)
#' @param databases the database the quesry against. NOTE: can be downloaded from meme site  (DEFAULT: c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'))
#' @param seq.type  (DEFAULT: names(mot.weights)[1])
#' @param e.value.cutoff  (DEFAULT: 1)
#' @param resid.cutoff  (DEFAULT: 0.8)
#' @param dist.meth  [allr|ed|kullback|pearson|sandelin] (DEFAULT: 'ed')
#' @param q.thresh  (DEFAULT: 0.5)
#' @param min.overlap  (DEFAULT: 4)
#' @param q.pseudo  (DEFAULT: 0)
#' @param t.pseudo  (DEFAULT: 0)
#' @param min.gene.overlap  (DEFAULT: NA)
#' @param verbose  (DEFAULT: T)
#' @return a list of TFs and scores
#' @export
#' @usage regList <- tomtom.clusters(ks=1:k.clust, seq.type=names(mot.weights)[1],databases = c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'),e.value.cutoff=10, resid.cutoff=0.8, dist.meth="ed", q.thresh=0.5, min.overlap=4, q.pseudo=0, t.pseudo=0, min.gene.overlap=NA, verbose=T)
tomtom.clusters <- function( ks=1:k.clust, seq.type=names(mot.weights)[1],
                                   databases = c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'),
                                   e.value.cutoff=1, resid.cutoff=0.8, dist.meth="ed", 
                                   q.thresh=0.5, min.overlap=4, q.pseudo=0, t.pseudo=0, 
                                   min.gene.overlap=NA, verbose=T ) {
	regList <- foreach (k = ks) %dopar% {
		tomtom.one.cluster(k, seq.type, databases, e.value.cutoff, resid.cutoff, dist.meth, q.thresh, min.overlap, q.pseudo, t.pseudo, min.gene.overlap, verbose)	
	}
	return(regList)
}

#' Run tomtom on a cMonkey cluster to get putative regulators
#'   tomtom requires meme4.6.1 because that's the version that the downloadable databases work on.
#'   That can be a real pain in the ass.  This will require a lot of programming to integrate into main cMonkey
#' 
#' @param k  The cluster number to run tomtom on
#' @param databases the database the quesry against. NOTE: can be downloaded from meme site  (DEFAULT: c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'))
#' @param seq.type  (DEFAULT: names(mot.weights)[1])
#' @param e.value.cutoff  (DEFAULT: 1)
#' @param resid.cutoff  (DEFAULT: 0.8)
#' @param dist.meth  [allr|ed|kullback|pearson|sandelin] (DEFAULT: 'ed')
#' @param q.thresh  (DEFAULT: 0.5)
#' @param min.overlap  (DEFAULT: 4)
#' @param q.pseudo  (DEFAULT: 0)
#' @param t.pseudo  (DEFAULT: 0)
#' @param min.gene.overlap  (DEFAULT: NA)
#' @param verbose  (DEFAULT: T)
#' @return a list of TFs and scores
#' @export
#' @usage hitTable <- tomtom.one.cluster(k, seq.type=names(mot.weights)[1],databases = c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'),e.value.cutoff=10, resid.cutoff=0.8, dist.meth="ed", q.thresh=0.5, min.overlap=4, q.pseudo=0, t.pseudo=0, min.gene.overlap=NA, verbose=T)
tomtom.one.cluster <- function( k, seq.type=names(mot.weights)[1],
                                   databases = c('./data/motif_databases/macisaac_yeast.v1.meme','./data/motif_databases/scpd_matrix.meme'),
                                   e.value.cutoff=1, resid.cutoff=0.8, dist.meth="ed", 
                                   q.thresh=0.5, min.overlap=4, q.pseudo=0, t.pseudo=0, 
                                   min.gene.overlap=NA, verbose=T ) {

  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  for (database in databases) {
      if (!file.exists(database)) {
        cat('Tomtom Database file:',database,'does not exist\n')
        return(NULL)
      }
  }

  cluster.motif.lines <- function( k, mot=NA ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    lines <- character()
    memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
    if ( is.null( memeOut ) ) return( lines )
    if ( all(clusterStack[[ k ]]$resid > resid.cutoff) ) return( lines )
    ##max.motifs <- max( max.motifs, length( memeOut ) )
    for ( i in 1:length( memeOut ) ) {
      if ( ! is.na( mot ) && i != mot ) next
      if ( memeOut[[ i ]]$e.value > e.value.cutoff ) next
      pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
      lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                         memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
      lines <- c( lines, apply( pssm, 1, function( i )
                               sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    }
    lines
  }

  memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
  hitTable <- NULL
  for (i in 1:length(memeOut)) {
    tlines <- c(lines, cluster.motif.lines(k,i))
    cmdStr <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f %s %s"
    tfile <- my.tempfile( "tomtom_in_", )
    
    if(is.null(memeOut[[i]]$pssm)) { next }
    
    memeSeq=paste(meme.let[apply(memeOut[[i]]$pssm,1,which.max)],collapse="")

    if (verbose) cat("CLUSTER",k,"TARGET MOTIF:", i,memeSeq,"\n")
    cat( tlines, file=tfile, sep="\n" ) ## Write out all target motifs
  
    for (database in databases) {
      cmd <- sprintf( cmdStr, progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile, database )
      tout <- system( cmd, intern=T, ignore.stderr = T )
      
      for (out in tout) {
      	 if(grepl('^#Query', out)) next
      	 out <- strsplit(out,'\t')[[1]]
      	 if (as.numeric(out[5]) <= e.value.cutoff) {
      	   hitTable <- rbind(hitTable, data.frame(TargetID=out[2],Pvalue=as.numeric(out[4]),Evalue=as.numeric(out[5]),Qvalue=as.numeric(out[6]),memeSeq=memeSeq,meme.Eval=memeOut[[i]]$e.value,database=database))
      	 }
      }
    } #for (database in databases) { 
    unlink( tfile )
  } #for (i in 1:length(memeOut)) {
  
  return(hitTable)
}