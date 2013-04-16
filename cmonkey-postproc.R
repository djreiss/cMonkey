###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## find # of genes per cluster with annoatated "long name" (e.g. "peroxisome") or short name (e.g. "PEX")
clusters.w.func <- function( func, ks=1:k.clust, short=F, max.rows=999, p.val=F ) {
  if ( p.val ) {
    long.names <- get.long.names( attr( ratios, "rnames" ), short=short )
    n2 <- length( grep( func, long.names, perl=T, ignore.case=T ) )
  }
  ##sapply( ks, function( i ) {  
  mc <- get.parallel( length( ks ) )
  unlist( mc$apply( ks, function( i ) { ##, ... ) {
    rows <- get.rows( i )
    if ( length( rows ) <= 1 ) return( NA )
    rows.l <- get.long.names( rows, short=short )
    if ( ! p.val ) {
      if ( length( rows ) >= max.rows ) NA else
      length( grep( func, rows.l, perl=T, ignore.case=T ) )
    } else {
      phyper( length( grep( func, rows.l, perl=T, ignore.case=T ) ),
             n2, attr( ratios, "nrow" ) - n2, length( rows ), lower=F ) * length( ks ) ## bonferroni baby!
    }
  } ) ) ##, mc.cores=mc$par ) )
}

## Find # of genes per cluster out of given list of genes
clusters.w.genes <- function( genes, ks=1:k.clust, p.val=F ) {
  ##sapply( ks, function( i ) {
  mc <- get.parallel( length( ks ) )
  unlist( mc$apply( ks, function( i ) { ##, ... ) {
    rows <- get.rows( i )
    if ( length( rows ) <= 1 ) return( NA )
    if ( ! p.val ) sum( rows %in% genes )
    else phyper( sum( rows %in% genes ), length( genes ), attr( ratios, "nrow" ) - length( genes ),
                length( rows ), lower=F ) * length( ks )
  } ) ) ##, mc.cores=mc$par ) )
}

## Find # of cols per cluster out of given list of cols
clusters.w.conds <- function( conds, ks=1:k.clust, p.val=F ) {
  ##sapply( ks, function( i ) {
  mc <- get.parallel( length( ks ) )
  unlist( mc$apply( ks, function( i ) { ##, ... ) {
    cols <- get.cols( i )
    if ( ! p.val ) sum( cols %in% conds )
    else phyper( sum( cols %in% conds ), length( conds ), attr( ratios, "ncol" ) - length( conds ),
                length( cols ), lower=F ) * length( ks )
  } ) ) ##, mc.cores=mc$par ) )
}

## clusters.w.cogs <- function( ks=1:k.clust, p.val=F ) {
##   cog.ns <- table( genome.info$cog.code )
##   sapply( ks, function( i ) {
##     if ( ! p.val ) {
##       out <- max( table( genome.info$cog.code[ get.rows( i ) ] ), na.rm=T )
##       tmp2 <- names( which.max( table( genome.info$cog.code[ get.rows( i ) ] ) ) )
##       names( out ) <- tmp2; out
##     } else {
##       tmp <- max( table( genome.info$cog.code[ get.rows( i ) ] ), na.rm=T )
##       if ( is.infinite( tmp ) ) return( NA )
##       tmp2 <- names( which.max( table( genome.info$cog.code[ get.rows( i ) ] ) ) )
##       out <- phyper( tmp, cog.ns[ tmp2 ], sum( cog.ns ) - cog.ns[ tmp2 ], length( get.rows( i ) ),
##                     lower=F ) / length( ks ) ## bonferroni, baby!
##       names( out ) <- tmp2; out
##     }
##   } )
## }

## cluster.groupings <- function( k, grp=names( grouping.weights )[ 1 ] ) {
##   if ( file.exists( grp ) ) {
##     sif <- load.sif.interactions( grp )
##   } else {
##     sif <- get( grp )
##     if ( ncol( sif ) == 2 ) sif <- cbind( sif, combined_score=rep( 1, nrow( sif ) ) ) ## Add a weights column
##   }
##   ## We will assume that column with fewer unique names is the "groups" column.
##   colnames( sif ) <- c( "group", "protein", "combined_score" )
##   if ( length( unique( as.character( sif$protein ) ) ) < length( unique( as.character( sif$group ) ) ) ) {
##     sif <- sif[ ,c( 2, 1, 3 ) ]
##     colnames( sif ) <- c( "group", "protein", "combined_score" )
##   }
##   sif <- sif[ as.character( sif$protein ) %in% attr( ratios, "rnames" ), ]
##   ##sif <- sif[ order( sif$group ), ]
##   ##tmp <- tapply( sif$protein, sif$group )
##   ##Bnames( tmp ) <- as.character( sif$protein )
##   ##tmp <- tmp[ get.rows( k ) ]; tmp <- tmp[ ! is.na( tmp ) ]
##   if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
##   else rows <- k
##   tmp <- as.character( sif[ as.character( sif$protein ) %in% rows, "group" ] )
##   table( tmp )
## }


cluster.summary <- function( e.cutoff=0.01, nrow.cutoff=5, seq.type=names( mot.weights )[ 1 ], plot=F,
                            sort=c("score.norm","score","resid","e.value1","e.value2","nrow") ) { ##"loglik",
  ms <- NULL
  if ( ! is.null( seq.type ) ) ms <- meme.scores[[ seq.type ]]
  if ( is.null( ms ) ) e.cutoff <- NA
  score <-
    sapply( 1:k.clust, function( k ) mean( row.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
      row.scaling[ iter ] + if ( ! is.null( mot.scores ) )
        sapply( 1:k.clust, function( k ) mean( mot.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
          mot.scaling[ iter ] else 0 + if ( ! is.null( net.scores ) )
            sapply( 1:k.clust, function( k ) mean( net.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
              net.scaling[ iter ] else 0
  nrow <- sapply( 1:k.clust, function( k ) length( get.rows( k ) ) )

  out <- data.frame( k=1:k.clust, nrow=nrow, score=score, ##score.norm=score.norm,
                    resid=sapply( 1:k.clust, cluster.resid, varNorm=F ), 
                    consensus1=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 1 ]]$pssm ) ),
                    e.value1=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) Inf else
                      ms[[ k ]]$meme.out[[ 1 ]]$e.value ),
                    consensus2=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) "" else
                      if ( length( ms[[ k ]]$meme.out ) == 1 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 2 ]]$pssm ) ),
                    e.value2=sapply( 1:k.clust,
                      function( k ) if ( is.null( ms ) || is.null( ms[[ k ]]$meme.out ) || length( ms[[ k ]] ) <= 3 ) Inf else
                      if ( length( ms[[ k ]]$meme.out ) <= 1 ) Inf else
                      ms[[ k ]]$meme.out[[ 2 ]]$e.value )
                    )
  if ( all( out$consensus2 == "" ) ) out$consensus2 <- out$e.value2 <- NULL
  if ( ! is.na( sort[ 1 ] ) && sort[ 1 ] %in% colnames( out ) ) out <- out[ order( out[[ sort[ 1 ] ]] ), ]
  out
}

## Code for reloading updated cMonkey package, updating existing env. obj. and restarting cmonkey run on that env:
##  unload("cMonkey");require(cMonkey);update.cmonkey.env(e);cmonkey(e,dont.init=T)
update.cmonkey.env <- function( object, ... ) { ## Update all funcs contained in env to latest from cmonkey package
  if ( file.exists( "cmonkey-funcs.R" ) ) {
    tmp.e <- new.env()
    sys.source( "cmonkey.R", envir=tmp.e ) ##cmonkey.env )
    ##sys.source( "cmonkey-postproc.R", envir=tmp.e ) ##cmonkey.env )
  } else {
    tmp.e <- environment( cMonkey:::cmonkey ) ## Packaged - get the env. that the "cmonkey" function is stored in
  }
  
  for ( i in ls( tmp.e ) ) {
    if ( i %in% c( "DATE", "VERSION" ) ) ##, "cm.version", "cmonkey.time.started", "cmonkey.init", "cmonkey.session.info",
##                   "cmonkey.re.init", "cog.org", "col.let", "date.run", "dlf", "extend.vec", "get.COG.code", 
##                   "get.condition.groups", "get.genome.info", "get.operon.predictions", "get.predictome.links",
##                   "get.prolinks.links", "getMastPValuesAndEValues", "getMemeMotifInfo", "getMemeMotifPssm",
##                   "load.ratios", "load.sif.interactions", "get.condition.groups", "get.STRING.links", 
##                   "pssm.to.string", "rnd.seed", "rev.comp", "viewPssm", "get.mast.pvals", "mkBgFile",
##                   "remove.low.complexity", "residual.submatrix", "runMast", "runMeme", "system.time.limit", 
##                   "rsat.species", "rsat.urls", "save.logfile", "seed.clusters", "seed.method",
##                   "set.param", "string.version", "taxon.id", "update.cmonkey.env" ) )
      next
    f <- try( get( i, envir=tmp.e ) ) ## Copy func from cmonkey package env.
    f2 <- try( get( paste( "super", i, sep="." ), envir=object ), silent=T ) ## Original (overridden) func that may exist in env
    if ( class( f ) == "function" ) {
      environment( f ) <- object
      if ( class( f2 ) != "function" ) assign( i, f ) else assign( paste( "super", i, sep="." ), f )
    }
  }
  rm( f, f2, tmp.e, i )
  for ( i in ls() ) {
    if ( i %in% c( "i", "object" ) ) next
    f <- get( i )
    if ( is.function( f ) ) assign( i, f, object )
  }

  ##invisible( env )
}
