###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss.isb@gmail.com.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

get.rows <- function( k ) clusterStack[[ k ]]$rows
         ##, rm=get("row.membership") ) { out <- unique( rownames( which( rm[] == k, arr=T ) ) );
         ##                                             if ( is.null( out ) ) out <- character(); out } 

get.cols <- function( k ) {
  out <- clusterStack[[ k ]]$cols
  if ( is.null( out ) || is.na( out ) || length( out ) <= 0 ) out <- attr( ratios, "cnames" )
  out }
         ##, cm=get("col.membership") ) { out <- unique( rownames( which( cm[] == k, arr=T ) ) );
         ##                                             if ( is.null( out ) ) out <- character(); out } 

cmonkey.one.iter <- function( env, dont.update=F, ... ) {
  env <- env$update.all.clusters( env, dont.update=F, ... )

  row.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "rnames" ) %in% get.rows( k ) )
  if ( is.vector( row.memb ) ) row.memb <- t( row.memb )
  rownames( row.memb ) <- attr( ratios, "rnames" )    
  col.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "cnames" ) %in% get.cols( k ) )
  if ( is.vector( col.memb ) ) col.memb <- t( col.memb )
  rownames( col.memb ) <- attr( ratios, "cnames" )    
  
  ## OUTPUT
  if ( iter %in% stats.iters ) { 
    ##env$clusterStack <- get.clusterStack( ks=1:k.clust )
    tmp.stats <- env$get.stats()
    stats <- env$stats
    if ( iter > 1 && nrow(stats) > 0 &&
        any( ! colnames( tmp.stats ) %in% colnames( stats ) ) ) { ## this is for multiple seq-types of motif-finder-types
      ## assume tmp.stats has more columns than stats ...
      for ( cn in colnames( tmp.stats )[ ! colnames( tmp.stats ) %in% colnames( stats ) ] )
        stats[[ cn ]] <- rep(NA, nrow(stats))
    }
    ##env$stats <- rbind( env$stats, env$get.stats() )
    env$stats <- rbind( stats, tmp.stats )
    cat( organism, as.matrix( env$stats[ nrow( env$stats ), ] ), '\n' )
  } else {
    cat( sprintf( "==> %04d %.3f %.3f %.3f\n", iter, ##%5d sum( row.memb != old.row.memb, na.rm=T ),
                 mean( env$row.scores[,][ row.memb ], na.rm=T ), ##mean( col.scores[ col.memb ], na.rm=T, trim=0.05 ),
                 if ( ! is.null( env$mot.scores ) )
                 mean( env$mot.scores[,][ row.memb & env$mot.scores[,] < 0 ], na.rm=T, trim=0.05 )
                 else NA,
                 if ( ! is.null( env$net.scores ) ) mean( env$net.scores[,][ row.memb ##& net.scores < 0
                                                                      ], na.rm=T, trim=0.05 ) else NA ) ) ##, "\n" )
  }
  
  ## PLOTTING
  if ( ! is.na( plot.iters ) && iter %in% plot.iters ) {
    ##try(
    env$plotStats( iter, plot.clust=env$favorite.cluster(), new.dev=T ) ##, silent=T ) ## Can be set for your given organism
  }
  
  if ( exists( "cm.func.each.iter" ) ) try( cm.func.each.iter(), silent=T ) ## User-defined func. to run each iteration
    
  ## Allow temp source file to be sourced (e.g. to change a param in the middle of a run, or print or plot or save
  ## some intermediate results). If first line of file is '## QUIET' then this is done quietly.
  if ( any( cm.script.each.iter != "" ) ) {
    for ( f in cm.script.each.iter ) {
      if ( file.exists( f ) && file.info( f )$size > 1 ) {
        tmp <- readLines( f )
        if ( all( substr( tmp, 1, 1 ) == "#" ) ) next ## All commented-out code
        if ( tmp[ 1 ] != "## QUIET" ) cat( "Sourcing the script '", f, "' ...\n", sep="" )
        try( source( f, echo=tmp[ 1 ] != "## QUIET", local=T ), silent=T )
        rm( tmp )
      }
    }
  }
  ##}

#ifndef PACKAGE
  ##if ( big.memory == TRUE || big.memory > 0 )
  ##ffify.env( env ) ## Big matrices and lists go to filebacked version
#endif
  
  ## Note: with the above code, when the env is saved via save.image(), all ff obj's are "closed" but their
  ##   filestores still exist, so you can "open.ff(x)" each of them after the env is re-loaded,
  ##   and that will reconnect them with their files.
  
  if ( get.parallel()$mc ) {
    if ( getDoParName() == "doMC" ) { ##require( parallel, quietly=T ) ) { ## Clean up any parallel spawned processes (as doc'ed in mclapply help)
      chld <- parallel:::children()
      if ( length( chld ) > 0 ) { try( { parallel::mckill( chld ); tmp <- parallel:::mccollect( chld ) }, silent=T ) }
    } else if ( getDoParName() == "doSNOW" && "data" %in% ls( pos=foreach:::.foreachGlobals ) ) {
      cl <- get( "data", pos=foreach:::.foreachGlobals ) ## Tricky, eh?
      if ( ! is.null( data ) ) stopCluster( cl )
    }
  }

  if ( ! dont.update ) env$iter <- env$iter + 1
  invisible( env )
}

## TODO: store/plot indiv. resids and network scores for each separate component
get.stats <- function( mean.func=mean ) { ##median
  changed <- NA

  if ( ! exists( "row.memb" ) ) {
    row.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "rnames" ) %in% get.rows( k ) ); if ( is.vector( row.memb ) ) row.memb <- t( row.memb )
    rownames( row.memb ) <- attr( ratios, "rnames" )    
    col.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "cnames" ) %in% get.cols( k ) ); if ( is.vector( col.memb ) ) col.memb <- t( col.memb )
    rownames( col.memb ) <- attr( ratios, "cnames" )    
  }

  if ( ! exists( "row.scores" ) || is.null( row.scores ) ) {
    if ( attr( get.all.scores, "version" ) == 1 ) {
      tmp <- get.all.scores()
      row.scores <- tmp$r; mot.scores <- tmp$m; net.scores <- tmp$n; col.scores <- tmp$c
    } else if ( attr( get.all.scores, "version" ) == 2 ) {
      tmp <- get.old.scores.matrices()
      row.scores <- tmp$r; mot.scores <- tmp$m; net.scores <- tmp$n; col.scores <- tmp$c
    }
  }

  cs <- as.list( clusterStack )
  resids <- sapply( cs, "[[", "resid" )
  if ( is.matrix( resids ) ) resids <- apply( resids, 1, function( r ) mean.func( r[ r != 1.0 ], na.rm=T ) )
  else resids <- mean.func( resids[ resids != 1.0 ], na.rm=T )
  p.clusts <- sapply( cs, "[[", "p.clust" )
  if ( is.list( p.clusts ) ) p.clusts <- sapply( cs[ ! sapply( p.clusts, is.null ) ], "[[", "p.clust" )
  if ( is.matrix( p.clusts ) ) p.clusts <- apply( p.clusts, 1, mean.func, na.rm=T ) 
  else p.clusts <- mean.func( p.clusts, na.rm=T )
  out <- data.frame( iter=iter, changed=changed,
             row.scores=mean.func( row.scores[,][ row.memb[,] ], na.rm=T ), ##, trim=0.05 ),
             col.scores=mean.func( col.scores[,][ col.memb[,] ], na.rm=T ), ##, trim=0.05 ),
             mot.scores=if ( ! is.null( mot.scores ) ) mean.func( mot.scores[,][ row.memb[,] ], na.rm=T ) else NA, 
             net.scores=if ( ! is.null( net.scores ) ) mean.func( net.scores[,][ row.memb[,] ], na.rm=T ) else NA,
             resid=weighted.mean( resids, row.weights, na.rm=T ),
             nrow=mean.func( sapply( cs, "[[", "nrows" ), na.rm=T ),
             ncol=mean.func( sapply( cs, "[[", "ncols" ), na.rm=T ),
             p.clust=if ( ! all( is.na( p.clusts ) ) ) weighted.mean( p.clusts, mot.weights, na.rm=T ) else NA
             )
  if ( length( resids ) > 1 ) for ( i in names( resids ) ) {
    out <- cbind( out, resids[ i ] )
    names( out )[ ncol( out ) ] <- paste( "resid", i, sep="." )
  }
  if ( length( p.clusts ) > 1 ) for ( i in names( p.clusts ) ) {
    out <- cbind( out, p.clusts[ i ] )
    names( out )[ ncol( out ) ] <- paste( "p.clust", i, sep="." )
  }
  ## net.scores weighted by n.rows of bicluster
  if ( length( networks ) > 1 ) {
    for ( i in names( net.weights ) ) {
      if ( exists( "cluster.net.scores" ) && i %in% colnames( cluster.net.scores ) )
        out <- cbind( out, weighted.mean( cluster.net.scores[ ,i ], sapply( cs, "[[", "nrows" ), na.rm=T ) )
      else out <- cbind( out, rep( NA, nrow( out ) ) )
      names( out )[ ncol( out ) ] <- paste( "net", i, sep="." )
    }
    if ( exists( "cluster.net.scores" ) && "net.scores" %in% colnames( cluster.net.scores ) ) {
      out[ ,"net.scores" ] <- weighted.mean( cluster.net.scores[ ,"net.scores" ],
                                            sapply( cs, "[[", "nrows" ), na.rm=T )
    } else {
      out[ ,"net.scores" ] <- mean.func( net.scores[,][ row.memb[,] ], na.rm=T ) ##, trim=0.05 )
    }
  }
  out
}

## Only override global parameters if they were not already set previously ----
##   but even if they were, copy them into the 'cmonkey.params' environment
set.param <- function( name, val, env=cmonkey.params, override=F, quiet=F ) {
  if ( ! exists( name, envir=env ) || override ) {
    if ( ! quiet ) try( { cat( name, "-> " ); str( val, digits.d=9, no.list=T ) } )
    assign( name, val, envir=env )
  } else {
    val <- get( name, envir=env )
    if ( ! quiet ) try( { cat( name, "= " ); str( val, digits.d=9, no.list=T ) } ) 
    assign( name, val, envir=env )
  }
  assign( name, val, envir=parent.frame() )
}

## Make sure all processes are killed via kill(children(),SIGKILL) ??
## Use "parallel.cores" param to determine if parallelization is desired - if it is FALSE or 0 or 1 or NA, then no.
## If it is >1 or TRUE then yes. If TRUE, then determine # of cores via parallel::detectCores()
## This should allow running on Windows boxes when there is no parallel package.
## Note that with the use of "foreach", this can now run on Windows if we use doSMP instead of doMC !!!
##    Set the "foreach.register.backend" function to a different backend if desired.

##is.windows.system <- function() FALSE

foreach.register.backend <- function( par ) { ##, force.snow=is.windows.system() ) {
  if ( ! require( foreach ) ) return( NULL )
  if ( par > 1 && require( doMC, quietly=T ) ) registerDoMC( cores=par )
  ##if ( par == 1 ) { registerDoSEQ(); return( NULL ) }
  ##if ( ! is.windows.system() && ! force.snow ) {
  ## else if require( doSNOW ) {
  ##  if ( ! "data" %in% ls( pos=foreach:::.foreachGlobals ) ) {
  ##    cl <- makeCluster( rep( 'localhost', par ), "SOCK" )
  ##    registerDoSNOW( cl )
  ##  }
  ##}
  else registerDoSEQ()
}

get.parallel <- function( X=k.clust, verbose=F, para.cores=get( "parallel.cores" ) ) {
  if (
#ifndef PACKAGE
      ( exists( "DEBUG" ) && ! is.function( DEBUG ) && DEBUG == TRUE ) ||
#endif
      is.na( para.cores ) || ( is.logical( para.cores ) && para.cores == FALSE ) ||
      ( is.numeric( para.cores ) && para.cores <= 1 ) ) {
    out <- list( mc=FALSE, par=para.cores, apply=lapply )
    if ( verbose ) cat( "NOT PARALLELIZING\n" )
  } else {
    try( has.multi <- require( parallel, quietly=T ), silent=T )
    if ( ! has.multi || ( has.multi && parallel:::isChild() ) ) {    
      out <- list( mc=FALSE, par=para.cores, apply=lapply )
      if ( verbose ) cat( "NOT PARALLELIZING\n" )
    } else {
      mc <- has.multi && ! parallel:::isChild() && X > 1 && ! is.na( para.cores ) &&
      ( is.numeric( para.cores ) && para.cores > 1 ) ||
      ( is.logical( para.cores ) && para.cores == TRUE )
      par <- para.cores
      out.apply <- lapply 
      if ( mc ) {
        if ( is.logical( par ) && par == TRUE ) par <- parallel:::detectCores() ## all.tests=TRUE )
        par <- min( c( X, par, parallel:::detectCores() ) ) ## all.tests=TRUE ) ) )
        if ( verbose ) cat( "PARALLELIZING:", par, ": " )
        ## if ( ! exists( "foreach.register.backend" ) || is.null( foreach.register.backend ) ||
        ##     is.null( foreach.register.backend( par ) ) ) { ##use.foreach ) {
        ##   out.apply <- mclapply
        ## } else {
        foreach.register.backend( par )
        if ( verbose ) cat( getDoParName(), getDoParWorkers(), "\n" )
        out.apply <- function( list, FUN, ... ) foreach( l=list ) %dopar% { FUN( l, ... ) }
        ## }
        ##if ( attr( ratios, "nrow" ) > big.run ) print( gc() ) ## gc() before we spawn new copied processes
      } else {
        par <- 1
        if ( verbose ) cat( "NOT PARALLELIZING:", par, "\n" )
      }
      out <- list( mc=mc, par=par, apply=out.apply )
    }
  }
  ##if ( attr( ratios, "nrow" ) > big.run ) print( gc() ) ## gc() before we spawn new copied processes
  if ( is.numeric( out$par ) && ! is.na( out$par ) ) options( mc.cores=out$par ) ## getOption("cores") is the default of how mclapply gets its mc.cores number
  else if ( is.na( out$par ) || ( is.logical( out$par ) && out$par == TRUE ) ) options( cores=NULL )
  else options( cores=1 )
  out
}

#' Create clusterStack-style object 
#'  SD Notes: force is a little mysterious to me.
#'
#' @params ks The clusters to include (DEFAULT: 1:k.clust)
#' @params force If FALSE, will return the existing clusterStack if possible (DEFAULT: F)
#' @params ...  Must include atleast "iter"
#' @export
#' @usage clusterStack <- get.clusterStack( ks=1:k.clust, force=F, ... ) 
get.clusterStack <- function( ks=1:k.clust, force=F, ... ) {
  if ( ! force && ! is.null( attr( clusterStack, "iter" ) ) && attr( clusterStack, "iter" ) == iter )
    return( clusterStack )
  mc <- get.parallel( length( ks ) )
  clusterStack <- mc$apply( ks, get.clust, ... ) ###, mc.cores=mc$par )
#ifndef PACKAGE
  tmp <- list.reference( clusterStack, file=sprintf( "%s/clusterStack", cmonkey.filename ), type="RDS" )
  clusterStack <- tmp
#endif
  attr( clusterStack, "iter" ) <- iter
  clusterStack
}

#' Get the unpreprocessed ratios matrix (i.e. ratios.raw)
#'  
#' @export
#' @usage ratios <- get.unpreprocessed.ratios()
get.unpreprocessed.ratios <- function( ... ) {
  return( ratios.raw )
}

residual.submatrix <- function( rats, rows, cols, varNorm=F, ... ) {
    rows <- rows[ rows %in% rownames( rats ) ]
    cols <- cols[ cols %in% colnames( rats ) ]
    if ( length( rows ) <= 1 || length( cols ) <= 1 ) return( 1 )
    maxRowVar <- attr( rats, "maxRowVar" )
    rats <- rats[ rows, cols ]
    if ( is.vector( rats ) || any( dim( rats ) <= 1 ) || mean( is.na( rats ) ) > 0.95 ) return( 1 )

    d.rows <- rowMeans( rats, na.rm=T )
    d.cols <- colMeans( rats, na.rm=T )
    d.all <- mean( d.rows, na.rm=T )

    ##rij <- rats + d.all
    ##rij[,] <- rij[,] - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
    ##rij[,] <- rij[,] - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
    ##rij[,] <- rij[,] - outer( d.rows, d.cols, '+' ) ## more elegant but slower!
    rats[,] <- rats[,] + d.all - outer( d.rows, d.cols, '+' ) ## more elegant but slower!

    ##average.r <- mean( abs( rij ), na.rm = TRUE )
    average.r <- mean( abs( rats ), na.rm = TRUE )
    if ( varNorm && ! is.null( maxRowVar ) ) {
        ##maxRowVar <- attr( rats, "maxRowVar" )
        row.var <- mean( apply( rats, 1, var, use = "pairwise.complete.obs" ), na.rm=T )
        if ( is.na( row.var ) || row.var > maxRowVar ) row.var <- maxRowVar
        average.r <- average.r / row.var
    }
    average.r
}

cluster.resid <- function( k, rats.inds="COMBINED", varNorm=F, in.cols=T, ... ) {
  ## FLOC residual number is a good statistic

  inds <- rats.inds
  if ( rats.inds[ 1 ] == "COMBINED" ) inds <- names( get( "row.weights" ) )
  rows <- get.rows( k ); cols <- get.cols( k )
  resids <- sapply( ratios[ inds ], function( rn ) {
#ifndef PACKAGE
    if ( row.score.func == "orig" ) { ## Original FLOC resid score works for co-expressed genes
#endif
      if ( in.cols ) residual.submatrix( rn, rows, cols, varNorm=varNorm )
      else residual.submatrix( rn, get.rows( k ), colnames( rn )[ ! colnames( rn ) %in% cols ], varNorm=varNorm )
#ifndef PACKAGE
    }
    else { ## using cor2 or mi will results in pos- and neg- correlated genes, FLOC score wont work!
      if ( in.cols ) mean( get.row.scores( k, for.rows=rows, ratios=rn, method=row.score.func ) )
      else mean( get.row.scores( k, cols=cols, for.rows=rows, ratios=rn, method=row.score.func ) )
    }
#endif
  } )
  if ( rats.inds[ 1 ] == "COMBINED" ) resids <- weighted.mean( resids, row.weights[ inds ], na.rm=T )
  if ( rats.inds[ 1 ] != "COMBINED" && length( resids ) < length( inds ) && all( is.na( resids ) ) ) {
    resids <- rep( NA, length( inds ) ); names( resids ) <- inds }

  resids
}

cluster.pclust <- function( k, mot.inds="COMBINED" ) { ## actually returns a list with e-vals and p-vals
  inds <- mot.inds
  if ( is.null( inds ) ) return( list( p.clusts=NA, e.vals=NA ) )
  if ( mot.inds[ 1 ] == "COMBINED" ) inds <- names( get( "mot.weights" ) )
  if ( is.null( inds ) ) return( list( p.clusts=NA, e.vals=NA ) )
  rows <- get.rows( k )
  p.clusts <- sapply( inds, function( n ) {
    ms <- meme.scores[[ n ]][[ k ]]
    out <- NA
    if ( is.null( ms ) || is.null( ms$pv.ev ) || length( ms$pv.ev ) <= 0 ) return( out )
    if ( length( rows ) > 0 && ! is.null( ms$pv.ev ) && ! is.null( ms$pv.ev[[ 1 ]] ) ) {
      if ( 'p.value' %in% colnames( ms$pv.ev[[ 1 ]] ) )
        out <- mean( log10( ms$pv.ev[[ 1 ]][ rownames( ms$pv.ev[[ 1 ]] ) %in% rows, "p.value" ] ), na.rm=T )
      else if ( 'pvals' %in% colnames( ms$pv.ev[[ 1 ]] ) )
        out <- mean( log10( ms$pv.ev[[ 1 ]][ rownames( ms$pv.ev[[ 1 ]] ) %in% rows, "pvals" ] ), na.rm=T )
    }
    ## pvs <- meme.scores[[ n ]]$all.pv
    ## if ( is.null( pvs ) ) return( NA )
    ## rows <- rows[ rows %in% rownames( pvs ) ]
    ## if ( length( rows ) > 0 ) mean( log10( pvs[ rows, k ] ), na.rm=T ) else NA
    out
  } )
  
  e.vals <- sapply( inds, function( n ) {
    ms <- meme.scores[[ n ]][[ k ]]
    sapply( 1:length( ms$meme.out ),
           function( i ) if ( length( rows ) > 0 && ! is.null( ms$meme.out ) && ! is.null( ms$meme.out[[ i ]] ) )
           ms$meme.out[[ i ]]$e.value else NA )
  } )
  if ( ! is.matrix( e.vals ) ) e.vals <- t( t( e.vals ) )

  if ( mot.inds[ 1 ] == "COMBINED" ) {
    p.clusts <- weighted.mean( p.clusts, mot.weights[ inds ], na.rm=T )
    e.vals <- apply( e.vals, 1, weighted.mean, mot.weights[ inds ], na.rm=T )
  } else if ( mot.inds[ 1 ] != "COMBINED" ) {
    if ( length( p.clusts ) < length( inds ) && all( is.na( p.clusts ) ) ) p.clusts <- rep( NA, length( inds ) )
    e.vals <- apply( e.vals, 1, function( i ) if ( length( i ) < length( inds ) && all( is.na( i ) ) )
                    rep( NA, length( inds ) ) else i )
  }
  if ( ! is.matrix( e.vals ) ) e.vals <- t( t( e.vals ) )
  if ( is.matrix( e.vals ) && ncol( e.vals ) != length( inds ) ) e.vals <- t( e.vals )
  if ( mot.inds[ 1 ] != "COMBINED" ) names( p.clusts ) <- colnames( e.vals ) <- inds
  else e.vals <- as.vector( e.vals )
  list( p.clusts=p.clusts, e.vals=e.vals )
}

## Create cluster-style cluster object 
get.clust <- function( k, fill=T, fill.motif=T, seq.type=names( mot.weights ), varNorm=F, ... ) { ##"upstream" ) {
  gen.clust <- function( rowNames, colNames=NA, fill=F, motif=F, n.motifs=3, ... ) {
    rowNames <- rowNames[ rowNames %in% attr( ratios, "rnames" ) ]
    if ( ! is.null( colNames ) && length( colNames ) > 1 && ! is.na( colNames ) )
      colNames <- colNames[ colNames %in% attr( ratios, "cnames" ) ]
    c.tmp <- list( nrows=length( rowNames ), ncols=length( colNames ), rows=rowNames, cols=colNames,
                  k=999, p.clust=1.0, e.val=rep( 999, n.motifs ), resid={out=rep(NA,length(row.weights));names(out)<-names(row.weights);resid=out} )
    if ( fill && c.tmp$nrows > 0 && c.tmp$ncols > 0 && ! all( is.na( colNames ) ) )
      c.tmp$resid <- cluster.resid( k, names( row.weights ), varNorm=varNorm, ... ) 
    return( c.tmp )
  }
  
  cols <- get.cols( k ); rows <- get.rows( k )
  if ( length( cols ) <= 0 ) cols <- NA
  clust <- gen.clust( rows, cols, fill=fill, motif=F, n.motifs=max( unlist( n.motifs ) ) ) 
  clust$k <- k

  if ( fill.motif ) {
    tmp <- cluster.pclust( k, seq.type )
    clust$e.val <- tmp$e.vals
    clust$p.clust <- tmp$p.clusts
  }  
  clust
}

## Allow weighting of columns (used in get.row.scores) (this can be "overridden")
get.col.weights <- function( rows, cols, ratios ) {
  NA ## Default: unweighted
##   if ( length( rows ) == 1 ) return( NA ) ## Below is example code weighting columns by their "index of dispersion" over the rows in the cluster
##   rats <- ratios[ rows, cols, drop=F ]
##   vars <- apply( rats, 2, var, na.rm=T )
##   vars <- vars / ( abs( apply( rats, 2, mean, na.rm=T ) ) + 0.05 )
##   vars[ is.na( vars ) | vars == 0 ] <- 1
##   weights <- 1 / vars
##   weights <- weights / sum( weights, na.rm=T ) * length( weights )
##   names( weights ) <- cols
##   weights
}

## Allow weighting of rows (used in get.col.scores) (this can be "overridden")
get.row.weights <- function( rows, cols, ratios ) NA

get.row.scores <- function( k, cols=get.cols( k ), for.rows="all", ratios=ratios[[ 1 ]],
#ifndef PACKAGE
                           method=c("cor2","abscor","cor","dist","orig")[1],
#endif
                           ... ) {
  ## Compute scores for ALL rows (over just the cols IN the cluster)
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( is.null( for.rows ) || for.rows[ 1 ] == "all" ) for.rows <- rownames( ratios ) ##attr( ratios, "rnames" )
  rows <- rows[ rows %in% rownames( ratios ) ]
  cols <- cols[ cols %in% colnames( ratios ) ]
  if ( length( rows ) < 1 || length( cols ) < 1 ) return( rep( NA, length( for.rows ) ) )

#ifndef PACKAGE
  if ( method == "orig" ) { ## ORIGINAL, GOOD, FAST SCORE FUNCTION: mean-squared deviation from mean profile
    ## (note this is practically the same as new method "dist". but it is a lot faster. go figure.)
#endif
    rats <- ratios[ for.rows, cols, drop=F ]
    rats.mn <- colMeans( rats[ rows, , drop=F ], na.rm=T )
    ##rats.mn <- matrix( rats.mn, nrow=nrow( rats ), ncol=ncol( rats ), byrow=T )
    ##rats[,] <- ( rats[,] - rats.mn )^2 ## abs(
    rats[,] <- t( t( rats ) - rats.mn )^2
    col.weights <- if ( exists( "get.col.weights" ) ) get.col.weights( rows, cols, ratios ) else NA
    if ( is.na( col.weights[ 1 ] ) ) rats <- rowMeans( rats, na.rm=T )
    else rats <- apply( rats, 1, weighted.mean, w=col.weights[ cols ], na.rm=T )
    rats <- log( rats + 1e-99 )
#ifndef PACKAGE
  }
  else if ( method == "pval" ) rats <- get.row.scores.pVals( k, cols, rows, ratios, method, ... ) ## Sam's new func.
  else if ( exists( "get.row.scores.NEW" ) ) rats <- get.row.scores.NEW( k, cols, rows, ratios, method, ... )
#endif
  return( rats )
}

#ifndef PACKAGE
## See minet and infotheo packages for information-based distances and implementations like ARACNE
## TODO: how do we include row.weights with the new scores?
get.row.scores.NEW <- function( k, cols=get.cols( k ), rows, ratios=ratios[[ 1 ]],
                           method=c("cor2","abscor","cor","dist")[1], ... ) {
  rats <- ratios[ ,cols, drop=F ]
  if ( method == "dist" || substr( method, 1, 5 ) == "dist." ) { ## Try a more flexible dist() based function -- mean of distances (slower)
    ## if ( method == "dist" ) method <- "dist.euc" ## Default - euclidean
    ## d <- as.matrix( dist( rats, method=substring( method, 6 ) ) ); diag( d ) <- NA
    ## rats <- log( apply( d[ ,rows, drop=F ], 1, mean, na.rm=T, ... ) + 1e-99 )
    require( proxy ) ## Better, faster!!! (10x)
    if ( method == "dist" ) method <- "dist.Euclidean" ## Default - euclidean    
    d <- proxy::dist( rats, rats[ rows, ,drop=F ], method=substring( method, 6 ) ); d[ cbind( rows, rows ) ] <- NA
    rats <- log( apply( d, 1, mean, na.rm=T, ... ) + 1e-99 )
  }
  else if ( method %in% c( "cor", "cor2", "abscor" ) ) { ## Try a cor^2 (allow anticorr) function -- median of cor^2 's
    rats <- t( rats )
    ##d <- cor( rats, use="pairwise", ... ); diag( d ) <- NA
    ##rats <- apply( log( 1 - d[ ,rows, drop=F ] + 1e-99 ), 1, mean, na.rm=T, ... ) ## + 1 + 1e-99 )
    d <- t( cor( rats[ ,rows, drop=F ], rats, use="pairwise", ... ) ); d[ cbind( rows, rows ) ] <- NA ## faster! (10x)
    if ( method == "cor2" ) d <- d^2
    else if ( method == "abscor" ) d <- abs( d )
    rats <- apply( log( 1 - d + 1e-99 ), 1, mean, na.rm=T, ... ) ## + 1 + 1e-99 )
  }
  ## Try mutual information!!! Options: estimator="mi.empirical", "mi.mm", "mi.shrink", "mi.sg" need a non-default
  ##    'discretizer': disc="none", "equalfreq", "equalwidth" or "globalequalwidth" (see infotheo:::discretize)
  ##    estimator="pearson","spearman","kendall" is not discretized (but use "cor2" option instead - same but faster)
  else if ( method == "mi" || substr( method, 1, 3 ) == "mi." ) { 
    require( minet )
    if ( method == "mi" ) method <- "mi.spearman" ## The default estimator
    d <- build.mim( t( rats ), estimator=substring( method, 4 ), disc="equalwidth" ); diag( d ) <- NA ## slow!!!
    rats <- -apply( d[ ,rows, drop=F ], 1, mean, na.rm=T, ... )
  }
  return( rats )
}
#endif

## get.row.scores.cpp <- function( rats, rats.mn ) { ##, rank=rank.ties ) {
##   if ( TRUE && ! is.loaded( "get_row_scores" ) ) {
##     so.file <- "cpp/get_row_scores"
##     if ( ! exists( "is.slave" ) || ! is.slave ) {
##       require( Rcpp ) ## Dont really need to load the package - but need to make sure it's installed
##       code <-
##         "RcppExport SEXP get_row_scores( SEXP rats, SEXP rats_mn ) {
##        SEXP rl = R_NilValue; // Use this when there is nothing to be returned.
##        char* exceptionMesg = NULL; 
##        try {
##          RcppMatrixView<double> r( rats );
##          RcppVectorView<double> rm( rats_mn );
##          int nr = r.dim1(), nc = r.dim2();
##          vector<double> out( nr );
##          vector<int> ni( nr ); //, nj( nc );

##          for ( int j = 0; j < nc; j ++ ) {
##            double rr = rm( j );
##            if ( ISNA( rr ) ) continue;
##            //nj[ j ] = nj[ j ] + 1;
##            for ( int i = 0; i < nr; i ++ ) {
##              double rrr = r( i, j );
##              if ( ISNA( rrr ) ) continue;
##              out[ i ] = out[ i ] + fabs( rrr - rr ); //pow( rrr - rr, 2 );
##              ni[ i ] = ni[ i ] + 1;
##            }
##          }
##          for ( int i = 0; i < nr; i ++ ) out[ i ] = out[ i ] / ni[ i ];

##          RcppResultSet rs; // Build result set to be returned as a list to R.
##          rs.add( \"out\", out );
##          rl = rs.getReturnList(); // Get the list to be returned to R.
##        } catch( std::exception& ex ) { 
##          exceptionMesg = copyMessageToR( ex.what() ); 
##        } catch( ... ) { 
##          exceptionMesg = copyMessageToR( \"unknown reason\" ); 
##        } 
##        if ( exceptionMesg != NULL ) error( exceptionMesg ); 
##        return rl;
##      }"

##       source( "Rcpp_compile3.R", local=T ) ## Loads compile.cpp.code() (requires Rcpp package)
##       compile.rcpp.code( code, fname=so.file, verbose=T )
##     } else {
##       dyn.load( paste( so.file, ".so", sep="" ) )
##     }
##   }
##   out <- .Call( "get_row_scores", rats, rats.mn )$out
##   out[ out == 0 ] <- NA
##   out
## }

## Default is to normalize (diff)^2 by mean expression level, similar to "index of dispersion"
##    http://en.wikipedia.org/wiki/Index_of_dispersion
get.col.scores <- function( k, for.cols="all", ratios=ratios[[ 1 ]], 
#ifndef PACKAGE
                           method=c("new","orig","ent")[1],                            
#endif
                           norm.method=c("mean","all.colVars","none")[1], ... ) {
  ## Compute scores for ALL cols (over just the rows IN the cluster)
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.cols[ 1 ] == "all" ) for.cols <- colnames( ratios ) ##attr( ratios, "cnames" )
  rows <- rows[ rows %in% rownames( ratios ) ]
  if ( length( rows ) <= 1 ) return( rep( NA, length( for.cols ) ) )

  rats <- ratios[ rows, for.cols, drop=F ]
  row.weights <- if ( exists( "get.row.weights" ) ) get.row.weights( rows, cols, ratios ) else NA

#ifndef PACKAGE
  if ( method == "orig" ) { ## OLD, GOOD EUCLIDEAN DISTANCE METRIC
#endif
    if ( is.na( row.weights[ 1 ] ) ) { ## Default
      rats.mn <- ##matrix(
        colMeans( rats, na.rm=T )##, nrow=nrow( rats ), ncol=ncol( rats ), byrow=T )
    } else { ## Custom row weights
      rats.mn <- ##matrix(
        apply( rats, 2, weighted.mean, w=row.weights[ rows ], na.rm=T )##, ncol=ncol( rats ), byrow=T )
    }
    
    rats[,] <- t( t( rats ) - rats.mn )^2 ## abs( Multiplying is faster than squaring
    rats <- colMeans( rats, na.rm=T )
#ifndef PACKAGE
  }

  ## else if ( method == "ent" ) {  ## ENTROPY BASED METHOD (do it on abs(rats) to allow for repression too
  ##   require( infotheo )
  ##   rats <- apply( discretize( rats, disc="globalequalwidth", nbins=ncol( ratios ) / 2 ), 2, entropy )
  ## }
  
  else if ( method == "new" ) { ## TRY -- score columns on amount by which they contribute to avg in-row score,
    ###  however that score is calculated! Note - this shouldn't have to worry about the varNorm anymore (???)
    mn <- mean( get.row.scores( k, ratios=ratios, method=row.score.func, ... ), na.rm=T )
    rows <- get.rows( k )
    cols <- get.cols( k )
    rats <- -sapply( for.cols, function( cc ) {
      ##if(row.score.func!="orig")cat(k,cc,"\n")
      if ( is.na( row.weights[ 1 ] ) ) {
        if ( cc %in% cols ) mn / mean( get.row.scores( rows, cols=cols[ cols != cc ], ratios=ratios,
                                                      method=row.score.func, ... ), na.rm=T ) ## - mn
        ##else mean( mn -
        else mean( get.row.scores( k, cols=c( cols, cc ), ratios=ratios, method=row.score.func, ... ), na.rm=T ) / mn
      } else {
        if ( cc %in% cols ) mn / weighted.mean( get.row.scores( rows, cols=cols[ cols != cc ], ratios=ratios,
                                                  method=row.score.func, ... ), w=row.weights, na.rm=T ) ## - mn
        ##else weighted.mean( mn -
        else weighted.mean( get.row.scores( rows, cols=c( cols, cc ), ratios=ratios, method=row.score.func, ... ),
                           w=row.weights, na.rm=T ) / mn
      }
    } )
    ##rats <- log( rats + 1e-99 )
    return( rats )
}
#endif
  
  var.norm <- 0.99
  if ( norm.method == "all.colVars" ) {
    all.colVars <- attr( ratios, "all.colVars" )
    if ( ! is.null( all.colVars ) ) var.norm <- all.colVars[ for.cols ]
  } else if ( norm.method == "mean" ) {
#ifndef PACKAGE
    if ( ! exists( "rats.mn" ) ) {
      row.weights <- if ( exists( "get.row.weights" ) ) get.row.weights( rows, cols, ratios ) else NA
      rats.tmp <- ratios[ rows, for.cols, drop=F ]
      if ( is.na( row.weights[ 1 ] ) ) { ## Default
        rats.mn <- ##matrix(
          colMeans( rats.tmp, na.rm=T )##, nrow=nrow( rats.tmp ), ncol=ncol( rats.tmp ), byrow=T )
      } else { ## Custom row weights
        rats.mn <- ##matrix(
          apply( rats.tmp, 2, weighted.mean, w=row.weights[ rows ], na.rm=T )##, ncol=ncol( rats.tmp ), byrow=T )
      }
    }
#endif
    var.norm <- abs( rats.mn ) ##[ 1, ] ) ##0.99 ## Use the mean expr. level (higher expressed expected to have higher noise)
  }
  
  ##col.weights <- get.col.weights( rows, cols )
  ##if ( is.na( col.weights ) )
  rats <- rats / ( var.norm + 0.01 ) ## default
  ##else rats <- colMeans( rats, na.rm=T ) / ( var.norm * col.weights[ cols ] + 0.01 ) ## customized col. weights
  ##return( log( rats + 1e-99 ) )
  rats
}

get.motif.scores <- function( k, meme.scores, seq.type="upstream meme", for.rows="all" ) { ##m=meme.scores$upstream[[ k ]], for.rows="all" ) { 
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- attr( ratios, "rnames" )
  if ( length( rows ) <= 1 || is.null( meme.scores[[ seq.type ]]$all.pv ) ) ##is.null( m ) || is.null( m$pv.ev ) )
    return( rep( NA, length( for.rows ) ) )
  ##m.scores <- log( m$pv.ev[[ 1 ]][ ,"p.value" ] ); names( m.scores ) <- rownames( m$pv.ev[[ 1 ]] )
  m.scores <- log( meme.scores[[ seq.type ]]$all.pv[ ,k ] )
  m.scores <- m.scores[ for.rows ] ##; names( m.scores ) <- for.rows
  return( m.scores )
}

## TODO: this could be a lot faster (?) if we used data.table rather than data.frame
get.network.scores <- function( k, net=networks$string, for.rows="all", p1.col="protein1", p2.col="protein2", 
                               score.col="combined_score", combine.func=sum ) {
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- attr( ratios, "rnames" )
  if ( length( rows ) < 1 ) return( rep( NA, length( for.rows ) ) )
  cons <- net[ as.character( net[[ p1.col ]] ) %in% rows, c( p2.col, score.col ), drop=F ] ## Note subset is slower!
  if ( is.null( cons ) || nrow( cons ) <= 0 ) return( rep( NA, length( for.rows ) ) )
  cons <- cons[ as.character( cons[[ p2.col ]] ) %in% for.rows, , drop=F ]
  ##cons <- get.cluster.network( rows, net, for.rows )
  if ( is.null( cons ) || nrow( cons ) <= 0 ) return( rep( NA, length( for.rows ) ) )
  tmp <- tapply( as.numeric( cons[[ score.col ]] ), as.character( cons[[ p2.col ]] ), combine.func, na.rm=T ) /
    length( rows )
  scores <- rep( NA, length( for.rows ) ); names( scores ) <- for.rows
  scores[ names( tmp ) ] <- tmp
  return( -log( scores + 1 ) )
}

pareto.adjust.weights <- function( kwin=51, diter=21, max.delta=0.05 ) {
  out.scaling <- c( row=row.scaling[ iter - 1 ], mot=mot.scaling[ iter - 1 ], net=net.scaling[ iter - 1 ] )
  orig.scaling <- c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ] )
  if ( iter < diter/2 ) return( out.scaling )
  diter <- min( diter, iter )
  kwin <- min( kwin, iter )
  avs <- NULL
  for ( i in names( out.scaling ) ) {##c( "row", "mot", "net" ) ) {
    if ( i == "row" ) col <- "row.scores" ##"resid"
    else if ( i == "mot" ) col <- "mot.scores" ##"p.clust"
    else col <- paste( i, "scores", sep="." )
    tmp <- stats[ ,col ]; tmp <- tmp[ ! is.na( tmp ) ]
    if ( length( tmp ) > 0 ) {
      av <- runmed( tmp, k=kwin )
      dy <- ( av[ length( av ) ] - av[ max( 1, length( av ) - diter ) ] ) / diff( range( av, na.rm=T ) )
      if ( ! is.na( dy ) ) out.scaling[ i ] <- out.scaling[ i ] + min( 0.1, dy )
      out.scaling[ i ] <- max( out.scaling[ i ], orig.scaling[ i ] )
    }
  }
  out.scaling
}

pareto.adjust.weights.OLD <- function( iter, delta.iter=200, delta.factor=1, n.avg=50, max.delta=0.05 ) { ##0.01
  if ( iter == 1 ) return( c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ] ) )
  out.scaling <- c( row=row.scaling[ iter - 1 ], mot=mot.scaling[ iter - 1 ], net=net.scaling[ iter - 1 ] )
  if ( iter < delta.iter + n.avg + 10 ) return( out.scaling )
  all.diffs <- numeric()
  for ( i in c( "row", "mot", "net" ) ) {
    if ( i == "row" ) col <- "resid"
    else if ( i == "mot" ) col <- "p.clust"
    else col <- paste( i, "scores", sep="." )
    tops <- stats[[ col ]][ (iter-n.avg+1):iter ]
    bots <- stats[[ col ]][ (iter-n.avg+1):iter - delta.iter ]
    all.diffs[ i ] <- mean( tops - bots, na.rm=T ) / abs( diff( range( stats[[ col ]], na.rm=T ) ) )
  }
  for ( i in which( all.diffs > 0 & ! is.na( all.diffs ) & ! is.na( out.scaling ) & out.scaling > 0 ) )
    out.scaling[ names( all.diffs )[ i ] ] <- out.scaling[ names( all.diffs )[ i ] ] +
      min( c( out.scaling[ i ] * max.delta, all.diffs[ i ] * delta.factor ) )
  out.scaling
}

seed.clusters <- function( k.clust, seed.method="rnd", col.method="rnd" ) {
  ## Allow it to be overridden by a custom function if desired (e.g. to seed from a gene list -- no that's a bad
  ## example -- there's the "list=" option below)
  if ( seed.method == "custom" && exists( "seed.clusters.custom" ) )
    return( seed.clusters.custom( k.clust, col.method ) )
  if ( substr( seed.method, 1, 3 ) == "net" && length( networks ) <= 0 ) {
    cat( "Seed method is", seed.method, ", but no networks -- using 'kmeans' instead.\n" )
    seed.method <- "kmeans"
  }
  if ( seed.method == "rnd" ) { ## RANDOM, ALL GENES ASSIGNED TO k.clust CLUSTERS
    rm <- t( sapply( 1:attr( ratios, "nrow" ),
                                function( i ) sample( 1:k.clust, n.clust.per.row[ 1 ],
                                                     replace=n.clust.per.row[ 1 ] > attr( ratios, "nrow" ) ) ) )
  } else if ( substr( seed.method, 1, 5 ) == "list=" ) { ## File with sets of genes - one set of genes per line
    rm <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] ) ## separated by space,tab,or comma
    rownames( rm ) <- attr( ratios, "rnames" )
    fname <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ]
    if ( exists( fname ) ) lists <- get( fname ) ## List of vectors already exists in memory
    else if ( file.exists( fname ) ) lists <- strsplit( readLines( fname ), split="[ \t,]", perl=T )
    for ( k in 1:min( c( k.clust, length( lists ) ) ) ) {
      probes <- unlist( lapply( get.synonyms( lists[[ k ]] ), function( i ) i %in% rownames( rm ) ) )
      rm[ probes[ rm[ probes, 1 ] == 0 ], 1 ] <- k
      rm[ probes[ rm[ probes, 1 ] != 0 ], 2 ] <- k
    }
    if ( length( lists ) < k.clust ) { ## Fill in additional clusters randomly
      for ( k in ( length( lists ) + 1 ):k.clust ) {
        rnames <- attr( ratios, "rnames" )[ ! attr( ratios, "rnames" ) %in% unlist( lists ) ]
        rows <- sample( rnames, 5 )
        rm[ rows[ rm[ rows, 1 ] == 0 ], 1 ] <- k
        rm[ rows[ rm[ rows, 1 ] != 0 & rm[ rows, 2 ] == 0 ], 2 ] <- k
      }
    }
  } else if ( substr( seed.method, 1, 4 ) == "rnd=" ) { ## RANDOM SAMPLED N per cluster
    n.samp <- as.integer( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    rm <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=n.clust.per.row[ 1 ] )
    rownames( rm ) <- attr( ratios, "rnames" )
    for ( i in 1:n.clust.per.row ) {
      sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
      for ( k in 1:k.clust ) {
        g <- sample( attr( ratios, "rnames" )[ ! sampled ], n.samp )
        rm[ g, 1 ] <- k
        sampled[ g ] <- TRUE
      }
    }
  } else if ( seed.method == "kmeans" ) { ## Kmeans seeded
    if ( ! exists( "ratios" ) ) stop( "kmeans seed method but no ratios" )
    tmp.rat <- get.cluster.matrix() ##ratios;
    tmp.rat[ is.na( tmp.rat ) ] <- 0
    ##cat(dim(tmp.rat),k.clust,"\n")
    km <- kmeans
#ifndef PACKAGE
    if ( ! is.na( parallel.cores ) && ( parallel.cores == TRUE || parallel.cores > 1 ) ) {
        ##&& require( biganalytics ) && require( doMC ) )
        get.parallel(); if ( require( biganalytics ) ) km <- bigkmeans }
#endif
    tmp.rat[ is.na( tmp.rat ) ] <- rnorm( sum( is.na( tmp.rat ) ) )
    rm <- km( tmp.rat, centers=k.clust, iter.max=20, nstart=2 )$cluster
    names( rm ) <- attr( ratios, "rnames" )
    if ( n.clust.per.row[ 1 ] > 1 ) rm <-
      cbind( rm, matrix( rep( 0, attr( ratios, "nrow" ) * ( n.clust.per.row[ 1 ] - 1 ) ),
                                         ncol=n.clust.per.row[ 1 ] - 1 ) )
#ifndef PACKAGE
  } else if ( substr( seed.method, 1, 11 ) == "trimkmeans=" ) { ## Trimmed k-means seeded
    if ( ! exists( "ratios" ) ) stop( "trimkmeans seed method but no ratios" )
    require( trimcluster )
    trim <- as.numeric( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    tmp.rat <- get.cluster.matrix() ##ratios;
    tmp.rat[ is.na( tmp.rat ) ] <- 0
    rm <- trimkmeans( tmp.rat, k.clust, trim=trim, maxit=20, runs=2 )$classification
    if ( n.clust.per.row[ 1 ] > 1 ) rm <-
      cbind( rm, matrix( rep( 0, attr( ratios, "nrow" ) * ( n.clust.per.row[ 1 ] - 1 ),
                                         ncol=n.clust.per.row[ 1 ] - 1 ) ) )
  } else if ( substr( seed.method, 1, 4 ) == "cor=" ) { ## RANDOM GENES AND THEIR N CORRELATED GENES
    if ( ! exists( "ratios" ) ) stop( "cor seed method but no ratios" )
    n.cor <- as.integer( strsplit( seed.method, "=" )[[ 1 ]][ 2 ] )
    rats <- get.cluster.matrix()
    cors <- if ( attr( ratios, "nrow" ) < 6000 ) cor( t( rats ), use="pairwise" ) else NULL
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) < n.cor ) sampled[ sample( 1:length( sampled ) ) ] <- FALSE
        rnames <- attr( ratios, "rnames" )[ ! sampled ]
        g <- sample( rnames, 1 )
        if ( ! is.null( cors ) ) g <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:n.cor ] ]
        else g <- rnames[ order( apply( rats[ ! sampled, ], 1, cor, rats[ g, ] ), decreasing=T )[ 1:n.cor ] ]
        rm[ g ] <- k
        if ( length( g ) == 1 ) { ## if option was cor=1, try to spread the wealth by preventing genes correlated 
          if ( ! is.null( cors ) ) tmp <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:10 ] ] ## with this one from getting
          else tmp <- rnames[ order( apply( rats[ ! sampled, ], 1, cor, rats[ g, ] ), decreasing=T )[ 1:10 ] ]
          sampled[ tmp ] <- TRUE ## added to other clusters
        }
        sampled[ g ] <- TRUE
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    rm <- do.call( cbind, tmp )
  } else if ( substr( seed.method, 1, 4 ) == "net=" ) { ## Network seeded; e.g. 'net=string:5'
    if ( ! exists( "networks" ) || length( networks ) <= 0 ) stop( "net seed method but no networks" )
    seed.method <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ] ## Name of network to use and # of seed nodes
    net.name <- strsplit( seed.method, ":" )[[ 1 ]][ 1 ]
    net <- networks[[ net.name ]]
    n.seed <- as.integer( strsplit( seed.method, ":" )[[ 1 ]][ 2 ] )
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, length( unique( as.character( net$protein1 ) ) ) )
    names( sampled ) <- unique( as.character( net$protein1 ) )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) <= 0 ) sampled[ 1:length( sampled ) ] <- FALSE
        rnames <- names( which( ! sampled ) )
        gs <- sample( rnames, 1 )
        qiter <- 0
        while( length( gs ) < n.seed && qiter < 20 ) {
          ns <- as.character( net$protein2[ as.character( net$protein1 ) %in% gs ] )
          if ( length( ns ) + length( gs ) >= n.seed ) ns <- sample( ns, size=n.seed-length( gs ),
                                      prob=net$combined_score[ as.character( net$protein1 ) %in% gs ] )
          gs <- unique( c( gs, ns ) )
          qiter <- qiter + 1
        }
        rm[ gs ] <- k
        sampled[ gs ] <- TRUE  ## try to spread the wealth by preventing genes connected to these from being added 
        if ( n.seed <= 2 ) sampled[ as.character( net$protein2[ as.character( net$protein1 ) %in% gs ] ) ] <- TRUE ## to other clusts
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    rm <- do.call( cbind, tmp )
  } else if ( substr( seed.method, 1, 7 ) == "netcor=" ) { ## Correld net-weighted seeds, e.g. 'netcor=string:5'
    if ( ! exists( "ratios" ) ) stop( "netcor seed method but no ratios" )
    if ( ! exists( "networks" ) || length( networks ) <= 0 ) stop( "netcor seed method but no networks" )
    seed.method <- strsplit( seed.method, "=" )[[ 1 ]][ 2 ] ## Name of network to use and # of seed nodes
    net.name <- strsplit( seed.method, ":" )[[ 1 ]][ 1 ]
    net <- networks[[ net.name ]]
    n.seed <- as.integer( strsplit( seed.method, ":" )[[ 1 ]][ 2 ] )
    rats <- get.cluster.matrix()
    cors <- cor( t( rats ), use="pairwise" )
    tmp.mat <- matrix( 0, nrow=nrow( cors ), ncol=ncol( cors ) ); dimnames( tmp.mat ) <- dimnames( cors )
    tmp.lookup <- 1:attr( ratios, "nrow" ); names( tmp.lookup ) <- attr( ratios, "rnames" )
    net <- net[ as.character( net$protein1 ) %in% attr( ratios, "rnames" ) &
               as.character( net$protein2 ) %in% attr( ratios, "rnames" ), ]
    tmp.mat[ cbind( tmp.lookup[ as.character( net$protein1 ) ], tmp.lookup[ as.character( net$protein2 ) ] ) ] <-
      net$combined_score / 1000
    cors <- cors + tmp.mat; rm( tmp.mat )
    rm <- rep( 0, attr( ratios, "nrow" ) ); names( rm ) <- attr( ratios, "rnames" )
    sampled <- rep( FALSE, attr( ratios, "nrow" ) ); names( sampled ) <- attr( ratios, "rnames" )
    mc <- get.parallel( n.clust.per.row )
    tmp <- mc$apply( 1:n.clust.per.row, function( i ) {
      for ( k in 1:k.clust ) {
        if ( sum( ! sampled ) < n.seed ) sampled[ sample( 1:length( sampled ) ) ] <- FALSE
        rnames <- attr( ratios, "rnames" )[ ! sampled ]
        g <- sample( rnames, 1 )
        g <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:n.seed ] ]
        rm[ g ] <- k
        if ( length( g ) == 1 ) { ## if option was cor=1, try to spread the wealth by preventing genes correlated 
          tmp <- rnames[ order( cors[ g, ! sampled ], decreasing=T )[ 1:10 ] ] ## with this one from getting 
          sampled[ tmp ] <- TRUE ## added to other clusters
        }
        sampled[ g ] <- TRUE
      }
      rm[ attr( ratios, "rnames" ) ]
    } ) 
    rm <- do.call( cbind, tmp )
#endif
  }
  
  if ( is.vector( rm ) ) rm <- t( rm ) ## probably n.clust.per.row == 1
  if ( nrow( rm ) == 1 ) rm <- t( rm ) ## probably n.clust.per.row == 1
 
  if ( col.method == "rnd" ) {
    cm <- t( sapply( 1:attr( ratios, "ncol" ), function( i )
                                sample( 1:k.clust, n.clust.per.col[ 1 ],
                                       replace=n.clust.per.col[ 1 ] > k.clust ) ) ) ##attr( ratios, "ncol" ) ) ) )
  } else if ( col.method == "best" ) {
    if ( ! exists( "ratios" ) ) stop( "best col seed method but no ratios" )
    all.rats <- get.cluster.matrix()
    attr( all.rats, "all.colVars" ) <- apply( all.rats, 2, var, use="pair", na.rm=T )
    col.scores <- -sapply( 1:k.clust, function( k )
                          if ( sum( rm == k, na.rm=T ) <= 0 ) rep( NA, attr( ratios, "ncol" ) ) else 
                          get.col.scores( k=rownames( which( rm == k, arr=T ) ),
                                         ratios=all.rats, method="orig" ) ) 
    cm <- t( apply( col.scores, 1, function( i ) order( i, decreasing=T )[ 1:n.clust.per.col[ 1 ] ] ) )
  }

  rownames( rm ) <- attr( ratios, "rnames" ) 
  rownames( cm ) <- attr( ratios, "cnames" )
  list( rm=rm, cm=cm )
}  

## SEEDING -- seed new clusters if no clusterStack exists.
##    ALSO re-set the random seed if env$cmonkey.params$rnd.seed doesn't exist (e.g. for multiple runs on the
##       same pre-initialized env.) -- need to this via:
##       rm(list="rnd.seed",envir=env$cmonkey.params);rm(list="rnd.seed",envir=env)
cmonkey.re.seed <- function( env ) {
  if ( ! exists( "rnd.seed", envir=env$cmonkey.params ) ) {
    ##tmp.rnd.seed <- as.integer( Sys.time() ) ## This will only actually change the seed if
    op <- options( digits.secs=10 )
    tmp.time <- as.character( Sys.time() )
    options( op ); rm( op )
    tmp.rnd.seed <- as.integer( substr( gsub( '[-:. ]', "", tmp.time ), 12, 20 ) )
    cat( "RESETTING RANDOM SEED: " ) ## env$rnd.seed and env$cmonkey.params$rnd.seed are GONE
    env$set.param( "date.run", tmp.time, env$cmonkey.params ) 
    env$date.run <- env$cmonkey.params$date.run
    env$set.param( "rnd.seed", tmp.rnd.seed, env$cmonkey.params ) 
    env$rnd.seed <- env$cmonkey.params$rnd.seed
    set.seed( env$rnd.seed )
    rm( tmp.rnd.seed )
    ##env$set.param( "cmonkey.filename", paste( "cmonkey", env$cmonkey.version, env$organism,
    ##                                         paste( sapply( env$ratios, dim ), collapse="x" ),
    ##                                         gsub( " ", "_", env$date.run ), sep="_" ), env$cmonkey.params )
    ##env$cmonkey.filename <- env$cmonkey.params$cmonkey.filename
  }
  if ( ! is.null( env$ratios ) && attr( env$ratios, "ncol" ) > 1 ) {
    cat( "Seeding all clusters using methods:", env$seed.method, "\n" ) ##seed.method.cols, "\n" )
    tmp <- env$seed.clusters( env$k.clust, seed.method=env$seed.method[ "rows" ], col.method=env$seed.method[ "cols" ] ) ##.cols )
  } else {
    cat( "Seeding all clusters using methods: rnd rnd\n" )
    tmp <- env$seed.clusters( env$k.clust, seed.method="rnd", col.method="rnd" )
  }
  env$clusterStack <- lapply( 1:env$k.clust, function( k )
                             list( rows=rownames( which( tmp$rm == k, arr=T ) ),
                                  cols=rownames( which( tmp$cm == k, arr=T ) ) ) )
  attr( env$clusterStack, "iter" ) <- env$iter - 1
  invisible( env )
}

get.cluster.matrix <- function( rows=NULL, cols=NULL, matrices=names( ratios ) ) {
  if ( is.null( rows ) ) rows <- attr( ratios, "rnames" )
  if ( is.null( cols ) ) cols <- attr( ratios, "cnames" )
  cols.b <- attr( ratios, "cnames" )[ attr( ratios, "cnames" ) %in% cols ]

  rats <- matrix( NA, nrow=length( rows ), ncol=length( cols.b ) )
  rownames( rats ) <- rows; colnames( rats ) <- cols.b
  cnames <- character()
  for ( n in matrices ) { ## names( ratios ) ) {
    r.tmp <- ratios[[ n ]][ rows[ rows %in% rownames( ratios[[ n ]] ) ],
                           cols.b[ cols.b %in% colnames( ratios[[ n ]] ) ], drop=F ]
    if ( is.null( r.tmp ) || all( is.na( r.tmp ) ) ) next
    if ( is.vector( r.tmp ) ) { r.tmp <- t( r.tmp ); rownames( r.tmp ) <- rows }
    cnames <- c( cnames, colnames( r.tmp ) )
    rats[ rownames( r.tmp ), colnames( r.tmp ) ] <- r.tmp; rm( r.tmp )
  }
  rats[ ,colnames( rats ) %in% cnames, drop=F ]
}

## Get all potential synonyms from feature.names table and possibly from "translation.tab" if it exists
get.synonyms <- function( gns, ft=genome.info$feature.names, ignore.case=T, verbose=F, fast=F, force=F ) { ##transl.table
  gns.input <- gns
  if ( exists( "no.genome.info" ) && no.genome.info ) { out <- as.list( gns ); names( out ) <- gns; return( out ) }
  out <- list()
  if ( ( ! force && exists( "genome.info" ) && ! is.null( genome.info$synonyms ) ) ) {
    gns.cached <- gns[ gns %in% names( genome.info$synonyms ) ]
    out <- genome.info$synonyms[ gns.cached ]
    gns <- gns[ ! gns %in% names( genome.info$synonyms ) ]
    if ( length( gns ) <= 0 || ( is.null( ft ) && ( ! exists( "translation.tab" ) || is.null( translation.tab ) ) ) )
      return( out )
  }
  ##cat( "Generating gene synonyms cache... this could take a little while.\n" )
  tmp.out <- as.list( gns ); names( tmp.out ) <- gns
  if ( is.null( ft ) && ( ! exists( "translation.tab" ) || is.null( translation.tab ) ) ) return( c( out, tmp.out ) )
  gns.orig <- gns
  gns <- gsub( "m$|_\\d$|\\-\\S$", "", gns, perl=T ) ## specific for Halo, also works for yeast
  gns <- gsub( "([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\\\\\1", gns, perl=T ) ## Escape the regex-y characters
  gns <- gns[ ! is.na( gns ) & gns != "" ]
  ft <- ft[ ,c( "id", "names" ) ]
  if ( exists( "translation.tab" ) && ! is.null( translation.tab ) )
    ft <- rbind( ft, data.frame( id=as.character( translation.tab$V1 ), names=as.character( translation.tab$V2 ) ) )
  ft <- subset( ft, names != "" )
  
  if ( verbose ) ggggg <- gns[ seq( 1, length( gns ), by=min(length(gns),100) ) ]

  mc <- get.parallel( length( gns ), verbose=F )
  tmp <- mc$apply( gns, function( g ) {
    if ( verbose && g %in% ggggg ) cat( " ...", g )
    greg <- paste( "^", g, "$", sep="" )
    tmp <- subset( ft, grepl( greg, id, perl=T, ignore=ignore.case ) |
                  grepl( greg, names, perl=T, ignore=ignore.case ) ) ## ) )
    if ( nrow( tmp ) <= 0 ) return( g ) ##tmp.out ) ##character() )
    tmp2 <- unique( c( g, as.character( tmp$id ), as.character( tmp$names ) ) )
    if ( ! fast ) {
      tmp2 <- subset( ft, id %in% tmp2 | names %in% tmp2 )
      tmp2 <- unique( c( g, as.character( tmp2[ ,1 ] ), as.character( tmp2[ ,2 ] ) ) )
    }
    tmp2 <- gsub( "\\\\([\\[\\]\\(\\)\\{\\}\\.\\+\\-'\"])", "\\1", tmp2, perl=T )
    unique( tmp2 )
  } ) 

  names( tmp ) <- gns.orig
  if ( verbose ) cat( "\n" )
  out <- c( tmp, out )
  out[ gns.input[ gns.input %in% names(out) ] ]
}

## Get id's from feature.names table possibly translating via "translation.tab" or other rules if necessary
## get.ids <- function( rows, ft=genome.info$feature.names, ... ) {
##   ids <- subset( ft, grepl( rows, names, ignore=T, perl=T ), c( id, names ) )
##   if ( exists( "translation.tab" ) && ! is.null( translation.tab ) &&
##       nrow( ids ) < length( rows ) ) { ## Use translation table first
##     rows2 <- as.character( subset( translation.tab, V1 %in% rows, V2 )[ ,1 ] )
##     ids <- rbind( ids, subset( ft, names %in% rows2, c( id, names ) ) )
##   }
##   unique( ids )
## }

get.gene.coords <- function( rows, op.shift=T, op.table=genome.info$operons, ... ) { ##, override=F ) {
  if ( is.null( rows ) ) {
    if ( organism == 'hal' ) rows <- grep( "^VNG", as.character( genome.info$feature.tab$canonical_Name ), perl=T, val=T )
    else if ( organism == 'sce' ) rows <- grep( e$genome.info$gene.regex, as.character( genome.info$feature.tab$id ), perl=T, val=T )
    else rows <- grep( "^NP_", as.character( genome.info$feature.tab$id ), perl=T, val=T )
  }
  rows <- unique( rows )
  syns <- get.synonyms( rows, ... )
  tab <- genome.info$feature.tab
  ##if ( "chrom_position" %in% colnames( tab ) ) tab <- tab[ ! duplicated( tab$chrom_position ), ]
  ids <- lapply( syns, function( s ) s[ s %in% tab$id ] )
  if ( all( sapply( ids, length ) < 1 ) ) {
    warning( "Could not find gene start/stop for any input genes", call.=F ); return( NULL ) }
  if ( any( sapply( ids, length ) < 1 ) ) warning( "Could not find gene start/stop for all input genes", call.=F )
  ids <- ids[ sapply( ids, length ) >= 1 ]
  ids <- sapply( ids, "[", 1 )
  ids <- data.frame( id=ids, names=names( ids ) ) ##, gene=rows )
  
  coos <- NULL
  if ( op.shift && exists( 'op.table' ) ) { ## Replace each gene's coords with the head of its operon, but keep the "locus_tag".
    if ( attr( op.table, "source" ) == "rsat" ) {
      ops <- merge( ids, op.table, by.x="id", by.y="query", all=F, sort=F )
      ops2 <- ops[ order( ops$lead ), ]
      coos <- merge( ops, tab, by.x="lead", by.y="name", all=F, sort=F )[ ,c( "id.x", "names", "contig",
                                                                               "strand", "start_pos", "end_pos" ) ]
      ## Due to not-unique mapping of short names in operons table to coords, we need to get rid of dupes
      ##    in an intelligent way:
##       if ( FALSE ) { ## BIG TODO!!!!
##         dupes <- as.character( coos$names[ which( duplicated( coos$names ) ) ] )
##         syns <- get.synonyms( dupes )
##         for ( s in names( syns ) ) {
##           tmp <- subset( genome.info$feature.tab, name %in% syns[[ s ]] )
##           tmp2 <- subset( coos, id %in% tmp$id )
##         }
##       }
    } else if ( attr( op.table, "source" ) == "microbes.online" ) {
      ops <- NULL
      if ( ! any( ids$names %in% op.table$gene ) ) {
        ids2 <- lapply( syns, function( s ) s[ s %in% op.table$gene ] )
        if ( all( sapply( ids2, length ) < 1 ) ) {
          warning( "Could not find operon info for any input genes", call.=F )
        } else {
          if ( any( sapply( ids2, length ) < 1 ) ) warning( "Could not find operon info for all input genes", call.=F )
          ids2 <- ids2[ sapply( ids2, length ) >= 1 ]
          ids2 <- sapply( ids2, "[", 1 )
          ids2 <- data.frame( id=ids2, names=names( ids2 ) ) ##, gene=rows )
          if ( any( ! rows %in% ids2$names ) ) {
            capture.output( prefix <- get.gene.regex( as.character( ids2$id ) )[ 1 ] )
            missing <- rows[ ! rows %in% ids2$names ]
            syns2 <- sapply( syns[ missing ], function( i ) grep( sprintf( '^%s', prefix ), i, val=T )[ 1 ] )
            syns2[ is.na( syns2 ) ] <- names( syns2 )[ is.na( syns2 ) ]; names( syns2 ) <- NULL
            ids2 <- rbind( ids2, data.frame( id=syns2, names=missing ) )
          }
          ops <- merge( ids2, op.table, by.x="id", by.y="gene", all.x=T, incomparables=NA, sort=F )
        }
      }
      if ( is.null( ops ) ) #ops <- merge( ids, op.table, by.x="names", by.y="gene", all.x=T, incomparables=NA, sort=F )
          ops <- op.table[ ids ]
      if ( any( is.na( ops$head ) ) ) {
        head <- as.character( ops$head ); head[ is.na( head ) ] <- as.character( ops$id[ is.na( head ) ] ) ##names[ is.na( head ) ] )
        ops$head <- as.factor( head )
      }
      head.genes <- unique( as.character( ops$head ) ); head.genes <- head.genes[ ! is.na( head.genes ) ]
      head.syns <- get.synonyms( head.genes )
      head.ids <- lapply( head.syns, function( s ) s[ s %in% tab$id ] )
      head.ids <- head.ids[ sapply( head.ids, length ) >= 1 ]
      head.ids <- data.frame( id=sapply( head.ids, "[", 1 ), names=names( head.ids ) )
      ops2 <- merge( ops, head.ids, by.x="head", by.y="names", all.x=T, sort=F )
      coos <- merge( ops2, tab, by.x="id.y", by.y="id", all.x=T, sort=F )[ ,c( "id.x", "names", "contig",
                                                                     "strand", "start_pos", "end_pos" ) ]
    }
  } else { ##if ( ! op.shift )
    coos <- merge( ids, tab, by="id", sort=F )[ ,c( "id", "names", "contig", "strand", "start_pos", "end_pos" ) ]
  }
  colnames( coos )[ 1 ] <- "id"
  if ( is.factor( coos$start_pos ) ) coos$start_pos <- as.numeric( levels( coos$start_pos ) )[ coos$start_pos ]
  if ( is.factor( coos$end_pos ) ) coos$end_pos <- as.numeric( levels( coos$end_pos ) )[ coos$end_pos ]
  coos[ ! duplicated( coos[ ,1:4 ] ), ]
}  

get.long.names <- function( k, shorter=F ) {
  if ( is.numeric( k[ 1 ] ) ) { rows <- get.rows( k ) }
  else { rows <- k }
  if ( is.null( genome.info$feature.tab ) ) {
    out <- rep( "", length( rows ) ); names( out ) <- rows; return( rows ) }
  ids <- get.synonyms( rows )
  mc <- list( apply=lapply )
  if ( ! shorter ) desc <- mc$apply( ids, function( i ) subset( genome.info$feature.tab, id %in% i, 
                                                               select=c( "id", "description" ) ) )
  else {
    desc <- mc$apply( ids, function( i ) subset( genome.info$feature.tab, id %in% i,
                                                select=c( "id", "name", "description" ) ) )
    for ( i in 1:length( desc ) ) if ( length( desc[[ i ]]$name ) > 0 && desc[[ i ]]$name %in% rows ) {
      if ( grepl( "(", desc[[ i ]]$description, fixed=T ) ) ## Try to parse out short name from description
        desc[[ i ]]$name <- strsplit( as.character( desc[[ i ]]$description ), "[()]", perl=T )[[ 1 ]][ 2 ]
    }
  }
  out <- sapply(desc, function(i) as.character(i[1, 2]))
  out <- out[ rows ]
  names( out ) <- rows
  out[ is.na( out ) | out == names( out ) ] <- "" 
  out
}

get.gene.regex <- function( names=NULL, verbose=F ) {
  if ( ! is.null( names ) ) tmp <- names
  else { ## Get common prefix from feature.names and use those genes (assume >40% of ORF names have this suffix)
    if ( exists( 'ratios' ) && ! is.null( ratios ) ) tmp <- toupper( attr( ratios, "rnames" ) )
    else if ( exists( 'genome.info' ) && ! is.null( genome.info$feature.names ) ) {
      tmp <- toupper( subset( genome.info$feature.names, type == "primary", select="names", drop=T ) )
      if ( exists( 'ratios' ) && ! is.null( ratios ) )
        tmp <- tmp[ toupper( tmp ) %in% toupper( attr( ratios, "rnames" ) ) ]
    }
  }
  
  qqq <- sapply( 1:4, function( nch ) max( table( substr( tmp, 1, nch ) ) ) / length( tmp ) ); nch <- 0
  if ( any( qqq > 0.9 ) ) { nch <- which( qqq > 0.9 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >60% of names
  else if ( any( qqq > 0.6 ) ) { nch <- which( qqq > 0.6 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >60% of names
  else if ( any( qqq > 0.4 ) ) { nch <- which( qqq > 0.4 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >40% of names
  ##nch <- 0; while( max( table( substr( tmp, 1, nch + 1 ) ) ) / length( tmp ) > 0.4 ) nch <- nch + 1
  prefix <- NA
  if ( nch > 0 ) {
    prefix <- names( which.max( table( substr( tmp, 1, nch ) ) ) )
    if ( verbose ) message( "Assuming gene/probe names have common prefix '", prefix, "'." )
    genome.info$gene.prefix <- prefix
  } else {
    if ( verbose ) message( "Could not find a common gene/probe identifier prefix. This only matters if there's no expression matrix." )
    prefix <- genome.info$gene.prefix <- NA
  }

  tmp2 <- tmp
  if ( length( unique( nchar( tmp2 ) ) ) > 1 ) {
    nc <- max( nchar( tmp2 ) )
    for ( i in 1:length( tmp2 ) )
      tmp2[ i ] <- paste( tmp2[ i ], rep( ' ', nc - nchar( tmp2[ i ] ) ), sep='', collapse='' )
  }
  tmp2 <- do.call( rbind, strsplit( tmp2, '' ) )
  regex <- apply( tmp2, 2, function( i ) sort( unique( i ) ) )
  for ( i in 1:length( regex ) ) {
    ii <- as.integer( regex[[ i ]] )
    if ( ! any( is.na( ii ) ) ) {
      if ( length( ii ) == length( ii[ 1 ]:ii[ length( ii ) ] ) && all( ii ) == ii[1]:ii[length(i)] )
        regex[[ i ]] <- paste( '[', paste( ii[ 1 ], ii[ length( ii ) ], sep='-' ), ']', sep='' )
    }
    if ( length( regex[[ i ]][ regex[[ i ]] != ' ' ] ) > 1 ) regex[[ i ]] <- c( '[', regex[[ i ]], ']' )
    if ( any( regex[[ i ]] == '' | regex[[ i ]] == ' ' | is.na( regex[[ i ]] ) ) )
      regex[[ i ]] <- c( regex[[ i ]][ regex[[ i ]] != ' ' ], '?' )
  }
  regex <- paste( unlist( lapply( regex, paste, sep='', collapse='' ) ), sep='', collapse='' )
  if ( verbose ) message( "Assuming gene/probe names have regex '", regex, "'." )
  c( prefix, regex )
}

extend.vec <- function( v, n=n.iter ) {
  if ( length( v ) < n ) v <- c( v, rep( v[ length( v ) ], n.iter - length( v ) ) ); v }

operon.list <- function() { ## convert operon table to list of operons
  out <- list()
  ops <- as.matrix( genome.info$operons )
  for ( i in 1:nrow( ops ) ) {
    if ( ! ops[ i, 1 ] %in% names( out ) ) out[[ ops[ i, 1 ] ]] <- ops[ i, 2 ]
    else out[[ ops[ i, 1 ] ]] <- c( out[[ ops[ i, 1 ] ]], ops[ i, 2 ] )
  }
  out
}

#ifndef PACKAGE
if ( file.exists( "cmonkey-init.R" ) ) {
  ##source( "cmonkey.R", local=T )
  source( "cmonkey-init.R", local=T )
  source( "cmonkey-data-load.R", local=T ) ## Functions for loading the data
  source( "cmonkey-optim1.R", local=T ) ## Pre-v5.0.0 optimization functions
  ##source( "cmonkey-optim2.R", local=T ) ## Post-v5.0.0 optimization functions
  source( "cmonkey-motif.R", local=T ) ## Functions for motif finding/scoring
  source( "cmonkey-plotting.R", local=T ) ## Functions for all cmonkey plotting
  source( "cmonkey-postproc.R", local=T ) ## Functions for all post-processing and analysis of cmonkey clusters
  source( "cmonkey-motif-other.R", local=T ) ## Functions for motif finding with other algorithms and blasting
  source( "cmonkey-bigmem.R", local=T ) ## Functions for using on-disk list and matrix storage for big organisms
}
#endif
