source( "cmonkey-motif-other.R" )

## Create an ensemble from a single run by setting, e.g.
##   cm.func.each.iter <- function() if ( iter %% 25 == 0 ) save.image( sprintf( "zzz_hpy_%04d.RData", iter ) )
## Then run this with rdata.glob="zzz_hpy_*.RData"
## This can also be used for RData output from multiple runs.
## Then, to get bicluster #24 in RData file from iter=1900, for example:
## which(names(e$fnames.to.cluster)=="./zzz_hpyu_1900.RData"&e$fnames.to.cluster==24)
cmonkey.ensemble <- function( rdata.glob, filter=NULL ) {
  files <- list.files( patt=glob2rx( rdata.glob ), full=T )
  all.ratios <- env <- NULL
  fnames.to.cluster <- numeric()
  for ( f in files ) {
    print( f )
    load( f )
    if ( ! exists( "e" ) ) next
    print( names( e$mot.weights ) )
    if ( ! is.null( filter ) ) {
      tmp <- filter( e )
      cat( "FILTER:", tmp, "\n" )
      if ( tmp == FALSE ) next
    }
    e$clusterStack <- e$get.clusterStack( force=T )
    if ( is.null( env ) ) {
      env <- e
      all.ratios <- e$get.cluster.matrix()
      tmp <- e$get.stats() ##stats[ nrow( e$stats ), ]
      env$stats <- tmp
      tmp.f <- 1:e$k.clust; names( tmp.f ) <- rep( f, length( tmp.f ) )
      fnames.to.cluster <- tmp.f
      env$from.files <- f
      next
    }

    env$from.files <- c( env$from.files, f )
    env$k.clust <- env$k.clust + e$k.clust
    rats <- e$get.cluster.matrix()
    if ( any( ! rownames( rats ) %in% rownames( all.ratios ) ) ||
        any( ! colnames( rats ) %in% colnames( all.ratios ) ) ||
        any( ! rownames( all.ratios ) %in% rownames( rats ) ) ||
        any( ! colnames( all.ratios ) %in% colnames( rats ) ) ) {
      tmp <- matrix( NA, nrow=length( unique( c( rownames( rats ), rownames( all.ratios ) ) ) ),
                    ncol=length( unique( c( colnames( rats ), colnames( all.ratios ) ) ) ) )
      rownames( tmp ) <- unique( c( rownames( rats ), rownames( all.ratios ) ) )
      colnames( tmp ) <- unique( c( colnames( rats ), colnames( all.ratios ) ) )
      tmp[ rownames( rats ), colnames( rats ) ] <- rats
      tmp[ rownames( all.ratios ), colnames( all.ratios ) ] <- all.ratios
      all.ratios <- tmp
      rm( tmp )
    }

    if ( as.integer( gsub( '.', '', env$cmonkey.version, fixed=T ) ) < 500 ) { ## No row.membership anymore!
      new.rm <- matrix( 0, nrow=nrow( all.ratios ), ncol=ncol( env$row.membership ) + ncol( e$row.membership ) )
      rownames( new.rm ) <- rownames( all.ratios )
      new.rm[ rownames( env$row.membership ), 1:ncol( env$row.membership ) ] <- env$row.membership
      tmp <- e$row.membership
      tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$row.membership )
      new.rm[ rownames( e$row.membership ), ( ncol( env$row.membership ) + 1 ):ncol( new.rm ) ] <- tmp
      env$row.membership <- new.rm; rm( new.rm, tmp )
      ## tmp <- tmp[ rownames( tmp ) %in% rownames( env$row.membership ), ,drop=F ]
      ## tmp2 <- matrix( 0, nrow=nrow( all.ratios ), ncol=ncol( tmp ) ); rownames( tmp2 ) <- rownames( all.ratios )
      ## tmp2[ rownames( tmp ), ] <- tmp
      ## env$row.membership <- cbind( env$row.membership, tmp2 )
    }
    
    tmp.f <- 1:e$k.clust + max( c( 0, fnames.to.cluster ) ) ##);
    names( tmp.f ) <- rep( f, length( tmp.f ) )
    fnames.to.cluster <- c( fnames.to.cluster, tmp.f )

    if ( as.integer( gsub( '.', '', env$cmonkey.version, fixed=T ) ) < 500 ) {
      new.cm <- matrix( 0, nrow=ncol( all.ratios ), ncol=ncol( env$col.membership ) + ncol( e$col.membership ) )
      rownames( new.cm ) <- colnames( all.ratios )
      new.cm[ rownames( env$col.membership ), 1:ncol( env$col.membership ) ] <- env$col.membership
      tmp <- e$col.membership
      tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$col.membership )
      new.cm[ rownames( e$col.membership ), ( ncol( env$col.membership ) + 1 ):ncol( new.cm ) ] <- tmp
      env$col.membership <- new.cm; rm( new.cm, tmp )
      
      ## tmp <- e$col.membership
      ## tmp[ tmp != 0 ] <- tmp[ tmp != 0 ] + max( env$col.membership )
      ## tmp <- tmp[ rownames( tmp ) %in% rownames( env$col.membership ), ]
      ## tmp2 <- matrix( 0, nrow=ncol( all.ratios ), ncol=ncol( tmp ) )
      ## rownames( tmp2 ) <- colnames( all.ratios )
      ## tmp2[ rownames( tmp ), ] <- tmp
      ## env$col.membership <- cbind( env$col.membership, tmp2 )
    }
    
    if ( length( e$clusterStack ) < e$k.clust ) e$clusterStack[[ e$k.clust ]] <- ""
    env$clusterStack <- c( env$clusterStack, e$clusterStack[ 1:e$k.clust ] )

    tmp <- e$get.stats() ##stats[ nrow( e$stats ), ]
    if ( nrow( e$stats ) > 0 && ! all( colnames( e$stats ) %in% colnames( tmp ) ) ) {
      tmp2 <- e$stats[ 1, ] * NA
      tmp2[ , colnames( tmp ) ] <- tmp; tmp <- tmp2
    } else if ( nrow( e$stats ) <= 0 || ! all( colnames( tmp ) %in% colnames( e$stats ) ) ) {
      tmp <- tmp[ ,colnames( e$stats ) ]
    }
    env$stats <- rbind( env$stats, tmp ); rm( tmp, tmp2 )
    rm( tmp, tmp2 )

    if ( ! is.null( e$meme.scores ) && ! is.null( e$meme.scores[[ 1 ]] ) ) {
      ##e$meme.scores[[ 1 ]] <- 
      if ( length( e$meme.scores[[ 1 ]] ) < e$k.clust ) e$meme.scores[[ 1 ]][[ e$k.clust ]] <- ""
      for ( i in names( e$meme.scores ) ) {
        if ( ! i %in% names( env$meme.scores ) ) as.list( rep( "", env$k.clust ) )
        env$meme.scores[[ i ]][ tmp.f ] <- e$meme.scores[[ 1 ]][ 1:e$k.clust ]
        names( env$meme.scores[[ i ]] ) <- NULL
      }
    }
    rm( e )
  }
  
  env$ratios <- list( ratios=all.ratios )
  attr( env$ratios, "rnames" ) <- sort( unique( unlist( lapply( env$ratios, rownames ) ) ) )
  attr( env$ratios, "cnames" ) <- sort( unique( unlist( lapply( env$ratios, colnames ) ) ) )
  attr( env$ratios, "nrow" ) <- length( attr( env$ratios, "rnames" ) )
  attr( env$ratios, "ncol" ) <- length( attr( env$ratios, "cnames" ) )
  attr( env$ratios$ratios, "maxRowVar" ) <- mean( apply( env$ratios$ratios, 1, var, use="pair" ), na.rm=T ) ##* 1.2
  attr( env$ratios$ratios, "all.colVars" ) <- apply( env$ratios$ratios, 2, var, use="pair", na.rm=T )
  ##env$meme.scores <- list( `upstream meme`=env$meme.scores )
  env$fnames.to.cluster <- fnames.to.cluster
  ##rm( n.scores, net.scores, row.scores, r.scores, rr.scores, mot.scores, m.scores, row.memb, col.memb, envir=env )
  rm( net.scores, row.scores, mot.scores, envir=env )
  env
}

cmonkey.ensemble.analysis <- function( e, cluster.motifs=F, ... ) { ## e from cmonkey.ensemble()
  ##if ( force.single ) mclapply <- lapply ## even 2 cores is too much for my mac
  ##!else
  require( multicore )
  ## Get a matrix of the # of times (total) that each pair of genes is in the same cluster
  ##all.rm <- lapply( 1:nrow( e$row.membership ), function( i ) { i <- e$row.membership[ i, ]; i[ i != 0 ] } )
  row.membership <- e$row.col.membership.from.clusterStack( e$clusterStack )$r
  all.rm <- lapply( 1:nrow( row.membership ), function( i ) { i <- row.membership[ i, ]; i[ i != 0 ] } )
  row.ov <- do.call( rbind, mclapply( all.rm, function( i ) sapply( all.rm, function( j ) sum( i %in% j ) ) ) )
  ## Matrix of the total # of (unique) clusters that genes i OR j are in (union)
  all.len <- lapply( all.rm, length )
  n.clust <- do.call( rbind, mclapply( all.len, function( i ) sapply( all.len, function( j ) min( c( i, j ) ) ) ) )
  rm( all.rm, all.len )

  row.ov[ lower.tri( row.ov, diag=T ) ] <- 0
  n.clust[ lower.tri( n.clust, diag=T ) ] <- 1
  rownames( row.ov ) <- rownames( n.clust ) <- rownames( row.membership )
  row.ov[ n.clust < 10 ] <- 0; n.clust[ n.clust < 10 ] <- Inf
  r.sif <- which( row.ov / n.clust > 0.2, arr=T )
  r.sif <- data.frame( g1=rownames( r.sif ), g2=rownames( row.ov )[ r.sif[ ,2 ] ], type=rep( "bic", nrow( r.sif ) ),
                      weight=row.ov[ r.sif ] / n.clust[ r.sif ] )
  r.sif <- r.sif[ order( r.sif$weight, decreasing=T ), ]; rownames( r.sif ) <- NULL

  gs <- unique( c( as.character( r.sif$g1 ), as.character( r.sif$g2 ) ) )
  r.na <- data.frame( g=gs, short=e$get.long.names( gs, short=T )[ gs ] )
  sh <- as.character( r.na$short ); sh[ which( sh == "" ) ] <- as.character( r.na$g[ which( sh == "" ) ] )
  r.na$short <- as.factor( sh ); rm( sh, gs )

  out <- new.env() ##list( sif=r.sif, na=r.na )
  out$sif <- r.sif
  out$na <- r.na

  if ( cluster.motifs ) {
    tt.out <- e$motif.similarities.tomtom( ... )
    e$parallel.cores <- 1
    tt.out2 <- e$cluster.tomtom.results( tt.out, ... )
    ## TODO: use mcl to do the clustering...? see cmPostProc2_mot_metaclustering2.R
    out$tt.out <- tt.out
    tt.out2 <- tt.out2[ sapply( tt.out2, function( i ) length( attr( i, "mot.names" ) ) ) > 0 ]
    out$tt.out2 <- tt.out2
  }

  out$e <- e
  sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out )
  
  out
}
