###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## This func. is the meat of cmonkey - run one iteration and update the scores/parameters.
## This is called from cmonkey() - see there for its usage.
## Warning - this relies on lots of global variables, that are updated inside the function and
##   the updated variables are returned. Again, see cmonkey() for its usage.

## BIG TODO: make it so all big.memory are updated in-place
update.all.clusters <- function( env, dont.update=F, ... ) {  
  tmp <- env$row.col.membership.from.clusterStack( env$clusterStack )
  row.membership <- tmp$r; col.membership <- tmp$c
  
  tmp <- get.all.scores( ... )
  env$row.scores <- tmp$r##[,];
  env$mot.scores <- tmp$m; env$net.scores <- tmp$n; env$col.scores <- tmp$c
  env$meme.scores <- tmp$ms
  if ( ! is.null( tmp$cns ) ) env$cluster.net.scores <- tmp$cns

  tmp <- get.combined.scores( quant=F ) ## quant=T )
  r.scores <- tmp$r; c.scores <- tmp$c ##; n.scores <- tmp$n; m.scores <- tmp$m

  if ( length( tmp$scalings ) > 0 ) {
    env$row.scaling[ iter ] <- tmp$scalings[ "row" ]
    env$mot.scaling[ iter ] <- tmp$scalings[ "mot" ]
    env$net.scaling[ iter ] <- tmp$scalings[ "net" ]
  }
  rm( tmp )

  row.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "rnames" ) %in% get.rows( k ) ); if ( is.vector( row.memb ) ) row.memb <- t( row.memb )
  rownames( row.memb ) <- attr( ratios, "rnames" )    
  col.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "cnames" ) %in% get.cols( k ) ); if ( is.vector( col.memb ) ) col.memb <- t( col.memb )
  rownames( col.memb ) <- attr( ratios, "cnames" )    
  
  ## Fuzzify scores a bit for stochasticity! (fuzz should be between 0.2 and 0 (decreasing with iter)
  if ( row.scaling[ iter ] > 0 && fuzzy.index[ iter ] > 1e-5 ) {
    r.scores[,] <- r.scores[,] +
      rnorm( length( r.scores[,] ), sd=sd( r.scores[,][ row.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
    if ( ! is.null( c.scores ) ) c.scores[,] <- c.scores[,] +
      rnorm( length( c.scores[,] ), sd=sd( c.scores[,][ col.memb[,] == 1 ], na.rm=T ) * fuzzy.index[ iter ] )
  }

  tmp <- get.density.scores( ks=1:k.clust, r.scores, col.scores )
  rr.scores <- tmp$r; cc.scores <- tmp$c; rm( tmp )
    
  ## NEW - will it work? -- help shrink big clusters, grow small clusters, both in rows and cols
  size.compensation.func.rows <- function( n ) exp( -n / ( attr( ratios, "nrow" ) * n.clust.per.row / k.clust ) )
  size.compensation.func.cols <- function( n ) exp( -n / ( attr( ratios, "ncol" ) * n.clust.per.col / k.clust ) )
  for ( k in 1:k.clust ) {
    tmp <- sum( row.memb[ ,k ] )
    rr.scores[ ,k ] <- rr.scores[ ,k ] * size.compensation.func.rows( max( tmp, cluster.rows.allowed[ 1 ] ) ) 
    if ( ! is.null( cc.scores ) ) {
      tmp <- sum( col.memb[ ,k ] )
      cc.scores[ ,k ] <- cc.scores[ ,k ] * size.compensation.func.cols( max( tmp, attr( ratios, "ncol" ) / 10 ) )
    }
  }
  
  ## Fuzzify it along the lines of fuzzy c-means clustering
  ##   -- see http://en.wikipedia.org/wiki/Data_clustering#Fuzzy_c-means_clustering
  ## No -- it doesnt affect things - same ordering (and updated memberships are based on ordering)
  ##   but should use these scores to weight the centroids that are selected in the next iteration.
  ##   for this, fuzz should vary between e.g. 10 and 2
  ##rr.scores <- ( rr.scores / sum( rr.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )
  ##cc.scores <- ( cc.scores / sum( cc.scores, na.rm=T ) )^( 2 / ( fuzz - 1 ) )

  tmp <- get.updated.memberships( row.membership, col.membership, rr.scores, cc.scores )
  row.membership <- tmp$r; col.membership <- tmp$c; rm( tmp )
    

  if ( env$post.adjust == TRUE && env$iter == env$n.iter ) { ## post-adjustment of all clusters???
    ##env$pre.adjusted.row.membership <- env$row.membership
    ##env2 <-
    env$adjust.all.clusters( env, expand.only=F )
    ##env <- env2; rm( env2 )
    gc()
##     row.membership.orig <- row.membership
##     adjust.all.clusters( expand=3 )
  }
  
  if ( ! dont.update ) {
    env$clusterStack <- lapply( 1:env$k.clust, function( k )
                               list( rows=rownames( which( row.membership == k, arr=T ) ),
                                    cols=rownames( which( col.membership == k, arr=T ) ) ) )
    env$clusterStack <- env$get.clusterStack( ks=1:k.clust )     
  }
  env
}

## WARNING: This function relies on A LOT of global variables!!!
get.combined.scores <- function( quantile.normalize=F ) {
  r.scores <- row.scores[,]
  r.scores <- matrix.reference( r.scores )
## qqqifndef 
  if ( ! quantile.normalize ) {
    row.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "rnames" ) %in% get.rows( k ) )
    rownames( row.memb ) <- attr( ratios, "rnames" )
    tmp <- r.scores[,] < -20; r.scores[,][ tmp ] <- min( r.scores[,][ ! tmp ], na.rm=T )
    rsm <- r.scores[,][ row.memb ]
    tmp <- mad( rsm, na.rm=T )
    if ( tmp != 0 ) r.scores[,] <- ( r.scores[,] - median( rsm, na.rm=T ) ) / tmp
    else { tmp <- sd( rsm, na.rm=T ); if ( tmp != 0 ) r.scores[,] <- ( r.scores[,] - median( rsm, na.rm=T ) ) / tmp }
    rm( tmp, rsm )
  }
## qqqendif
  
  tmp <- r.scores[,] < -20; r.scores[,][ tmp ] <- min( r.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## -220
  r.scores[,][ is.infinite( r.scores[,] ) ] <- NA
  r.scores[,][ is.na( r.scores[,] ) ] <- max( r.scores[,], na.rm=T )
  ## }
  ##cat( "HERE: row", r.scores[1,1], "\n" )
  
##qqq ifndef 
  if ( ! quantile.normalize && ! is.null( mot.scores ) || ! is.null( net.scores ) )
    rs.quant <- quantile( r.scores[,], 0.01, na.rm=T )
##qqq endif

  if ( ! is.null( mot.scores ) ) {
    m.scores <- mot.scores[,]
    ##m.scores <- matrix.reference( m.scores ) ## not useful as m.scores (and n.scores) are not saved in the env.
  } else m.scores <- NULL
  if ( ! is.null( mot.scores ) && ! is.null( m.scores ) ) {
    tmp <- m.scores[,] < -20; m.scores[,][ tmp ] <- min( m.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## effective zero is 10^-20
##qqq ifndef 
    if ( ! quantile.normalize ) {
      m.scores[,] <- m.scores[,] - quantile( m.scores[,], 0.99, na.rm=T ) ## Make it so no mot -> zero score -> no penalization
      m.scores[,] <- m.scores[,] / abs( quantile( m.scores[,], 0.01, na.rm=T ) ) * abs( rs.quant ) ## Make it so min(mot.scores) == min(row.scores) so no mot.score is too dominant
    }
##qqq endif
  } ##!else m.scores <- NULL
  
  if ( ! is.null( net.scores ) ) {
    n.scores <- net.scores[,]
    n.scores <- matrix.reference( n.scores ) ## not useful as m.scores (and n.scores) are not saved in the env.
  } else n.scores <- NULL
  if ( ! is.null( net.scores ) && ! is.null( n.scores ) ) {
    n.scores[,] <- n.scores[,] - quantile( n.scores[,], 0.99, na.rm=T ) ## Make it so no edge -> zero score -> no penalization
##qqq ifndef 
    if ( ! quantile.normalize ) {
      qqq <- abs( quantile( n.scores[,], 0.01, na.rm=T ) )
      if ( qqq == 0 ) qqq <- sort( n.scores[,] )[ 10 ] ## For really sparse networks
      if ( qqq == 0 ) qqq <- min( n.scores[,], na.rm=T ) ## For really sparse networks
      if ( qqq != 0 ) n.scores[,] <- n.scores[,] / qqq * abs( rs.quant ) ## Make it so min(net.scores) == min(row.scores) so no net.score is too dominant
      rm( qqq )
    }
##qqq endif
  }

  if ( ! is.null( col.scores ) ) {
    ## c.scores <- NULL
    c.scores <- col.scores[,] * 0
    c.scores <- matrix.reference( c.scores )
    tmp <- c.scores[,] < -20; c.scores[,][ tmp ] <- min( c.scores[,][ ! tmp ], na.rm=T ); rm( tmp ) ## effective zero ##-220
    ##}
  } else c.scores <- NULL

  new.weights <- c( row=row.scaling[ iter ], mot=mot.scaling[ iter ], net=net.scaling[ iter ] ) ##numeric()
  
  ## Hey, instead of just making each distrib. have the same mean and variance, let's make them have the same
  ##   exact distribution!
  if ( quantile.normalize ) {
    tmp <- list( row=r.scores, mot=m.scores, net=n.scores )
    if ( sum( sapply( tmp, function( i ) ! is.null( i ) ) ) > 1 ) {
      wts <- new.weights[ ! sapply( tmp, is.null ) ]
      tmp <- quantile.normalize.scores( tmp, weights=wts )
      if ( ! is.null( r.scores ) ) r.scores[,] <- tmp$row[,]
      if ( ! is.null( m.scores ) ) m.scores[,] <- tmp$mot[,]
      if ( ! is.null( n.scores ) ) n.scores[,] <- tmp$net[,]
    }
    rm( tmp )
  }
  
  if ( new.weights[ "row" ] != 1 ) r.scores[,] <- r.scores[,] * new.weights[ "row" ] 
  if ( ! is.null( m.scores ) ) {
    tmp <- ! is.na( m.scores[,] )
    if ( is.null( r.scores ) ) r.scores <- m.scores[,] * 0
    r.scores[,][ tmp ] <- r.scores[,][ tmp ] + m.scores[,][ tmp ] * new.weights[ "mot" ]
  }
  if ( ! is.null( n.scores ) ) {
    tmp <- ! is.na( n.scores[,] )
    if ( is.null( r.scores ) ) r.scores <- n.scores[,] * 0
    r.scores[,][ tmp ] <- r.scores[,][ tmp ] + n.scores[,][ tmp ] * new.weights[ "net" ]
  }
  r.scores <- matrix.reference( r.scores )
  c.scores <- matrix.reference( c.scores )
  ##if ( ! is.null( n.scores ) ) n.scores <- matrix.reference( n.scores )
  ##if ( ! is.null( m.scores ) ) m.scores <- matrix.reference( m.scores )
  invisible( list( r=r.scores, c=c.scores, scalings=new.weights ) ) ##n=net.scores, m=mot.scores, 
}

## Quick and dirty quantile normalization of scores matrices against each other!!!
## Optionally weight the scores matrices!!!
quantile.normalize.scores <- function( scores, weights=NULL, keep.nas=F ) {
  ## scores is a list of matrices of same dimension
  if ( ! is.list( scores ) || sum( ! sapply( scores, is.null ) ) <= 1 ) return( scores )
  scores <- scores[ ! sapply( scores, is.null ) ]
  ##ns <- names( scores )
  d <- dim( scores[[ 1 ]] )
  dn <- dimnames( scores[[ 1 ]] )
  ##scores <- sapply( scores, function( s ) as.vector( s[,] ) )
  ##tmp <- apply( scores, 2, rank, ties="min", na="keep" )
  ##tmp2 <- apply( scores, 2, sort, na.last=T ) ##sapply( 1:ncol( scores ), function( i ) scores[ tmp[ ,i ], i ] )
  tmp2 <- sapply( scores, function( i ) sort( i[,], na.last=T ) ) ##sapply( 1:ncol( scores ), function( i ) scores[ tmp[ ,i ], i ] )
  if ( is.null( weights ) ) tmp2.mn <- rowMeans( tmp2, na.rm=T )
  else { ##tmp2.mn <- apply( tmp2, 1, weighted.mean, weights, na.rm=T )
    for ( i in 1:ncol( tmp2 ) ) tmp2[ ,i ] <- tmp2[ ,i ] * weights[ i ]
    tmp2.mn <- rowMeans( tmp2, na.rm=T ) / sum( weights, na.rm=T )
  }
  rm( tmp2 )
  ##out <- lapply( 1:ncol( scores ), function( i ) tmp2.mn[ tmp[ ,i ] ] )
  ##out <- lapply( out, function( i ) { z <- matrix( i, nrow=d[ 1 ], ncol=d[ 2 ] ); dimnames( z ) <- dn; z } ) 
  ##names( out ) <- ns
  ##out
  ##tmp <- sapply( scores, function( i ) rank( i[,], ties="min", na="keep" ) )
  out <- list()
  for ( n in names( scores ) ) {
    tmp <- rank( scores[[ n ]][,], ties="min", na="keep" )
    ##scores[[ n ]][,] <-
    z <- matrix( tmp2.mn[ tmp ], nrow=d[ 1 ], ncol=d[ 2 ] ) ##tmp[ ,n ] ]
    dimnames( z ) <- dn
    if ( keep.nas ) z[ is.na( scores[[ n ]][,] ) ] <- NA
    out[[ n ]] <- z
  }
  out ##scores
}

get.all.scores <- function( ks=1:k.clust, force.row=F, force.col=F, force.motif=F, force.net=F,
                           quantile.normalize=F ) {
  mc <- get.parallel( length( ks ) )

  if ( is.null( row.scores ) ) {
    row.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) )
    rownames( row.scores ) <- attr( ratios, "rnames" )
    row.scores <- matrix.reference( row.scores )
  }
  
  ## Compute row.scores (microarray data)
  if ( force.row || ( row.scaling[ iter ] > 0 && ! is.na( row.iters[ 1 ] ) && iter %in% row.iters ) ) {
    ##if ( ! exists( "row.scores" ) || is.null( row.scores ) || nrow( row.scores ) != attr( ratios, "nrow" ) ||
    ##    ncol( row.scores ) != max( ks ) ) {
    if ( is.null( row.scores ) ) {
      row.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) )
      rownames( row.scores ) <- attr( ratios, "rnames" )
      row.scores <- matrix.reference( row.scores )
    } else row.scores[ ,ks ] <- 0
    ## TODO: Try to quantile normalize the ratios scores too (as is done w/ mot.scores, net.scores below) but
    ##    this is harder because each different "tmp.row" matrix may have different nrows.
    ##rs.func <- function() { ## for profiling
    row.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "rnames" ) %in% get.rows( k ) )
    rownames( row.memb ) <- attr( ratios, "rnames" )
    for ( i in names( ratios ) ) { 
      if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
      tmp.row <- do.call( cbind, mc$apply( ks, get.row.scores, ratios=ratios[[ i ]]
                                          ) )
      tmp <- is.infinite( tmp.row ) | is.na( tmp.row )
      if ( any( tmp ) ) tmp.row[ tmp ] <-
        quantile( tmp.row[ row.memb[ rownames( tmp.row ), ] & ! tmp ], 0.95 ) ##0 ## No measurement is NOT an NA - it is an observation
      tmp <- rownames( row.scores )[ rownames( row.scores ) %in% rownames( tmp.row ) ]
      row.scores[ tmp, ks ] <- row.scores[ tmp, ks ] + tmp.row[ tmp, ] * row.weights[ i ]
      rm( tmp.row, tmp )
    }
    ##row.scores }; row.scores <- rs.func()    
  }
  
  ## Compute col.scores (microarray data)
  ##if ( ! exists( "col.scores" ) || is.null( col.scores ) || nrow( col.scores ) != attr( ratios, "ncol" ) ||
  ##    ncol( col.scores ) != max( ks ) ) {
  if ( n.clust.per.col < k.clust && ## If using all columns, can skip this step...  
      ( force.col || ( row.scaling[ iter ] > 0 && ! is.na( col.iters[ 1 ] ) && iter %in% col.iters ) ) ) {
    if ( is.null( col.scores ) ) {
      col.scores <- matrix( 0, nrow=attr( ratios, "ncol" ), ncol=max( ks ) )
      rownames( col.scores ) <- attr( ratios, "cnames" )
      col.scores <- matrix.reference( col.scores )
    } else col.scores[ ,ks ] <- 0
    ##col.scores <- matrix( 0, nrow=attr( ratios, "ncol" ), ncol=max( ks ) )
    ## cs.func <- function() { ## for profiling
    col.memb <- sapply( 1:k.clust, function( k ) attr( ratios, "cnames" ) %in% get.cols( k ) )
    rownames( col.memb ) <- attr( ratios, "cnames" )
    for ( i in names( row.weights ) ) { 
      if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
      tmp.col <- do.call( cbind, mc$apply( ks, get.col.scores, ratios=ratios[[ i ]]
                                          ) )
      tmp <- is.infinite( tmp.col ) | is.na( tmp.col )
      if ( any( tmp ) ) tmp.col[ tmp ] <-
        quantile( tmp.col[ col.memb[,][ rownames( tmp.col ), ] == 1 & ! tmp ], 0.95 ) ##0 ## No measurement is NOT an NA - it is an observation
      tmp <- rownames( col.scores )[ rownames( col.scores ) %in% rownames( tmp.col ) ]
      col.scores[ tmp, ks ] <- col.scores[ tmp, ks ] + tmp.col[ tmp, ] * row.weights[ i ]
      rm( tmp.col, tmp )
    }
    ## col.scores }; col.scores <- cs.func()
  }

  ## Run meme on each cluster (every meme.iters iterations)
  for ( i in names( mot.weights ) ) {
    if ( force.motif == "run.meme" || ( mot.scaling[ iter ] > 0 && ! is.na( meme.iters[[ i ]][ 1 ] ) &&
           iter %in% meme.iters[[ i ]] && exists( "genome.info" ) && ! no.genome.info ) ) {
      if ( mot.weights[ i ] == 0 || is.na( mot.weights[ i ] ) ) next
      tmp <- motif.all.clusters( ks, seq.type=i, verbose=T ) ##strsplit( i, " " )[[ 1 ]][ 1 ],
                                               ##algo=strsplit( i, " " )[[ 1 ]][ 2 ] )
      meme.scores[[ i ]] <- tmp
    }
  }
  
  ## Compute mot.scores from meme output
  if ( force.motif == TRUE || force.motif == "run.meme" || ( mot.scaling[ iter ] > 0 && ! is.na( mot.iters[ 1 ] ) &&
                         ##iter %in% c( meme.iters, mot.iters ) && exists( "genome.info" ) && ! no.genome.info ) ) {
                                iter %in% mot.iters && exists( "genome.info" ) && ! no.genome.info ) ) {
    if ( is.null( mot.scores ) ) {
      mot.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) ) 
      rownames( mot.scores ) <- attr( ratios, "rnames" )
      mot.scores <- matrix.reference( mot.scores )
    } else mot.scores[ ,ks ] <- 0
    tmp.mots <- list()
    for ( i in names( mot.weights ) ) {
      if ( mot.weights[ i ] == 0 || is.na( mot.weights[ i ] ) ) next
      tmp.mot <- do.call( cbind, mc$apply( ks, get.motif.scores, seq.type=i ) )
      tmp.mot[ is.infinite( tmp.mot ) | is.na( tmp.mot ) ] <- 0 ## No measurement is NOT an NA - it is an observation
      if ( quantile.normalize && sum( mot.weights > 0 & ! is.na( mot.weights ) ) > 1 ) tmp.mots[[ i ]] <- tmp.mot
      else mot.scores[ ,ks ] <- mot.scores[ ,ks ] + tmp.mot[,] * mot.weights[ i ]
      rm( tmp.mot )
    }
    if ( quantile.normalize && length( tmp.mots ) > 1 ) {
      ##combine <- "sum" ## TODO: allow one to take eg. the "min" of each element in the matrices (instead of the sum)
      tmp.mots <- quantile.normalize.scores( tmp.mots, weights=mot.weights[ mot.weights != 0 ] )
      for ( i in names( tmp.mots ) ) mot.scores[ ,ks ] <- mot.scores[ ,ks ] + tmp.mots[[ i ]][,] * mot.weights[ i ]
      rm( tmp.mots )
    }
  }

  cluster.ns <- NULL
  ## Compute net.scores from STRING ... add weighted scores for other networks (if exist)
  if ( force.net || ( net.scaling[ iter ] > 0 && ! is.na( net.iters[ 1 ] ) && exists( "genome.info" ) &&
                     iter %in% net.iters ) ) {
    ##net.scores <- NULL
    if ( is.null( net.scores ) ) {
      net.scores <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=max( ks ) ) 
      rownames( net.scores ) <- attr( ratios, "rnames" )
      net.scores <- matrix.reference( net.scores )
    } else net.scores[ ,ks ] <- 0
    tmp.nets <- list()
    ## ns.func <- function() { ## for profiling
    for ( i in names( networks ) ) { 
      if ( net.weights[ i ] == 0 || is.na( net.weights[ i ] ) ) next
      if ( nrow( subset( networks[[ i ]], protein1 %in% attr( ratios, "rnames" ) & protein2 %in%
                        attr( ratios, "rnames" ) ) ) <= 0 ) next
      tmp.net <- do.call( cbind, mc$apply( ks, get.network.scores, net=networks[[ i ]] ) )
      if ( all( is.na( tmp.net ) ) || all( is.character( tmp.net ) ) ) next
      ##tmp.net[,] <- tmp.net[,] - max( tmp.net[ ! is.infinite( tmp.net ) ], na.rm=T ) -
      ##  abs( diff( range( tmp.net[ ! is.infinite( tmp.net ) ], na.rm=T ) ) ) / 10
      tmp.net[ is.infinite( tmp.net ) | is.na( tmp.net ) ] <- 0 ## No edge is NOT an NA - it is an observation
      if ( quantile.normalize && sum( net.weights > 0 & ! is.na( net.weights ) ) > 1 ) tmp.nets[[ i ]] <- tmp.net
      else net.scores[ ,ks ] <- net.scores[ ,ks ] + tmp.net[,] * net.weights[ i ]
      cluster.ns <- cbind( cluster.ns, do.call( c, mc$apply( ks, function( k ) mean( tmp.net[ get.rows( k ), k ],
                                                                                    na.rm=T, trim=0.05 ) ) ) )
      colnames( cluster.ns )[ ncol( cluster.ns ) ] <- i
      rm( tmp.net )
    }    
    if ( quantile.normalize && length( tmp.nets ) > 1 ) {
      tmp.nets <- quantile.normalize.scores( tmp.nets, weights=net.weights[ net.weights != 0 ] )
      ##combine <- "sum" ## TODO: allow one to take eg. the "min" of each element in the matrices (instead of the sum)
      for ( i in names( tmp.nets ) ) net.scores[ ,ks ] <- net.scores[ ,ks ] + tmp.nets[[ i ]][,] * net.weights[ i ]
      rm( tmp.nets )
    }

    cluster.ns <- cbind( cluster.ns, do.call( c, mc$apply( ks, function( k ) mean( net.scores[ get.rows( k ), k ],
                                                                                  na.rm=T, trim=0.05 ) ) ) )
    colnames( cluster.ns )[ ncol( cluster.ns ) ] <- "net.scores"
    ##list( net.scores, cluster.ns ) }; tmp <- ns.func(); net.scores <- tmp[[ 1 ]]; cluster.ns <- tmp[[ 2 ]]
  }
  list( r=row.scores, m=mot.scores, ms=meme.scores, n=net.scores, c=col.scores, cns=cluster.ns )  ##r=row.scores, 
}

attr( get.all.scores, "version" ) <- 1

## Hack to get more genes added to really small (<10 rows) clusters
get.density.scores <- function( ks=1:k.clust, r.scores, c.scores, plot="none",
                               bw.scale=function( nr ) exp( -nr / 10 ) * 10 ) { ##r.scores, col.scores, 
  rr <- attr( ratios, "rnames" )
  rs <- r.scores ##[,]
  ## TODO: adjust bw to be function of iter (is this necessary with fuzzy.index turned on???)
  bw.r <- max( diff( range( rs[,], na.rm=T ) ) / 100, 0.001 ) ##"nrd0"
  get.rr.scores <- function( k ) {
    rows <- get.rows( k ) 
    cols <- get.cols( k ) 
    rsk <- rs[ ,k, drop=T ]
    if ( length( rows ) > 0 && length( cols ) > 0 ) {
      bw <- bw.r * bw.scale( length( rows ) )
      d <- density( rsk[ rows ], na.rm=T, bw=bw, adjust=2,
                   from=min( rsk, na.rm=T ) - 1, to=max( rsk, na.rm=T ) + 1, n=256 )
      p <- approx( d$x, rev( cumsum( rev( d$y ) ) ), rsk )$y
      if ( "rows" %in% plot ) { h=hist(rsk,breaks=50,main=NULL,xlab="Combined scores");tmp.scale <- round( attr( ratios, "nrow" ) / length( rows )/4 );hist(rep(rsk[rows],tmp.scale),breaks=h$breaks,col="red",border="red",add=T);hist(rsk,breaks=h$breaks,add=T);lines(d$x,d$y/max(d$y,na.rm=T)*attr(ratios,"nrow")/50,col="blue");lines(sort(rsk),p[order(rsk)]/max(p,na.rm=T)*attr(ratios,"nrow")/50,col="green") }
    } else {
      p <- rep( 1, attr( ratios, "nrow" ) )
    }
    return( p / sum( p, na.rm=T ) ) ##* length( rows )
  }

  rr.scores <- NULL
  mc <- get.parallel( length( ks ) )

  if ( ! is.null( row.scores ) ) {
    rr.scores <- row.scores[,] * 0
    rr.scores <- matrix.reference( rr.scores )
    rr.scores[,] <- do.call( cbind, mc$apply( ks, get.rr.scores ) )
    rr.scores[,][ is.infinite( rr.scores[,] ) ] <- NA
    ##rownames( rr.scores ) <- attr( ratios, "rnames" )
  }
  
  cc.scores <- NULL ## this is required for n.clust.per.col == k.clust !!! dont delete it!
  if ( ! is.null( col.scores ) ) {
    cs <- c.scores ##[,]    
    bw.c <- max( diff( range( cs[,], na.rm=T ) ) / 100, 0.001 )
    get.cc.scores <- function( k ) {
      cols <- get.cols( k ) 
      rows <- get.rows( k ) 
      csk <- cs[ ,k, drop=T ]
      if ( length( cols ) > 0 && length( rows ) > 0 &&
          ! all( is.na( csk[ cols ] ) ) && ! all( is.infinite( csk[ cols ] ) ) &
          ! all( csk[ cols ][ ! is.na( csk[ cols ] ) ] == csk[ cols[ ! is.na( csk[ cols ] ) ][ 1 ] ] ) ) {
        d <- density( csk[ cols ], na.rm=T, from=min( csk, na.rm=T ) - 1, to=max( csk, na.rm=T ) + 1,
                     bw=bw.c, adjust=2, n=256 )
        p <- approx( d$x, rev( cumsum( rev( d$y ) ) ), csk )$y
        if ( "cols" %in% plot ) { h=hist(csk,breaks=50,main=NULL,xlab="Combined scores");tmp.scale <- round( attr( ratios, "ncol" ) / length( cols )/4 ) + 1;hist(rep(csk[cols],tmp.scale),breaks=h$breaks,col="red",border="red",add=T);hist(csk,breaks=h$breaks,add=T);lines(d$x,d$y/max(d$y,na.rm=T)*attr(ratios,"ncol")/50,col="blue");lines(sort(csk),p[order(csk)]/max(p,na.rm=T)*attr(ratios,"ncol")/50,col="green") }
      } else {
        p <- rep( 1, attr( ratios, "ncol" ) )
      }
      return( p / sum( p, na.rm=T ) ) ##* length( cols )
    }

    cc.scores <- col.scores[,] * 0
    cc.scores <- matrix.reference( cc.scores )

    if ( ! is.null( c.scores ) ) { ##&& ! is.na( c.scores ) ) {
      cc.scores[,] <- do.call( cbind, mc$apply( ks, get.cc.scores ) )
      cc.scores[,][ is.infinite( cc.scores ) ] <- NA
      ##rownames( cc.scores ) <- attr( ratios, "cnames" )
    }
  }
  
  invisible( list( r=rr.scores, c=cc.scores ) )
}  

motif.all.clusters <- function( ks=1:k.clust, seq.type=names( mot.weights )[ 1 ], verbose=T, debug=F, ... ) { 
  out.ms <- meme.scores[[ seq.type ]]
  mc <- get.parallel( length( ks ), verbose=T, para.cores=get( "parallel.cores.motif" ) )
  ## Make sure random seed is different for each process (for tempfile generation!)
  if ( any( grepl( "foreach", deparse( mc$apply ) ) ) && getDoParName() == "doMC" )
    mc$apply <- function( list, FUN, ... ) 
      foreach( l=list, .options.multicore=list( preschedule=F, set.seed=T ) ) %dopar% { FUN( l, ... ) }
  if ( ! debug ) {
    out.ms <- mc$apply( ks, FUN=function( k ) try( motif.one.cluster( k, seq.type=seq.type, verbose=F, ... ) ) )
  } else {
    message( "DEBUG MODE: NOT PARALLELIZING!\n" )
    out.ms <- lapply( ks, FUN=function( k ) motif.one.cluster( k, seq.type=seq.type, verbose=T, ... ) )
  }
  out.ms[[ k.clust + 1 ]] <- ""

  for ( k in ks ) {
    if ( length( out.ms ) < k || is.null( out.ms[[ k ]] ) || class( out.ms[[ k ]] ) == "try-error" ||
        out.ms[[ k ]]$k != k || ( ! is.null( out.ms[[ k ]]$iter ) && out.ms[[ k ]]$iter != iter ) ) {
      out <- try( motif.one.cluster( k, seq.type=seq.type, verbose=T, ... ) ) 
    } else {
      out <- out.ms[[ k ]]
    }
    if ( class( out ) == "try-error" ) ## try one last time
      out <- try( motif.one.cluster( k, seq.type=seq.type, verbose=T, ... ) ) 
    if ( class( out ) == "try-error" || is.null( out ) || out$k != k ) { ## fail (skip)
      message( "ERROR on cluster ", k )
      out <- list() 
    } else if ( verbose ) {
      cat( iter, k, length( get.rows( k ) ), seq.type, '\t' ) 
    }      
    if ( verbose ) {
      if ( is.null( out ) || is.null( out$meme.out ) ) cat( 'Inf \n' )
      else {
        ind <- 1
        if ( ! is.null( out$pv.ev ) ) {
          if ( "p.value" %in% colnames( out$pv.ev[[ 1 ]] ) )
            mn <- mean( log10( out$pv.ev[[ ind ]][ rownames( out$pv.ev[[ ind ]] ) %in% get.rows( k ), "p.value" ] ), na.rm=T )
          else mn <- mean( log10( out$pv.ev[[ ind ]]$pvals ), na.rm=T )
        } else {
          mn <- 'Inf'
        }
        cat( k, if ( attr( out$meme.out, "is.pal" ) ) "pal" else "non", 
            sapply( out$meme.out[ 1:min( 3, length( out$meme.out ) ) ], "[[", "e.value" ), mn,
            ##if ( ! is.null( out$pv.ev ) )
            ##mean( log10( out$pv.ev[[ 1 ]][ rownames( out$pv.ev[[ 1 ]] ) %in% get.rows( k ), "p.value" ] ), na.rm=T )
            ##!else 'Inf',
            '\t', pssm.to.string( out$meme.out[[ 1 ]]$pssm ), "\n" )
      }
    }
    out$iter <- iter
    out$k <- k
    out.ms[[ k ]] <- out
  }

  ## mot.rows <- character()
  ## for ( k in 1:k.clust ) {
  ##   if ( is.null( out.ms[[ k ]]$pv.ev ) ) next
  ##   mot.rows <- unique( c( mot.rows, rownames( out.ms[[ k ]]$pv.ev[[ 1 ]] ) ) )
  ## }
  ## mot.rows <- sort( mot.rows )
  ## out.pv <- out.ev <- NULL
  ## for ( k in 1:k.clust ) {
  ##   m <- out.ms[[ k ]]
  ##   if ( is.null( m ) || is.null( m$pv.ev ) ) {
  ##     out.pv <- cbind( out.pv, rep( NA, length( mot.rows ) ) )
  ##     out.ev <- cbind( out.ev, rep( NA, length( mot.rows ) ) )
  ##   } else {
  ##     m.scores <- numeric( length=length( mot.rows ) )
  ##     tmp <- m$pv.ev[[ 1 ]][ ,"p.value" ]; names( tmp ) <- rownames( m$pv.ev[[ 1 ]] )
  ##     m.scores <- tmp[ mot.rows ]
  ##     out.pv <- cbind( out.pv, m.scores ); colnames( out.pv ) <- NULL
  ##     m.scores <- numeric( length=length( mot.rows ) )
  ##     tmp <- m$pv.ev[[ 1 ]][ ,"e.value" ]; names( tmp ) <- rownames( m$pv.ev[[ 1 ]] )
  ##     m.scores <- tmp[ mot.rows ]
  ##     out.ev <- cbind( out.ev, m.scores ); colnames( out.ev ) <- NULL
  ##     out.ms[[ k ]]$pv.ev[[ 1 ]] <- NULL ## Consolidate all pv.ev[[1]] data frames into 2 matrices
  ##   }
  ## }
  ## rownames( out.pv ) <- mot.rows
  ## if ( ! is.null( out.pv ) ) rownames( out.pv ) <- mot.rows
  
  ##out.ms[[ k.clust + 1 ]] <- ""
  ##if ( require( ff ) && object.size( out.pv ) / 1048600 > big.memory ) { ## Matrices > 50MB get filebacked
  ##  dir.create( cmonkey.filename, recursive=T, show=F )
  ##  out.ms$all.pv <- as.ff( out.pv, filename=paste( cmonkey.filename, "/all.pv.", seq.type, sep="" ), overwrite=T )
  ##} else
  ## Consolidate all pv.ev[[1]] data frames into 2 matrices
  out.ms$all.pv <- make.pv.ev.matrix( out.ms ) ##out.pv
  
  if ( FALSE ) { ##big.memory == TRUE ) {
    for ( k in 1:k.clust ) { ## This consolidates space but breaks things if meme was not run again (i.e.
      m <- out.ms[[ k ]]     ## cluster didn't change), so only use it in case of 'big.memory==TRUE'
      if ( ! is.null( m ) && ! is.null( m$pv.ev ) ) out.ms[[ k ]]$pv.ev[[ 1 ]] <- NULL
    }
  }
  
  ##if ( require( ff ) && object.size( out.ev ) / 1048600 > big.memory ) { ## Matrices > 50MB get filebacked
  ##  dir.create( cmonkey.filename, recursive=T, show=F )
  ##  out.ms$all.ev <- as.ff( out.ev, filename=paste( cmonkey.filename, "/all.ev.", seq.type, sep="" ), overwrite=T )
  ##} else out.ms$all.ev <- out.ev
  
  attr( out.ms, "seq.type" ) <- seq.type
  invisible( out.ms )
}  

make.pv.ev.matrix <- function( out.ms, make.ev=F ) {
  mot.rows <- character()
  for ( k in 1:k.clust ) {
    if ( is.null( out.ms[[ k ]]$pv.ev ) ) next
    mot.rows <- unique( c( mot.rows, rownames( out.ms[[ k ]]$pv.ev[[ 1 ]] ) ) )
  }
  mot.rows <- sort( mot.rows )
  out.pv <- out.ev <- NULL
  for ( k in 1:k.clust ) {
    m <- out.ms[[ k ]]
    if ( is.null( m ) || is.null( m$pv.ev ) ) {
      out.pv <- cbind( out.pv, rep( NA, length( mot.rows ) ) )
      if ( make.ev ) out.ev <- cbind( out.ev, rep( NA, length( mot.rows ) ) )
    } else {
      m.scores <- numeric( length=length( mot.rows ) )
      tmp <- m$pv.ev[[ 1 ]][ ,"p.value" ]; names( tmp ) <- rownames( m$pv.ev[[ 1 ]] )
      m.scores <- tmp[ mot.rows ]
      out.pv <- cbind( out.pv, m.scores ); colnames( out.pv ) <- NULL
      if ( make.ev ) {
        m.scores <- numeric( length=length( mot.rows ) )
        tmp <- m$pv.ev[[ 1 ]][ ,"e.value" ]; names( tmp ) <- rownames( m$pv.ev[[ 1 ]] )
        m.scores <- tmp[ mot.rows ]
        out.ev <- cbind( out.ev, m.scores ); colnames( out.ev ) <- NULL
      }
      out.ms[[ k ]]$pv.ev[[ 1 ]] <- NULL ## Consolidate all pv.ev[[1]] data frames into 2 matrices
    }
  }
  rownames( out.pv ) <- mot.rows
  if ( ! is.null( out.pv ) ) rownames( out.pv ) <- mot.rows
  out.pv
}

## This reassigning of all rows over all clusters is too much ... problem is crappy clusters end up losing
##     ALL rows, and okay clusters get TOO MANY. IDEA: allow only 1 change per gene per iter
## TODO: When updating, don't let number of genes in a cluster get > max.cluster.rows or < min.cluster.rows
get.updated.memberships <- function( row.membership, col.membership, rr.scores, cc.scores ) {
  ##if ( ! is.null( rr.scores ) ) {
  ##n.rows <- tabulate( row.membership )
  rm <- t( apply( rr.scores, 1, order, decreasing=T )[ 1:n.clust.per.row, ,drop=F ] ) ##[ iter ]
  rm <- t( apply( rm, 1, sort ) ); if ( n.clust.per.row == 1 ) rm <- t( rm ) ##[ iter ]
  ##cra <- cluster.rows.allowed
  if ( ncol( rm ) < ncol( row.membership ) ) rm <- cbind( rm, matrix( 0, nrow=nrow( rm ),
                                                                     ncol=ncol( row.membership ) - ncol( rm ) ) )

  for ( i in 1:nrow( rm ) ) {
    if ( all( rm[ i, ] %in% row.membership[ i, ] ) ) next
    mc <- max.changes[ "rows" ] ## If 0<max.changes<1, it's a "prob of seeing a change" and choose whether to change randomly
    if ( mc < 1 && mc > 0 && runif( 1 ) > mc ) next ##!else mc <- 1
    for ( ii in 1:mc ) { 
    if ( sum( ! rm[ i, ] %in% row.membership[ i, ] ) >= mc ) {
      if ( any( row.membership[ i, ] == 0 ) ) {
        col.change <- which( row.membership[ i, ] == 0 )[ 1 ] ## If any are zero (no cluster assigned) use that index
      } else {
        ttmp <- tabulate( row.membership[ i, ] ) ## Any clusters have this gene assigned more than once?
        if ( any( ttmp > 1 ) ) { ## If so, use that index
          col.change <- which( row.membership[ i, ] %in% which( ttmp > 1 ) )[ 1 ]
        } else { ## Only make change if it improves the rr.score value 
          delta <- rr.scores[ i, rm[ i, ] ] - rr.scores[ i, row.membership[ i, ] ]
          if ( any( row.membership[ i, ] %in% rm[ i, ] ) ) delta[ row.membership[ i, ] %in% rm[ i, ] ] <- 0
          if ( all( is.na( delta ) | delta <= 0 ) ) next
          col.change <- which.max( delta ) 
        }
      }
      if ( exists( "maintain.seed" ) && ! is.null( maintain.seed ) && ! is.null( maintain.seed$rows ) && ## Don't remove row i from cluster
          ! is.null( maintain.seed$rows[[ as.character( row.membership[ i, col.change ] ) ]] ) &&
          rownames( row.membership )[ i ] %in%
          maintain.seed$rows[[ as.character( row.membership[ i, col.change ] ) ]] ) next 
      if ( ! rm[ i, col.change ] %in% row.membership[ i, ] ) row.membership[ i, col.change ] <- rm[ i, col.change ]
    }
  }
  }
  ##} 
  
  if ( ! is.null( cc.scores ) ) {
  ##n.cols <- tabulate( col.membership )
  cm <- t( apply( cc.scores, 1, order, decreasing=T )[ 1:n.clust.per.col, ,drop=F ] ) ##[ iter ]
  if ( ncol( cm ) < ncol( col.membership ) ) cm <- cbind( cm, matrix( 0, nrow=nrow( cm ),
                                                                     ncol=ncol( col.membership ) - ncol( cm ) ) )

  for ( i in 1:nrow( cm ) ) {
    mc <- max.changes[ "cols" ] ## If 0<max.changes<1, it's a "prob of seeing a change" and choose whether to change randomly
    if ( mc < 1 && mc > 0 && runif( 1 ) > mc ) next ##!else mc <- 1
    for ( ii in 1:mc ) {
    if ( sum( ! cm[ i, ] %in% col.membership[ i, ] ) >= mc ) {
      if ( any( col.membership[ i, ] == 0 ) ) {
        col.change <- which( col.membership[ i, ] == 0 )[ 1 ]
      } else {
        ttmp <- tabulate( col.membership[ i, ] )
        if ( any( ttmp > 1 ) ) {
          col.change <- which( col.membership[ i, ] %in% which( ttmp > 1 ) )[ 1 ]
        } else { ## Only make change if it improves the cc.score value 
          delta <- cc.scores[ i, cm[ i, ] ] - cc.scores[ i, col.membership[ i, ] ]
          if ( all( is.na( delta ) | delta <= 0 ) ) next
          col.change <- which.max( delta ) 
        }
      }
      if ( exists( "maintain.seed" ) && ! is.null( maintain.seed ) && ! is.null( maintain.seed$cols ) && ## Don't remove col i from cluster
          ! is.null( maintain.seed$cols[[ as.character( col.membership[ i, col.change ] ) ]] ) &&
          colnames( col.membership )[ i ] %in%
          maintain.seed$cols[[ as.character( col.membership[ i, col.change ] ) ]] ) next 
      col.membership[ i, col.change ] <- cm[ i, col.change ]
    }
  } 
  }
  }
  
  ##cat("HERE:",sum(sapply(1:nrow(rm),function(i)sum(!rm[i,]%in%row.membership[i,]))),
  ##    sum(sapply(1:nrow(cm),function(i)sum(!cm[i,]%in%col.membership[i,]))),"\n")
  invisible( list( r=row.membership, c=col.membership ) )
}


## Remove clusters (set their indexes in row/col.membership to zero) if they have too many rows
## remove.clusters.toobig <- function( toobig=cluster.rows.allowed[ 2 ] ) {
##   if ( any( tabulate( row.membership, k.clust ) >= toobig ) ) { 
##     cat( "These", sum( tabulate( row.membership, k.clust ) >= toobig, na.rm=T ),
##       "clusters have TOO MANY members: ", which( tabulate( row.membership, k.clust ) >= toobig ), "\n" )
##     has.too.many <- which( tabulate( row.membership, k.clust ) >= toobig )
##     row.membership[ row.membership %in% has.too.many ] <- 0
##     ##col.membership[ col.membership %in% has.too.many ] <- 0
##     for ( k in has.too.many ) meme.scores[[ k ]] <- list( iter=iter )
##     ##rows.changed[ has.too.many ] <- rows.changed.motif[ has.too.many ] <- rows.changed.net[ has.too.many ] <- TRUE
##   }
##   invisible( list( r=row.membership, ms=meme.scores ) )
## }  

row.col.membership.from.clusterStack <- function( cs ) {
  row.memb <- col.memb <- NULL ##col.membership * 0
  for ( k in 1:length( cs ) ) {
    row.memb <- cbind( row.memb, rep( 0, attr( ratios, "nrow" ) ) )
    if ( ncol( row.memb ) == 1 ) rownames( row.memb ) <- attr( ratios, "rnames" )
    rows <- cs[[ k ]]$rows; rows <- rows[ ! is.na( rows ) ]
    row.memb[ rows, k ] <- k
    col.memb <- cbind( col.memb, rep( 0, attr( ratios, "ncol" ) ) )
    if ( ncol( col.memb ) == 1 ) rownames( col.memb ) <- attr( ratios, "cnames" )
    cols <- cs[[ k ]]$cols; cols <- cols[ ! is.na( cols ) ]
    col.memb[ cols, k ] <- k
  }
  row.memb <- t( apply( row.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
  row.memb <- row.memb[ ,apply( row.memb, 2, sum ) != 0, drop=F ]
  colnames( row.memb ) <- NULL
  col.memb <- t( apply( col.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
  col.memb <- col.memb[ ,apply( col.memb, 2, sum ) != 0, drop=F ]
  colnames( col.memb ) <- NULL

  ## note: need to do this afterwards, too:
  # e$row.memb <- t( apply( e$row.membership, 1, function( i ) 1:e$k.clust %in% i ) )
  # e$col.memb <- t( apply( e$col.membership, 1, function( i ) 1:e$k.clust %in% i ) )
  # e$clusterStack <- e$get.clusterStack( force=T )
  ##

  if ( ncol( row.memb ) < n.clust.per.row ) row.memb <- cbind( row.memb, rep( 0, nrow( row.memb ) ) )
  if ( ncol( col.memb ) < n.clust.per.col ) col.memb <- cbind( col.memb, rep( 0, nrow( col.memb ) ) )
  list( r=row.memb, c=col.memb )
}

## TODO: add another function to randomly seed clusters with no rows or no cols
re.seed.empty.clusters <- function( row.membership, col.membership,
                                   toosmall.r=cluster.rows.allowed[ 1 ], toosmall.c=0,
                                   toobig.r=cluster.rows.allowed[ 2 ], 
                                   n.r=cluster.rows.allowed[ 1 ] * 2, n.c=5 ) {
  ## TODO: for zero-row clusters, take a random gene(s) and assign it to this cluster.
  rm <- row.membership
  rats <- get.cluster.matrix()
  if ( any( tabulate( unlist( apply( rm, 1, unique ) ), k.clust ) <= toosmall.r ) ) {
    which.zero <- which( tabulate( unlist( apply( rm, 1, unique ) ), k.clust ) <= toosmall.r )
    cat( "These", length( which.zero ), "clusters have TOO FEW rows: ", which.zero, "\n" )
    which.toobig <- which( tabulate( unlist( apply( rm, 1, unique ) ), k.clust ) >= toobig.r )
    cat( "These", length( which.toobig ), "clusters have TOO MANY rows: ", which.toobig, "\n" )
    which.zero <- c( which.zero, which.toobig )
    for ( k in which.zero ) {
      all.zero <- names( which( apply( rm, 1, function( i ) all( i <= toosmall.r ) ) ) )
      if ( length( all.zero ) < n.r ) {
        all.zero <- unique( c( all.zero, rownames( which( rm == 0, arr=T ) ) ) )
        all.zero <- unique( c( all.zero, names( which( apply( rm, 1, function( i ) all( i == i[ 1 ] ) ) ) ) ) )
      }
      if ( length( all.zero ) <= 1 ) break
      gs <- sample( all.zero, 1 ) ##min( length( all.zero ), n.r ) )
      cors <- apply( rats[ all.zero, ], 1, cor, rats[ gs, ], use="pairwise" )
      gs <- names( cors[ order( cors, decreasing=T )[ 1:n.r ] ] ); gs <- gs[ ! is.na( gs ) ]
      ##cat(all.zero,"\t",gs,"\n")
      ##cat(k,gs,"\n")
      for ( g in gs ) {
        if ( any( rm[ g, ] == 0 ) ) rm[ g, which( rm[ g, ] == 0 )[ 1 ] ] <- k
        else rm[ g, 1 ] <- k
      }
    }
    for ( tt in names( mot.weights ) ) for ( k in which.zero ) meme.scores[[ tt ]][[ k ]] <- list( iter=iter )
    ##rows.changed[ which.zero ] <- rows.changed.motif[ which.zero ] <- rows.changed.net[ which.zero ] <- TRUE
  }
  
  ## TODO: for zero-col clusters, take a random col(s) and assign it to this cluster.
  cm <- col.membership
  if ( any( tabulate( cm, k.clust ) <= toosmall.c ) ) {
    which.zero <- which( tabulate( cm, k.clust ) <= toosmall.c )
    cat( "These", length( which.zero ), "clusters have TOO FEW columns: ", which.zero, "\n" )
    for ( k in which.zero ) {
      all.zero <- names( which( apply( cm, 1, function( i ) all( i <= toosmall.c ) ) ) )
      if ( length( all.zero ) <= n.c )
        all.zero <- unique( c( all.zero, rownames( which( cm == 0, arr=T ) ) ) )
      if ( length( all.zero ) <= 1 ) break
      cs <- unique( sample( all.zero, min( length( all.zero ), n.c ) ) ); cs <- cs[ ! is.na( cs ) ]
      ##cat( "\tSetting", k, "<-", cs, "\n" )
      for ( cc in cs ) cm[ cc, which( cm[ cc, ] == 0 )[ 1 ] ] <- k
    }
    ##rows.changed[ which.zero ] <- TRUE
  }

  invisible( list( r=rm, c=cm, ms=meme.scores ) )
}

id.duplicate.clusters <- function( scores=r.scores, cor.cutoff=0.9 ) {
  cors <- cor( scores[,], use="pairwise", method="pearson" )
  cors[ lower.tri( cors, diag=T ) ] <- NA
  tmp <- which( cors >= cor.cutoff, arr=T )
  cbind( tmp, cors[ tmp ] )
}

## Consolidate highly-correlated pairs of clusters into one (set row/col.membership of other one to zero)
consolidate.duplicate.clusters <- function( row.membership, col.membership, scores=r.scores, cor.cutoff=0.9,
                                           n.cutoff=5, motif=F, seq.type="upstream meme" ) {
  row.m <- row.membership; ms <- meme.scores ##$upstream ##col.m <- col.membership
  cors <- id.duplicate.clusters( scores, cor.cutoff )
  if ( nrow( cors ) <= 0 ) return( invisible( list( r=row.m, ms=meme.scores, scores=scores ) ) )
  cr <- max( cors[ ,3 ], na.rm=T )
  n.cut <- 1
  while( cr > cor.cutoff && ! is.infinite( cr ) && n.cut <= n.cutoff ) {
    tmp <- cors[ which( cors[ ,3 ] == cr ), 1:2 ]
    if ( any( get.rows( tmp[ 1 ] ) %in% get.rows( tmp[ 2 ] ) ) ) {
      ev1 <- if ( is.null( meme.scores[[ seq.type ]][[ tmp[ 1 ] ]]$meme.out ) ) Inf else
             min( sapply( meme.scores[[ seq.type ]][[ tmp[ 1 ] ]]$meme.out, "[[", "e.value" ), na.rm=T )
      ev2 <- if ( is.null( meme.scores[[ seq.type ]][[ tmp[ 2 ] ]]$meme.out ) ) Inf else
             min( sapply( meme.scores[[ seq.type ]][[ tmp[ 2 ] ]]$meme.out, "[[", "e.value" ), na.rm=T )
      ## Reorder - keep cluster that has best upstream motif or else keep cluster with more rows
      if ( ! ( is.infinite( ev1 ) && is.infinite( ev2 ) ) && ev2 < ev1 ) tmp <- tmp[ c( 2, 1 ) ]
      else if ( length( get.rows( tmp[ 1 ] ) ) < length( get.rows( tmp[ 2 ] ) ) ) tmp <- tmp[ c( 2, 1 ) ]

      row.m[ row.m == tmp[ 2 ] ] <- tmp[ 1 ]
      cat( "MERGING:", tmp, "\t", length( get.rows( tmp[ 1 ] ) ), length( get.rows( tmp[ 2 ] ) ), "\t", 
          length( unique( c( get.rows( tmp[ 1 ] ), get.rows( tmp[ 2 ] ) ) ) ), "\t", cr, "\n" )
      scores[ ,tmp[ 2 ] ] <- NA
      for ( tt in names( mot.weights ) ) {
        ms[[ tt ]][[ tmp[ 2 ] ]] <- list( iter=iter )
        if ( motif && sum( ! get.rows( tmp[ 1 ] ) %in% get.rows( tmp[ 2 ] ) ) > 0 )
          ms[[ tt ]][[ tmp[ 1 ] ]] <- try( meme.one.cluster( tmp[ 1 ], verbose=T, consens=meme.consensus, seq.type=tt ) ) ##, run.mast=run.mast ) )
      }
      n.cut <- n.cut + 1
    }
    ##cors[ ,tmp[ 2 ] ] <- cors[ tmp[ 2 ], ] <- NA
    cors[ which( cors[ ,3 ] == cr ), ] <- NA
    cr <- max( cors[ ,3 ], na.rm=T )
  }
  invisible( list( r=row.m, ms=ms, scores=scores ) )
}

## Can use this func to decide if any memberships shouln't exist
## this can be tweaked so that not all genes end up in exactly 2 biclusters
## Right now (by default) it does nothing (if rows=0, cols=0).
filter.updated.memberships <- function( row.membership, col.membership, rr.scores, cc.scores,
                                       quant.cutoff=c( rows=0, cols=0 ) ) {
  rm <- row.membership
  if ( quant.cutoff[ "rows" ] > 0 ) {
    ##if ( ! exists( "row.memb" ) ) row.memb <- t( apply( row.membership, 1, function( i ) 1:k.clust %in% i ) )
    qc <- quantile( rr.scores[,][ row.memb[,] == 1 ], prob=quant.cutoff[ "rows" ] )
    for ( i in 1:nrow( rm ) ) {
      tmp <- which( rm[ i, ] != 0 )
      rm[ i, tmp[ rr.scores[ i, rm[ i, tmp ] ] < qc ] ] <- 0
    }
  }  

  cm <- col.membership
  if ( quant.cutoff[ "cols" ] > 0 ) {
    ##if ( ! exists( "col.memb" ) ) col.memb <- t( apply( col.membership, 1, function( i ) 1:k.clust %in% i ) )
    qc <- quantile( cc.scores[,][ col.memb[,] == 1 ], prob=quant.cutoff[ "cols" ] )
    for ( i in 1:nrow( cm ) ) {
      tmp <- which( cm[ i, ] != 0 )
      cm[ i, tmp[ cc.scores[ i, cm[ i, tmp ] ] < qc ] ] <- 0
    }
  }

  ##invisible( list( r=rm, c=cm ) )
  NULL
}

## Hacky way to improve cluster in one swoop - add the best outside gene with a better score than the worst gene
##   already in, then remove that worst gene. Repeat until there are no outside genes better than any inside genes.
## Meme the cluster (TODO: during each iteration?); TODO: update row/mot/net/col scores too?
adjust.clust <- function( k, row.memb=get("row.membership"), expand.only=T, limit=100, ##motif=F, plot=F, 
                         ##scores="rr.scores", quant.cutoff=0.1, force.expand=0 ) { ##0.25 ) {
                         scores="r.scores", quant.cutoff=0.33, force.expand=0 ) {
  if ( scores == "rr.scores" || scores == "r.scores" ) {
    tmp <- get.combined.scores( quant=T )
    r.scores <- tmp$r
    if ( scores == "rr.scores" ) {
      scores <- get.density.scores( ks=1:k.clust )$r
      scores <- 1 - scores[,]
    } else {
      scores <- r.scores
    }
    rm( r.scores )
  } else {
    scores <- get( scores )
  }

  get.rows2 <- function( k, rm ) rownames( which( rm == k, arr=T ) )
  
  scores <- scores[,] ## In case it's an ff
  old.rows <- get.rows( k )
  if ( force.expand == 0 ) {
    wh <- names( which( scores[ which( ! attr( ratios, "rnames" ) %in% old.rows ), k ] <
                       quantile( scores[ old.rows, k ], quant.cutoff, na.rm=T ) ) )
  } else {
    expand.only <- TRUE
    wh <- names( sort( scores[ ! attr( ratios, "rnames" ) %in% old.rows, k ], decreasing=F )[ 1:force.expand ] )
  }
  if ( length( wh ) > limit ) { warning( "Surpassing limit." ); return( invisible( list( r=row.memb ) ) ) }
  else if ( length( wh ) <= 0 ) return( invisible( list( r=row.memb ) ) )
  tries <- 0
  while( length( wh ) > 0 && tries < 50 ) {
    wh2 <- names( which.max( scores[ wh, k ] ) )
    wh2.scores <- scores[ wh2, row.memb[ wh2, ] ]
    wh2a <- names( which.max( scores[ get.rows2( k, rm=row.memb ), k ] ) )
    for ( col in 1:ncol( row.memb ) ) if ( all( row.memb[ wh2, col ] == 0 ) ) break
    if ( col == ncol( row.memb ) && any( row.memb[ wh2, col ] != 0 ) ) {
      row.memb <- cbind( row.memb, rep( 0, nrow( row.memb ) ) ); col <- col + 1 }
    row.memb[ wh2, col ] <- k ##which.min( wh2.scores ) ] <- k
    if ( ! expand.only ) row.memb[ wh2a, row.memb[ wh2a, ] == k ] <- 0
    if ( force.expand == 0 ) {
      wh <- names( which( scores[ which( ! attr( ratios, "rnames" ) %in% get.rows2( k, rm=row.memb ) ), k ] <
                         quantile( scores[ get.rows2( k, rm=row.memb ), k ], quant.cutoff, na.rm=T ) ) )
    } else {
      wh <- wh[ ! wh %in% wh2 ]
    }
    if ( length( get.rows2( k, rm=row.memb ) ) > cluster.rows.allowed[ 2 ] ) break
    tries <- tries + 1
  }
  new.rows <- get.rows2( k, rm=row.memb )
  if ( any( ! new.rows %in% old.rows ) || any( ! old.rows %in% new.rows ) )
    cat( "ADJUSTED CLUSTER:", k, length( old.rows ), length( new.rows ), sum( ! old.rows %in% new.rows ), "\n" )
  ## if ( motif ) {
  ##   ms <- list()
  ##   for ( seq.types in names( meme.scores ) ) {
  ##     ms[[ seq.type ]] <- meme.one.cluster( new.rows, ms=meme.scores[[ seq.type ]][[ k ]], verbose=T,
  ##                                          consens=meme.consensus, seq.type=seq.type )
  ##   }
  ##   return( invisible( list( r=row.memb, ms=ms ) ) )
  ## }
  row.memb <- t( apply( row.memb, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
  row.memb <- row.memb[ ,apply( row.memb, 2, sum ) != 0, drop=F ]
  colnames( row.memb ) <- NULL
  invisible( list( r=row.memb ) )
}

## #!ifndef 
## ## Now we either:
## ##   add outside genes that are better than in-genes with scores at the 66% quantile (of in-gene scores) level.
## ##       (quant.cutoff > 0) or add N best outside genes (quant.cutoff > 1)
## ##   or remove inside genes that are worse than out-genes with scores at the 0.1% quantile (of out-gene scores) level.
## ##       (quant.cutoff < 0) or remove N worst inside genes (quant.cutoff < 1)
## adjust.clust.2 <- function( k, row.memb=get("row.membership"), limit=10,
##                          scores="r.scores", quant.cutoff=0.33 ) {
##   if ( scores == "rr.scores" ) {
##     if ( ! exists( "rr.scores" ) ) scores <- get.density.scores( ks=1:k.clust )$r
##     else scores <- get( scores )
##     scores <- 1 - scores[,]
##   } else {
##     scores <- get( scores )
##   }
##   scores <- scores[,] ## In case it's an ff
##   old.rows <- get.rows( k )
##   if ( quant.cutoff > 0 ) { ## Find genes to add to bicluster
##     if ( quant.cutoff < 1 ) wh <- names( which( scores[ which( ! attr( ratios, "rnames" ) %in% old.rows ), k ] <
##                        quantile( scores[ old.rows, k ], quant.cutoff, na.rm=T ) ) )
##     else wh <- names( sort( scores[ ! attr( ratios, "rnames" ) %in% old.rows, k ], decreasing=F )[ 1:quant.cutoff ] )
##   } else { ## Find genes to remove from bicluster
##     if ( quant.cutoff > -1 ) wh <- names( which( scores[ old.rows, k ] >
##                                                 quantile( scores[ which( ! attr( ratios, "rnames" ) %in%
##                                                                 old.rows ), k ], abs( quant.cutoff ), na.rm=T ) ) )
##     else wh <- names( sort( scores[ old.rows, k ], decreasing=T )[ 1:abs( quant.cutoff ) ] )
##   }
##   if ( length( wh ) > limit ) { warning( "Surpassing limit." ); return( invisible( list( r=row.memb ) ) ) }
##   else if ( length( wh ) <= 0 ) return( invisible( list( r=row.memb ) ) )
##   tries <- 0
##   rm <- row.memb + 0 ## make a copy
##   rm <- cbind( rm, rep( 0, nrow( rm ) ) )
##   while( length( wh ) > 0 && tries < 50 ) {
##     if ( quant.cutoff > 0 ) { ## add outside genes to this bicluster
##       rm[ wh, ncol( rm ) ] <- k
##     } else {
##       for ( w in wh ) rm[ w, which( rm[ w, ] == k ) ] <- 0
##     }
##     if ( length( get.rows( k, rm=rm ) ) > cluster.rows.allowed[ 2 ] ) break
##     tries <- tries + 1
##   }
##   new.rows <- get.rows( k, rm=rm )
##   if ( any( ! new.rows %in% old.rows ) || any( ! old.rows %in% new.rows ) )
##     cat( "ADJUSTED CLUSTER:", k, length( old.rows ), length( new.rows ), "\n" )
##   rm <- t( apply( rm, 1, function( i ) c( i[ i != 0 ], i[ i == 0 ] ) ) )
##   rm <- rm[ ,apply( rm, 2, sum ) != 0, drop=F ]
##   colnames( rm ) <- NULL
##   invisible( list( r=rm ) )
## }
## #!endif

adjust.all.clusters <- function( env, ks=1:env$k.clust, force.motif=T, ... ) {
  old.stats <- env$stats

  tmp <- env$row.col.membership.from.clusterStack( env$clusterStack )
  row.membership <- tmp$r; col.membership <- tmp$c
  
  mc <- env$get.parallel( length( ks ) )
  new.rm <- mc$apply( ks, function( k ) env$adjust.clust( k, row.membership, ... )$r )
  rm <- ##cbind( env$row.membership[,],
    do.call( cbind, new.rm ) ##)

  ## Consolidate all new row.membership matrices into one... note this only allows for cluster expansion!
  ##   TODO: need to look at each new.rm, see which genes' clusters were set to zero, and set it to zero in
  ##         the new rm.
  for ( i in 1:nrow( rm ) ) {
    tmp <- unique( rm[ i, rm[ i, ] != 0 ] )
    rm[ i, ] <- c( tmp, rep( 0, ncol( rm ) - length( tmp ) ) )
  }
  rm <- rm[ ,apply( rm, 2, sum ) != 0, drop=F ]
  colnames( rm ) <- NULL
  
  ##if ( any( dim( row.membership ) != dim( rm ) ) || any( row.membership != rm ) ) {
    ##env$row.membership <- rm
    ##attr( env$clusterStack, "iter" ) <- NULL ## force it to update
    env$clusterStack <- lapply( 1:env$k.clust, function( k )
                               list( rows=rownames( which( rm == k, arr=T ) ),
                                    cols=env$clusterStack[[ k ]]$cols ) )
    env$clusterStack <- env$get.clusterStack( ks=1:k.clust )
  ##env$iter <- env$n.iter + 1
  env$post.adjust <- FALSE ## HACK to prevent it from getting done again.
    env$cmonkey.one.iter( env, dont.update=T, force.row=T, force.col=T,
                         force.motif=if ( force.motif & ! no.genome.info ) "run.meme", force.net=T )
  ##}
  
##   tmp <- env$get.all.scores( force.row=T, force.col=T, force.motif=force.motif & ! no.genome.info, force.net=T )
##   env$row.scores <- tmp$r[,]; env$mot.scores <- tmp$m[,]; env$net.scores <- tmp$n[,]; env$col.scores <- tmp$c[,]
##   env$meme.scores <- tmp$ms
## #!ifndef 
##   for ( i in names( env$meme.scores ) ) {
##     for ( j in c( "all.pv", "all.ev" ) ) {
##       if ( ! is.null( env$meme.scores[[ i ]][[ j ]] ) && "ff" %in% class( env$meme.scores[[ i ]][[ j ]] ) ) {
##         if ( is.null( tmp[[ i ]] ) ) tmp[[ i ]] <- list()
##         tmp[[ i ]][[ j ]] <- env$meme.scores[[ i ]][[ j ]]
##         env$meme.scores[[ i ]][[ j ]] <- env$meme.scores[[ i ]][[ j ]][,]
##       }
##     }
##   }
## #!endif
##   env$row.memb <- t( apply( env$row.membership, 1, function( i ) 1:k.clust %in% i ) )
##   env$col.memb <- t( apply( env$col.membership, 1, function( i ) 1:k.clust %in% i ) )
  
##   tmp <- env$get.combined.scores()
##   env$r.scores <- tmp$r[,]; env$c.scores <- tmp$c[,]
##   rm( tmp )
##   env$clusterStack <- env$get.clusterStack( ks=1:k.clust, force=T )
  
##   env$stats <- rbind( old.stats, env$get.stats() )
  print( rbind( OLD=old.stats[ nrow( old.stats ), ], NEW=env$stats[ nrow( env$stats ), ] ) )
  invisible( env )
}

## ## TODO: include motif comparison via "motif.similarities.tomtom"
## compare.clusters <- function( k1, k2, scores=r.scores ) {
##   plot( scores[ ,k1 ], scores[ ,k2 ], pch=20, cex=0.5 ) ##, ## + 0.5 * attr( ratios, "rnames" ) %in% get.rows( k1 ),
##   points( scores[ get.rows( k1 ), k1 ], scores[ get.rows( k1 ), k2 ], col="red", cex=0.5, pch=20 )
##   points( scores[ get.rows( k2 ), k1 ], scores[ get.rows( k2 ), k2 ], col="green", cex=0.5, pch=20 )
##   points( scores[ get.rows( k1 )[ get.rows( k1 ) %in% get.rows( k2 ) ], k1 ],
##          scores[ get.rows( k2 )[ get.rows( k2 ) %in% get.rows( k1 ) ], k2 ], col="blue", cex=0.5, pch=20 )
##   cat( length( get.rows( k1 ) ), length( get.rows( k2 ) ), sum( get.rows( k1 ) %in% get.rows( k2 ) ), "\t",
##       cor( scores[ ,k1 ], scores[ ,k2 ], use="pairwise", method="pearson" ), "\n" )
## }
