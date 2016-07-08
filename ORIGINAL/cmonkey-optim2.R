###################################################################################
## cMonkey - version 5, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss.isb@gmail.com.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

cm.version <- "5.0.0"

update.all.clusters <- function( env, dont.update=F, ... ) {
  mc <- get.parallel( k.clust )
  all.scores <- mc$apply( 1:k.clust, get.all.scores, return.scores=F, densities=T, verbose=F, members=T,
                         force.motif=F, ... )

  if ( ! dont.update ) {
    for ( k in 1:length( all.scores ) ) {
      if ( ! is.null( all.scores[[ k ]]$ms ) ) {
        for ( i in names( all.scores[[ k ]]$ms ) ) if ( ! is.null( all.scores[[ k ]]$ms[[ i ]] ) )
          env$meme.scores[[ i ]][[ k ]] <- all.scores[[ k ]]$ms[[ i ]]
      }
      env$clusterStack[[ k ]]$rows <- all.scores[[ k ]]$members$r
      env$clusterStack[[ k ]]$nrows <- length( all.scores[[ k ]]$members$r )
      env$clusterStack[[ k ]]$cols <- all.scores[[ k ]]$members$c
      env$clusterStack[[ k ]]$ncols <- length( all.scores[[ k ]]$members$c )
    }
    env$clusterStack <- env$get.clusterStack( ks=1:k.clust, force=T )
    ##env$iter <- env$iter + 1
  }
  ## print(env$clusterStack[[68]]$rows)
  ## ##dev.set(2);env$plotClust(68)
  ## env$stats <- rbind( env$stats, env$get.stats() )
  ## cat( organism, as.matrix( env$stats[ nrow( env$stats ), ] ), '\n' )
  ## if(env$iter%%10==9){dev.set(2);env$plotClust(68)}
  ## if(env$iter%%5==0){dev.set(3);env$plotStats()}
  env
}

## Same as previous but get it straight from the bicluster's pv.ev
get.motif.scores <- function( k, out.ms, for.rows="all" ) { ##m=meme.scores$upstream[[ k ]], for.rows="all" ) { 
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k ) 
  else rows <- k
  if ( for.rows[ 1 ] == "all" ) for.rows <- attr( ratios, "rnames" )
  if ( length( rows ) <= 1 || is.null( out.ms$pv.ev ) ) ##is.null( m ) || is.null( m$pv.ev ) )
    return( rep( NA, length( for.rows ) ) )
  m.scores <- rep( NA, attr( ratios, "nrow" ) ); names( m.scores ) <- attr( ratios, "rnames" )
  tmp <- out.ms$pv.ev[[ 1 ]][ ,"p.value" ]; names( tmp ) <- rownames( out.ms$pv.ev[[ 1 ]] )
  m.scores[ names( tmp ) ] <- tmp
  m.scores <- log( m.scores )
  m.scores <- m.scores[ for.rows ] ##; names( m.scores ) <- for.rows
  return( m.scores )
}

get.all.scores <- function( k, return.scores=F, densities=F, verbose=F, force.motif=F, allow.motif=T,
                           members=F, remove.nas=T, plot=F, ... ) {
  if ( verbose ) cat( "Computing scores for cluster:", k, "\n" )

  row.memb <- attr( ratios, "rnames" ) %in% get.rows( k )
  names( row.memb ) <- attr( ratios, "rnames" )
  col.memb <- attr( ratios, "cnames" ) %in% get.cols( k )
  names( col.memb ) <- attr( ratios, "cnames" )

  x <- NULL
  for ( i in names( ratios ) ) { 
    if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
    r <- get.row.scores( k, ratios=ratios[[ i ]], method=row.score.func, ... )
    if ( remove.nas ) { ## No measurement is NOT an NA - it is an observation
      tmp <- is.infinite( r ) | is.na( r )
      if ( any( tmp ) ) r[ tmp ] <- quantile( r[ row.memb[ names( r ) ] & ! tmp ], 0.95 ) ##0
    }
    x <- cbind( x, r ); colnames( x )[ ncol( x ) ] <- paste( "ratios", i, sep="." )
  }
  
  y <- NULL
  for ( i in names( ratios ) ) { 
    if ( row.weights[ i ] == 0 || is.na( row.weights[ i ] ) ) next
    c <- get.col.scores( k, ratios=ratios[[ i ]], method=col.score.func, ... )
    if ( all( is.na( c ) ) ) names( c ) <- colnames( ratios[[ i ]] )
    if ( remove.nas ) { ## No measurement is NOT an NA - it is an observation
      tmp <- is.infinite( c ) | is.na( c )
      if ( any( tmp ) ) c[ tmp ] <- quantile( c[ col.memb[ names( c ) ] & ! tmp ], 0.95 ) ##0
    }
    y <- cbind( y, c ); colnames( y )[ ncol( y ) ] <- paste( "ratios", i, sep="." )
  }
  
  out.ms <- list()
  for ( i in names( mot.weights ) ) {
    if ( mot.weights[ i ] == 0 || is.na( mot.weights[ i ] ) ) next
    out <- meme.scores[[ i ]][[ k ]]
    if ( allow.motif == TRUE && ( force.motif == TRUE || force.motif == "run.meme" ||
           ( mot.scaling[ iter ] > 0 && ! is.na( meme.iters[[ i ]][ 1 ] ) &&
            iter %in% meme.iters[[ i ]] && exists( "genome.info" ) && ! no.genome.info ) ) ) {
      out <- motif.one.cluster( k, seq.type=i, verbose=verbose, force=force.motif, ... )
      if ( is.null( out ) || class( out ) == "try-error" ||
          out$k != k || ( ! is.null( out$iter ) && out$iter != iter ) ) {
        out <- try( motif.one.cluster( k, seq.type=i, verbose=verbose, force=force.motif, ... ) ) 

        if ( class( out ) == "try-error" ) ## try one last time
          out <- try( motif.one.cluster( k, seq.type=i, verbose=verbose, force=force.motif, ... ) ) 
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
                mn <- mean( log10( out$pv.ev[[ ind ]][ rownames( out$pv.ev[[ ind ]] ) %in% get.rows( k ),
                                                      "p.value" ] ), na.rm=T )
              else mn <- mean( log10( out$pv.ev[[ ind ]]$pvals ), na.rm=T )
            } else {
              mn <- 'Inf'
            }
            cat( k, if ( attr( out$meme.out, "is.pal" ) ) "pal" else "non", 
                sapply( out$meme.out[ 1:min( 3, length( out$meme.out ) ) ], "[[", "e.value" ), mn,
                '\t', pssm.to.string( out$meme.out[[ 1 ]]$pssm ), "\n" )
          }
        }
      }
      out$iter <- iter
      out$k <- k
    }

    m <- get.motif.scores( k, out, ... )
    if ( remove.nas ) m[ is.infinite( m ) | is.na( m ) ] <- 0 ## No measurement is NOT an NA - it is an observation
    if ( ! is.null( out ) && ( is.null( out$iter ) || out$iter == iter ) ) out.ms[[ i ]] <- out
    x <- cbind( x, m ); colnames( x )[ ncol( x ) ] <- paste( "motif", i, sep="." )
  }

  for ( i in names( networks ) ) { 
    if ( net.weights[ i ] == 0 || is.na( net.weights[ i ] ) ) next
    if ( nrow( subset( networks[[ i ]], protein1 %in% attr( ratios, "rnames" ) & protein2 %in%
                      attr( ratios, "rnames" ) ) ) <= 0 ) next
    n <- get.network.scores( k, net=networks[[ i ]] )
    if ( remove.nas ) n[ is.infinite( n ) | is.na( n ) ] <- 0 ## No measurement is NOT an NA - it is an observation
    x <- cbind( x, n ); colnames( x )[ ncol( x ) ] <- paste( "network", i, sep="." )
  }

  rm( r, c, m, n )

  x.d <- y.d <- NULL
  ## Compute the density scores (ecdf) for each separate input score, used for sampling new gene/condition members
  if ( densities || members ) {
    p <- rep( 1 / nrow( x ), nrow( x ) )
    x.d <- apply( x, 2, function( xx ) {
      if ( ! all( is.na( xx ) ) && ! all( xx == xx[ 1 ] ) ) {
        fun <- ecdf( c( xx[ row.memb ], max( xx[ row.memb ] ) + 1e-5 ) )
        p <- 1 - fun( xx )
        ## fun2 <- ecdf( c( xx[ ! row.memb ], max( xx[ ! row.memb ] ) + 1e-5 ) )
        ## p <- fun( xx ) / fun2( xx )
        ## if ( any( is.infinite( p ) ) )
        ##   p[ is.infinite( p ) | is.na( p ) ] <- max( p[ ! is.infinite( p ) & ! is.na( p ) ], na.rm=T )
        ## p <- p / max( p, na.rm=T )
      }
      p
    } )
    rownames( x.d ) <- rownames( x )

    p <- rep( 1 / nrow( y ), nrow( y ) )
    y.d <- apply( y, 2, function( yy ) {
      if ( ! all( is.na( yy ) ) && ! all( yy == yy[ 1 ] ) ) {
        fun <- ecdf( c( yy[ col.memb ], max( yy[ col.memb ] ) + 1e-5 ) )
        p <- 1 - fun( yy )
        ## fun2 <- ecdf( c( yy[ ! col.memb ], max( yy[ ! col.memb ] ) + 1e-5 ) )
        ## p <- fun( yy ) / fun2( yy )
        ## if ( any( is.infinite( p ) ) )
        ##   p[ is.infinite( p ) | is.na( p ) ] <- max( p[ ! is.infinite( p ) & ! is.na( p ) ], na.rm=T )
        ## p <- p / max( p, na.rm=T )
      }
      p
    } )
    rownames( y.d ) <- rownames( y )
  }

  scores <- list( r=x, c=y, r.d=x.d, c.d=y.d, ms=out.ms, k=k )
  if ( plot ) plot.cluster.scores( scores )

  ## Compute new members here... so it can be done in parallelization with rest of scores
  if ( members ) scores$members <- get.cluster.members( scores, ... )
  if ( ! return.scores ) scores$r <- scores$c <- scores$r.d <- scores$c.d <- NULL
  return( scores )
}

plot.cluster.scores <- function( scores ) {
  row.memb <- attr( ratios, "rnames" ) %in% get.rows( scores$k )
  names( row.memb ) <- attr( ratios, "rnames" )

  par( mfrow=c( 2, 2 ) )
  for ( i in colnames( scores$r )[ 1:4 ] ) {
    h <- hist( scores$r[ ,i ], breaks=200, main=i )
    hist( rep( scores$r[ row.memb, i ], 10 ), breaks=h$breaks, add=T, col='red', border='red' )
    p <- scores$r.d[ ,i ]
    lines( sort( scores$r[ ,i ] ), p[ order( scores$r[ ,i ] ) ] / max( p ) * max( h$counts ) / 2, col='green' )
  }
}

attr( get.all.scores, "version" ) <- 2

get.cluster.members <- function( scores, weights="calc", pseudo=0.01, TEMP="calc", max.change=c(rows=5,cols=10),
                                count.power=c(rows=8,cols=2) ) {

  if ( TEMP == "calc" ) TEMP <- seq( 0.15, 0.05, length=n.iter )[ iter ]

  wts <- rep( 0, ncol( scores$r.d ) )
  if ( ! is.na( weights ) && weights == "calc" ) {
    cn <- sapply( strsplit( colnames( scores$r.d ), ".", fixed=T ), "[", 1 )
    cn2 <- sapply( lapply( strsplit( colnames( scores$r.d ), ".", fixed=T ), "[", -1 ), paste, collapse="." )
    wts[ cn == "ratios" ] <- row.scaling[ iter ] * row.weights[ cn2[ cn == "ratios" ] ]
    wts[ cn == "motif" ] <- mot.scaling[ iter ] * mot.weights[ cn2[ cn == "motif" ] ]
    wts[ cn == "network" ] <- net.scaling[ iter ] * net.weights[ cn2[ cn == "network" ] ]
    wts[ is.na( wts ) ] <- 0
    wts <- wts / sum( wts, na.rm=T )
  }
  if ( all( wts == 0 | is.na( wts ) ) ) wts[] <- 1
  
  probs <- apply( scores$r.d, 1, weighted.mean, w=wts, na.rm=T ) ## + pseudo
  probs <- probs / max( probs, na.rm=T )
  rows <- get.rows( scores$k )

  ## Adjust prob based on poisson func. based on how many clusters each row is already in (counts)
  if ( length( rows ) > 0 && ! all( is.na( rows ) ) && ! is.na( count.power ) && ! is.na( count.power[ "rows" ] ) ) {
    counts <- table( unlist( lapply( clusterStack, "[[", "rows" ) ) )[ names( probs ) ]
    counts[ is.na( counts ) ] <- 0; names( counts ) <- names( probs )
    count.prob <- ppois( counts, n.clust.per.row, lower=F ) / ppois( 1, n.clust.per.row, lower=F )
    count.prob[ counts <= 1 ] <- 1.1 ## prob of adding; don't make too high
    count.prob[ rows ] <- ppois( counts[ rows ], n.clust.per.row, lower=T ) / ppois( 2, n.clust.per.row, lower=T )
    ##count.prob <- count.prob / max( count.prob, na.rm=T ) ##[ count.prob > 1 ] <- 1
    count.prob <- count.prob^count.power[ "rows" ]
  }
  
  if ( FALSE ) { ##iter >= n.iter ) {
    rows <- unique( c( rows, names( which( probs >= 0.5 ) ) ) )
  } else {
    probs.add <-  exp( -( 1 - probs ) / TEMP )
    probs.drop <- exp( -      probs   / TEMP )

    row.memb <- attr( ratios, "rnames" ) %in% rows
    names( row.memb ) <- attr( ratios, "rnames" )
    sample.probs <- ifelse( row.memb, probs.drop, probs.add ) 
    if ( exists( "count.prob" ) ) sample.probs <- sample.probs * count.prob
    ##sample.probs <- sample.probs / max( sample.probs )
    n.row.per.clust <- n.clust.per.row / k.clust * length( row.memb ) ## balance probs to get right # of rows
    balance <- sum( sample.probs[ row.memb ] ) / sum( sample.probs[ ! row.memb ] ) * n.row.per.clust / sum( row.memb )
    if ( is.infinite( balance ) || is.na( balance ) ) balance <- 1
    if ( balance >= 1 ) sample.probs[ ! row.memb ] <- sample.probs[ ! row.memb ] * balance
    else sample.probs[ row.memb ] <- sample.probs[ row.memb ] / balance
    ##sample.probs <- sample.probs / max( sample.probs )
    allowed.moves <- sample.probs[ sample.probs > runif( length( sample.probs ) ) ]
    if ( length( allowed.moves ) > max.change[ 'rows' ] ) allowed.moves <- ##sort( allowed.moves, decreasing=T )[ 1:max.change ]
      sample( allowed.moves, max.change[ 'rows' ], prob=sample.probs[ names( allowed.moves ) ] )
    rows <- c( rows, names( row.memb[ names( allowed.moves ) ][ row.memb[ names( allowed.moves ) ] == FALSE ] ) )
    rows <- rows[ ! rows %in% names( row.memb[ names( allowed.moves ) ][ row.memb[ names( allowed.moves ) ] == TRUE ] ) ]
  }
  
  wts <- rep( 0, ncol( scores$c.d ) )
  if ( ! is.na( weights ) && weights == "calc" ) {
    cn <- sapply( strsplit( colnames( scores$c.d ), ".", fixed=T ), "[", 1 )
    cn2 <- sapply( lapply( strsplit( colnames( scores$c.d ), ".", fixed=T ), "[", -1 ), paste, collapse="." )
    wts[ cn == "ratios" ] <- row.scaling[ iter ] * row.weights[ cn2[ cn == "ratios" ] ]
    wts[ cn == "motif" ] <- mot.scaling[ iter ] * mot.weights[ cn2[ cn == "motif" ] ]
    wts[ cn == "network" ] <- net.scaling[ iter ] * net.weights[ cn2[ cn == "network" ] ]
    wts[ is.na( wts ) ] <- 0
    wts <- wts / sum( wts, na.rm=T )
  }
  if ( all( wts == 0 | is.na( wts ) ) ) wts[] <- 1
  
  probs <- apply( scores$c.d, 1, weighted.mean, w=wts, na.rm=T ) ## + pseudo
  probs <- probs / max( probs, na.rm=T )
  cols <- get.cols( scores$k )

  ## Adjust prob based on poisson func. based on how many clusters each col is already in (counts)
  if ( length( cols ) > 0 && ! all( is.na( cols ) ) && ! is.na( count.power ) && ! is.na( count.power[ "cols" ] ) ) {
    counts <- table( unlist( lapply( clusterStack, "[[", "cols" ) ) )[ names( probs ) ]
    counts[ is.na( counts ) ] <- 0; names( counts ) <- names( probs )
    count.prob <- ppois( counts, n.clust.per.col, lower=F ) / ppois( 2, n.clust.per.col, lower=F )
    count.prob[ counts <= 1 ] <- 1.1
    count.prob[ cols ] <- ppois( counts[ cols ], n.clust.per.col, lower=T ) /
      ppois( n.clust.per.col, n.clust.per.col, lower=T )
    ##count.prob <- count.prob / max( count.prob, na.rm=T ) ##[ count.prob > 1 ] <- 1
    count.prob <- count.prob^count.power[ "cols" ]
  }
  
  if ( iter >= n.iter ) {
    cols <- unique( c( cols, names( which( probs >= 0.5 ) ) ) )
  } else {
    probs.add  <- exp( -( 1 - probs ) / TEMP )
    probs.drop <- exp( -      probs   / TEMP )

    col.memb <- attr( ratios, "cnames" ) %in% cols
    names( col.memb ) <- attr( ratios, "cnames" )
    sample.probs <- ifelse( col.memb, probs.drop, probs.add )
    if ( exists( "count.prob" ) ) sample.probs <- sample.probs * count.prob
    ##sample.probs <- sample.probs / max( sample.probs )
    n.col.per.clust <- n.clust.per.col / k.clust * length( col.memb ) ## balance probs to get right # of cols
    balance <- sum( sample.probs[ col.memb ] ) / sum( sample.probs[ ! col.memb ] ) * n.col.per.clust / sum( col.memb )
    if ( is.infinite( balance ) || is.na( balance ) ) balance <- 1
    if ( balance >= 1 ) sample.probs[ ! col.memb ] <- sample.probs[ ! col.memb ] * balance
    else sample.probs[ col.memb ] <- sample.probs[ col.memb ] / balance
    ##sample.probs <- sample.probs / max( sample.probs )
    allowed.moves <- sample.probs[ sample.probs > runif( length( sample.probs ) ) ]
    if ( length( allowed.moves ) > max.change[ 'cols' ] ) allowed.moves <- ##sort( allowed.moves, decreasing=T )[ 1:max.change ]
      sample( allowed.moves, max.change[ 'cols' ], prob=sample.probs[ names( allowed.moves ) ] )
    cols <- c( cols, names( col.memb[ names( allowed.moves ) ][ col.memb[ names( allowed.moves ) ] == FALSE ] ) )
    cols <- cols[ ! cols %in% names( col.memb[ names( allowed.moves ) ][ col.memb[ names( allowed.moves ) ] == TRUE ] ) ]
  }  
  list( r=rows, c=cols )
}

## Get scores for every bicluster and format into the old 'row.scores', etc. matrices that were used
## in cmonkey-optim1.R (pre-5.0 version). This is for backward compatibility, mostly with get.stats() and
## plotStats().
get.old.scores.matrices <- function( ks=1:k.clust ) {
  mc <- get.parallel( length( ks ) )
  all.scores <- mc$apply( ks, get.all.scores, return.scores=T, densities=F, verbose=F, members=F,
                         force.motif=F, allow.motif=F )
  row.scores <- sapply( all.scores, function( i ) i$r[ ,"ratios.ratios" ] )
  mot.scores <- lapply( all.scores, function( i ) i$r[ ,grep( "^motif\\.", colnames( i$r ) ), drop=F ] )
  if ( ! is.null( unlist( mot.scores ) ) ) {
    mot.scores <- lapply( colnames( mot.scores[[ 1 ]] ), function( i )
                         sapply( mot.scores, function( j ) j[ ,i ] ) )
    if ( length( mot.scores ) > 1 )
      for ( m in 2:length( mot.scores ) ) mot.scores[[ 1 ]] <- mot.scores[[ 1 ]] + mot.scores[[ m ]]
    mot.scores <- mot.scores[[ 1 ]]
  } else rm( mot.scores )
  net.scores <- lapply( all.scores, function( i ) i$r[ ,grep( "^network\\.", colnames( i$r ) ), drop=F ] )
  if ( ! is.null( unlist( net.scores ) ) ) {
    net.scores <- lapply( colnames( net.scores[[ 1 ]] ), function( i )
                         sapply( net.scores, function( j ) j[ ,i ] ) )
    if ( length( net.scores ) > 1 )
      for ( m in 2:length( net.scores ) ) net.scores[[ 1 ]] <- net.scores[[ 1 ]] + net.scores[[ m ]]
    net.scores <- net.scores[[ 1 ]]
  } else rm( net.scores )
  col.scores <- sapply( all.scores, function( i ) i$c[ ,"ratios.ratios" ] )
  list( r=row.scores, m=mot.scores, n=net.scores, c=col.scores )
}

## Some code taken from ~/Documents/scratch/halo/biclust-v3/R_scripts/old/fit-cluster.R
## Remember to source cmonkey-optim2.R into the env to override the new get.all.scores(k) function !!!
test.fit.cluster <- function( k, verbose=F, plot=F, ... ) {
  scores <- get.all.scores( k, verbose, ... )
  x <- scores$r
  
  row.memb <- attr( ratios, "rnames" ) %in% get.rows( k )
  names( row.memb ) <- attr( ratios, "rnames" )

  y <- as.factor( row.memb )
  xx <- x[ ,apply( x[ y == TRUE, ], 2, sum ) != 0 ]

  if ( plot ) {
    len <- 25
    xp <- expand.grid( as.data.frame( apply( x, 2, function( i ) seq( min(i), max(i), length=len ) ) ) )
    xxp <- expand.grid( as.data.frame( apply( xx, 2, function( i ) seq( min(i), max(i), length=len ) ) ) )
  }
  
  out0 <- glm( y ~ . - 1, data=as.data.frame( x ), family='binomial', ... )
  if ( verbose ) print( summary( out0 ) )
  prob0 <- predict( out0, type="response" ); names( prob0 ) <- row.names( row.scores )
  if ( plot ) prob0p <- predict( out0, newdata=as.data.frame( xp ), type="response" )
  require( varSelRF )
  out1 <- varSelRF( x, y, ntree=5000, ntreeIterat=2000, vars.drop.frac=0.5, whole.range=F, keep.forest=T, ... )
  if ( verbose ) print( out1 )
  prob1 <- predict( out1$rf.model, type="prob", newdata=subset( x, select=out1$selected.vars ) )[ ,2 ]
  if ( plot ) prob1p <- predict( out1$rf.model, type="prob", newdata=subset( xp, select=out1$selected.vars ) )[ ,2 ]
  require( brglm )
  out2 <- brglm( y ~ . - 1, data=as.data.frame( x ), ... )
  if ( verbose ) print( summary( out2 ) )
  prob2 <- predict( out2, type="response" ); names( prob2 ) <- row.names( row.scores )
  if ( plot ) prob2p <- predict( out2, newdata=as.data.frame( xp ), type="response" )
  require( nnet )
  out3 <- nnet( y ~ . - 1, data=as.data.frame( x ), skip=T, softmax=F, size=3, decay=0.1, maxit=1000, tra=F, ... )
  prob3 <- predict( out3 )[ ,1 ]
  if ( plot ) prob3p <- predict( out3, newdata=as.data.frame( xp ) )
  require( MASS )
  prob4 <- prob5 <- rep( NA, length( prob1 ) )
  out4 <- try( qda( y ~ . - 1, data=as.data.frame( xx ), method="mle", ... ) )
  if ( ! "try-error" %in% class( out4 ) ) {
    prob4 <- predict( out4 )$posterior[ ,2 ]
    if ( plot ) prob4p <- predict( out4, newdata=as.data.frame( xp ) )$posterior[ ,2 ]
  }
  out5 <- try( lda( y ~ . - 1, data=as.data.frame( xx ), method="mle", ... ) )
  if ( ! "try-error" %in% class( out5 ) ) {
    prob5 <- predict( out5 )$posterior[ ,2 ]
    if ( plot ) prob5p <- predict( out5, newdata=as.data.frame( xp ) )$posterior[ ,2 ]
  }
  require( class )
  out6 <- knn.cv( x, as.integer( y ), k=attr(ratios,"nrow")/k.clust*n.clust.per.row, prob=T, use.all=T, ... )
  prob6 <- as.numeric( attr( out6, "prob" ) )
  if ( plot ) {
    prob6p <- knn( x, xp, as.integer( y ), k=attr(ratios,"nrow")/k.clust*n.clust.per.row, prob=T, use.all=T, ... )
    prob6p <- as.numeric( attr( prob6p, "prob" ) )
  }
  require( e1071 )
  out7 <- svm( y ~ . - 1, data=as.data.frame( cbind( x, y=as.integer( y ) - 1 ) ), probability=T, scale=F, ... )
  prob7 <- predict( out7, as.data.frame( x ), prob=T )
  if ( plot ) prob7p <- predict( out7, as.data.frame( xp ), prob=T )
  require( mda )
  out8 <- fda( y ~ . - 1, data=as.data.frame( cbind( x, y=as.integer( y ) - 1 ) ), ... )
  prob8 <- predict( out8, type="posterior" )[ ,2 ]
  if ( plot ) prob8p <- predict( out8, as.data.frame( xp ), type="posterior" )[ ,2 ]
  out9 <- mda( y ~ . - 1, data=as.data.frame( cbind( x, y=as.integer( y ) - 1 ) ), ... )
  prob9 <- predict( out9, newdata=as.data.frame( x ), type="posterior" )[ ,2 ]
  if ( plot ) prob9p <- predict( out9, as.data.frame( xp ), type="posterior" )[ ,2 ]
  
  probs <- cbind( prob0, prob1, prob2, prob3, prob4, prob5, prob6, prob7, prob8, prob9 )

  if ( plot ) {
    ## plot( x[ ,1 ], apply( probs, 1, sum ), col=as.integer( y ) )
    xp2 <- apply( x, 2, function( i ) seq( min(i), max(i), length=len ) )
    probsp <- cbind( prob0p, prob1p, prob2p, prob3p, prob4p, prob5p, prob6p, prob7p, prob8p, prob9p )
    probsp.tot <- array( apply( probsp, 1, mean, na.rm=T ), dim=rep( len, ncol( xp ) ) )
    comb <- t( combn( 1:ncol( x ), 2 ) )
    par( mfrow=c( ceiling( sqrt( nrow( comb ) ) ), floor( sqrt( nrow( comb ) ) ) ) )
    for ( i in 1:nrow( comb ) ) {
      i1 <- comb[ i, 1 ]; i2 <- comb[ i, 2 ]
      plot( x[ ,i1 ], x[ ,i2 ], xlab=colnames( x )[ i1 ], ylab=colnames( x )[ i2 ], col=as.integer( y ),
           pch=20, cex=0.5 )
      tmp <- aperm( probsp.tot, c( i1, i2, ( 1:ncol( x ) )[ ! 1:ncol( x ) %in% comb[ i, ] ] ) )
      contour( xp2[ ,i1 ], xp2[ ,i2 ], tmp[ ,,1,1 ], add=T, labex=0 ) ##levels=c( 0.7, 0.8, 0.9 ), 
    }
  }

  ## good.methods = c(
  ##  "fit.cluster.lda.qda( cl, which='qda', plot=F, predict=T )",
  ##  "fit.cluster.lda.qda( cl, which='lda', plot=F, predict=T )",
  ##  "knn.cluster( cl, plot=F, predict=T, k=cl$nrows/2 )",
  ##  ##"knn.cluster( cl, plot=F, predict=T, k=cl$nrows/3 )",
  ##  ##"knn.cluster( cl, plot=F, predict=T, k=cl$nrows/4 )",
  ## "fit.cluster.nnet( cl, plot=F, predict=T, repeats=10, skip=T, softmax=F, size=3, decay=0.02, maxit=1000, tra=F )",
  ## "fit.cluster.nnet( cl, plot=F, predict=T, repeats=10, skip=T, softmax=F, size=2, decay=0.02, maxit=1000, tra=F )",
  ## "fit.cluster.nnet( cl, plot=F, predict=T, repeats=10, skip=T, softmax=F, size=4, decay=0.05, maxit=1000, tra=F )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='polynomial', net.names='none', degree=1, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='polynomial', net.names='none', degree=2, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='polynomial', net.names='none', degree=3, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='polynomial', net.names='none', degree=3, weight=0.5, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='radial', net.names='none', gamma=0.1, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='radial', net.names='none', gamma=0.5, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='radial', net.names='none', gamma=1, weight=0.1, cutoff=0.95 )",
  ##  "fit.cluster.svm( cl, plot=F, predict=T, kernel='radial', net.names='none', gamma=0.5, weight=0.01, cutoff=0.95 )",
  ##  "fit.cluster.fda( cl, plot=F, predict=T, cutoff=0.95, weight=0.1, method=mars )",
  ##  "fit.cluster.fda( cl, plot=F, predict=T, cutoff=0.95, weight=0.1, method=polyreg )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight=0.1, method=mars, subclasses=c(2,1) )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight=0.1, method=mars, subclasses=c(5,2) )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight=0.1, method=mars, subclasses=c(15,3) )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight=0.5, method=mars, subclasses=c(2,1) )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight=0.5, method=mars, subclasses=c(5,2) )",
  ##  "fit.cluster.mda( cl, plot=F, predict=T, cutoff=0.95, weight='relative', method=mars, subclasses=c(5,2) )",
  ##  "fit.cluster.logist.new( cl, plot=F, predict=T, joint=T, row.cutoff=0.8 )",
  ##  "fit.cluster.logist.new( cl, plot=F, predict=T, joint=F, row.cutoff=0.5 )"
  ##  )
  invisible( probs )  
}
