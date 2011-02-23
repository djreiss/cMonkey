###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## Find # of genes per cluster with annoatated "long name" (e.g. "peroxisome") or short name (e.g. "PEX")
clusters.w.func <- function( func, ks=1:k.clust, short=F, max.rows=999, p.val=F ) {
  if ( p.val ) {
    long.names <- get.long.names( attr( ratios, "rnames" ), short=short )
    n2 <- length( grep( func, long.names, perl=T, ignore.case=T ) )
  }
  ##sapply( ks, function( i ) {  
  mc <- get.parallel( length( ks ) )
  unlist( mc$apply( ks, function( i ) { ##, ... ) {
    rows <- get.long.names( get.rows( i ), short=short )
    if ( ! p.val ) {
      if ( length( get.rows( i ) ) >= max.rows ) NA else
      length( grep( func, rows, perl=T, ignore.case=T ) )
    } else {
      phyper( length( grep( func, rows, perl=T, ignore.case=T ) ),
             n2, attr( ratios, "nrow" ) - n2, length( get.rows( i ) ), lower=F ) * length( ks ) ## bonferroni baby!
    }
  } ) ) ##, mc.cores=mc$par ) )
}

## Find # of genes per cluster out of given list of genes
clusters.w.genes <- function( genes, ks=1:k.clust, p.val=F ) {
  ##sapply( ks, function( i ) {
  mc <- get.parallel( length( ks ) )
  unlist( mc$apply( ks, function( i ) { ##, ... ) {
    rows <- get.rows( i )
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
                            sort=c("score.norm","score","resid","e.value1","e.value2","nrow") ) {
  ms <- meme.scores[[ seq.type ]]
  score <-
    sapply( 1:k.clust, function( k ) mean( r.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
      row.scaling[ iter ] + if ( ! is.null( mot.scores ) )
        sapply( 1:k.clust, function( k ) mean( mot.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
          mot.scaling[ iter ] else 0 + if ( ! is.null( net.scores ) )
            sapply( 1:k.clust, function( k ) mean( net.scores[ get.rows( k ), k ], na.rm=T, trim=0.01 ) ) *
              net.scaling[ iter ] else 0
  nrow <- tabulate( unlist( apply( row.membership, 1, unique ) ), k.clust )

  out <- data.frame( k=1:k.clust, nrow=nrow, score=score, ##score.norm=score.norm,
                    resid=sapply( 1:k.clust, cluster.resid, varNorm=F ), 
                    consensus1=sapply( 1:k.clust,
                      function( k ) if ( length( ms[[ k ]] ) <= 2 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 1 ]]$pssm ) ),
                    e.value1=sapply( 1:k.clust,
                      function( k ) if ( length( ms[[ k ]] ) <= 2 ) Inf else
                      ms[[ k ]]$meme.out[[ 1 ]]$e.value ),
                    consensus2=sapply( 1:k.clust,
                      function( k ) if ( length( ms[[ k ]] ) <= 2 ) "" else
                      if ( length( ms[[ k ]]$meme.out ) == 1 ) "" else
                      pssm.to.string( ms[[ k ]]$meme.out[[ 2 ]]$pssm ) ),
                    e.value2=sapply( 1:k.clust,
                      function( k ) if ( length( ms[[ k ]] ) <= 2 ) Inf else
                      if ( length( ms[[ k ]]$meme.out ) <= 1 ) Inf else
                      ms[[ k ]]$meme.out[[ 2 ]]$e.value )
                    )
  if ( all( out$consensus2 == "" ) ) out$consensus2 <- out$e.value2 <- NULL
  if ( ! is.na( sort[ 1 ] ) && sort[ 1 ] %in% colnames( out ) ) out <- out[ order( out[[ sort[ 1 ] ]] ), ]
  out
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

  ## note: need to do this afterwards, too:
  # e$row.memb <- t( apply( e$row.membership, 1, function( i ) 1:e$k.clust %in% i ) )
  # e$col.memb <- t( apply( e$col.membership, 1, function( i ) 1:e$k.clust %in% i ) )
  # e$clusterStack <- e$get.clusterStack( force=T )
  ##
  
  list( r=row.memb, c=col.memb )
}

## TODO: add another function to randomly seed clusters with no rows or no cols
re.seed.empty.clusters <- function( toosmall.r=cluster.rows.allowed[ 1 ], toosmall.c=0,
                                   n.r=cluster.rows.allowed[ 1 ] * 2, n.c=5 ) {
  ## TODO: for zero-row clusters, take a random gene(s) and assign it to this cluster.
  rm <- row.membership
  rats <- get.cluster.matrix()
  if ( any( tabulate( unlist( apply( rm, 1, unique ) ), k.clust ) <= toosmall.r ) ) {
    which.zero <- which( tabulate( unlist( apply( rm, 1, unique ) ), k.clust ) <= toosmall.r )
    cat( "These", length( which.zero ), "clusters have TOO FEW rows: ", which.zero, "\n" )
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
consolidate.duplicate.clusters <- function( scores=r.scores, cor.cutoff=0.9, n.cutoff=5, motif=F,
                                           seq.type="upstream meme" ) {
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
filter.updated.memberships <- function( quant.cutoff=c( rows=0, cols=0 ) ) {
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

## Code for reloading updated cMonkey package, updating existing env. obj. and restarting cmonkey run on that env:
##  unload("cMonkey");require(cMonkey);update.cmonkey.env(e);cmonkey(e,dont.init=T)
update.cmonkey.env <- function( object, ... ) { ## Update all funcs contained in env to latest from cmonkey package
  if ( file.exists( "cmonkey-funcs.R" ) ) {
    tmp.e <- new.env()
    sys.source( "cmonkey-funcs.R", envir=tmp.e ) ##cmonkey.env )
    sys.source( "cmonkey-postproc.R", envir=tmp.e ) ##cmonkey.env )
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

adjust.clust.2 <- function( k, row.memb=get("row.membership"), expand.only=T, plot=F, limit=100, ##motif=F, 
                           scores="r.scores", quant.cutoff=0.33, force.expand=0 ) {
  if ( scores == "rr.scores" ) {
    if ( ! exists( "rr.scores" ) ) scores <- get.density.scores( ks=1:k.clust )$r
    else scores <- get( scores )
    scores <- 1 - scores[,]
  } else if ( scores == "r.scores" ) {
    scores <- get.combined.scores()$r
  } else {
    scores <- get( scores )
  }
  scores <- scores[,] ## In case it's an ff
  old.rows <- get.rows( k )
  sc.in <- scores[ old.rows, k ]
  sc.out <- scores[ ! rownames( scores ) %in% old.rows, k ]
  pv <- t.test( sc.out, sc.in, alt='g' )$p.value
  to.add <- to.remove <- character()
  new.rows <- orig.rows <- old.rows
  for ( g in 1:10000 ) {
    g <- sample( rownames( scores ), 1 )
    ##new.rows <- orig.rows
    if ( g %in% old.rows ) new.rows <- new.rows[ old.rows != g ]
    else new.rows <- unique( c( old.rows, g ) )
    sc.in1 <- scores[ new.rows, k ]
    sc.out1 <- scores[ ! rownames( scores ) %in% new.rows, k ]
    pv1 <- t.test( sc.out1, sc.in1, alt='g' )$p.value
    if ( pv1 <= pv ) {
      cat( g, pv, pv1, g %in% old.rows, "\n" )
      if ( g %in% orig.rows ) to.remove <- unique( c( to.remove, g ) )
      else to.add <- unique( c( to.add, g ) )
      pv <- pv1
    } else {
      new.rows <- old.rows
    }
    old.rows <- new.rows
  }
  list( add=to.add, remove=to.remove )
}

## Hacky way to improve cluster in one swoop - add the best outside gene with a better score than the worst gene
##   already in, then remove that worst gene. Repeat until there are no outside genes better than any inside genes.
## Meme the cluster (TODO: during each iteration?); TODO: update row/mot/net/col scores too?
## Now we add outside genes that are better than in-genes with scores at the 66% quantile level.
adjust.clust <- function( k, row.memb=get("row.membership"), expand.only=T, plot=F, limit=100, ##motif=F, 
                         ##scores="rr.scores", quant.cutoff=0.1, force.expand=0 ) { ##0.25 ) {
                         scores="r.scores", quant.cutoff=0.33, force.expand=0 ) {
  if ( scores == "rr.scores" ) {
    if ( ! exists( "rr.scores" ) ) scores <- get.density.scores( ks=1:k.clust )$r
    else scores <- get( scores )
    scores <- 1 - scores[,]
  } else {
    scores <- get( scores )
  }
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
    wh2a <- names( which.max( scores[ get.rows( k, rm=row.memb ), k ] ) )
    for ( col in 1:ncol( row.memb ) ) if ( all( row.memb[ wh2, col ] == 0 ) ) break
    if ( col == ncol( row.memb ) && any( row.memb[ wh2, col ] != 0 ) ) {
      row.memb <- cbind( row.memb, rep( 0, nrow( row.memb ) ) ); col <- col + 1 }
    row.memb[ wh2, col ] <- k ##which.min( wh2.scores ) ] <- k
    if ( ! expand.only ) row.memb[ wh2a, row.memb[ wh2a, ] == k ] <- 0
    if ( force.expand == 0 ) {
      wh <- names( which( scores[ which( ! attr( ratios, "rnames" ) %in% get.rows( k, rm=row.memb ) ), k ] <
                         quantile( scores[ get.rows( k, rm=row.memb ), k ], quant.cutoff, na.rm=T ) ) )
    } else {
      wh <- wh[ ! wh %in% wh2 ]
    }
    if ( length( get.rows( k, rm=row.memb ) ) > cluster.rows.allowed[ 2 ] ) break
    tries <- tries + 1
  }
  new.rows <- get.rows( k, rm=row.memb )
  if ( any( ! new.rows %in% old.rows ) || any( ! old.rows %in% new.rows ) )
    cat( "ADJUSTED CLUSTER:", k, length( old.rows ), length( new.rows ), "\n" )
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

adjust.all.clusters <- function( env, ks=1:env$k.clust, force.motif=T, ... ) {
  old.stats <- env$stats
  
  mc <- env$get.parallel( length( ks ) )
  new.rm <- mc$apply( ks, function( k ) env$adjust.clust( k, env$row.membership, ... )$r )
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
  
  if ( any( dim( env$row.membership ) != dim( rm ) ) || any( env$row.membership != rm ) ) {
    env$row.membership <- rm
    attr( env$clusterStack, "iter" ) <- NULL ## force it to update
    env$cmonkey.one.iter( env, dont.update=T, force.row=T, force.col=T,
                         force.motif=if ( force.motif & ! no.genome.info ) "run.meme", force.net=T )
  }
  
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

## TODO: include motif comparison via "motif.similarities.tomtom"
compare.clusters <- function( k1, k2, scores=r.scores ) {
  plot( scores[ ,k1 ], scores[ ,k2 ], pch=20, cex=0.5 ) ##, ## + 0.5 * attr( ratios, "rnames" ) %in% get.rows( k1 ),
  points( scores[ get.rows( k1 ), k1 ], scores[ get.rows( k1 ), k2 ], col="red", cex=0.5, pch=20 )
  points( scores[ get.rows( k2 ), k1 ], scores[ get.rows( k2 ), k2 ], col="green", cex=0.5, pch=20 )
  points( scores[ get.rows( k1 )[ get.rows( k1 ) %in% get.rows( k2 ) ], k1 ],
         scores[ get.rows( k2 )[ get.rows( k2 ) %in% get.rows( k1 ) ], k2 ], col="blue", cex=0.5, pch=20 )
  cat( length( get.rows( k1 ) ), length( get.rows( k2 ) ), sum( get.rows( k1 ) %in% get.rows( k2 ) ), "\t",
      cor( scores[ ,k1 ], scores[ ,k2 ], use="pairwise", method="pearson" ), "\n" )
}
