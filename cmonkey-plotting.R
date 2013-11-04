###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## resids.in.out <- function( cluster ) {
##   out <- rep( 1, 2 ); names( out ) <- c( "in", "out" )
##   if ( length( cluster$nrows ) <= 0 || is.null( cluster$nrows ) || cluster$nrows <= 1 ) return( out )
##   if ( cluster$ncol > 1 )
##     out[ "in" ] <- residual.submatrix( ratios[ cluster$rows, cluster$cols ] )
##   if ( sum( ! attr( ratios, "cnames" ) %in% cluster$cols ) > 1 )
##     out[ "out" ] <- residual.submatrix( ratios[ cluster$rows, attr( ratios, "cnames" )[ ! attr( ratios, "cnames" ) %in% cluster$cols ] ] )
##   return( out )
## }##,

###############################
plotCluster <- function( cluster, imag=F, cond.labels=F, o.genes=NULL, col.func=if(imag) topo.colors else rainbow,
                        rats.names=names( ratios ), main=NULL, range.r=NULL, no.par=F, sort=F, box.plot=F, ... ) { ##, title=NULL, for.figure=F ) {##regulators=NULL,
  if ( length( cluster$rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  k <- cluster$k
  if ( is.null( main ) ) main <- paste( sprintf( "Cluster: %04d %s; resid: %s; r/c: %d/%d", k, organism,
                                                paste( sprintf( "%.2f", cluster$resid[ rats.names ] ), collapse=" " ),
                                                length( cluster$rows ), length( cluster$cols ) ) )
  rats <- get.cluster.matrix( unique( c( cluster$rows, o.genes ) ), cluster$cols, matrices=rats.names )
  cols.b <- colnames( rats )[ colnames( rats ) %in% cluster$cols ]
  
  if ( sort ) {
    o1 <- order( apply( rats[ cluster$rows, cols.b, drop=F ], 2, mean, na.rm=T ) )
    cols.b <- cols.b[ o1 ]
    rats <- rats[ ,cols.b, drop=F ]
  }
  
  if ( all( is.na( rats ) ) ) { plot( 0, 0, typ="n", min=main, ... ); return() }
  if ( is.vector( rats ) ) { rats <- t( rats ); rownames( rats ) <- cluster$rows }

  if ( imag ) {
    grey.image <- function( mat, n.gray=32, x=1:nrow(mat), y=1:ncol( mat ), col=gray((0:n.gray)/n.gray), ... )
      image( x, y, mat, col=col, ... )
    grey.image( t( rats ), col=col.func( 256 ) )
    return()
  }
  
  ##if ( renorm ) rats <- standardize.row( rats, NULL, 0 )
  ##if ( "ylim" %in% names( list( ... ) ) ) range.r <- ylim
  ##!else
  if ( is.null( range.r ) ) ##range.r <- range( rats, na.rm=T )
    range.r <- range( rats[ rats != min( rats, na.rm=T ) & rats != max( rats, na.rm=T ) ], na.rm=T )
  if ( cond.labels && cluster$ncols < 100 ) range.r[ 1 ] <- range.r[ 1 ] * 1.5
  
  if ( ! no.par ) ##&& ! cond.labels )
    par( mar=rep(2.0,4), mgp=c(3,1,0)*0.5 )##; on.exit( par( old.pars ) )

  plot( 1:length( cols.b ), ylim=range.r, xlab=NA, ylab=NA, main=main, typ="n", xaxs="i", ... )

  if ( length( rats.names ) > 1 ) { ## Shade different ratios matrices backgrounds
    ind <- 0.5
    rts <- NULL
    for ( i in 1:length( rats.names ) ) {
      col <- sapply( col2rgb( i+1 ) / 255 + 0.9, function( cc ) min( cc, 1 ) )
      col <- rgb( col[ 1 ], col[ 2 ], col[ 3 ] )
      rect( ind, range.r[1]+0.05, ind+sum( colnames( rats ) %in% colnames( ratios[[ rats.names[ i ] ]] ) ),
           range.r[2]-0.05, col=col, dens=NA )
      ind <- ind + sum( colnames( rats ) %in% colnames( ratios[[ rats.names[ i ] ]] ) )
      rts <- cbind( rts, rats[ ,cols.b[ cols.b %in% colnames( ratios[[ rats.names[ i ] ]] ), drop=F ] ] )
    }
    rats <- rts; rm( rts )
    cols.b <- colnames( rats )
  }

  if ( exists( "col.rug" ) ) {
    if ( is.integer( col.rug ) ) colmap <- col.func( max( col.rug ) )[ col.rug[ cols.b ] ]
    else colmap <- col.rug[ cols.b ]
  } else if ( all( deparse( col.func ) == deparse( rainbow ) ) ) {
    colmap <- col.func( length( cols.b ) ) ##rep( 'black', length( cols.b ) )
  } else {
    colmap <- col.func( cols.b )
  }
  ##colmap <- rug1[ cols.b, 2 ]
  ##!else 

  if ( box.plot ) {
    colMeans <- apply( rats[ cluster$rows, ,drop=F ], 2, mean, na.rm=T )
    colSd <- apply( rats[cluster$rows, ,drop=F ], 2, sd, na.rm=T )
    matlines( 1:length( cols.b ), cbind( colMeans - 2 * colSd, colMeans + 2 * colSd ), lty=1, col='lightgrey' )
    boxplot( as.data.frame( rats[ cluster$rows, ,drop=F ] ), ylim=range.r, 
            names=NA, main=main, col=colmap, outline=FALSE, border=FALSE, add=T, xaxs="i", xaxt="n", ... )
    if ( sort ) lines( 1:length( cols.b ), colMeans, lty=1, lwd=1, col='red' )
  } else {
    cmap <- col.func( cluster$nrows )
    matlines( 1:length( cols.b ), t( rats[ cluster$rows, ,drop=F ] ), ylim=range.r, xlab=NA, ylab=NA, main=main,
             col=cmap, lty=1, ... ) ##typ="l",
    if ( exists( "col.rug" ) ) for ( i in unique( col.rug ) )
      rug( which( cols.b %in% names( which( col.rug == i ) ) ), col=colmap[ which( col.rug == i )[ 1 ] ] )
  }
  
  if ( cond.labels ) { 
    tmp.y <- rep( range.r[1] * 0.85, cluster$ncols ) 
    cols <- if ( box.plot ) colmap else 'black'
    text( 1:cluster$ncols, tmp.y, cols.b, srt=90, col=cols, ... ) 
  }

  if ( names( dev.cur() ) == "devSVG" ) {
    par( family="Arial" )
    for ( c in 1:length( cols.b ) ) {
      setSVGShapeToolTip( cols.b[ c ] )
      rect( c, range.r[1], c+1, range.r[2], col=NA, border=NA ) 
    }
  }

  if ( ! is.null( o.genes ) ) {
    matlines( 1:length( cols.b ), t( rats[ o.genes, , drop=F ] ), lty=1, lwd=3, col=2:6 )
    legend( "bottomright", legend=o.genes, lty=1, lwd=3, col=2:6, cex=0.7, bty="n" )
  }
}

###############################
plotCluster.all.conds <- function( cluster, imag=F, cond.labels=F, o.genes=NULL, rats.names=names( ratios ), 
                                  range.r=NULL, sort=F, box.plot=F, col.func=if( imag ) topo.colors else rainbow,
                                  only.in.conds=F, ... ) {
  if ( length( cluster$rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  k <- cluster$k
  main <- paste( sprintf( "Cluster: %04d %s; resid: %s; r/c: %d/%d", k, organism,
                         ##paste( sprintf( "%.2f", cluster$resid[ rats.names ] ), collapse=" " ),
                         sprintf( "%.2f", weighted.mean( cluster$resid[ rats.names ], row.weights[ rats.names ], na.rm=T ) ),
                         length( cluster$rows ), length( cluster$cols ) ) )

  rats <- get.cluster.matrix( unique( c( cluster$rows, o.genes ) ), NULL, matrices=rats.names ) ##sort( cluster$cols )
  cols.b <- c( colnames( rats )[ colnames( rats ) %in% cluster$cols ],
              colnames( rats )[ ! colnames( rats ) %in% cluster$cols ] )
  if ( only.in.conds ) cols.b <- colnames( rats )[ colnames( rats ) %in% cluster$cols ]
  
  ##cols.b <- attr( ratios, "cnames" )[ attr( ratios, "cnames" ) %in% cluster$cols ]
  ##cols.b <- c( cols.b, attr( ratios, "cnames" )[ ! attr( ratios, "cnames" ) %in% cluster$cols ] )
  ##rats <- get.cluster.matrix( cluster$rows, cols.b, matrices=rats.names )

  if ( sort ) {
    inClust <- colnames( rats )[ colnames( rats ) %in% cluster$cols ]
    o1 <- order( apply( rats[ cluster$rows, inClust ,drop=F ], 2, mean, na.rm=T ) )
    outClust <- colnames( rats )[ ! colnames( rats ) %in% cluster$cols ]
    o2 <- order( apply( rats[ cluster$rows, outClust ,drop=F ], 2, mean, na.rm=T ) )
    cols.b <- c( inClust[ o1 ], outClust[ o2 ] )
    if ( only.in.conds ) cols.b <- inClust[ o1 ]
  }
  
  len.b <- length( cols.b )
  rats <- rats[ ,cols.b ,drop=F ]

  par( mar=rep(2.0,4), mgp=c(3,1,0)*0.5 )

  if ( all( is.na( rats ) ) ) { plot( 0, 0, typ="n", main=main, ... ); return() }
  if ( is.vector( rats ) ) { rats <- t( rats ); rownames( rats ) <- cluster$rows }

  if ( imag ) {
    grey.image( t( rats ), col=col.func( 256 ) )
    lines( rep( cluster$ncols + 0.5, 2 ), c( -999, 9999 ), col=2, lwd=3, lty=2 )
    return()
  }

  if ( is.null( range.r ) )
    range.r <- range( rats[ rats != min( rats, na.rm=T ) & rats != max( rats, na.rm=T ) ], na.rm=T )
  if ( cond.labels && len.b < 100 ) range.r[ 1 ] <- range.r[ 1 ] * 1.5

  plot( 1:len.b, xlim=c( 0.95, len.b+0.05 ), ylim=range.r, xlab=NA, ylab=NA, main=main, typ="n", xaxs="i", ... )
  
  if ( length( ratios ) > 1 ) { ## Shade different ratios matrices backgrounds
    ind <- 0.5
    rts.in <- rts.out <- NULL
    for ( in.out in 1:2 ) {
      if ( in.out == 1 ) cols <- cols.b[ cols.b %in% cluster$cols ]
      else if ( in.out == 2 ) cols <- cols.b[ ! cols.b %in% cluster$cols ]
      for ( i in 1:length( ratios ) ) {
        #col <- sapply( col2rgb( i+1 ) / 255 + 0.9, function( cc ) min( cc, 1 ) )
        #col <- rgb( col[ 1 ], col[ 2 ], col[ 3 ] )
        col <- col.func( length( ratios ), s=0.2 )[ i ]
        rect( ind, range.r[1]+0.01, ind+sum( cols %in% colnames( ratios[[ i ]] ) ), range.r[2]-0.01, col=col,
             dens=NA )
        ind <- ind + sum( cols %in% colnames( ratios[[ i ]] ) )
        if ( in.out == 1 ) rts.in <- cbind( rts.in, rats[ ,cols[ cols %in% colnames( ratios[[ i ]] ) ] ,drop=F ] )
        else if ( in.out == 2 ) rts.out <- cbind( rts.out, rats[ ,cols[ cols %in% colnames( ratios[[ i ]] ) ] ,drop=F ] )
      }
    }
    rats <- cbind( rts.in, rts.out )
    cols.b <- colnames( rats )
    len.b <- length( cols.b )
    rm( rts.in, rts.out )
  }

  if ( exists( "col.rug" ) ) {
    if ( is.integer( col.rug ) ) colmap <- col.func( max( col.rug ) )[ col.rug[ cols.b ] ]
    else colmap <- col.rug[ cols.b ]
  } else if ( length( rats.names ) > 1 ) {  ## use a different color for each ratio
    colmap <- sapply( cols.b, function( col ) which( sapply( ratios[rats.names], function(i) col %in% colnames(i) ) )[ 1 ] )
    if ( is.list( colmap ) ) { colmap <- unlist( colmap ); names( colmap ) <- cols.b }
    colmap <- col.func( max( colmap ) )[ colmap ]
  } else if ( all( deparse( col.func ) == deparse( rainbow ) ) ) {
    colmap <- col.func( length( cols.b ) ) ##rep( 'black', length( cols.b ) )
  } else {
    colmap <- col.func( cols.b )
  }
  ##colmap <- rug1[ cols.b, 2 ]
  ##!else

  if ( box.plot ) {
    ##if ( exists( "rug1" ) ) colmap <- rug1[ cols.b, 2 ]
    ##!else colmap <- rep( 'black', length( cols.b ) )
##     if ( exists( "col.rug" ) ) colmap <- col.func( max( col.rug ) )[ col.rug[ cols.b ] ]
##     else if ( all( deparse( col.func ) == deparse( rainbow ) ) ) colmap <- col.func( len.b ) ##rep( 'black', length( cols.b ) )
##     else colmap <- col.func( cols.b )
    colMeans <- apply( rats[ cluster$rows, ,drop=F ], 2, mean, na.rm=T )
    colSd <- apply( rats[cluster$rows, ,drop=F ], 2, sd, na.rm=T )
    matlines( 1:length( cols.b ), cbind( colMeans - 2 * colSd, colMeans + 2 * colSd ), lty=1, col='lightgrey' )
    boxplot( as.data.frame( rats[ cluster$rows, ,drop=F ] ), ylim=range.r, 
            names=NA, main=main, col=colmap, outline=FALSE, border=FALSE, add=T, xaxs="i", xaxt="n", ... )
    if ( sort ) lines( 1:length( cols.b ), colMeans, lty=1, lwd=1, col='red' ) ## This hides the boxes if not sorted
  } else {
    cmap <- col.func( cluster$nrows )
    matlines( 1:len.b, t( rats[ cluster$rows, ,drop=F ] ), ylim=range.r, col=cmap, main=main, ##type="l", 
             xlab=NA, ylab=NA, lty=1, ... )
    if ( exists( "col.rug" ) ) for ( i in unique( col.rug ) )
      rug( which( cols.b %in% names( which( col.rug == i ) ) ), col=colmap[ which( col.rug == i )[ 1 ] ] )
  }
  
  ##lines( rep( cluster$ncols + 0.5, 2 ), range.r, col=2, lwd=3, lty=2 )
  cols.in <- colnames( rats )[ colnames( rats ) %in% cluster$cols ]
  if ( ! only.in.conds ) lines( rep( length( cols.in ) + 0.5, 2 ), range.r, col=2, lwd=3, lty=2 )

  if ( ! is.null( o.genes ) ) {
    matlines( 1:len.b, t( rats[ o.genes, , drop=F ] ), lty=1, lwd=3, col=2:6 )
    legend( "bottomright", legend=o.genes, lty=1, lwd=3, col=2:6, cex=0.7, bty="n" )
  }
  
  if ( cond.labels ) { ## || len.b < 100 ) {
    tmp.y <- rep( range.r[1] * 0.85, len.b ) ##(1:len.b / 1:len.b ) * range.r[1] * 0.85
    cols <- if ( box.plot ) colmap else 'black'
    text( 1:len.b, tmp.y, cols.b, srt=90, col=cols, ... ) ##cex=0.7 )
  }

  if ( names( dev.cur() ) == "devSVG" ) {
    par( family="Arial" )
    for ( c in 1:length( cols.b ) ) {
      setSVGShapeToolTip( cols.b[ c ] )
      rect( c, range.r[1], c+1, range.r[2], col=NA, border=NA ) ##, xpd=NA )
    }
  }
}

plotCluster.motif <- function( cluster, seqs=cluster$seqs, layout=NULL, colors=NULL, motif.e.cutoff=Inf,
                              no.plotCluster=T, addl.text=NULL, ... ) { 
  if ( length( cluster$rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  if ( names( dev.cur() ) == "devSVG" ) par( family="Arial" )
  
  if ( any( ! cluster$rows %in% attr( ratios, "rnames" ) ) ) {
    cluster$rows <- cluster$rows[ cluster$rows %in% attr( ratios, "rnames" ) ]
    cluster$nrows <- length( cluster$rows )
    warning( cluster$k, ": Some cluster rows are not in the ratios. Will plot without these rows.\n" )
  }
  if ( any( ! cluster$cols %in% attr( ratios, "cnames" ) ) ) {
    cluster$cols <- cluster$cols[ cluster$cols %in% attr( ratios, "cnames" ) ]
    cluster$ncols <- length( cluster$cols )
    warning( cluster$k, ": Some cluster cols are not in the ratios. Will plot without these cols.\n" )
  }

  seq.types <- cluster$seq.type
  if ( length( seq.types ) == 1 ) {
    n.pssm.plot <- 3
  } else {
    n.pssm.plot <- 6
  }

  if ( is.null( layout ) ) {
    if ( length( seq.types ) == 1 ) {
      layout <- matrix( c( 1, 1, 1, 1, 1, 1, 1, 1, 8, 2, 2, 2, 2, 2, 2, 2, 2, 
                          1, 1, 1, 1, 1, 1, 1, 1, 8, 2, 2, 2, 2, 2, 2, 2, 2, 
                          3, 3, 3, 3, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 
                          5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7 ), ncol=17, byrow=T )
    } else { ## Allow for up to 3 seq.type motifs (2 each type)
      layout <- matrix( c( 1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2,
                          1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2,
                          1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2,
                          1, 1, 1, 1, 1, 1, 1, 1, 11, 2, 2, 2, 2, 2, 2, 2, 2,
                          3, 3, 3, 3, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                          5, 5, 5, 5, 6, 6, 6, 6, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                          7, 7, 7, 7, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                          8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10 ), ncol = 17, byrow = T )
    }
    if ( no.plotCluster ) { ## Dont include panel with plot of cluster with only in-cols (plotCluster())
      layout[ layout == 1 | layout == max( layout ) ] <- 2
      layout <- layout - 1
      layout[ ,1 ][ layout[ ,1 ] == 1 ] <- max( layout ) + 1
    }
  }
  layout( layout )
    
  k <- cluster$k
  
  if ( ! is.null( ratios ) ) {
    args <- list( ... )
    args <- args[ names( args ) != "p.val.shade.cutoff" ] ## To prevent the warning: 
    args$cluster <- cluster ## 'In plot.window(...) : "p.val.shade.cutoff" is not a graphical parameter'
    do.call( plotCluster.all.conds, args ) ##o.genes=o.genes, ... )
    if ( ! no.plotCluster ) do.call( plotCluster, args ) ##o.genes=o.genes, ... )
  }

  rows <- cluster$rows
  if ( is.null( colors ) ##&& exists( "gene.coords" )
      || ! is.null( cluster$cog.code ) ) {
    tmp.lett <- 1:26
    names( tmp.lett ) <- LETTERS
    ##if ( exists( "gene.coords" ) && ! is.null( gene.coords$gene.code ) ) coo <- gene.coords$gene.code[ rows ]
    ##!else
    if ( ! is.null( cluster$cog.code ) ) coo <- cluster$cog.code[ rows ]
    else coo <- 1:length( rows )
    tmp <- unique( tmp.lett[ coo ] )
    names( tmp ) <- names( tmp.lett[ coo ][ ! duplicated( tmp.lett[ coo ] ) ] )
    cols <- rainbow( length( tmp ) )
    names( cols ) <- names( tmp )
    cols <- cols[ names( tmp.lett[ coo ] ) ]
    cols[ is.na( names( cols ) ) ] <- "darkgrey"
    names( cols ) <- rows
  } else {
    cols <- rainbow( length( rows ) )
    names( cols ) <- rows
  }
  colors <- cols

  ##if ( length( seq.types ) > 1 ) seq.types <- seq.types[ 1:3 ]
  n.plotted <- 1
  for ( seq.type in seq.types ) {
    if ( n.plotted > n.pssm.plot ) break
    ##if ( is.na( seq.type ) ||
    if ( is.null( cluster[[ seq.type ]]$e.val ) ||
        all( is.na( cluster[[ seq.type ]]$e.val ) ) || is.null( cluster[[ seq.type ]]$motif.out ) || ##min( cluster$e.val, na.rm=T ) > min.e.val ||
        ##all( is.na( cluster[[ seq.type ]]$p.clust ) ) ||
        is.null( cluster[[ seq.type ]]$motif.out$pssms ) ) next
##       for ( i in 1:2 ) plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" )
##     } else {
    pssm <- cluster[[ seq.type ]]$motif.out$pssms
    for ( ppp in 1:min( floor( n.pssm.plot / length( seq.types ) ), length( pssm ) ) ) {
      if ( n.plotted > n.pssm.plot ) break
      ##opar <- par( mar=rep(0.5,4), mgp=c(3,1,0)*0.5 )##; on.exit( par( old.pars ) )
      if ( cluster[[ seq.type ]]$motif.out$e.values[ ppp ] > motif.e.cutoff ) next
      viewPssm( pssm[[ ppp ]], mot.ind=ppp, main.title=sprintf( "%s PSSM #%d; e=%.3g", seq.type, ppp, 
                                              cluster[[ seq.type ]]$motif.out$e.values[ ppp ] ), cex.main=0.9 )
      n.plotted <- n.plotted + 1
      ##par( opar )
    }
    ##  if ( length( pssm ) < 2 ) for ( i in ( length( pssm ) + 1 ):2 )
    ##    plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" )
    ##}
    ##if ( length( seq.types ) == 1 ) plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" )
  }
  while ( n.plotted <= n.pssm.plot ) {
    plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" )
    n.plotted <- n.plotted + 1
  }

  suppressWarnings( cluster <- plotCluster.network( cluster, ... ) ) ## suppress printing of confusing internal igraph warnings

  ##seqs <- cluster$seqs
  if ( is.null( seqs ) ) { seqs <- rep( "", length( cluster$rows ) ); names( seqs ) <- cluster$rows }
  if ( ! is.null( seq.type ) && ! is.null( seqs ) && length( seqs ) > 0 )
    plotClusterMotifPositions( cluster, seqs, colors=colors, ... ) ##o.genes=o.genes, 
  ##p.val.shade.cutoff=p.val.shade.cutoff, ... )
  else plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" )

  try ( {
    par( mar=rep(0.5,4), mgp=c(2,1,0)*0.5 )##; on.exit( par( old.pars ) )
    plot( c( 0.5, 2.5 ), c( -1, 1 ), type="n", tck=0.01, cex.lab=0.2, cex.sub=0.2, cex.axis=0.2, axes=F )
    if ( names( dev.cur() ) == "devSVG" ) {
      par( family="Arial" )
      setSVGShapeToolTip( title=paste( "Cluster:", sprintf( "%04d", k ), organism, cmonkey.version ),
                         desc1=sprintf( "resid = %s; genes = %d; conds = %d",
                           paste( sprintf( "%.2f", cluster$resid ), collapse=" " ),
                           length( cluster$rows ), length( cluster$cols ) ) ) 
      setSVGShapeURL( paste( "http://www.genome.ad.jp/dbget-bin/www_bget?", 
                            paste( organism, ":", cluster$rows, sep="", collapse="+" ), sep="" ) )
      rect( 0.5, -1, 3.25, +1, col="lightgreen", border=NA )
    }
  } ) 

  ##if ( names( dev.cur() ) != "devSVG" ) { ## TODO: this dies on the SVG setting - fix it!?
  if ( ! is.null( cluster$name ) ) text( 0.65, 0, cluster$name, srt=90, xpd=NA, cex=1 )
  else if ( ! is.null( addl.text ) ) text( 0.65, 0, addl.text, srt=90, xpd=NA, cex=1 )
  text( 1.5, 0, sprintf( "%s iter=%d", date.run, iter ), srt=90, xpd=NA, cex=1 )
  text( 2.35, 0, paste( "cMonkey Version", cmonkey.version, organism ), srt=90, xpd=NA, cex=1 )
  ##}

  invisible( cluster )
  ##if ( ! is.null( layout ) && exists( "plotCluster.custom" ) ) plotCluster.custom( cluster, seqs=cluster$seqs, ... )
}

plotCluster.network <- function( cluster, network="all", o.genes=NULL, colors=NULL, cex=0.7, no.legend=F, ... ) {
  if ( length( cluster$rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  require( igraph0 )
  rows <- cluster$rows
  if ( is.null( cluster$network ) ) {
    if ( network == "all" ) network <- names( networks )
    for ( i in network ) { 
      if ( ! i %in% names( networks ) ) next
      tmp.net <- networks[[ i ]][ networks[[ i ]]$protein1 %in% rows & networks[[ i ]]$protein2 %in% rows, ]
      tmp.net <- cbind( tmp.net, net=rep( i, nrow( tmp.net ) ) )
      cluster$network <- if ( ! is.null( cluster$network ) ) rbind( cluster$network, tmp.net ) else tmp.net
    }
  }
  network <- cluster$network
  nrows <- rows; if ( ! is.null( o.genes ) ) nrows <- unique( c( nrows, o.genes ) )
  if ( is.null( cluster$cog.code ) && "cog.code" %in% names( genome.info ) )
    cluster$cog.code <- genome.info$cog.code[ rows ]
  
  if ( is.null( cluster$colors ) ) {
    if ( is.null( colors ) ##&& exists( "gene.coords" )
        || "cog.code" %in% names( genome.info ) ) {
      tmp.lett <- 1:26
      names( tmp.lett ) <- LETTERS
      ##if ( exists( "gene.coords" ) && ! is.vector( gene.coords ) && ! is.null( gene.coords$gene.code ) )
      ##  coo <- gene.coords$gene.code[ rows ]
      ##!else
      if ( ! is.null( cluster$cog.code ) ) coo <- cluster$cog.code[ rows ]
      else coo <- 1:length( rows )
      tmp <- unique( tmp.lett[ coo ] )
      names( tmp ) <- names( tmp.lett[ coo ][ ! duplicated( tmp.lett[ coo ] ) ] )
      cols <- rainbow( length( tmp ) )
      names( cols ) <- names( tmp )
      cols <- cols[ names( tmp.lett[ coo ] ) ]
      cols[ is.na( names( cols ) ) ] <- "darkgrey"
      names( cols ) <- rows
      cluster$colors <- cols
    } else {
      cols <- rainbow( length( rows ) )
      names( cols ) <- rows
      cluster$colors <- cols
    }
  }
  colors <- cluster$colors
  
  if ( is.null( network ) || nrow( network ) <= 0 )
    network <- data.frame( protein1=nrows, protein2=nrows, combined_score=jitter( rep( 1/50, length( nrows ) ) ),
                          net=rep( "none", length( nrows ) ) )
  not.in <- nrows[ ! nrows %in% network$protein1 & ! nrows %in% network$protein2 ]
  for ( i in not.in ) network <- rbind( network, data.frame( protein1=i, protein2=i, combined_score=0, net="none" ) )
  gr <- graph.edgelist( as.matrix( network[ ,1:2 ] ), directed=F )
  net.wts <- as.numeric( network$combined_score ); names( net.wts ) <- as.character( network$net )
  for ( n in unique( names( net.wts ) ) ) {
    if ( n == "none" ) next
    net.wts[ names( net.wts ) == n ] <- net.wts[ names( net.wts ) == n ] /
      max( net.wts[ names( net.wts ) == n ], na.rm=T )
  }
  gr.layout <- layout.fruchterman.reingold( gr, niter=3000, weights=net.wts / 5 ) ##network$combined_score / 5000 )
  gr.layout <- layout.norm( gr.layout, -1, 1, -1, 1 )

  edge.colors <- character(); curves <- rep( 0, nrow( network ) )
  nets <- unique( as.character( network$net ) )
  if ( "none" %in% nets ) {
    nets <- unique( c( "none", nets ) )
    inds <- c( 1, 1:( length( nets ) - 1 ) ); inds <- inds[ inds != 0 ]; inds <- inds[ 1:length( nets ) ]
    net.colors <- t( col2rgb( inds, T ) ) / 255
    net.colors[ 1, 4 ] <- 0
  } else {
    net.colors <- t( col2rgb( 1:length( nets ), T ) ) / 255
  }
  rownames( net.colors ) <- nets
  for ( i in 1:nrow( network ) ) {
    net <- as.character( network$net )[ i ]
    shade <- net.wts[ i ] ##network$combined_score[ i ] / 1000
    nodes <- c( as.character( network$protein1 )[ i ], as.character( network$protein2 )[ i ] )
    sub.net <- subset( network,protein1 %in% nodes & protein2 %in% nodes )
    sub.nets <- unique( as.character( sub.net$net ) )
    curve.it <- max( 0, nrow( sub.net ) - 2 ) / 2 * 0.33
    col <- net.colors[ net, ]
    col2 <- col; col2[ col2 == 0 ] <- 1 - shade
    edge.colors[ i ] <- if ( names( dev.cur() ) != "X11" ) rgb( col[ 1 ], col[ 2 ], col[ 3 ], shade ) else
                   rgb( col2[ 1 ], col2[ 2 ], col2[ 3 ] )
    curves[ i ] <- curve.it * floor( which( sub.nets == net ) / 2 ) *
      ( if ( which( sub.nets == net ) %% 2 == 0 ) -1 else 1 )
  }
  if ( all( curves == 0 ) ) curves <- FALSE ## same but faster to draw
  if ( ! no.legend ) {
    labels <- try( get.long.names( get.vertex.attribute( gr, "name" ), short=T ), silent=T )
    if ( class( labels ) == "try-error" ) labels <- get.vertex.attribute( gr, "name" )
    labels[ is.na( labels ) | labels == "" ] <- get.vertex.attribute( gr, "name" )[ is.na( labels ) | labels == "" ]
  } else labels <- NA ## Assume it's for the thumbnail (so no vertex labels needed either)
  plot( gr, layout=gr.layout, margin=0, rescale=F, edge.curved=curves,
       vertex.color=colors[ get.vertex.attribute( gr, "name" ) ],
       vertex.frame.color=colors[ get.vertex.attribute( gr, "name" ) ], vertex.label.cex=cex, vertex.size=7,
       vertex.label=labels,
       vertex.label.family=if ( names( dev.cur() ) != "pdf" ) "Arial" else "sans",
       edge.color=edge.colors, edge.width=round( net.wts ) + 1 ) 
  if ( ! no.legend && length( nets[ nets != "none" ] ) > 0 )
    legend( "bottomright", legend=nets[ nets != "none" ], col=1:length( nets[ nets != "none" ] ),
                                                            lty=1, lwd=2, bty="n", cex=0.5 )

  if ( names( dev.cur() ) == "devSVG" ) {
    names <- cluster$gene.coords ##get.long.names( cluster$k )
    for ( i in 1:nrow( gr.layout ) ) {
      gene <- get.vertex.attribute( gr, "name" )[ i ]
      setSVGShapeToolTip( title=gene, desc1=ifelse( is.na( names[ gene ] ) || is.null( names[ gene ] ), "",
                                        names[ gene ] ) )
      setSVGShapeURL( paste( "http://www.genome.ad.jp/dbget-bin/www_bget?", organism, ":", gene, sep="" ) )
      points( gr.layout[ i, 1 ], gr.layout[ i, 2 ], col="#FF000001", cex=10/3 )
    }
  }
  cluster
}  

col.let <- c( "A", "C", "G", "T" )
                   
viewPssm <- function( pssm, e.val=NA, mot.ind=NA, use.char=T, main.title=NA, no.par=F, scale.e=NA, boxes=F, new=T,
                     xoff=0, yoff=0, no.axis.labels=F, min.height.drawn=1e-5, ... ) { 
  if ( is.null( pssm ) ) return()
  getEntropy <- function( pssm ) {
    pssm[ pssm == 0 ] <- 0.00001
    entropy <- apply( pssm, 1, function( i ) -sum( i * log2( i ) ) )
    return( entropy )
  }

  char.coords = list( T=list( x=c( 0.45, 0.55, 0.55, 1, 1, 0, 0, 0.45 ), y=c( 0, 0, 0.9, 0.9, 1, 1, 0.9, 0.9 ), color=2 ),
    A=list( x=c( 0, 0.1, 0.28, 0.72, 0.68, 0.32, 0.5, 0.9, 1, 0.55, 0.45, 0 ), y=c( 0, 0, 0.4, 0.4, 0.5, 0.5, 0.9, 0, 0, 1, 1, 0 ), color=3 ),
    C=list( x=c( 1, 1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 0.85, 1.0, 1.0, 0.9, 0.9, 0.8, 0.55, 0.45, 0.2, 0.1, 0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9 ), y=c( 0.6, 0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 0.35, 0.4, 0.4, 0.35, 0.2, 0.1, 0.1, 0.2, 0.42, 0.58, 0.8, 0.9, 0.9, 0.8, 0.65, 0.6 ), color=4 ),
    G=list( x=c( 1, 1, 0.85, 0.55, 0.45, 0.15, 0, 0, 0.15, 0.45, 0.55, 0.85, 1.0, 1.0, 0.7, 0.7, 0.9, 0.8, 0.55, 0.45, 0.2, 0.1, 0.1, 0.2, 0.45, 0.55, 0.8, 0.9, 0.9 ), y=c( 0.6, 0.7, 0.9, 1, 1, 0.9, 0.65, 0.35, 0.1, 0, 0, 0.1, 0.35, 0.5, 0.5, 0.4, 0.4, 0.2, 0.1, 0.1, 0.2, 0.42, 0.58, 0.8, 0.9, 0.9, 0.8, 0.65, 0.6 ), color="orange" ) )##,
  
  draw.char <- function( char=col.let, rect=c( 0, 0, 1, 1 ), border=NULL, ... ) {
    ## rect is (x,y,width,height)
    if ( rect[ 4 ] <= min.height.drawn ) return()
    x <- char.coords[[ char ]]$x * rect[ 3 ] + rect[ 1 ]
    y <- char.coords[[ char ]]$y * rect[ 4 ] + rect[ 2 ]
    color <- char.coords[[ char ]]$color
    if ( is.null( border ) ) border <- color
    if ( ! boxes ) {
      polygon( x + xoff, y + yoff, col=color, border=border, density=NA, ... )
    } else {
      rect[ 1 ] <- rect[ 1 ] + xoff
      rect[ 2 ] <- rect[ 2 ] + yoff
      rect( rect[ 1 ], rect[ 2 ], rect[ 1 ] + rect[ 3 ], rect[ 2 ] + rect[ 4 ], col=color, border=border,
           density=NA, ... )
    }
  }

  win.size <- nrow( pssm )
  if ( ! no.par ) par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )##; on.exit( par( old.pars ) )
  if ( any( pssm <= 0 ) ) pssm <- pssm + 1e-10
  if ( any( pssm > 1 ) ) pssm <- t( apply( pssm, 1, function( i ) i / ( sum( i ) + 1e-10 ) ) )

  entr <- getEntropy( pssm )
  if ( is.na( scale.e ) ) scale.e <- (2 - entr) / 2
  ##scale.e[ scale.e < 0.05] <- 0.05
  x.range <- c(0.5, win.size + 0.5 )
  y.range <- c(0,max(scale.e))
  if ( new ) plot( x.range, y.range, type="n", tck=0.01, cex.lab=0.2, cex.sub=0.2, cex.axis=0.2, axes=F )
  if ( ! is.na( main.title[ 1 ] ) ) {
    if ( ! is.na( mot.ind ) ) title( main.title, col.main=mot.ind+1, xpd=NA, ... )
    else title( main.title, xpd=NA, ... )
  } else if ( ! is.na( mot.ind ) || ! is.na( e.val ) ) {
    if ( ! is.na( e.val ) ) tmp.tit <- sprintf( "PSSM #%d; E=%.3g", mot.ind, e.val )
    else tmp.tit <- tmp.tit <- sprintf( "PSSM #%d", mot.ind )
    title( tmp.tit, col.main=mot.ind+1, xpd=NA, ... )
    ##} else title( tmp.tit, xpd=NA, ... )
  }
  pssm.sc <- scale.e * pssm
  for (j in 1:win.size) {
    inds <- sort( pssm.sc[ j, ], index=T )$ix
    for (i in 1:4) {
      ind <- inds[ i ]
      if ( i == 1 ) {
        if ( ! use.char ) {
          rect( (j-0.5), 0, (j+0.5), pssm.sc[j,ind], col=colMap[ind])
          if (pssm[j,ind] > 0.05) text(j, 0 + pssm.sc[j,ind]/2, colLet[ind])
        } else {
          draw.char( col.let[ ind ], c( (j-0.4), 0, 0.9, pssm.sc[j,ind]-0.001 ), ... )
        }
        prev.h <-  pssm.sc[j,ind]
      } else {
        if ( ! use.char ) {
          rect( (j-0.5), prev.h, (j+0.5), (pssm.sc[j,ind]+prev.h), col=colMap[ind])
          if (pssm.sc[j,ind] > 0.05 ) {
            if (i == 2) text(j,  prev.h + 0.5 * pssm.sc[j,ind], colLet[ind], col=8)
            else text(j,  prev.h + 0.5 * pssm.sc[j,ind], colLet[ind])
          }
        } else {
          draw.char( col.let[ ind ], c( (j-0.4), prev.h, 0.9, pssm.sc[j,ind]-0.001 ), ... )
        }
        prev.h <- prev.h + pssm.sc[j,ind]
      }
    }
  }
  if ( ! no.axis.labels ) {
    if ( win.size < 10 )
      text( 1:win.size, rep( -0.01, win.size ), as.character( 1:win.size ), cex=0.7, adj=c(0.5,1), xpd=NA )
    else if ( win.size < 20 ) text( seq( 1, win.size, 2 ), rep( -0.01, win.size ),
                                   as.character( seq( 1, win.size, 2 ) ), cex=0.7, adj=c(0.5,1), xpd=NA )
    else if ( win.size < 50 ) text( seq( 1, win.size, 5 ), rep( -0.01, win.size ),
                                   as.character( seq( 1, win.size, 5 ) ), cex=0.7, adj=c(0.5,1), xpd=NA )
    else text( seq( 1, win.size, 25 ), rep( -0.01, win.size ),
              as.character( seq( 1, win.size, 25 ) ), cex=0.7, adj=c(0.5,1), xpd=NA )
  }
  invisible( y.range )
}

## Note: to get list of all positions of motif #1 in SEARCHED sequences, this now works:
## posns=unlist(sapply(e$meme.scores$upstream[1:e$k.clust],function(i)i$meme.out[[1]]$posns$start))
## And this will give the posns for motif #1 in all SCANNED sequences (in-clust):
## posns=unlist(sapply(e$meme.scores$upstream[1:e$k.clust],function(i)i$pv.ev[[1]]$posns))
## So the 'out.posns' returned here is now commented out! And I guess we don't need 'no.plot' either.
## This is all moved to new func in cmonkey-postproc.R -- "plot.all.clusterMotifPositions"
###############################
plotClusterMotifPositions <- function( cluster, seqs=cluster$seqs, long.names=T, shade=T, p.val.shade.cutoff=999,
                                      colors=NULL, sort.by="p.value", o.genes=NULL, no.key=F,  ##no.plot=F, 
                                      short.names=organism == "sce", seq.type=cluster$seq.type[ 1 ], ... ) { 
  if ( length( cluster$rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  k <- cluster$k
  rows <- cluster$rows
  if ( ! is.null( o.genes ) ) rows <- unique( c( rows, o.genes ) )

  motif.out <- NULL
  if ( ! is.null( seq.type ) && seq.type %in% names( cluster ) ) motif.out <- cluster[[ seq.type ]]$motif.out
  is.dup.seq <- get.dup.seqs( cluster$seqs )
  p.clust <- cluster$p.clust ##[ seq.type ]
  e.clust <- cluster$e.val ##[ seq.type ]
  motif.info <- NULL
  if ( ( ! all( is.na( p.clust ) ) || ! all( is.na( e.clust ) ) ) && ! is.null( motif.out ) &&
      ! is.null( motif.out$pv.ev ) && length( motif.out$pv.ev ) > 1 )
    motif.info <- subset( motif.out$pv.ev[[ 2 ]], gene %in% rows )

  ##if ( no.plot ) {
  ## out.posns <- list()
  ## out.posns$starts <- out.posns$ends <- out.posns$mots <- integer()
  ## out.posns$p.vals <- out.posns$e.vals <- numeric()
  ##}

  if ( is.null( colors ) ##&& exists( "gene.coords" )
      || "cog.code" %in% names( genome.info ) ) {
    tmp.lett <- 1:26
    names( tmp.lett ) <- LETTERS
    ##if ( exists( "gene.coords" ) && ! is.vector( gene.coords ) && ! is.null( gene.coords$gene.code ) )
    ##  coo <- gene.coords$gene.code[ rows ]
    ##!else
    if ( ! is.null( cluster$cog.code ) ) coo <- cluster$cog.code[ rows ]
    else coo <- 1:length( rows )
    tmp <- unique( tmp.lett[ coo ] )
    names( tmp ) <- names( tmp.lett[ coo ][ ! duplicated( tmp.lett[ coo ] ) ] )
    cols <- rainbow( length( tmp ) )
    names( cols ) <- names( tmp )
    cols <- cols[ names( tmp.lett[ coo ] ) ]
    cols[ is.na( names( cols ) ) ] <- "darkgrey"
    names( cols ) <- rows
  } else {
    cols <- rainbow( length( rows ) )
    names( cols ) <- rows
  }
  cluster$colors <- colors <- cols

  no.motif <- FALSE
  p.values <- motif.widths <- pssm <- NULL ##diagrams <- NULL
  if ( ! is.null( motif.out ) && ! is.null( motif.info ) && nrow( motif.info ) > 0 &&
      ( ! all( is.na( p.clust ) ) || ! all( is.na( e.clust ) ) ) ) { 
    p.values <- motif.out$p.values[ rows ]
    motif.widths <- sapply( motif.out$pssms, nrow, simplify=T )
    pssm <- motif.out$pssms
  } else {
    no.motif <- TRUE
    p.values <- numeric( length( rows ) )
    motif.widths <- 0
  }
  seqs <- seqs[ rows ]
  names( seqs ) <- names( p.values ) <- rows ##names( diagrams ) <-

  seq.lengths <- nchar( seqs ); seq.lengths[ seq.lengths == 2 ] <- NA ## NA seqs become nchar==2 for some reason!
  if ( any( seq.lengths[ ! is.na( seq.lengths ) ] > median( seq.lengths, na.rm=T ) ) ) {
    seqs <- substr( seqs, 1, median( seq.lengths, na.rm=T ) ) 
    seq.lengths <- nchar( seqs )
  }

  maxlen <- max( seq.lengths, na.rm=T )
  ##print(k);print(rows);print(seqs);print(is.dup.seq);print(seq.lengths);print(maxlen)
  if ( ! is.null( seq.type ) && ( maxlen == 0 || is.infinite( maxlen ) ) )
    maxlen <- diff( motif.upstream.search[[ seq.type ]] )
  inds <- integer()
  if ( no.motif && ( sort.by == "p.value" || sort.by == TRUE ) ) sort.by <- "gene.name"
  if ( sort.by == "gene.name" ) inds <- sort( rows, decreasing=T, index=T )$ix
  else if ( sort.by == "p.value" || sort.by == TRUE ) inds <- order( p.values[ rows ], decreasing=T, na.last=F )
  else if ( sort.by == "resid" ) inds <- order( row.scores[ rows, k ], decreasing=T )
  ##!else if ( sort.by == "total" ) inds <- order( rr.scores[ rows, k ], decreasing=T )
  if ( length( inds ) < length( rows ) ) inds <- c( (1:length( rows ))[ ! 1:length( rows ) %in% inds ], inds )
  x.range <- c( -maxlen*0.08, maxlen*1.15 ) ##c( -maxlen*0.15, maxlen*1.08 )
  y.range <- c( 0.5, length( rows ) + 1 )
  ##if ( ! no.plot ) {
  plot( x.range, y.range, type="n", axes=F, xlab="sequence position", ylab="" )
  cexes <- 1.0
  if ( ! no.key ) axis( side=1, pos=0.6, tck=0.01, mgp=c(0.1,0.1,0.1), 
                       labels=c( -1, seq( -100, -maxlen, -100 ) ),
                       at=seq( maxlen, 0, -100 ) + motif.upstream.scan[[ seq.type ]][ 1 ], ... )
  if ( max( seq.lengths, na.rm=T ) > 0 )
    sapply( maxlen - c( 0, motif.upstream.search[[ seq.type ]] ) + motif.upstream.scan[[ seq.type ]][ 1 ],
           function( i ) lines( rep( i, 2 ), c( -999, 999 ), col="lightgray", lty=2 ) )
  colmap <- rainbow( length( rows ) )
  ##}
  mots.used <- numeric()

  if ( is.list( motif.widths ) ) {
    if ( length( motif.widths ) <= 0 ) motif.widths <- 0
    else {
      for ( i in 1:length( motif.widths ) ) if ( is.null( motif.widths[[ i ]] ) ) motif.widths[[ i ]] <- 0
      motif.widths <- unlist( motif.widths )
    }
  }

  lwd <- 3
  if ( length( rows ) > 20 ) lwd <- 1
  else if ( length( rows ) > 10 ) lwd <- 2
  if ( no.key ) lwd <- 1

  if ( ! no.motif ) {
    tmp.mot.info <- subset( motif.info, gene %in% rows )
    tmp.mot.info <- subset( tmp.mot.info, posns <= diff( motif.upstream.scan[[ seq.type ]] ) )
    p.min <- quantile( log10( tmp.mot.info$pvals ), 0.1, na.rm=T )
    if ( is.na( p.min ) ) p.min <- -5
    p.max <- quantile( log10( tmp.mot.info$pvals ), na.rm=T, 0.9 ) ##+ 1
    if ( is.na( p.max ) ) p.max <- log10( p.val.shade.cutoff )
  }
  
  for ( j in 1:length( rows ) ) {
    jj <- inds[ j ]
    cur.gene <- rows[ jj ]
    seq.len <- seq.lengths[ jj ]

    if ( ! is.null( rows ) ) {
      label <- rows[ jj ]
      if ( ! is.null( colors ) ) rect( maxlen+5, j-0.18, maxlen*1.195, j+0.18, col=colors[ label ],
                                      border=colors[ label ], lwd=3 )
    }

    if ( ! no.motif ) {
      rects <- NULL
      mot.info <- subset( tmp.mot.info, gene == cur.gene ) ##motif.out$mast.info[[ cur.gene ]]
      if ( nrow( mot.info ) > 0 ) {
        mots <- mot.info$mots
        starts <- mot.info$posns ##is still f-ed up (bugs in get.mast.pvals())
        widths <- motif.widths[ abs( mots ) ]
        for ( i in 1:length( mots ) ) {
          mot <- mots[ i ]
          if ( is.na( mot ) || is.na( seq.len ) ) next
          start <- starts[ i ]
          if ( start > seq.len ) next
          end <- start + widths[ i ]
          ##if ( no.plot ) {
          ## out.posns$starts <- c( out.posns$starts, start )
          ## out.posns$ends <- c( out.posns$ends, end )
          ## out.posns$mots <- c( out.posns$mots, abs( mot ) )
          ## out.posns$p.vals <- c( out.posns$p.vals, mot.info$pvals[ i ] )
          ## out.posns$e.vals <- c( out.posns$e.vals, motif.out$e.values[ abs( mot ) ] )
          ## if ( no.plot ) next
          ##}
          
          mots.used <- unique( c( mots.used, abs( mot ) ) )
          col <- abs( mot ) + 1

          if ( shade ) {
            if ( ! is.null( mot.info ) ) p.val <- mot.info$pvals[ i ]
            else p.val <- 1e-5
            if ( is.na( p.val ) || p.val > p.val.shade.cutoff || p.val > 1 ) next ##p.val <- 10
            else if ( p.val <= 0 ) p.val <- 1e-5
            ## Map range p.val of -5 to log10(p.val.shade.cutoff) (dark to light)
            ## e.g. for red: -5 -> (1.0,0,0)     log10(p.val.shade.cutoff) -> (1.0,1.0,1.0)
            p.val <- log10( p.val )
            col <- col2rgb( palette()[ col ] ) / 255
            col[ col > 0 ] <- 1
            tmp <- if ( p.val < 10 ) min( 1, max( 0, ( p.val - p.min ) / ( p.max - p.min ) ) ) else 0.99
            if ( names( dev.cur() ) != "X11" ) {
              alpha <- tmp
            } else {
              col[ col == 0 ] <- tmp
              alpha <- 0
            }
            col[ col < 0 ] <- 0
            col[ col > 1 ] <- 1 #0.9
            col <- rgb( col[ "red", 1 ], col[ "green", 1 ], col[ "blue", 1 ], 1-alpha )
          }
          
          start.1 <- start + maxlen - seq.len
          end.1 <- end + maxlen - seq.len
          if ( names( dev.cur() ) == "devSVG" ) {
            par( family="Arial" )
            setSVGShapeToolTip( title=sprintf( "Motif # %2d", abs( mot ) ),
                               desc1=paste( ifelse ( mot < 0, "Rev.", "For." ), "strand,",
                                 start - maxlen, "to", end-maxlen ),
                               desc2=sprintf( "p-value = %.2g", 10^p.val ) )
          }
          if ( ! is.null( mot.info ) ) {
            if ( names( dev.cur() ) == "devSVG" ) {
              if ( mot > 0 ) rect( start.1, j+0.01, end.1, j+0.3, col=col, border=col )
              else if ( mot < 0 ) rect( start.1, j-0.3, end.1, j-0.01, col=col, border=col )
            } else {
              if ( mot > 0 ) rects <- rbind( rects, c( start.1, j+0.01, end.1, j+0.3, col=col, border=col ) )
              else if ( mot < 0 ) rects <- rbind( rects, c( start.1, j-0.3, end.1, j-0.01, col=col, border=col ) )
            }
          } else {
            if ( names( dev.cur() ) == "devSVG" ) {
              if ( mot > 0 ) rect( start.1, j+0.01, end.1, j+0.3, border=col )
              else if ( mot < 0 ) rect( start.1, j-0.3, end.1, j-0.01, border=col )
            } else {
              if ( mot > 0 ) rects <- rbind( rects, c( start.1, j+0.01, end.1, j+0.3, col=NA, border=col ) )
              else if ( mot < 0 ) rects <- rbind( rects, c( start.1, j-0.3, end.1, j-0.01, col=NA, border=col ) )
            }
          }
        }
      }
      if ( ! is.null( rects ) ) rect( rects[ ,1 ], rects[ ,2 ], rects[ ,3 ], rects[ ,4 ], col=rects[ ,5 ],
                                     border=rects[ ,6 ] ) ## Faster over remote x11 than doing each individually
    }
    
    ##if ( no.plot ) next

    ## Make gene start be at the right
    slen <- seq.lengths[ jj ]; if ( all( seq.lengths[ ! is.na( seq.lengths ) ] ) == 0 ) slen <- 50
    lines( c( maxlen - slen, maxlen ), c( j, j ), lwd=lwd + as.integer( rows[ jj ] %in% o.genes ), col=colmap[jj] )

    ## New: if there is a masked region of sequence ("N"s), then grey it out!
    if ( grepl( "N", seqs[ jj ] ) ) {
      locs <- which( strsplit( seqs[ jj ], "" )[[ 1 ]] == "N" )
      diff.locs <- c( diff( locs ), 999 )
      ##locs <- unique( c( locs[ 1 ], locs[ locs %% 5 == 0 ], locs[ length( locs ) ] ) )
      ##sapply( locs, function( jjj ) lines( c( jjj, jjj + 1 ) + maxlen - slen, c( j, j ),
      ##                                  lwd=lwd + as.integer( rows[ jj ] %in% o.genes ), col="gray" ) )
      for ( i in 1:sum( diff.locs > 1 ) ) {
        if ( i == 1 ) lines( c( locs[ 1 ], locs[ diff.locs > 1 ][ i ] ) + maxlen - slen, c( j, j ),
               lwd=lwd + as.integer( rows[ jj ] %in% o.genes ), col="gray", lty=2 )
        else lines( c( locs[ which( diff.locs > 1 ) + 1 ][ i - 1 ], locs[ diff.locs > 1 ][ i ] ) + maxlen - slen, c( j, j ),
               lwd=lwd + as.integer( rows[ jj ] %in% o.genes ), col="gray", lty=2 )
      }
    }
    
    if ( ! is.null( rows ) ) {
      label <- rows[ jj ]
      col <- "black"
      if ( exists( "all.tfs" ) && label %in% all.tfs ) col <- "tomato3"
      if ( names( dev.cur() ) == "devSVG" || ! long.names && ! no.key ) {
        label <- substr( label, 1, 80 )
        text( maxlen*1.2, j, labels=label, adj=c(1,0.5), col=col, xpd=NA, ... ) ##cex=0.7,
      }
      if ( long.names || names( dev.cur() ) == "devSVG" ) {
        g.name <- toupper( label )
        ##gene.coords <- cluster$gene.coords
        ##if ( exists( "gene.coords" ) && is.vector( gene.coords ) && is.character( gene.coords ) ) {
        if ( ! is.null( cluster$gene.coords ) ) g.name <- cluster$gene.coords[ label ]
        ##} else {
          ##if ( exists( "gene.coords" ) && ! is.null( gene.coords$gene.func ) )
          ##  g.name <- gene.coords$gene.func[ g.name ]
        ##   if ( exists( "gene.coords" ) && ! is.null( gene.coords$gene.name ) &&
        ##       gene.coords$gene.name[ label ] != "-" && ! is.na( gene.coords$gene.name[ label ] ) &&
        ##       toupper( gene.coords$gene.name[ label ] ) != toupper( g.name ) ) {
        ##     ll <- gene.coords$gene.name[ label ]
        ##     ##if ( nchar( ll ) > 60 ) ll <- substr( ll, nchar( ll ), nchar( ll ) - 60 )
        ##     g.name <- paste( ll, ": ", g.name, sep="" )
        ##   }
        ## }

        g.name[ is.na( g.name ) ] <- label[ is.na( g.name ) ]
        if ( is.na( g.name ) ) g.name <- label
        if ( names( dev.cur() ) == "devSVG" ) {
          par( family="Arial" )
          setSVGShapeURL( paste( "http://www.genome.ad.jp/dbget-bin/www_bget?", organism, ":", label, sep="" ) )
          if ( ! is.na( g.name ) ) setSVGShapeToolTip( label, g.name )
          else setSVGShapeToolTip( label )
          rect( maxlen*1.2, j-0.18, maxlen, j+0.18, col=NA, border=NA, xpd=NA )
        } else if ( ! is.na( g.name ) ) { ##|| for.figure ) {
          lab <- label
          if ( toupper( g.name ) != toupper( label ) && g.name != "" ) {
            g.name <- gsub( "^[:\\s+]+", "", gsub( "\\s+$", "", g.name, perl=T ), perl=T ) ##maxlen/8
            if ( names( dev.cur() ) == "X11" ) g.name <- strtrim( g.name, 40 ) else gname <- strtrim( g.name, 60 ) 
            if ( label != "" ) lab <- paste( g.name, ": ", label, sep="" ) ##label, " (", g.name, ")", sep="" )
            else lab <- g.name
          }
          if ( nchar( lab ) > 60 ) lab <- substr( lab, nchar( lab ) - 60, nchar( lab ) )
          if ( ! no.key ) text( maxlen*1.2, j, labels=lab, adj=c(1,0.5), col=col, xpd=NA, ... ) ##cex=0.9,
        }
      }

      if ( ( ! all( is.na( p.clust ) ) || ! all( is.na( e.clust ) ) ) && ! no.key )
        text( -maxlen*0.07, j, labels=sprintf( "%.2f", p.values[ label ] ), 
             xpd=NA, col=if ( label %in% names( which( ! is.dup.seq ) ) ) "black" else "blue", ... ) 
    }
  }

  if ( ! no.key && ( ! all( is.na( p.clust ) ) && ! all( is.na( e.clust ) ) ) ) {##&& ! for.figure ) ##! no.plot && 
    text( -maxlen*0.15, length(rows)+0.9, labels=sprintf( "log10(P) %s", seq.type ), pos=4, ... ) ##cex=0.7 )
    mots.used <- sort( unique( mots.used ) )
    if ( length( mots.used ) > 1 ) {
      ##text( maxlen*0.03, length(rows)+0.9, "Motif legend:", xpd=NA, adj=c(0,0.5), ... )
      sapply( 1:length( mots.used ), function( j )
             text( maxlen*0.24 + (j+0)*maxlen*0.03, length(rows)+0.9, as.character( mots.used[ j ] ),
                  col=mots.used[ j ] + 1, xpd=NA, adj=c(0,0.5), ... ) )
    }
    n.unique.seqs <- sum( ! is.dup.seq )
    text( maxlen*1.2, length(rows)+0.9, sprintf( "log10(P.clust)=%.2f; %d seqs; %d uniq", p.clust[ seq.type ],
                                                length( seqs ), n.unique.seqs ), xpd=NA, adj=c(1,0.5), ... )
  }
  ##invisible( out.posns )
}

plotClust <- function( k, cluster=NULL, w.motifs=T, all.conds=T, title=NULL, o.genes=NULL, dont.plot=F,
                      network="all", short.names=organism == "sce", seq.type=names( mot.weights ), ... ) {
  if ( ! dont.plot && names( dev.cur() ) != "devSVG" ) opar <- par( no.readonly=T )
  if ( ! is.null( cluster ) ) {
    if ( ! dont.plot ) plotCluster.motif( cluster, seqs=cluster$seqs, p.val.shade.cutoff=1, o.genes=o.genes,
                                         no.plotCluster=all.conds, ... )
    return( invisible( cluster ) )
  }
  c <- get.clust( k, varNorm=F )
  rows <- get.rows( k ); if ( ! is.null( o.genes ) ) rows <- unique( c( rows, o.genes ) )
  if ( length( rows ) <= 0 ) { warning( "Trying to plot a cluster with no rows!" ); return() }
  if ( ! w.motifs && ! dont.plot ) { 
    if ( all.conds ) plotCluster.all.conds( c, o.genes=o.genes, ... )
    else plotCluster( c, o.genes=o.genes, ... )
  } else {
    c$seq.type <- seq.type
    for ( st in seq.type ) {
      c[[ st ]] <- list()
      c[[ st ]]$motif.out <- meme.scores[[ st ]][[ k ]]
      tmp <- cluster.pclust( k, st )
      c[[ st ]]$e.val <- tmp$e.vals
      c[[ st ]]$p.clust <- tmp$p.clusts
      c[[ st ]]$motif.out$pssms <- lapply( c[[ st ]]$motif.out$meme.out, "[[", "pssm" )
      c[[ st ]]$motif.out$e.values <- c[[ st ]]$e.val
      if ( ! is.null( c[[ st ]]$motif.out$pv.ev ) ) { ##&& ! is.null( meme.scores[[ st ]]$all.pv ) ) { ## Replicate old pv.ev list of 2 data frames
        if ( "gene" %in% colnames( c[[ st ]]$motif.out$pv.ev[[ 1 ]] ) ) c[[ st ]]$motif.out$pv.ev[[ 2 ]] <- c[[ st ]]$motif.out$pv.ev[[ 1 ]]
        ##!else if ( "gene" %in% colnames( c[[ st ]]$motif.out$pv.ev[[ 2 ]] ) ) c[[ st ]]$motif.out$pv.ev[[ 2 ]] <- c[[ st ]]$motif.out$pv.ev[[ 2 ]]
        if ( ! is.null( meme.scores[[ st ]]$all.pv ) ) {
          tmp <- cbind( p.value=meme.scores[[ st ]]$all.pv[ ,k ],
                 e.value=if ( "all.ev" %in% names( meme.scores[[ st ]] ) ) meme.scores[[ st ]]$all.ev[ ,k ] else NA )
        } else {
          pv.ev <- meme.scores[[ st ]][[ k ]]$pv.ev[[ 1 ]]
          if ( ncol( pv.ev ) <= 2 ) pv.ev <- meme.scores[[ st ]][[ k ]]$pv.ev[[ 2 ]]
          tmp <- NULL
          if ( ncol( pv.ev ) > 0 ) {
            tmp <- as.matrix( pv.ev[ ,2:ncol( pv.ev ) ] )
            rownames( tmp ) <- pv.ev[ ,1 ]; colnames( tmp ) <- c( "p.value", "posns", "mots" )
          }
        }         
        c[[ st ]]$motif.out$pv.ev[[ 1 ]] <- tmp
        if ( ! is.null( tmp ) ) {
          c[[ st ]]$motif.out$p.values <- log10( c[[ st ]]$motif.out$pv.ev[[ 1 ]][ ,"p.value" ] )
          names( c[[ st ]]$motif.out$p.values ) <- rownames( c[[ st ]]$motif.out$pv.ev[[ 1 ]] )
        }
      }
    }
  }

  if ( ! is.na( mot.iters[ 1 ] ) && ! no.genome.info ) {
    ## Filtering is used for plotting masked subseqs (which are not input to MEME) as dashed lines
    c$seqs <- get.sequences( rows, distance=motif.upstream.scan[[ seq.type[ 1 ] ]],
                                seq.type=seq.type[ 1 ], filter=T, uniq=F ) ##, remove.repeats=T, mask=T, blast=T ) ##[ rows ]
    tmp <- c$seqs[ rows ]
    if ( ! is.null( tmp ) ) names( tmp ) <- rows
    attr( tmp, "start.stops" ) <- attr( c$seqs, "start.stops" ); c$seqs <- tmp; rm( tmp )
  } else c$seqs <- NULL
  if ( ! is.na( net.iters[ 1 ] ) ) {
    if ( network == "all" ) network <- names( networks )
    for ( i in network ) {
      if ( ! i %in% names( networks ) ) next
      tmp.net <- networks[[ i ]][ networks[[ i ]]$protein1 %in% rows & networks[[ i ]]$protein2 %in% rows, ]
      tmp.net <- cbind( tmp.net, net=rep( i, nrow( tmp.net ) ) )
      c$network <- if ( ! is.null( c$network ) ) rbind( c$network, tmp.net ) else tmp.net
    }
  }
  c$gene.coords <- get.long.names( rows, short=short.names )
  if ( "cog.code" %in% names( genome.info ) ) c$cog.code <- genome.info$cog.code[ rows ] 
  if ( ! is.null( title ) ) c$name <- title
  ##cluster<<-c
  if ( ! dont.plot ) {
    plotCluster.motif( c, seqs=c$seqs, p.val.shade.cutoff=1, o.genes=o.genes, no.plotCluster=all.conds, ... )
    if ( names( dev.cur() ) != "devSVG" && ! "layout" %in% names( list( ... ) ) ) par( opar )
  }
  invisible( c )
}

plotScores <- function( k, o.genes=NULL, b.genes=NULL, recompute=F ) { 
  opar <- par( no.readonly=T )
  rows <- get.rows( k ) 

  if ( recompute || ! exists( "row.scores" ) || is.null( row.scores ) ) {
    if ( attr( get.all.scores, "version" ) == 1 ) {
      tmp <- get.all.scores( k, force.row=T, force.col=T, force.motif=T, force.net=T )
      rs <- tmp$r; ms <- tmp$m; ns <- tmp$n; cs <- tmp$c
    } else if ( attr( get.all.scores, "version" ) == 2 ) {
      tmp <- get.old.scores.matrices( k )
      rs <- tmp$r[ ,1 ]; ms <- tmp$m[ ,1 ]; ns <- tmp$n[ ,1 ]
      if ( all( is.na( ms ) ) || all( ms == ms[ 1 ] ) ) rm( ms )
      if ( all( is.na( ns ) ) || all( ns == ns[ 1 ] ) ) rm( ns )
    }
  }

  tmp.scale <- round( attr( ratios, "nrow" ) / length( rows ) / 4 )
  layout( matrix( c( 1, 2, 3, 4, 4, 5, 4, 4, 6 ), 3, 3, byrow=T ) )
  if ( ! exists( "rs" ) ) rs <- row.scores[ ,k, drop=T ]
  rs[ rs < -220 ] <- min( rs[ rs > -220 ], na.rm=T )
  h <- try( hist( rs, breaks=20, main=paste( "Cluster", k ), xlab="Ratios scores" ), silent=T )
  if ( class( h ) != "try-error" ) {
    try( hist( rep( rs[ rows ], tmp.scale ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
    try( hist( rs, breaks=h$breaks, add=T ), silent=T )
  }

  if ( exists( "ms" ) || ( ! is.null( mot.scores ) && ! all( is.na( mot.scores[ ,k ] ) ) && ! no.genome.info ) ) { 
    if ( ! exists( "ms" ) ) ms <- mot.scores[ ,k, drop=T ]; ms[ ms < -20 ] <- min( ms[ ms > -20 ], na.rm=T )
    h <- try( hist( ms, breaks=20, main=NULL, xlab="Motif scores" ), silent=T ) 
    if ( class( h ) != "try-error" ) {
      try( hist( rep( ms[ rows ], tmp.scale * 3 ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( ms, breaks=h$breaks, add=T ), silent=T )
    }
  } else { ms <- NULL; plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" ) }

  if ( exists( "ns" ) || ( ! is.null( net.scores ) && ! all( net.scores[ ,k ] == 0 ) ) ) { 
    if ( ! exists( "ns" ) ) { ns <- net.scores[ ,k, drop=T ]; ns[ ns < -20 ] <- min( ns[ ns > -20 ], na.rm=T ); ns <- -log10( -ns ) }
    ns[ is.infinite( ns ) ] <- max( ns[ ! is.infinite( ns ) ], na.rm=T ) + 0.1
    h <- try( hist( ns, breaks=20, main=NULL, xlab="-log10(-Network scores)" ), silent=T ) 
    if ( class( h ) != "try-error" ) {
      try( hist( rep( ns[ rows ], tmp.scale / 3 ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( ns, breaks=h$breaks, add=T ), silent=T )
    }
  } else { ns <= NULL; plot( 1, 1, typ="n", axes=F, xaxt="n", yaxt="n", xlab="", ylab="" ) }
  
  row.memb <- attr( ratios, "rnames" ) %in% rows
  if ( ! is.null( ms ) && ! all( is.na( ms ) ) && ! all( ns == 0 ) ) {
    plot( rs, ms, typ="n", main=paste( "Cluster", k ), xlab="Ratios scores", ylab="Mot scores" ) ##,log="xy")
    text( rs, ms, label=1:length(rs), col=row.memb+1, cex=0.5 )
  } else if ( ! is.null( ns ) && ! all( ns == 0 ) ) {
    plot( rs, ns, typ="n", main=paste( "Cluster", k ), xlab="Ratios scores", ylab="Net scores" )
    text( rs, ns, label=1:length(rs), col=row.memb+1, cex=0.5 )
  } else {
    plot( rs, jitter(rep(0,length(rs))), typ="n", main=paste( "Cluster", k ), xlab="Ratios scores", ylab="" ) 
    text( rs, jitter(rep(0,length(rs))), label=1:length(rs), col=row.memb+1, cex=0.5 )
  }
  if ( ! is.null( o.genes ) ) text( rs[ o.genes ], ms[ o.genes ],
                                   label=which( attr( ratios, "rnames" ) %in% o.genes ), col="green", cex=0.5 )
  if ( ! is.null( b.genes ) ) text( rs[ b.genes ], ms[ b.genes ],
                                   label=which( attr( ratios, "rnames" ) %in% b.genes ), col="blue", cex=0.5 )

  try( {
    tmp <- get.combined.scores( quant=F )
    r.scores <- tmp$r; c.scores <- tmp$c ##; n.scores <- tmp$n; m.scores <- tmp$m
    rr <- get.density.scores( ks=k, r.scores, c.scores, plot="rows" )$r ##ks=k, plot=T )$r
    rr <- rr[ ,k, drop=T ] ##; names( rr ) <- attr( ratios, "rnames" ) ## hack!
    h <- try( hist( log10( rr ), breaks=50, main=NULL, xlab="Density (membership) scores" ), silent=T ) 
    if ( class( h ) != "try-error" ) {
      try( hist( rep( log10( rr[ rows ] ), tmp.scale ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( log10( rr ), breaks=h$breaks, add=T ), silent=T )
    }
  }, silent=T )
  par( opar )
}

## Note can do pdf("xxx.pdf");plotStats(new.dev=F,plot.clust=1);dev.off() to redirect all 3 plots to a single pdf.
plotStats <- function( iter=stats$iter[ nrow( stats ) ], plot.clust=NA, new.dev=F, ... ) { ##new=F,
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
      tmp <- get.combined.scores( quant=T )
      r.scores <- tmp$r; c.scores <- tmp$c ##; n.scores <- tmp$n; m.scores <- tmp$m
    } else if ( attr( get.all.scores, "version" ) == 2 ) {
      tmp <- get.old.scores.matrices()
      row.scores <- tmp$r; mot.scores <- tmp$m; net.scores <- tmp$n; col.scores <- tmp$c
    }
  }  
  
  opar <- par( no.readonly=T )
  tmp.scale <- round( 1 / mean( row.memb, na.rm=T ) / 4 )
  if ( new.dev ) { if ( length( dev.list() ) < 1 ) dev.new(); dev.set( 2 ) } ##;
  layout( matrix( c( 1,2,3,
                     1,2,3,
                     4,5,6,
                     4,5,6,
                     7,9,11,
                     8,10,12 ), byrow=T, ncol=3 ) )
  par( mar=c(3,3,2,0.1), mgp=c(3,1,0)*0.5 )
  stats <- stats[ stats[, "iter" ] <= iter, ,drop=F ]
  try( matplot( stats[ ,"iter" ], stats[ ,grep( "resid", colnames( stats ), val=T ) ], typ="l", xlab="iter",
               ylab="Mean resid", main=sprintf( "Iter: %d", iter ), lty=1 ), silent=T )
  sapply( c( 51, 101, 21 ), function( kwin )
         try( matlines( stats[ ,"iter" ], apply( stats[ ,grep( "resid", colnames( stats ), val=T ), drop=F ], 2,
                           function( i ) runmed( i, k=min( length( i ), kwin ) ) ), lty=2, lwd=0.6 ), silent=T ) )
  if ( ( nn <- length( grep( "resid", colnames( stats ) ) ) ) > 1 )
    legend( "bottomleft", legend=gsub( "resid.", "", grep( "resid", colnames( stats ), val=T ) ), lwd=1, bty="n", 
           col=1:nn, lty=1:nn, cex=0.5 )
  ##try( matplot( stats[ ,"iter", drop=F ], stats[ ,c("nrow","ncol"), drop=F ], typ="l", xlab="iter",
  ##             ylab="Mean nrow, ncol" ), silent=T ) ## This stat is pretty useless
  if ( exists( "row.scores" ) && ! is.null( mot.scores ) && ! all( is.na( mot.scores[,] ) ) ) {
    rs <- row.scores[]; rs[ rs < -20 ] <- min( rs[ rs > -20 ], na.rm=T )
    h <- try( hist( rs, breaks=50, main=NULL, xlab="Ratios scores" ), silent=T ) 
    if ( class( h ) != "try-error" ) {
      try( hist( rep( rs[ row.memb == 1 ], tmp.scale ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( rs, breaks=h$breaks, add=T ), silent=T )
    }
  }
  if ( exists( "mot.scores" ) && ! is.null( mot.scores ) && ! all( is.na( mot.scores[,] ) ) ) {
    ms <- mot.scores[,]; ms[ ms < -20 ] <- min( ms[ ms > -20 ], na.rm=T ); ms[ ms >= 0 ] <- NA
    h <- try( hist( ms, breaks=50, main=NULL, xlab="Motif scores" ), silent=T )
    if ( class( h ) != "try-error" ) {
      try( hist( rep( ms[ row.memb == 1 ], tmp.scale * 3 ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( ms, breaks=h$breaks, add=T ), silent=T ) ##main="Mot scores",
    }
  }
  tmp <- stats[ ,grep( "p.clust", colnames( stats ), val=T ), drop=F ]
  if ( ! all( is.na( tmp ) ) ) {
    try( matplot( stats[ ,"iter" ], tmp, typ="l", xlab="iter", ylab="Mean motif p-value",
                 main=sprintf( "Motif scaling: %.3f", mot.scaling[ max( 1, iter - 1 ) ] ), lty=1 ), silent=T )
    sapply( c( 51, 101, 21 ), function( kwin ) try( matlines( stats[ ! is.na( tmp ), "iter" ],
                      apply( tmp, 2, function( i ) runmed( i[ ! is.na( i ) ], k=min( sum( ! is.na( i ) ), kwin ) ) ),
                         lty=2, lwd=0.6 ), silent=T ) )
    if ( ( nn <- length( grep( "p.clust", colnames( stats ) ) ) ) > 1 )
      legend( "bottomleft", legend=gsub( "p.clust.", "", grep( "p.clust", colnames( stats ), val=T ) ), lwd=1,
             bty="n", col=1:nn, lty=1:nn, cex=0.5 )
  }
  if ( exists( "net.scores" ) && ! is.null( net.scores ) && ! all( net.scores[,] == 0 ) ) {
    ns <- net.scores[,]; ns[ ns < -20 ] <- min( ns[ ns > -20 ], na.rm=T ); ns[ ns >= 0 ] <- NA; ns[,] <- -log10( -ns )
    tmp.scale <- ceiling( tmp.scale * mean( ! is.na( ns ), na.rm=T ) )
    h <- try( hist( ns, breaks=50, main=NULL, xlab="-log10(-Network scores)" ), silent=T ) 
    if ( class( h ) != "try-error" ) {
      try( hist( rep( ns[ row.memb == 1 ], tmp.scale ), breaks=h$breaks, col="red", border="red", add=T ), silent=T )
      try( hist( ns, breaks=h$breaks, add=T ), silent=T ) 
    }
    tmp <- stats[ ,grep( "net.", colnames( stats ), val=T, fixed=T ), drop=F ]
    try( matplot( stats[ ,"iter" ], tmp, typ="l", xlab="iter", ylab="Mean net-score",
                 main=sprintf( "Net scaling: %.3f", net.scaling[ max( 1, iter - 1 ) ] ), lty=1 ), silent=T )
    try( matlines( stats[ ,"iter" ],
                  apply( tmp, 2, function( i ) runmed( i[ ! is.na( i ) ], k=min( sum( ! is.na( i ) ), 51 ) ) ),
                  lty=2, lwd=0.6 ), silent=T )
    if ( ( nn <- length( grep( "net.", colnames( stats ) ) ) ) > 1 )
      try( legend( "bottomleft", legend=gsub( "net.", "", grep( "net.", colnames( stats ), val=T ) ), lwd=1, bty="n",
                  col=1:nn, lty=1:nn, cex=0.5 ), silent=T )
  }
  
  clusterStack <- get.clusterStack( ks=1:k.clust )
  resids <- sapply( as.list( clusterStack ), "[[", "resid" )
  try( hist( resids[ resids <= 1.5 ], main=NULL, xlab="Cluster Residuals",
            xlim=if ( all( resids > 0, na.rm=T ) ) c( 0, 1.5 ) else range( resids, na.rm=T ), breaks=k.clust/4 ),
      silent=T )
  if ( ! exists( "mot.scores" ) || is.null( mot.scores ) || all( is.na( mot.scores[,] ) ) ) {
    pclusts <- sapply( as.list( clusterStack ), "[[", "p.clust" )
    try( hist( pclusts[ pclusts <= 1 ], main=NULL, xlab="Cluster Motif P-values",
              xlim=if ( all( pclusts > 0, na.rm=T ) ) c( 0, 1 ) else range( pclusts, na.rm=T ), breaks=k.clust/4 ),
        silent=T )
  }
  
  if ( mot.scaling[ iter ] > 0 ) { ##! is.null( mot.scores ) && ! all( is.na( mot.scores[,] ) ) && ! all( mot.scores[,] == 0 ) ) {
    plot.all.clusterMotifPositions <- function( ks=1:k.clust, mots=1, e.cutoff=1, p.cutoff=0.05,
                                               seq.type=names( mot.weights )[ 1 ], breaks=100, ... ) {
      if ( seq.type == "ALL" ) seq.type <- names( mot.weights )
      df <- NULL
      for ( st in seq.type ) {
        ms <- meme.scores[[ st ]][ ks ]
        ind <- 1
        if ( ! "posns" %in% colnames( ms[[ 1 ]]$pv.ev[[ 1 ]] ) ) ind <- 2
        posns <- as.vector( unlist( sapply( ms, function( i ) i$pv.ev[[ ind ]]$posns ) ) )
        pvals <- as.vector( unlist( sapply( ms, function( i ) i$pv.ev[[ ind ]]$pvals ) ) )
        imots <- as.vector( unlist( sapply( ms, function( i ) i$pv.ev[[ ind ]]$mots ) ) )
        genes <- as.vector( unlist( sapply( ms, function( i ) as.character( i$pv.ev[[ ind ]]$gene ) ) ) )
        slens <- nchar( genome.info$all.upstream.seqs[ st ][[ st ]][ genes ] )
        clusts <- as.vector( unlist( sapply( ms, function( i )
                                 rep( i$k, if ( is.null( i$pv.ev[[ ind ]] ) ) 0 else nrow( i$pv.ev[[ ind ]] ) ) ) ) )
        ms <- meme.scores[[ st ]]
        evals <- sapply( 1:length( imots ), function( i )
                        ms[[ clusts[ i ] ]]$meme.out[[ abs( imots[ i ] ) ]]$e.value )
        df <- rbind( df, data.frame( clusts, posns, pvals, imots, evals, genes, slens,
                                    seq.type=rep( st, length( clusts ) ) ) )
      }
      df2 <- subset( df, evals < e.cutoff & pvals < p.cutoff & abs( imots ) %in% mots )
      psns <- df2$posns - df2$slens - do.call( rbind, motif.upstream.scan[ df2$seq.type ] )[ ,1 ]
      ##for ( i in 1:length( psns ) ) psns[ i ] <- psns[ i ] - motif.upstream.scan[[ df2$seq.type[ i ] ]][ 2 ]
      ##stop()
      ##psns <- psns - do.call( rbind, motif.upstream.scan[ df2$seq.type ] )[ ,2 ]
      h <- hist( psns, breaks=breaks,
                xlab=sprintf( "%s %s", paste( seq.type, collapse=" " ), paste( mots, collapse=" " ) ), ... )
      dd <- density( psns, bw=5 )
      lines( dd$x, dd$y * max( h$counts ) / max( dd$y ) * 0.9, col="red" )
      invisible( df ) ##data.frame( clusts, posns, pvals, imots, evals ) )
    }
    
    try( plot.all.clusterMotifPositions( main="Positions of motif #1", ... ), silent=T ) ##xlab="Position upstream", 
    if ( any( sapply( meme.scores$upstream[ 1:k.clust ], function( i ) length( i$meme.out ) ) == 2 ) )
      try( plot.all.clusterMotifPositions( mots=2, main="Positions of motif #2", ... ), ##xlab="Position upstream", 
          silent=T )
  }
  
  n.rows <- sapply( 1:k.clust, function( k ) length( get.rows( k ) ) )
  try( hist( n.rows, main=NULL, xlab="Cluster Nrows", breaks=k.clust/4,
            xlim=c( -5, max( n.rows, na.rm=T ) ) ), silent=T )
  n.cols <- sapply( 1:k.clust, function( k ) length( get.cols( k ) ) )
  try( hist( n.cols, main=NULL, xlab="Cluster Ncols", breaks=k.clust/4,
            xlim=c( -5, attr( ratios, "ncol" ) ) ), silent=T )
  nr <- table( unlist( sapply( clusterStack, "[[", "rows" ) ) )
  if ( length( nr ) < attr( ratios, "nrow" ) ) nr <- c( nr, rep( 0, attr( ratios, "nrow" ) - length( nr ) + 1 ) )
  try( hist( nr, breaks=seq(-0.5,max(nr,na.rm=T)+0.5,by=1), main=NULL, xlab="NClust per gene" ), silent=T )
  ##}, silent=T )
  
  if ( ! is.na( plot.clust ) ) { 
    if ( new.dev ) { if ( length( dev.list() ) < 2 ) dev.new(); dev.set( 3 ) }
    try( plotClust( plot.clust, w.motifs=T, cex=0.7 ), silent=T )
    if ( new.dev ) { if ( length( dev.list() ) < 3 ) dev.new(); dev.set( 4 ) }
    try( plotScores( plot.clust ), silent=T )
  }
  par( opar )
}

write.project <- function( ks=sapply( as.list( clusterStack ), "[[", "k" ), para.cores=1, ##save.session=T, ##pdfs=T, ##dev="SVG", 
                          out.dir=NULL, gaggle=T, seq.type=names( mot.weights )[ 1 ], gzip=T,
                          output=c("svg","pdf","png","html","main","rdata"), ... ) { ##network=F, 
  if ( is.null( out.dir ) ) {
    out.dir <- cmonkey.filename
    if ( iter != n.iter ) out.dir <- sprintf( "%s_%04d", out.dir, iter )
  }
  cat( "Outputing to", out.dir, "\n" )
  ##if ( ! is.null( dev ) && dev == "SVG" )
  ##require( RSVGTipsDevice )
  if ( ! file.exists( out.dir ) ) ##&& ! is.null( dev ) && dev == "SVG" )
    dir.create( out.dir, recursive=T, showWarnings=F )
  ##require( fUtilities )

  ##clusterStack <- get.clusterStack( ks=1:k.clust )
  ##ks <- sapply( clusterStack, "[[", "k" )
  clusterStack <- clusterStack[ ks ]
  mc <- get.parallel( length( ks ), para.cores=para.cores )
  
  if ( ! file.exists( paste( out.dir, "/svgs", sep="" ) ) )
    dir.create( paste( out.dir, "/svgs", sep="" ), showWarnings=F )
  if ( ! file.exists( paste( out.dir, "/pdfs", sep="" ) ) ) ##pdfs && 
    dir.create( paste( out.dir, "/pdfs", sep="" ), showWarnings=F )
  if ( ! file.exists( paste( out.dir, "/htmls", sep="" ) ) )
    dir.create( paste( out.dir, "/htmls", sep="" ), showWarnings=F )
  
  ##if ( ! is.null( dev ) && dev == "SVG" ) {
  if ( "svg" %in% output ) {
    require( RSVGTipsDevice )
    if ( ! file.exists( sprintf( "%s/svgs/stats.svg", out.dir ) ) &&
        ! file.exists( sprintf( "%s/svgs/stats.svgz", out.dir ) ) ) {
      ##cat( out.dir, "/svgs/stats.svg", "\n", sep="" )
      cat( "STATS...\n" )
      devSVGTips( sprintf( "%s/svgs/stats.svg", out.dir ), toolTipMode=2, title="Biclustering statistics",
                 xmlHeader=T )
      par( family="Arial" )
      plotStats( new.dev=F )
      dev.off()
    }

    require( igraph0 )
    cat( "SVGS: " )
    for ( qqq in 1:3 ) {
    ##mc$
      ##lapply( ks, function( k ) { ## will this work in parallel? seems to. No, it doesn't.
      for ( k in ks ) {
        ##k <- ks[ i ]
        if ( k %% 25 == 0 ) cat( k ) else cat( "." )
        if ( file.exists( sprintf( "%s/svgs/cluster%04d.svg", out.dir, k ) ) ||
            file.exists( sprintf( "%s/svgs/cluster%04d.svgz", out.dir, k ) ) ) next ##return( NULL )
        ## Note this is a HACK to get correctly formatted svgs for some reason
        tmp.cl <- try( plotClust( k, w.motifs=T, seq.type=seq.type, dont.plot=T, ... ) )
        devSVGTips( sprintf( "%s/svgs/cluster%04d.svg", out.dir, k ), toolTipMode=2,
                   title=sprintf( "Bicluster %04d", k ), xmlHeader=T )
        try( plotClust( cluster=tmp.cl, w.motifs=T, seq.type=seq.type, ... ) )
        dev.off()
      } ##)
    }
    cat( "\n" )
  }

  ##if ( pdfs ) {
  if ( "pdf" %in% output ) {
    require( igraph0 ) ## load it here
    cat( "PDFS: " )
    ##mc$
    lapply( ks, function( k ) { ## will this work in parallel? seems to.
      ##k <- ks[ i ]
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/pdfs/cluster%04d.pdf", out.dir, k ) ) ) next ##return( NULL )
        pdf( sprintf( "%s/pdfs/cluster%04d.pdf", out.dir, k ) ) ## Will be compressed later
      #try(
      plotClust( k, w.motifs=T, seq.type=seq.type, ... )#, silent=T )
      dev.off()
    } )
    cat( "\n" )
  }
  
  if ( gaggle && "html" %in% output ) { ## Embed the SVG in some gaggleized html for use w/ firegoose --
    require( hwriter ) ## see http://www.ebi.ac.uk/~gpau/hwriter/
    ## see http://gaggle.systemsbiology.net/docs/geese/firegoose/microformat/
    ## Note that "embed" might not be the best way to embed the SVG --
    ##    see http://www.w3schools.com/svg/svg_inhtml.asp for other options.
    ##'%^%' <- function( a, b ) paste( a, b, sep="\n" )
    cat( "HTMLS: " )
    ##mc$
    lapply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d.html", out.dir, k ) ) ) return()
      rows <- sort( get.rows( k ) )
      if ( length( rows ) <= 0 ) return()
      short.names <- try( get.long.names( rows, short=T ) )
      if ( class( short.names ) == 'try-error' ) short.names <- rep( '', length( rows ) )
      short.names <- cbind( rows, short.names )
      rownames( short.names ) <- colnames( short.names ) <- NULL
      long.names <- try( get.long.names( rows, short=F ) )
      if ( class( long.names ) == 'try-error' ) long.names <- rep( '', length( rows ) )
      long.names <- cbind( rows, long.names )
      rownames( long.names ) <- colnames( long.names ) <- NULL
      refseq.names <- unique( unlist( get.synonyms( rows ) ) )
      refseq.names <- grep( "^NP_", refseq.names, val=T )
      upstream.seqs <- try( get.sequences( k, filter=F, uniq=F ), silent=T ) ##uniquify=F, silent=T )
      if ( class( upstream.seqs ) == "try-error" || is.null( upstream.seqs ) || length( upstream.seqs ) == 0 ) {
        upstream.seqs <- rep( "", length( rows ) ); names( upstream.seqs ) <- rows }
      upstream.seqs <- cbind( names( upstream.seqs ), upstream.seqs )
      rownames( upstream.seqs ) <- colnames( upstream.seqs ) <- NULL

      htmltext <-
        paste( c( "<html><head><title>Bicluster %K (%FILE)</title>", ##, k, cmonkey.filename ),
                 "<style type=\"text/css\">",
                 "  .hidden {", "     display: none;", "   }",
                 "  .gaggle-data {", "     color: green;", "     font-size: xx-small;", "   }",
                 "   p {", "     color: red;", "     font-size: x-small;", "   }",
                 "</style>",
                 "<script type=\"text/javascript\">",
                 "   function toggleVisible(id){",
                 "      if (document.getElementById){",
                 "         obj = document.getElementById(id);",
                 "         if (obj) {",
                 "            if (obj.style.display == 'none'){",
                 "               obj.style.display = 'block';",
                 "            } else {",
                 "               obj.style.display = 'none';",
                 "            }",
                 "         }",
                 "      }",
                 "   }",
                 "</script>",
                 "</head>",
                 "<table><tr><td>",
                 "<iframe src=\"../svgs/cluster%K03d%K.svg\" width=\"600\" height=\"520\" frameborder=\"0\"></iframe>",
                 "</td><td>",
                 
                 "<p><a href=\"#bicluster%K03d%K\" onclick=\"toggleVisible('bicluster%K03d%K'); return false;\">[+]</a>",
                 "Show/hide bicluster #%K rows and columns.</p>", 
                 "<div id=\"bicluster%K03d%K\" style=\"display:none;\" class=\"gaggle-data bicluster\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>", length( rows ), length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-cluster\">",
                 "      <ol class=\"gaggle-rowNames\">",
                 paste( "<li>", sort( rows ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   <ol class=\"gaggle-columnNames\">",
                 paste( "<li>", sort( get.cols( k ) ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_genes\" onclick=\"toggleVisible('bicluster%K03d%K_genes'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (genes).</p>", 
                 "<div id=\"bicluster%K03d%K_genes\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K genes</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( rows ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( rows ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_short_names\" onclick=\"toggleVisible('bicluster%K03d%K_short_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (short gene names).</p>", 
                 "<div id=\"bicluster%K03d%K_short_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K short names</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( short.names ) ),
                 "   <span class=\"gaggle-namelist-tag hidden\">short_name</span>",
                 hwrite( short.names, table.class="toc", col.class=list( NA, "short_name" ), border=1,
                        table.style="font-family: monospace; font-size: xx-small; color: green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_long_names\" onclick=\"toggleVisible('bicluster%K03d%K_long_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (long gene names).</p>", 
                 "<div id=\"bicluster%K03d%K_long_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K long names</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( long.names ) ),
                 "   <span class=\"gaggle-namelist-tag hidden\">long_name</span>",
                 hwrite( long.names, table.class="toc", col.class=list( NA, "long_name" ), border=1,
                        table.style="font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_refseq_names\" onclick=\"toggleVisible('bicluster%K03d%K_refseq_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (NCBI RefSeq gene IDs).</p>", 
                 "<div id=\"bicluster%K03d%K_refseq_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K NCBI RefSeq IDs</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( refseq.names ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( refseq.names ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_upstream_seqs\" onclick=\"toggleVisible('bicluster%K03d%K_upstream_seqs'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K gene upstream sequences.</p>", 
                 "<div id=\"bicluster%K03d%K_upstream_seqs\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K upstream sequences</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( upstream.seqs ) ), 
                 "   <span class=\"gaggle-namelist-tag hidden\">upstream</span>",
                 hwrite( upstream.seqs, table.class="toc", col.class=list( NA, "upstream" ), border=1, 
                        table.style="font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_arrays\" onclick=\"toggleVisible('bicluster%K03d%K_arrays'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K columns (arrays; conditions).</p>", 
                 "<div id=\"bicluster%K03d%K_arrays\" style=\"display:none;\" class=\"gaggle-data arrays\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K arrays</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( get.cols( k ) ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_ratios\" onclick=\"toggleVisible('bicluster%K03d%K_ratios'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K ratios.</p>", 
                 "<div id=\"bicluster%K03d%K_ratios\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K ratios</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>", length( rows ), length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-matrix-tsv\">",
                 "        RATIOS",
                 "   </div>",
                 "</div>",

                 if ( ! is.null( seq.type ) && ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out ) &&
                     ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]] ) )
                 paste( "<p><a href=\"#bicluster%K03d%K_pssm1\" onclick=\"toggleVisible('bicluster%K03d%K_pssm1'); return false;\">[+]</a>", 
                       "Show/hide bicluster #%K motif PSSM #1.</p>", 
                       "<div id=\"bicluster%K03d%K_pssm1\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                       "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #1</span>", 
                       "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                       sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>",
                               nrow( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm ),
                               ncol( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm ) ),
                       "   <div class=\"gaggle-matrix-tsv\">",
                       "           MOTIF1",
                       "   </div>",
                       "</div>" ) else "",

                 if ( ! is.null( seq.type ) && length( meme.scores[[ seq.type ]][[ k ]]$meme.out ) >= 2 &&
                     ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]] ) )
                 paste( "<p><a href=\"#bicluster%K03d%K_pssm2\" onclick=\"toggleVisible('bicluster%K03d%K_pssm2'); return false;\">[+]</a>", 
                       "Show/hide bicluster #%K motif PSSM #2.</p>", 
                       "<div id=\"bicluster%K03d%K_pssm2\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                       "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #2</span>", 
                       "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                       sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>",
                               nrow( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm ),
                               ncol( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm ) ), 
                       "   <div class=\"gaggle-matrix-tsv\">",
                       "           MOTIF2",
                       "   </div>",
                       "</div>" ) else "",
                 
                 "</td></table>",
                 if ( "pdf" %in% output ) sprintf( "<a href=\"../pdfs/cluster%04d.pdf\">View PDF version</a>", k ) else "",
                 "</html>" ), collapse="\n" )

      rm( short.names, long.names, refseq.names, upstream.seqs )

      htmltext <- gsub( "%K03d%K", sprintf( "%04d", k ), htmltext )
      htmltext <- gsub( "%K", k, htmltext )
      htmltext <- gsub( "%FILE", cmonkey.filename, htmltext )
      htmltext <- gsub( "%SPECIES", gsub( "_", " ", rsat.species ), htmltext )

      tmp <- as.data.frame( get.cluster.matrix( rows, get.cols( k ) ) )
      tmp <- cbind( GENES=rownames( tmp ), tmp )
      tf <- tempfile()
      write.table( tmp, file=tf, sep="\t", quote=F, row.names=F ); rm( tmp )
      htmltext <- sub( "RATIOS", paste( readLines( tf ), collapse="\n" ), htmltext )
      unlink( tf )

      if ( ! is.null( seq.type ) && ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out ) ) {
        if ( ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]] ) ) {
          tmp <- as.data.frame( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm )
          if ( ! is.null( tmp ) && nrow( tmp ) > 0 ) {
            tmp <- cbind( 1:nrow( tmp ), tmp )
            colnames( tmp ) <- c( "POSITION", "A", "C", "G", "T" )
            write.table( tmp, file=tf, sep="\t", quote=F, row.names=F )
            htmltext <- sub( "MOTIF1", paste( readLines( tf ), collapse="\n" ), htmltext )
            unlink( tf )
          }
          rm( tmp )
        }
        
        if ( length( meme.scores[[ seq.type ]][[ k ]]$meme.out ) >= 2 &&
            ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]] ) ) {
          tmp <- as.data.frame( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm )
          if ( ! is.null( tmp ) && nrow( tmp ) > 0 ) {
            tmp <- cbind( 1:nrow( tmp ), tmp )
            colnames( tmp ) <- c( "POSITION", "A", "C", "G", "T" )
            write.table( tmp, file=tf, sep="\t", quote=F, row.names=F )
            htmltext <- sub( "MOTIF2", paste( readLines( tf ), collapse="\n" ), htmltext )
            unlink( tf )
          }
          rm( tmp )
        }
      }
      rm( tf )

      cat( htmltext, file=sprintf( "%s/htmls/cluster%04d.html", out.dir, k ), sep="\n" )
      rm( htmltext )
    } ) ##, mc.preschedule=F )
    cat( "\n" )
  }

  ##require( Cairo )
  if ( "png" %in% output ) {
    ##parallel.cores <- 1
    mc <- get.parallel( length( ks ), para=1 )
    cat( "PROFILES: " )
    for ( k in ks ) {
    ##mc$
    ##lapply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d_profile.png", out.dir, k ) ) ) return() ##next
      try( {
        c <- get.clust( k )
        png( sprintf( "%s/htmls/cluster%04d_profile.png", out.dir, k ), width=128, height=64, antialias="subpixel" )
        par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
        try( plotCluster( c, main="", no.par=T, ... ) )
        dev.off() }, silent=T )
    } ##)
    cat( "\n" )
    
    cat( "NETWORKS: " )
    require( igraph0 )
    ##mc$
    ##lapply( ks, function( k, ... ) {
    for ( k in ks ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d_network.png", out.dir, k ) ) ) return() ##next
      try( {
        png( sprintf( "%s/htmls/cluster%04d_network.png", out.dir, k ), width=64, height=64, antialias="subpixel" )
        par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
        c <- get.clust( k )
        try( plotCluster.network( c, cex=0.3, no.legend=T, ... ) )
        dev.off() }, silent=T )
    } ##)
    cat( "\n" )

    if ( ! is.null( seq.type ) ) {
      cat( "MOTIFS: " )
      for ( k in ks ) {
      ##mc$
      ##lapply( ks, function( k, ... ) {
        if ( k %% 25 == 0 ) cat( k ) else cat( "." )
        e.vals <- lapply( meme.scores[[ seq.type ]][[ k ]]$meme.out, "[[", "e.value" )
        pssms <- lapply( meme.scores[[ seq.type ]][[ k ]]$meme.out, "[[", "pssm" )
        if ( length( pssms ) < 2 ) {
          for ( i in ( length( pssms ) + 1 ):2 ) {
            pssms[[ i ]] <- matrix( 0.25, nrow=6, ncol=4 )
            e.vals[[ i ]] <- Inf
          }
        }
        for ( pp in 1:length( pssms ) ) {
          if ( file.exists( sprintf( "%s/htmls/cluster%04d_pssm%d.png", out.dir, k, pp ) ) ) return() ##next
          try( { 
            png( sprintf( "%s/htmls/cluster%04d_pssm%d.png", out.dir, k, pp ), width=128, height=64,
                antialias="subpixel" )
            if ( is.matrix( pssms[[ pp ]] ) )
              try( viewPssm( pssms[[ pp ]], e.val=NA, mot.ind=pp, main.title=sprintf( "e=%.3g", e.vals[[ pp ]] ),
                            cex.main=0.7 ), silent=T ) ##e.val=NA, mot.ind=pp
            dev.off() }, silent=T )
        }
      } ##)
      cat( "\n" )
      
      cat( "MOTIF POSITIONS: " )
      for ( k in ks ) {
      ##mc$
      ##lapply( ks, function( k, ... ) {
        if ( k %% 25 == 0 ) cat( k ) else cat( "." )
        if ( file.exists( sprintf( "%s/htmls/cluster%04d_mot_posns.png", out.dir, k ) ) ) return() ##next
        try( {
          png( sprintf( "%s/htmls/cluster%04d_mot_posns.png", out.dir, k ), width=128, height=12+6*length(get.rows(k)),
              antialias="subpixel" )
          par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
          c <- plotClust( k, dont.plot=T, ... ) ##seq.type=seq.type, 
          try( plotClusterMotifPositions( c, cex=0.4, no.key=T, ... ) )
          dev.off() }, silent=F )
      } ##)
      cat( "\n" )
    }
  }
  
  if ( "main" %in% output ) {
    mc <- get.parallel( length( ks ) )
    cat( "WRITING MAIN HTML TABLE..." )
    require( hwriter ) ## see http://www.ebi.ac.uk/~gpau/hwriter/
    dlf( paste( out.dir, "hwriter.css", sep="/" ), "http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css" )
    dlf( paste( out.dir, "sorttable.js", sep="/" ), "http://www.kryogenix.org/code/browser/sorttable/sorttable.js" )
    cat( "..." )
    ##cluster.summ <- cluster.summary( ... )
    ##if ( nrow( cluster.summ ) <= 0 )
    cluster.summ <- cluster.summary( e.cutoff=NA, nrow.cutoff=2 )
    write.table( cluster.summ, file=paste( out.dir, "/cluster.summary.tsv", sep="" ), quote=F, sep="\t" )
    cat( "..." )
    html <- openPage( paste( out.dir, "/index.html", sep="" ), link.javascript="sorttable.js",
                     title=paste( "cMonkey bicluster summary for run", cmonkey.filename ), link.css='hwriter.css' )
    hwrite( paste( "<h2>cMonkey bicluster summary for run", cmonkey.filename, "</h2>" ), html )
    hwrite( '<ul><li>Download a tab-delimited version of this table', html, link="cluster.summary.tsv",
           style='font-size:75%' )
    hwrite( "<li>Download a list of each bicluster's gene members", html, link="cluster.members.genes.txt",
           style='font-size:75%' )
    hwrite( "<li>Download a list of each bicluster's array/condition members", html,
           link="cluster.members.arrays.txt", style='font-size:75%' )
    cat( "..." )
    ##hwrite( "<li>Network summary of biclusters", html, link="svgs/main.svg", style='font-size:75%' )
    ##if ( "svg" %in% output )
    hwrite( "<li>Plots of summary statistics of biclustering run", html,
                                    link="svgs/stats.svg", style='font-size:75%' )
    ##if ( "rdata" %in% output )
    hwrite( "<li>Saved cMonkey R session file", html, link="cm_session.RData",
                                      style='font-size:75%' )
    hwrite( "<li>Summary of cMonkey input parameters</ul>", html, link="cm.params.txt", style='font-size:75%' )
    hwrite( "<br><center><b>Bicluster summary</b></center><br>", html )
    hwrite( "<br><center><b>Sort the table by a given column by clicking on the column's header.<br>Click on bicluster link in first column for more info.</b></center><br>", html,
           style='font-size:60%' )
    cat( "..." )
    himg0 <- hwriteImage( sprintf( "htmls/cluster%04d_profile.png", as.integer( rownames( cluster.summ ) ) ),
                         table=F )
    himg0 <- hwrite( paste( himg0, sprintf( "Residual = %.3f", cluster.summ$resid ), sep="<br>" ), center=TRUE,
                    table=F )
    himg0a <- hwriteImage( sprintf( "htmls/cluster%04d_network.png", as.integer( rownames( cluster.summ ) ) ),
                          table=F )
    if ( ! is.null( seq.type ) ) {
      e.val.1 <- lapply( meme.scores[[ seq.type ]][ as.integer( rownames( cluster.summ ) ) ],
                      function( i ) { mo <- i$meme.out; if ( length( mo ) >= 1 ) mo[[ 1 ]]$e.value else Inf } ) ##i$meme.out[[ 1 ]]$e.value )
      for ( i in 1:length( e.val.1 ) ) if ( is.null( e.val.1[[ i ]] ) ) e.val.1[[ i ]] <- NA
      himg1 <- hwriteImage( sprintf( "htmls/cluster%04d_pssm1.png", as.integer( rownames( cluster.summ ) ) ),
                           table=F, title=sprintf( "E-val = %.3g", unlist( e.val.1 ) ) )
      ##himg1a <- hwrite( as.character( cluster.summ$consensus1 ), table=T )
      himg1 <- hwrite( paste( himg1, as.character( cluster.summ$consensus1 ), sep="<br>" ), center=TRUE, table=F )
      if ( ! is.null( seq.type ) )
        e.val.2 <- lapply( meme.scores[[ seq.type ]][ as.integer( rownames( cluster.summ ) ) ],
                          function( i ) { mo <- i$meme.out; if ( length( mo ) > 1 ) mo[[ 2 ]]$e.value else Inf } ) ##i$meme.out[[ 2 ]]$e.value )
      else e.val.2 <- as.list( rep( NA, k.clust ) )
      for ( i in 1:length( e.val.2 ) ) if ( is.null( e.val.2[[ i ]] ) ) e.val.2[[ i ]] <- NA
      himg2 <- hwriteImage( sprintf( "htmls/cluster%04d_pssm2.png", as.integer( rownames( cluster.summ ) ) ),
                           table=F, title=sprintf( "E-val = %.3g", unlist( e.val.2 ) ) )
      himg2 <- hwrite( paste( himg2, as.character( cluster.summ$consensus2 ), sep="<br>" ), center=TRUE, table=F )
      himg2a <- hwriteImage( sprintf( "htmls/cluster%04d_mot_posns.png", as.integer( rownames( cluster.summ ) ) ),
                            table=F )
      e.val.1[ is.na( e.val.1 ) ] <- 9e9; e.val.2[ is.na( e.val.2 ) ] <- 9e9
    } else {
      e.val.1 <- e.val.2 <- as.list( rep( NA, k.clust ) )
      himg1 <- himg2 <- himg2a <- NULL
    }
    cluster.summ$score <- sprintf( "%.3f", cluster.summ$score )
    ##cluster.summ$resid <- sprintf( "%.3f", cluster.summ$resid ) 
    rn <- rownames( cluster.summ )
    cat( "..." )
    cluster.summ.orig <- cluster.summ
    cluster.summ <- cbind( bicluster=cluster.summ$k, n.genes=cluster.summ$nrow,
                          n.arrays=sapply( as.integer( rownames( cluster.summ ) ),
                            function( i ) length( get.cols( i ) ) ), score=cluster.summ$score,
                          residual=sprintf( "%.3f", cluster.summ$resid ) )
    if ( "score.norm" %in% colnames( cluster.summ.orig ) )
      cluster.summ <- cbind( cluster.summ, score.norm=sprintf( "%.3f", cluster.summ.orig$score.norm ) ) ##
    rownames( cluster.summ ) <- rn
    rows <- list(); for ( k in as.integer( rn ) ) rows[[ k ]] <- sort( get.rows( k ) )
    himg3 <- hwrite( sapply( as.integer( rn ), function( k ) paste( rows[[ k ]], collapse=" " ) ), table=F )
    cat( "...\n" )
    if ( ! no.genome.info ) {
      himg4 <- hwrite( unlist( mc$apply( as.integer( rn ), function( k ) {
        if ( k %% 25 == 0 ) cat( k ) else cat( "." ); if ( length( rows[[ k ]] ) <= 0 ) return();
        tmp <- get.long.names( rows[[ k ]], short=T ); tmp <- unique( tmp[ ! tmp %in% rows[[ k ]] & tmp != "" ] )
        paste( tmp, collapse=" " ) } ) ), table=F ); cat( "\n" )
      himg5 <- hwrite( unlist( mc$apply( as.integer( rn ), function( k ) {
        if ( k %% 25 == 0 ) cat( k ) else cat( "." ); if ( length( rows[[ k ]] ) <= 0 ) return();
        tmp <- get.long.names( rows[[ k ]], short=F ); tmp <- unique( tmp[ ! tmp %in% rows[[ k ]] & tmp != "" ] )
        paste( tmp, collapse=" | " ) } ) ), table=F ); cat( "\n" )
    } else {
      himg4 <- himg5 <- NULL
    }
    ## Let's make the table sortable using code from http://www.kryogenix.org/code/browser/sorttable/
    nas <- rep( NA, nrow( cluster.summ ) )
    hwrite( cbind( cluster.summ[ ,1:min( ncol( cluster.summ ), 6 ) ], profile=himg0, network=himg0a,
                  motif1=himg1, motif2=himg2, motif.posns=himg2a, probe.names=himg3, short.names=himg4, 
                  long.names=himg5 ), html, row.names=F, table.style='text-align:center;font-size:70%;font-family:Arial', table.class='sortable',
           row.style=list( 'font-weight:bold;text-align:center;font-size:70' ), 
           col.style=list( probe.names='font-size:70%', orf.names='font-size:50%', short.names='font-size:50%', long.names='font-size:50%', motif1='font-size:50%', motif2='font-size:50%' ),
           col.sorttable_customkey=list( residual=sprintf( "%.3f", cluster.summ.orig$residual ),
             score.norm=if ( "score.norm" %in% colnames( cluster.summ.orig ) )
             sprintf( "%.3f", cluster.summ.orig$score.norm ) else NULL,
             profile=sprintf( "%.3f", cluster.summ.orig$resid ),
             motif1=sprintf( "%.30f", unlist( e.val.1 ) ), e.val1=sprintf( "%.30f", unlist( e.val.1 ) ),
             motif2=sprintf( "%.30f", unlist( e.val.2 ) ), e.val2=sprintf( "%.30f", unlist( e.val.2 ) ) ),
           col.class=list( network=c( "sorttable_nosort", nas ), ##motif1=c( "sorttable_nosort", nas ), motif2=c( "sorttable_nosort", nas ),
             motif.posns=c( "sorttable_nosort", nas ) ),
           col.link=list( sprintf( "htmls/cluster%04d.html", as.integer( rownames( cluster.summ ) ) ) ) )
    closePage( html, splash=F )
    ##}
    
    for ( i in sapply( 1:k.clust, function( k ) c( k, sort( get.rows( k ) ) ) ) )
      cat( i, "\n", file=paste( out.dir, "/cluster.members.genes.txt", sep="" ), append=T )
    for ( i in sapply( 1:k.clust, function( k ) c( k, sort( get.cols( k ) ) ) ) )
      cat( i, "\n", file=paste( out.dir, "/cluster.members.arrays.txt", sep="" ), append=T )

    tmp <- capture.output( for ( name in ls( cmonkey.params ) ) {
      cat( name, "= " ); str( get( name, envir=cmonkey.params ), no.list=T ) } ) 
    cat( tmp, file=paste( out.dir, "/cm.params.txt", sep="" ), sep="\n", collapse="\n" )

    ## This is to make it display from the ISB central web server (but set xmlHeader=T works instead)
    ##if ( file.exists( "/sw/bin/rpl" ) || file.exists( "/usr/bin/rpl" ) )
    ##  system( paste( "rpl '<svg version' '<?xml version=\"1.0\"?><svg version' ", out.dir, "/svgs/*.svg" ) )
  }
    
  ## TODO: compress the svg's to svgz's (using gzip) and replace '.svg' with '.svgz' in all htmls.
  ## Optional because some web servers (ahem, Microsoft?) can't handle svgz's
  ## TODO: use pdftk to compress the pdfs, too
  if ( gzip ) { ##&& length( system( "which rpl", intern=T ) ) > 0 ) {
    rpl <- function( find, replace, file, ... ) {
      f <- readLines( file ); f <- gsub( find, replace, f, ... ); writeLines( f, con=file ) }
    ##if ( "svg" %in% output ) {
    system( sprintf( "gzip -v %s/svgs/*.svg", out.dir ) )
    for ( f in list.files( paste( out.dir, "/svgs", sep="" ), full=T ) ) if ( grepl( ".svg.gz", f, fixed=T ) )
      system( sprintf( "mv -v %s %s", f, sub( ".svg.gz", ".svgz", f, fixed=T ) ) )
    ##system( sprintf( "rpl .svg .svgz %s/*.html %s/htmls/*.html", out.dir, out.dir ) )
    ##for ( f in
    ##mc$
    lapply( c( list.files( sprintf( "%s/htmls", out.dir ), pattern=glob2rx( "*.html" ), full=T ),
                list.files( out.dir, pattern=glob2rx( "*.html" ), full=T ) ), function( f ) {
                  cat( f, "\n" ); rpl( '.svg"', '.svgz"', f, fixed=T ) } )
    ##for ( f in list.files( out.dir, pattern=glob2rx( "*.html" ), full=T ) ) {
    ##  cat( f, "\n" ); rpl( '.svg"', '.svgz"', f, fixed=T ) }
    ##}
  }
  if ( "rdata" %in% output ) save.cmonkey.env( file=paste( out.dir, "/cm_session.RData", sep="" ) )
  out.dir
}

write.bicluster.network <- function( out.dir=NULL, ks=1:k.clust, seq.type=names( mot.weights )[ 1 ], tomtom=T,
                                    tt.filter=function( tt ) subset( tt, overlap >= 4 & q.value <= 0.05 ),
                                    m.filter=function( m ) subset( m, e.value <= Inf ),
                                    gene.url=function( g ) sprintf( "http://microbesonline.org/cgi-bin/keywordSearch.cgi?taxId=%d&keyword=%s", taxon.id, g ),
                                    image.urls=T,
                                    ... ) {
  if ( is.null( out.dir ) ) {
    out.dir <- paste( cmonkey.filename, "network", sep="/" )
    if ( iter != n.iter ) out.dir <- sprintf( "%s_%04d/network", cmonkey.filename, iter )
  }
  if ( ! file.exists( out.dir ) ) dir.create( out.dir, recursive=T, showWarnings=F )
  cat( "Outputing to", out.dir, "\n" )

  ## Make r.sif: network of genes to bicluster nodes
  r.sif <- do.call( rbind, lapply( ks, function( k ) data.frame( get.rows( k ), sprintf( "bicluster_%04d", k ) ) ) )
  r.sif <- data.frame( r.sif[ ,1 ], "gene_member", r.sif[ ,2 ] ); colnames( r.sif ) <- c( "node1", "int", "node2" ) 

  ## Make c.sif: network of conditions to bicluster nodes
  ##c.sif <- do.call( rbind, lapply( ks, function( k ) data.frame( get.cols( k ), as.character( k ) ) ) )
  ##c.sif <- data.frame( c.sif[ ,1 ], "cond_member", c.sif[ ,2 ] ); colnames( c.sif ) <- c( "node1", "int", "node2" )

  ## Make m.sif: network of motif nodes to bicluster nodes. Include eda's with motif stats
  ms <- meme.scores[[ seq.type ]]
  m.sif <- do.call( rbind, lapply( ks, function( k ) if ( length( ms[[ k ]]$meme.out ) <= 0 ) NULL else
    data.frame( sprintf( "motif_%04d_%d", k, 1:length( ms[[ k ]]$meme.out ) ),
               sprintf( "bicluster_%04d", k ), sapply( ms[[ k ]]$meme.out, "[[", "width" ),
               sapply( ms[[ k ]]$meme.out, "[[", "sites" ), sapply( ms[[ k ]]$meme.out, "[[", "llr" ),
               sapply( ms[[ k ]]$meme.out, "[[", "e.value" ) ) ) )
  m.sif <- data.frame( m.sif[ ,1 ], "motif", m.sif[ ,2:ncol( m.sif ) ] )
  colnames( m.sif ) <- c( "node1", "int", "node2", "width", "nsites", "llr", "e.value" ) ## include eda's
  if ( ! is.null( m.filter ) ) m.sif <- m.filter( m.sif )
  
  ## Run tomtom to get motif similarities, add edges between motif nodes; include eda's with stats of motif match
  
  ## Make r.sif, m.sif and tt.sif have same number of columns (eda's from m.sif and tt.sif) 
##   for ( i in 1:( ncol( m.sif ) + ncol( tt.sif ) - 6 ) ) r.sif <- cbind( r.sif, rep( NA, nrow( r.sif ) ) )
##   colnames( r.sif ) <- c( colnames( r.sif )[ 1:3 ], colnames( m.sif )[ 4:ncol( m.sif ) ],
##                          colnames( tt.sif )[ 4:ncol( tt.sif ) ] )
##   nc <- ncol( m.sif )
##   for ( i in 1:( ncol( tt.sif ) - 3 ) ) m.sif <- cbind( m.sif, rep( NA, nrow( m.sif ) ) )
##   colnames( m.sif ) <- c( colnames( m.sif )[ 1:nc ], colnames( tt.sif )[ 4:ncol( tt.sif ) ] )
##   tmp <- tt.sif[ ,1:3 ]; for ( i in 1:( nc - 3 ) ) tmp <- cbind( tmp, rep( NA, nrow( tmp ) ) )
##   tmp <- cbind( tmp, tt.sif[ ,4:ncol( tt.sif ) ] ); colnames( tmp ) <- colnames( m.sif )
##   tt.sif <- tmp; rm( tmp )

  out.sif <- rbind( r.sif[ ,1:3 ], m.sif[ ,1:3 ], tt.sif[ ,1:3 ] )
  m.eda <- data.frame( edge=paste( m.sif[ ,1 ], " (", m.sif[ ,2 ], ") ", m.sif[ ,3 ], sep="" ),
                     m.sif[ ,4:ncol( m.sif ) ] )

  node.type <- as.data.frame( rbind( as.matrix( data.frame( attr( ratios, "rnames" ), "gene" ) ),
                                    as.matrix( data.frame( attr( ratios, "cnames" ), "condition" ) ),
                                    as.matrix( data.frame( sprintf( "bicluster_%04d", ks ), "bicluster" ) ),
                                    as.matrix( data.frame( as.character( m.sif$node1 ), "motif" ) ) ) )
  colnames( node.type ) <- c( "node", "type" )
  syn.names <- do.call( rbind, lapply( attr( ratios, "rnames" ),
             function( g ) data.frame( g, paste( get.synonyms( g )[[ g ]], collapse="::" ) ) ) )
  colnames( syn.names ) <- c( "gene", "synonyms" )
  l.names <- do.call( rbind, lapply( attr( ratios, "rnames" ),
             function( g ) data.frame( g, paste( get.long.names( g )[[ g ]], collapse="::" ) ) ) )
  colnames( l.names ) <- c( "gene", "long.name" )
  s.names <- do.call( rbind, lapply( attr( ratios, "rnames" ),
             function( g ) data.frame( g, paste( get.long.names( g, short=T )[[ g ]], collapse="::" ) ) ) )
  colnames( s.names ) <- c( "gene", "short.name" )
  g.url <- NULL
  if ( ! is.null( gene.url ) ) {
    g.url <- do.call( rbind, lapply( attr( ratios, "rnames" ), function( g ) data.frame( g, gene.url( g ) ) ) )
    colnames( g.url ) <- c( "gene", "url" )
  }
  mot.info <- do.call( rbind, lapply( ks, function( k ) if ( length( ms[[ k ]]$meme ) <= 0 ) NULL else
                                     data.frame( sprintf( "motif_%04d_%d", k, 1:length( ms[[ k ]]$meme.out ) ),
                                                sapply( ms[[ k ]]$meme.out, function( i ) pssm.to.string( i$pssm ) ),
                                                sapply( ms[[ k ]]$meme.out, "[[", "width" ),
                                                sapply( ms[[ k ]]$meme.out, "[[", "sites" ),
                                                sapply( ms[[ k ]]$meme.out, "[[", "llr" ),
                                                sapply( ms[[ k ]]$meme.out, "[[", "e.value" ),
                                                sprintf( "file://%s/%s/htmls/cluster%04d_pssm%d.png",
                                                        getwd(), out.dir, k, 1:length( ms[[ k ]]$meme.out ) )
                                                ) ) )
  colnames( mot.info ) <- c( "motif", "consensus", "width", "n.sites", "llr", "e.value", "imgURL" )
  clust.info <- lapply( c( "resid", "p.clust", "e.val", "nrows", "ncols" ),
                       function( i ) sapply( ks, function( j ) clusterStack[[ j ]][[ i ]] ) )
  names( clust.info ) <- c( "resid", "p.clust", "e.val", "nrows", "ncols" )
  clust.info[[ names( unlist( sapply( clust.info, nrow ) ) ) ]] <-
    t( clust.info[[ names( unlist( sapply( clust.info, nrow ) ) ) ]] )
  for ( n in names( clust.info ) ) if ( ! is.null( ncol( clust.info[[ n ]] ) ) &&
                                       is.null( colnames( clust.info[[ n ]] ) ) )
    colnames( clust.info[[ n ]] ) <- paste( n, 1:ncol( clust.info[[ n ]] ), sep="." )
  clust.info <- cbind( bicluster=sprintf( "bicluster_%04d", ks ), do.call( cbind, clust.info ),
                      url=sprintf( "file://%s/%s/htmls/cluster%04d.html", getwd(), out.dir, ks ),
                      imgURL=sprintf( "file://%s/%s/htmls/cluster%04d_profile.png", getwd(), out.dir, ks ) )
  rownames( clust.info ) <- NULL
  ## NOW MERGE ALL NOAs INTO ONE TAB-DELIMITED TABLE FOR LOADING INTO CYTOSCAPE
  noa <- merge( node.type, syn.names, by.x="node", by.y="gene", all=T )
  noa <- merge( noa, l.names, by.x="node", by.y="gene", all=T )
  noa <- merge( noa, s.names, by.x="node", by.y="gene", all=T )
  noa <- merge( noa, mot.info, by.x="node", by.y="motif", all=T )
  noa <- merge( noa, clust.info, by.x="node", by.y="bicluster", all=T )
  if ( ! is.null( g.url ) ) noa <- merge( noa, g.url, by.x="node", by.y="gene", all=T )
  noa$imgURL <- as.character( noa$imgURL.x )
  noa$imgURL[ ! is.na( noa$imgURL.y ) ] <- as.character( noa$imgURL.y[ ! is.na( noa$imgURL.y ) ] )
  noa <- noa[ , ! colnames( noa ) %in% c( "imgURL.x", "imgURL.y" ) ]
  write.table( out.sif, quote=F, sep="\t", col.names=T, row.names=F, file=paste( out.dir, "all.sif", sep="/" ) )
  ##write.table( c.sif, quote=F, sep="\t", col.names=T, row.names=F, file=paste( out.dir, "c.sif", sep="/" ) )
  write.table( noa, quote=F, sep="\t", col.names=T, row.names=F, file=paste( out.dir, "all.noa", sep="/" ) )
  write.table( m.eda, quote=F, sep="\t", col.names=T, row.names=F, file=paste( out.dir, "m.eda", sep="/" ) )
  write.table( tt.eda, quote=F, sep="\t", col.names=T, row.names=F, file=paste( out.dir, "tt.eda", sep="/" ) )
  cat( "Wrote", nrow( noa ), "nodes and", nrow( out.sif ), "edges to", out.dir, "\n" )
  invisible( list( sif=out.sif, noa=noa, m.eda=m.eda, tt.eda=tt.eda ) )
}
