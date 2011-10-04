get.biclusters <- function( genes, conditions, motifs, tfs ) {
  if ( ! missing( genes ) ) {
    out <- lapply( genes, function( g ) { out <- which( sapply( e$clusterStack, function( i ) g %in% i$rows ) )
                                          paste( "BIC", out, sep="_" ) } )
    names( out ) <- genes
  }

  if ( ! missing( conditions ) ) {
    out <- lapply( conditions, function( cc ) { out <- which( sapply( e$clusterStack, function( i ) cc %in% i$cols ) )
                                          paste( "BIC", out, sep="_" ) } )
    names( out ) <- conditions
  }

  if ( ! missing( motifs ) ) {
    out <- lapply( motifs, function( m ) paste( "BIC", sapply( strsplit( m, "_" ), "[", 2 ), sep="_" ) )
    names( out ) <- motifs
  }

  if ( ! missing( tfs ) ) {
    out <- lapply( tfs, function( tf ) lapply( e.coeffs, function( k ) {
      tmp <- k[[ 1 ]]$coeffs; tmp[ grep( tf, names( tmp ) ) ] } ) )
    for ( i in 1:length( out ) ) {
      out[[ i ]] <- out[[ i ]][ sapply( out[[ i ]], length ) > 0 ]
      names( out[[ i ]] ) <- paste( "BIC", as.integer( names( out[[ i ]] ) ), sep="_" )
    }
    names( out ) <- tfs
  }

  out
}

get.genes <- function( biclusters, motifs ) {
  if ( ! missing( biclusters ) ) {
    inds <- as.integer( sapply( strsplit( biclusters, "_" ), "[", 2 ) )
    out <- lapply( inds, function( i ) e$clusterStack[[ i ]]$rows )
    names( out ) <- biclusters
  } else if ( ! missing( motifs ) ) {
    out <- lapply( motifs, function( m ) {
      tmp <- get.motif.info( motif=m )[[ 1 ]]
      if ( is.null( tmp ) ) return( NULL )
      if ( 'pvals' %in% colnames( tmp$mast ) ) unique( as.character( subset( tmp$mast, pvals <= 0.05 &
                                   abs( mots ) == as.integer( strsplit( m, "_" )[[ 1 ]][ 3 ] ), gene, drop=T ) ) )
      else unique( as.character( subset( tmp$mast, p.value <= 0.05, gene, drop=T ) ) )
    } )
    names( out ) <- motifs
  }
  out
}

get.conditions <- function( biclusters ) {
  if ( ! missing( biclusters ) ) {
    inds <- as.integer( sapply( strsplit( biclusters, "_" ), "[", 2 ) )
    out <- lapply( inds, function( i ) e$clusterStack[[ i ]]$cols )
    names( out ) <- biclusters
  }
  out
}


get.motifs <- function( biclusters, motif.clusters, genes, positions, window=10 ) {
  if ( ! missing( biclusters ) ) {
    out <- lapply( biclusters, function( bic ) paste( gsub( "BIC_", "MOT_", bic ),
                       1:length( e$clusterStack[[ as.integer( gsub( "BIC_", "", bic ) ) ]]$e.val ), sep="_" ) )
    names( out ) <- biclusters
  }

  if ( ! missing( motif.clusters ) ) {
    tmp <- as.integer( sapply( strsplit( motif.clusters, "_" ), "[", 2 ) )
    out <- lapply( tmp, function( i ) paste( "MOT", attr( tt.out2[[ i ]], "mot.names" ), sep="_" ) )
    names( out ) <- motif.clusters
  }

  if ( ! missing( genes ) ) {
    out <- lapply( genes, function( g ) {
      unlist( lapply( paste( "MOT", c( paste( 1:e$k.clust, 1, sep="_" ),
                                      paste( 1:e$k.clust, 2, sep="_" ) ), sep="_" ), function( m ) {
        tmp <- get.motif.info( motif=m )[[ 1 ]]
        if ( is.null( tmp ) ) return( NULL )
        if ( 'pvals' %in% colnames( tmp$mast ) ) tmp <- unique( as.character( subset( tmp$mast, pvals <= 0.05 &
                                   abs( mots ) == as.integer( strsplit( m, "_" )[[ 1 ]][ 3 ] ), gene, drop=T ) ) )
        else tmp <- unique( as.character( subset( tmp$mast, p.value <= 0.05, gene, drop=T ) ) )
        if ( g %in% tmp ) return( m )
        return( NULL )
      } ) )
    } )
    names( out ) <- genes
  }

  if ( ! missing( positions ) ) {
    ## Positions is a 2-vector (start,stop) w/ names(positions)[1]=Chr, e.g.
    ##  positions=c(Chr=10000,11000)
    ##  assume default is Chr and if positions is 1-vector then it's position +/- 5
    if ( ! exists( "pssm.scans" ) ) pssm.scans <- get.motif.positions( NULL )
    if ( positions[ 1 ] == "choose" ) { ## get it from a plot (probably output from plot.promoter.architecture())
      pn <- names( positions )
      positions <- locator( 1, 'p' )$x
      names( positions ) <- pn
    } else if ( grepl( '[, ]', positions[ 1 ], perl=T ) ) {
      positions <- strsplit( positions, '[, ]', perl=T )[[ 1 ]]
    }
    if ( length( positions ) == 3 ) {
      tmp <- as.integer( positions[ 2:3 ] ); names( tmp )[ 1 ] <- positions[ 1 ]; positions <- tmp
    } else if ( is.character( positions ) && length( positions ) == 2 ) {
      positions <- as.integer( positions )
    }
    if ( is.null( names( positions )[ 1 ] ) ) names( positions )[ 1 ] <- 'Chr'
    if ( length( positions ) == 1 ) {
      tmp <- c( positions[1] - window, positions[1] + window ); names( tmp )[ 1 ] <- names( positions )[ 1 ]
      positions <- tmp }
    print( positions )
    chr <- names( positions )[ 1 ]
    ##require( data.table )
    ##scans <- pssm.scans[ gene == chr & posns %between% positions, ]
    scans <- subset( pssm.scans, gene == chr & posns %betw% positions )
    out <- paste( "MOT", scans$bic, scans$mots, sep="_" )
  }
  
  out
}

get.motif.clusters <- function( motifs ) {
  if ( ! missing( motifs ) ) {
    tmp <- get.motifs( motif.clust=paste( "MOTC", 1:length( tt.out2 ), sep="_" ) )
    out <- lapply( motifs, function( m ) names( which( sapply( tmp, function( i ) m %in% i ) ) ) )
    names( out ) <- motifs
  }    
  
  out
}

get.expression <- function( genes, biclusters ) {
  if ( ! missing( genes ) ) {
    out <- lapply( genes, function( g ) e$ratios$ratios[ g, ] )
    names( out ) <- genes
  }
  
  if ( ! missing( biclusters ) ) {
    out <- lapply( biclusters, function( bic ){print(bic);
                  e$ratios$ratios[ get.genes( biclusters=bic )[[ 1 ]], get.conditions( biclusters=bic )[[ 1 ]] ]} )
    names( out ) <- biclusters
  }

  out
}

get.expressions <- get.expression ## for agglom func

get.motif.info <- function( motifs ) {
  if ( ! missing( motifs ) ) {
    inds <- strsplit( motifs, "_" )
    out <- lapply( inds, function( i ) {
      ms <- e$meme.scores[[ 1 ]][[ as.integer( i[ 2 ] ) ]]
      if ( is.null( ms$meme.out ) || length( ms$meme.out ) < as.integer( i[ 3 ] ) ) return( NULL )
      out <- ms$meme.out[[ as.integer( i[ 3 ] ) ]]
      tmp <- ms$pv.ev[[ 2 ]]
      if ( ! is.null( tmp ) ) tmp <- subset( tmp, abs( mots ) == as.integer( i[ 3 ] ) )
      else tmp <- out$posns
      out$mast <- tmp
      out
    } )
    names( out ) <- motifs
  }
  
  out
}

get.motif.cluster.info <- function( motif.clusters ) {
  if ( ! missing( motif.clusters ) ) {
    out <- tt.out2[ as.integer( gsub( "MOTC_", "", motif.clusters ) ) ]
    names( out ) <- motif.clusters
  }
  out
}

## default: agglomerate all genes/conditions/tfs that are in all biclusters containing VNG0826C
##  source is a gene ID / condition name / motif.cluster / tf
##  srcType is 'gene' / 'condition' / 'motif.cluster' / 'tf'
##  targetType is 'gene' or 'condition' or 'motif.cluster' or 'tf'
##  path (to take from source to target) is 'bicluster', 'bicluster,motif', etc. <-- the latter can be used
##             to get motif.clusters for a given gene, for example

## examples:
## egrin2.agglom() -- the default, counts of genes co-expressed with VNG0826C
## egrin2.agglom( target='condition' ) -- counts of conditions in biclusters containing VNG0826C
## egrin2.agglom( target='tf' ) -- counts of TFs regulating biclusters containing VNG0826C
## egrin2.agglom( target='motif.cluster', path='bicluster,motif' ) -- counts of motif clusters in biclusters containing
##                  VNG0826C -- agglomeration is cumulative, first gene->bicluster->motif->motif.cluster
##     NOTE - p-value calculation does not work right now for paths that don't end in 'bicluster', so try this:
## egrin2.agglom( srcType='motif.cluster', src='MOTC_5', path='motif,bicluster' ) -- counts of genes in biclusters
##                  containing a motif that is in motif cluster MOTC_5
## egrin2.agglom( srcType='motif.cluster', src='MOTC_5', target='condition', path='motif,bicluster' ) -- counts of
##                  conditions in biclusters containing a motif that is in motif.cluster MOTC_5
## egrin2.agglom( src='circadian_drk8_cycling_1440min_vs_NRC-1e.sig', srcType='condition', target='condition' ) --
##                  counts of conditions co-occurring with this cond (should be normalized by total.cond.count)
## egrin2.agglom( src=c(Chr=716450,716550), srcType="position", target="motif.cluster", path="motif" )
##                  counts of motif.clusters that have motif instances that overlap the given chromosome coords
##                  NOTE: can use src='choose' to choose a position by clicking on a plot; see get.motifs()
## NEW: if 'path' contains 'bicluster', then 'cond.filter' can be used to filter out only biclusters that
##      contain all of the given conditions (in cond.filter); e.g. cond.filter=colnames(all.ratios)[1:2]
##      if 'p.val'==TRUE compute hypergeometric p-values
agglom <- function( src='MOTC_1', srcType='motif.cluster', targetType='gene', path='motif,bicluster', p.val=T, 
                   q.val="BH", cond.filter=NULL, cond.filter.frac=0.5, verbose=F ) {
  by <- unlist( strsplit( path, ',' ) )
  orig.src <- src; orig.srcType <- srcType
  for ( b in by ) {
    text <- paste( 'bys <- unlist( get.', b,   's( ', srcType, "=src ) )", sep="" )
    if ( verbose ) print( text )
    eval( parse( text=text ) )
    if ( ! is.null( cond.filter ) && b == 'bicluster' ) {
      if ( verbose ) print( 'cc <- get.conditions( biclust=bys )' )
      cc <- get.conditions( biclust=bys )
      tmp <- sapply( cc, function( i ) mean( cond.filter %in% i ) >= cond.filter.frac ) ## 50% of cond.filter in the biclusters
      if ( ! any( tmp ) ) warning( "No biclusters pass the cond.filter! Perhaps too many conditions?" )
      else bys <- bys[ tmp ]
    }
    src <- bys
    srcType <- b
  }
  text <- paste( 'out <- get.', targetType, 's( ', by[ length( by ) ],      '=bys )', sep="" )
  if ( verbose ) print( text )
  eval( parse( text=text ) )
  if ( targetType == 'tf' ) out <- lapply( out, lapply, names )
  tab <- sort( table( unlist( out ) ), decreasing=T )
  out2 <- t( t( tab ) )
  colnames( out2 ) <- "count"
  tab2 <- tab[ tab > 1 ]
  pvals <- qvals <- numeric()
  if ( p.val && length( tab2 ) > 0 && by[ length( by ) ] == 'bicluster' ) { ## Right now this is only valid for path='..,..,bicluster'
    bics1 <- as.list( names( tab ) ); names( bics1 ) <- bics1
    if ( targetType != "bicluster" ) {
      text <- paste( 'bics1 <- get.biclusters( ', targetType,      '=names( tab2 ) )', sep="" )
      if ( verbose ) print( text )
      eval( parse( text=text ) )
    }
    n.bics <- e$k.clust
    if ( ! is.null( cond.filter ) ) {
      bics <- table( unlist( get.biclusters( cond=cond.filter ) ) )
      bics <- bics[ bics >= length( cond.filter ) * cond.filter.frac ] ##as.integer( gsub( "BIC_", "",
      bics1 <- lapply( bics1, function( i ) i[ i %in% names( bics ) ] )
      n.bics <- length( bics )
    }
    pvals <- phyper( tab2, length( out ), n.bics, sapply( bics1, length )[ names( tab2 ) ], lower=F )
    if ( ! is.na( q.val ) ) {
      if ( q.val == TRUE ) q.val <- "BH"
      pv <- c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) )
      qvals <- p.adjust( pv, q.val )
    }
  }
  if ( p.val ) {
    out2 <- cbind( out2, p.value=c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) ) )
    if ( ! is.na( q.val ) ) out2 <- cbind( out2, q.value=c( qvals, rep( 1, nrow( out2 ) - length( qvals ) ) ) )
  }
  out2 <- as.data.frame( out2 )
  attr( out2, "total" ) <- length( bys )
  out2
}

## Translate name of obj into type (NOTE doesn't catch TFs!)
get.type <- function( obj ) {
  prefix <- sapply( strsplit( obj, "_" ), "[", 1 )
  out <- ifelse( prefix == "BIC", "bicluster",
                ifelse( prefix == "MOT", "motif",
                       ifelse( prefix == "MOTC", "motif.cluster",
                              ifelse( obj %in% attr( e$ratios, "rnames" ), "gene", "condition" ) ) ) )
  out
}

## Agglomerate from two "sources", to biclusters, and compute the p-value of the bicluster overlap, e.g.
##   agglom2( "MOTC_1", "motif,bicluster", "MOTC_170", "motif,bicluster" )
agglom2 <- function( src="MOTC_1", srcType=get.type( src ), srcPath="motif,bicluster",
                    target="MOTC_170", targetType=get.type( target ), targetPath="motif,bicluster",
                    cond.filter=NULL, cond.filter.frac=0.5, verbose=F ) {
  by <- unlist( strsplit( srcPath, ',' ) )
  orig.src <- src; orig.srcType <- srcType
  for ( b in by ) {
    text <- paste( 'bys1 <- unlist( get.', b,   's( ', srcType, "=src ) )", sep="" )
    if ( verbose ) print( text )
    eval( parse( text=text ) )
    if ( ! is.null( cond.filter ) && b == 'bicluster' ) {
      if ( verbose ) print( 'cc <- get.conditions( biclust=bys1 )' )
      cc <- get.conditions( biclust=bys1 )
      tmp <- sapply( cc, function( i ) mean( cond.filter %in% i ) >= cond.filter.frac ) ## 50% of cond.filter in the biclusters
      if ( ! any( tmp ) ) warning( "No biclusters pass the cond.filter! Perhaps too many conditions?" )
      else bys1 <- bys1[ tmp ]
    }
    src <- bys1
    srcType <- b
  }

  by <- unlist( strsplit( targetPath, ',' ) )
  orig.target <- target; orig.targetType <- targetType
  for ( b in by ) {
    text <- paste( 'bys2 <- unlist( get.', b,   's( ', targetType, "=target ) )", sep="" )
    if ( verbose ) print( text )
    eval( parse( text=text ) )
    if ( ! is.null( cond.filter ) && b == 'bicluster' ) {
      if ( verbose ) print( 'cc <- get.conditions( biclust=bys2 )' )
      cc <- get.conditions( biclust=bys2 )
      tmp <- sapply( cc, function( i ) mean( cond.filter %in% i ) >= cond.filter.frac ) ## 50% of cond.filter in the biclusters
      if ( ! any( tmp ) ) warning( "No biclusters pass the cond.filter! Perhaps too many conditions?" )
      else bys2 <- bys2[ tmp ]
    }
    target <- bys2
    targetType <- b
  }

  out <- data.frame( src=orig.src, target=orig.target, nSrc=length( bys1 ), nTarget=length( bys2 ),
                    nBoth=sum( bys1 %in% bys2 ), p.value=phyper( sum( bys1 %in% bys2 ),
                                                   length( bys1 ), e$k.clust - length( bys1 ), length( bys2 ),
                                                   lower=F ) )
  out
}

## Compute full gene-centric network of genes,motif.clusters,TFs connected to each gene
##   optionally over a subset of genes and conditions
network <- function( genes, conditions, q.val=T ) {
  if ( missing( genes ) ) genes <- attr( e$ratios, "rnames" )
  if ( missing( conditions ) ) conditions <- NULL

  all.rm <- get.biclusters( genes ) ##lapply( genes, function( i ) { i <- e$row.membership[ i, ]; unique( i[ i != 0 ] ) } )
  names( all.rm ) <- genes ##rownames( e$row.membership )
  if ( ! is.null( conditions ) ) {
    all.cm <- get.biclusters( conditions ) ##lapply( conditions, function( i ) { i <- e$col.membership[ i, ]; unique( i[ i != 0 ] ) } )
    names( all.cm ) <- conditions ##rownames( e$col.membership )
  }
  
  require( multicore ); apply.func <- mclapply
  out <- apply.func( genes, function( g ) {
    print( g )
    tmp <- agglom( src=g, srcType='gene', targetType='gene', path='bicluster', q.val=q.val,
                         cond.filter=conditions )
    done <- genes[ 1:which( genes == g ) ]
    tmp <- subset( tmp, p.value < 1 & count > 2 & ! rownames( tmp ) %in% done )
    if ( nrow( tmp ) <= 0 ) return( NULL ) ##next
    out <- data.frame( gene1=g, int="gene_gene", gene2=rownames( tmp ), count=tmp[ ,1 ],
                                  p.value=tmp[ ,2 ], q.value=tmp[ ,3 ] ) ##)
    rownames( out ) <- NULL
    print( dim( out ) )
    out
  } )
  out <- do.call( rbind, out )

  out2 <- apply.func( 1:length( tt.out2 ), function( mc ) {
    mc <- paste( "MOTC", mc, sep="_" )
    print( mc )
    tmp <- agglom( src=mc, srcType='motif.cluster', targetType='gene', path='motif',
                         q.val=q.val, cond.filter=conditions )
    tmp <- subset( tmp, count > 2 ) ## p.value < 1 & 
    if ( nrow( tmp ) <= 0 ) return( NULL ) ##next
    out2 <- data.frame( gene1=mc, int="motc_gene", gene2=rownames( tmp ), count=tmp[ ,1 ],
                                  p.value=tmp[ ,2 ], q.value=tmp[ ,3 ] ) ##)
    rownames( out2 ) <- NULL
    print( dim( out2 ) )
    out2
  } )
  out2 <- do.call( rbind, out2 )

  noa <- data.frame()
  genes <- unique( c( as.character( subset( out, int == "gene_gene" )$gene1 ),
                     as.character( subset( out, int == "gene_gene" )$gene2 ) ) )
  noa <- rbind( noa, data.frame( name=genes, type="gene" ) )
  motcs <- unique( as.character( subset( out2, int == "motc_gene" )$gene1 ) )
  noa <- rbind( noa, data.frame( name=motcs, type="motif_cluster" ) )
  
  list( out=out, out2=out2, noa=noa )
}

get.motif.positions <- function( motifs="ALL", seqs=e$genome.info$all.upstream.seqs[[ 1 ]],
                                seq.type=names( e$meme.scores )[ 1 ], counts=T, verbose=F ) {
  if ( ! exists( "pssm.scans" ) ) pssm.scans <- get.pssm.scans( motifs, seqs, seq.type )
  p.scans <- pssm.scans
  p.scans$mots <- abs( p.scans$mots ); gc() ## Don't use strand or p-values for these funcs, so remove em
  p.scans <- p.scans[ ,c( "bic", "mots", "gene", "posns" ) ]; gc() ##, "fwd", "pvals" ) ]; gc()
  require( data.table )
  p.scans <- as.data.table( p.scans ); gc()
  setkey( p.scans, bic, mots, gene, posns )

  if ( motifs != "ALL" ) {
    mots <- strsplit( gsub( "MOT_", "", motifs ), "_" )
    bi <- as.integer( sapply( mots, "[", 1 ) )
    mo <- as.integer( sapply( mots, "[", 2 ) )
    scans <- p.scans[ J( bi, mo, names( seqs ) ) ]
  } else {
    scans <- p.scans
  }
  scans <- scans[ ! is.na( posns ) ]
  if ( ! counts ) return( invisible( scans ) )

  ms <- e$meme.scores[[ seq.type ]]

  counts <- integer()
  if ( nrow( scans ) > 0 ) {
    for ( i in 1:nrow( scans ) ) {
      if ( verbose && i %% 1000 == 0 ) cat( i, nrow( scans ), length( counts ), "\n" )
      k <- scans$bic[ i ]
      mot <- scans$mots[ i ]
      width <- ms[[ k ]]$meme.out[[ mot ]]$width
      if ( width > 0 ) {
        posn <- scans$posns[ i ]
        inds <- ( posn - 1 ):( posn - 2 + width )
        if ( length( counts ) < inds[ length( inds ) ] ) {
          counts[ inds[ length( inds ) ] ] <- 0; counts[ is.na( counts ) ] <- 0 }
        counts[ inds ] <- counts[ inds ] + 1
      }
    }
  }
  invisible( counts )
}

get.pssm.scans <- function( motifs="ALL", seqs=e$genome.info$all.upstream.seqs[[ 1 ]],
                           seq.type=names( e$meme.scores )[ 1 ], p.cutoff='0.0001' ) {
  mast.cmd <- './progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.0001 -ev 9999 -comp' ## Different from mast cmd used for cMonkey; e.g. doesn't have -revcomp, also set the p-value cutoff to 1e-4 (-mt 0.0001)
  mast.cmd <- sprintf( './progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt %s -ev 9999 -comp', p.cutoff )

  ks <- 1:e$k.clust
  if ( motifs != "ALL" && ! is.null( motifs ) )
    ks <- as.integer( gsub( "BIC_", "", unlist( get.biclusters( motif=motifs ) ) ) )
  out <- do.call( rbind, mclapply( ks, function( k ) {
    print( k )
    if ( is.null( e$meme.scores[[ seq.type ]][[ k ]] ) || e$meme.scores[[ seq.type ]][[ k ]] == "" ||
        is.null( e$meme.scores[[ seq.type ]][[ k ]]$meme.out ) ) return( NULL )
    tmp1 <- e$cluster.meme.motif.lines( k, logodds=T )
    if ( length( tmp1 ) < 10 ) return( NULL ) ## No motifs
    tmp3 <- e$runMast( tmp1, mast.cmd, names( seqs ), seqs, verbose=F, seq.type=seq.type,
                      bg.list=NULL, bg.fname=NULL, unlink=T )
    if ( length( tmp3 ) <= 0 ) return( NULL )
    tmp4 <- e$getMastPValuesAndEValues( tmp3, names( seqs ) )[[ 2 ]]
    if ( nrow( tmp4 ) <= 0 ) return( NULL )
    ##print( dim( tmp4 ) )
    tmp4 <- cbind( tmp4, bic=rep( k, nrow( tmp4 ) ) )
    tmp4
  } ) )

  out  
}
