require( data.table ) ## %chin% is ~5x fas

get.biclusters <- function( genes, conditions, motifs, tfs, tf.weight.cutoff=0, ... ) {
  if ( ! missing( genes ) ) {
    if ( exists( "genes.to.biclust" ) ) {
      tmp <- genes.to.biclust[ genes ]
      out <- lapply( tmp, function( i ) paste( "BIC", i, sep="_" ) )
    } else {
      out <- lapply( genes, function( g ) {
        out <- which( unlist( lapply( e$clusterStack, function( i ) g %chin% i$rows ) ) )
        paste( "BIC", out, sep="_" ) } )
    }
    names( out ) <- genes
  }

  if ( ! missing( conditions ) ) {
    if ( exists( "conds.to.biclust" ) ) {
      tmp <- conds.to.biclust[ conditions ]
      out <- lapply( tmp, function( i ) paste( "BIC", i, sep="_" ) )
    } else {
      out <- lapply( conditions, function( cc ) {
        out <- which( unlist( lapply( e$clusterStack, function( i ) cc %chin% i$cols ) ) )
        paste( "BIC", out, sep="_" ) } )
    }
    names( out ) <- conditions
  }

  if ( ! missing( motifs ) ) {
    out <- lapply( motifs, function( m ) paste( "BIC", sapply( strsplit( m, "_" ), "[", 2 ), sep="_" ) )
    names( out ) <- motifs
  }

  if ( ! missing( tfs ) ) {
    if ( exists( "coeffs.to.biclust" ) && tf.weight.cutoff == 0 ) out <- lapply( coeffs.to.biclust[ tfs ],
                                                       function( i ) paste( "BIC", i, sep="_" ) )
    else {
      out <- lapply( tfs, function( tf ) lapply( e.coeffs, function( k ) { ##cat(class(k),k[[1]]$k,"\n");
        if ( is.null( k ) ) return( NULL ); tmp <- k[[ 1 ]]$coeffs; tmp <- tmp[ grep( tf, names( tmp ) ) ];
                                                                      tmp[ abs( tmp ) >= tf.weight.cutoff ] } ) )
      for ( i in 1:length( out ) ) {
        out[[ i ]] <- out[[ i ]][ sapply( out[[ i ]], length ) > 0 ]
        names( out[[ i ]] ) <- paste( "BIC", as.integer( names( out[[ i ]] ) ), sep="_" )
      }
    }
    names( out ) <- tfs
  }

  out
}

## TODO: get.genes from positions (i.e. if motif falls in gene's promoter region)
get.genes <- function( biclusters, motifs, ... ) {
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

get.conditions <- function( biclusters, ... ) {
  if ( ! missing( biclusters ) ) {
    inds <- as.integer( sapply( strsplit( biclusters, "_" ), "[", 2 ) )
    out <- lapply( inds, function( i ) e$clusterStack[[ i ]]$cols )
    names( out ) <- biclusters
  }
  out
}

get.motifs <- function( biclusters, motif.clusters, genes, positions, window=10,
                       expand.motif.cluster.clusters=T, op.shift=F, distance=e$motif.upstream.search[[ 1 ]], ... ) {
  if ( ! missing( biclusters ) ) {
    out <- lapply( biclusters, function( bic ) paste( gsub( "BIC_", "MOT_", bic ),
                       1:length( e$clusterStack[[ as.integer( gsub( "BIC_", "", bic ) ) ]]$e.val ), sep="_" ) )
    names( out ) <- biclusters
  }

  if ( ! missing( motif.clusters ) && exists( 'motif.clusts' ) ) {
    if ( expand.motif.cluster.clusters && exists( "motif.cluster.clusters" ) ) {
      mc.tmp <- as.integer( gsub( "MOTC_", "", motif.clusters ) )
      motif.clusters <- c( motif.clusters[ motif.cluster.clusters[ mc.tmp ] == 0 ],
                      paste( "MOTC", which( motif.cluster.clusters != 0 &
                                           motif.cluster.clusters %in% motif.cluster.clusters[ mc.tmp ] ), sep="_" ) )
      motif.clusters <- motif.clusters[ motif.clusters != "MOTC_" ]
    }
    tmp <- as.integer( sapply( strsplit( motif.clusters, "_" ), "[", 2 ) )
    ##out <- lapply( tmp, function( i ) paste( "MOT", attr( tt.out2[[ i ]], "mot.names" ), sep="_" ) )
    out <- motif.clusts[ tmp ]
    names( out ) <- motif.clusters
  }

  if ( ! missing( genes ) ) {
    out <- lapply( genes, function( g ) {
      if ( exists( 'genes.to.motifs' ) ) {
        return( genes.to.motifs[[ g ]] )
      } else if ( exists( 'meme.hits' ) ) {
        mh <- meme.hits[ J( g ) ]
        return( paste( 'MOT', mh$bic, mh$mot, sep='_' ) )
      } else {
        all.motifs <- unlist( lapply( 1:nrow( motif.widths ),
                                     function( i ) { ii <- which( motif.widths[ i, ] > 0 )
                                                     if ( length( ii ) <= 0 ) return( NULL )
                                                     paste( 'MOT', i, ii, sep='_' ) } ) ) ## This is faster than above!
        unlist( lapply( all.motifs, function( m ) { ##paste( "MOT", c( paste( 1:e$k.clust, 1, sep="_" ),
          ##                paste( 1:e$k.clust, 2, sep="_" ) ), sep="_" ), function( m ) {
          tmp <- get.motif.info( motif=m )[[ 1 ]]
          if ( is.null( tmp ) ) return( NULL )
          if ( 'pvals' %in% colnames( tmp$mast ) ) tmp <- unique( as.character( subset( tmp$mast, pvals <= 0.05 &
                                   abs( mots ) == as.integer( strsplit( m, "_" )[[ 1 ]][ 3 ] ), gene, drop=T ) ) )
          else tmp <- unique( as.character( subset( tmp$mast, p.value <= 0.05, gene, drop=T ) ) )
          if ( g %in% tmp ) return( m )
          return( NULL )
        } ) )
      }
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
    } else if ( guess.type( positions[ 1 ] ) == 'gene' ) { ## use gene's promoter region
      ##gene <- positions[ 1 ]
      posns <- matrix( NA, ncol=2 )
      for ( gene in positions ) {
        tmp <- e$get.sequences( gene, seq.type='upstream meme', op.shift=op.shift,
                               distance=distance, filter=F, ... )
        if ( is.null( tmp ) ) next
        posns <- rbind( posns, as.matrix( attr( tmp, 'start.stops' )[ 1:2 ] )[ 1, ] )
        rownames( posns )[ nrow( posns ) ] <- as.matrix( attr( tmp, 'start.stops' )[ 4 ] )[ 1, ]
      }
      positions <- posns[ -1, ,drop=F ]
    }
    if ( is.vector( positions ) ) {
      if ( length( positions ) == 3 ) {
        tmp <- as.integer( positions[ 2:3 ] ); names( tmp )[ 1 ] <- positions[ 1 ]; positions <- tmp
      } else if ( is.character( positions ) && length( positions ) == 2 ) {
        positions <- as.integer( positions )
      }
      if ( is.null( names( positions )[ 1 ] ) ) names( positions )[ 1 ] <- names( e$genome.info$genome.seqs )[ 1 ]
      if ( length( positions ) == 1 ) {
        tmp <- c( positions[1] - window, positions[1] + window ); names( tmp )[ 1 ] <- names( positions )[ 1 ]
        positions <- tmp }
    }
    positions <- unique( positions )
    chr <- names( positions )[ 1 ]
    if ( is.null( chr ) ) chr <- names( e$genome.info$genome.seqs )[ 1 ]
    if ( is.vector( positions ) ) {
      positions <- matrix( positions, ncol=2 )
      rownames( positions )[ 1 ] <- chr
    }
    print( positions )
    ##require( data.table )
    scans <- NULL
    for ( i in 1:nrow( positions ) ) {
      if ( key(pssm.scans)[1] == 'gene' && key(pssm.scans)[2] == 'posns' ) { ## faster if set up right
        scans <- rbind( scans, pssm.scans[ CJ( gene, positions[ i, 1 ]:positions[ i, 2 ] ) ] )
      } else {
        scans <- rbind( scans, pssm.scans[ gene == rownames( positions )[ i ] & posns %in% positions[ i, 1 ]:positions[ i, 2 ] ] ) ##posns %betw% positions ]
      }
    }
    ## scans <- subset( pssm.scans, gene == chr & posns %betw% positions )
    out <- paste( "MOT", scans$bic, abs( scans$mots ), sep="_" )
  }
  
  out
}

get.motif.clusters <- function( motifs, expand.motif.cluster.clusters=T, filter.bad=T, ... ) {
  if ( ! missing( motifs ) ) {
    if ( exists( 'motifs.to.motif.clusts' ) ) {
      out <- motifs.to.motif.clusts[ motifs ]
    } else {
      tmp.mc <- paste( "MOTC", 1:mc.length, sep='_' ) ##length( motif.clusts ), sep="_" )
      if ( filter.bad && exists( "bad.clusts" ) ) tmp.mc <- tmp.mc[ ! tmp.mc %in% bad.clusts ]
      tmp <- get.motifs( motif.clust=tmp.mc, expand=F )
      out <- lapply( motifs, function( m ) names( which( sapply( tmp, function( i ) m %chin% i ) ) ) )
    }
    names( out ) <- motifs

    if ( exists( "motif.cluster.clusters" ) && expand.motif.cluster.clusters ) {
      for ( i in 1:length( out ) ) {
        if ( length( out[[ i ]] ) <= 0 ) next
        tmp <- as.integer( gsub( "MOTC_", "", out[[ i ]] ) )
        if ( motif.cluster.clusters[ tmp ] == 0 ) next
        tmp2 <- which( motif.cluster.clusters == motif.cluster.clusters[ tmp ] )
        out[[ i ]] <- paste( "MOTC", tmp2[ 1 ], sep="_" ) ## Use lowest-indexed mot clust in this mot-cluster-cluster.
      }
      out <- out[ ! sapply( out, function( i ) is.null( i ) || i == '' ||
                as.integer( gsub( "MOTC_", "", i ) ) > mc.length ) ] ## Don't count unclustered small motif clusters
    }    
  }    
  
  out
}

get.expression <- function( genes, biclusters, ... ) {
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

get.motif.info <- function( motifs, ... ) {
  if ( ! missing( motifs ) ) {
    inds <- strsplit( motifs, "_" )
    out <- lapply( inds, function( i ) {
      ms <- e$meme.scores[[ 1 ]][[ as.integer( i[ 2 ] ) ]]
      if ( is.null( ms$meme.out ) || length( ms$meme.out ) < as.integer( i[ 3 ] ) ) return( NULL )
      out <- ms$meme.out[[ as.integer( i[ 3 ] ) ]]
      tmp <- NULL
      if ( ! is.null( ms$pv.ev ) ) {
        if ( length( ms$pv.ev ) > 1 ) tmp <- ms$pv.ev[[ 2 ]]
        else if ( length( ms$pv.ev ) == 1 && 'mots' %chin% names( ms$pv.ev[[ 1 ]] ) ) tmp <- ms$pv.ev[[ 1 ]] ## catchall for Halo ensemble run
      }
      if ( ! is.null( tmp ) && nrow( tmp ) > 0 ) tmp <- subset( tmp, abs( mots ) == as.integer( i[ 3 ] ) )
      else tmp <- out$posns
      out$mast <- tmp
      out
    } )
    names( out ) <- motifs
  }
  
  out
}

get.bicluster.info <- function( biclusters, ... ) {
  bs <- as.integer( sapply( strsplit( biclusters, '_' ), '[', 2 ) )
  e$clusterStack[ bs ]
}

get.motif.cluster.info <- function( motif.clusters, ... ) {
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
## out$agglom() -- the default, counts of genes co-expressed with VNG0826C
## out$agglom( target='condition' ) -- counts of conditions in biclusters containing VNG0826C
## out$agglom( target='tf' ) -- counts of TFs regulating biclusters containing VNG0826C
## out$agglom( target='motif.cluster', path='bicluster,motif' ) -- counts of motif clusters in biclusters containing
##                  VNG0826C -- agglomeration is cumulative, first gene->bicluster->motif->motif.cluster
##     NOTE - p-value calculation does not work right now for paths that don't end in 'bicluster', so try this:
## out$agglom( srcType='motif.cluster', src='MOTC_5', path='motif,bicluster' ) -- counts of genes in biclusters
##                  containing a motif that is in motif cluster MOTC_5
## out$agglom( srcType='motif.cluster', src='MOTC_5', target='condition', path='motif,bicluster' ) -- counts of
##                  conditions in biclusters containing a motif that is in motif.cluster MOTC_5
## out$agglom( src='circadian_drk8_cycling_1440min_vs_NRC-1e.sig', srcType='condition', target='condition' ) --
##                  counts of conditions co-occurring with this cond (should be normalized by total.cond.count)
## out$agglom( src=c(Chr=716450,716550), srcType="position", target="motif.cluster", path="motif" )
##                  counts of motif.clusters that have motif instances that overlap the given chromosome coords
##                  NOTE: can use src='choose' to choose a position by clicking on a plot; see get.motifs()
## out$agglom( "MPN332", srcType='position', 'motif.cluster', 'motif' ) ## to get motif.clusters hitting gene's
##                  promoter region
## out$agglom( src=set.of.biclusters, 'bicluster', target='gene', path=NULL ) ## If it's only a one-step agglom.
## NEW: if 'path' contains 'bicluster', then 'cond.filter' can be used to filter out only biclusters that
##      contain all of the given conditions (in cond.filter); e.g. cond.filter=colnames(all.ratios)[1:2]
##      if 'p.val'==TRUE compute hypergeometric p-values
## NEW: now filter param is a list, can be used to filter out any of the steps in 'path', e.g.
##      filter=list(bicluster=c('BIC_1','BIC_2','BIC_3'), motif=c('MOT_1_1','MOT_2_1','MOT_3_1'))

## Translate name of obj into type (NOTE doesn't catch TFs!)
guess.type <- function( chr ) {
  if ( grepl( 'BIC_', chr ) ) return( 'bicluster' )
  if ( grepl( 'MOT_', chr ) ) return( 'motif' )
  if ( grepl( 'MOTC_', chr ) ) return( 'motif.cluster' )
  if ( ! is.null( e$genome.info$gene.regex ) && grepl( e$genome.info$gene.regex, chr, ignore=T ) ) return( 'gene' )
  if ( ! is.null( e$genome.info$gene.prefix ) && grepl( e$genome.info$gene.prefix, chr, ignore=T ) ) return( 'gene' )
  if ( is.integer( chr ) || is.numeric( chr ) ) return( 'position' )
}

## get.type <- function( obj ) {
##   prefix <- sapply( strsplit( obj, "_" ), "[", 1 )
##   out <- ifelse( prefix == "BIC", "bicluster",
##                 ifelse( prefix == "MOT", "motif",
##                        ifelse( prefix == "MOTC", "motif.cluster",
##                               ifelse( obj %in% attr( e$ratios, "rnames" ), "gene", "condition" ) ) ) )
##   out
## }

agglom <- function( src='VNG0826C', targetType='gene', path='bicluster', srcType=guess.type( src ), ##'gene',
                   p.val=T, q.val="BH", ##cond.filter=NULL, cond.filter.frac=0.5, cond.filter.die=F, resid.filter=NA,
                   ##biclust.filter=NULL, resid.filter.die=F,
                   filter=list(), filter.die=F, verbose=F, tf.weights=T, ... ) {
  orig.src <- src; orig.srcType <- srcType
  if ( ! is.null( path ) ) {
  by <- unlist( strsplit( path, ',' ) )
  for ( b in by ) {
    text <- paste( 'bys <- unlist( get.', b,   's( ', srcType, "=src, ... ) )", sep="" )
    if ( verbose ) print( text )
    eval( parse( text=text ) )
    if ( length( bys ) <= 0 ) {
      cat( "Oops!\n" )
      return( NULL )
    }
    if ( verbose ) print( length( bys ) )
    if ( b %in% names( filter ) ) {
      bys <- bys[ bys %chin% filter[[ b ]] ]
      if ( verbose ) print( length( bys ) )
    }
    if ( length( bys ) <= 0 && filter.die ) {
      cat( "No", b, "s pass the filter! Perhaps too stringent?\n" )
      return( NULL )
    }
    bys <- bys[ ! is.na( bys ) ]
    src <- bys
    srcType <- b
  }
  } else {
    bys <- src; by=srcType
  }
  text <- paste( 'out <- get.', targetType, 's( ', by[ length( by ) ],      '=bys, ... )', sep="" )
  if ( verbose ) print( text )
  eval( parse( text=text ) )
  ##if ( verbose ) print( length( out ) )
  if ( targetType == 'tf' ) {
    out.orig <- out
    out <- lapply( out, lapply, names )
  }
  tab <- sort( table( unlist( out ) ), decreasing=T )
  if ( verbose ) print( sum( tab ) )
  out2 <- t( t( tab ) )
  colnames( out2 ) <- "count"
  tab2 <- tab[ tab > 1 ]
  pvals <- qvals <- numeric()
  ## Right now p-values can only be computed for for path='..,..,bicluster' but this can be fixed!
  if ( p.val && length( tab2 ) > 0 ) {
    n.bics <- e$k.clust
    bics1 <- as.list( names( tab ) ); names( bics1 ) <- bics1
    if ( by[ length( by ) ] == 'bicluster' ) {
      ##bics1 <- as.list( names( tab ) ); names( bics1 ) <- bics1
      if ( targetType != "bicluster" ) {
        text <- paste( 'bics1 <- get.biclusters( ', targetType, 's=names( tab2 ), ... )', sep="" )
        if ( verbose ) print( text )
        eval( parse( text=text ) )
        ##if ( verbose ) print( length( unlist( bics1 ) ) )
      }
      if ( 'bicluster' %in% names( filter ) ) {
        bys <- bys[ bys %chin% filter$bicluster ]
        bics1 <- lapply( bics1, function( b ) b[ b %chin% filter$bicluster ] )
        n.bics <- length( bics1 )
      }
    } else if ( by[ length( by ) ] == 'motif' ) {
      bics1 <- as.list( names( tab ) ); names( bics1 ) <- bics1
      if ( targetType != "motif" ) {
        text <- paste( 'bics1 <- get.motifs( ', targetType, '=names( tab2 ), ... )', sep="" )
        if ( verbose ) print( text )
        eval( parse( text=text ) )
        if ( verbose ) print( length( unlist( bics1 ) ) )
      }
      n.bics <- e$k.clust * 2
      if ( targetType == 'motif.cluster' ) n.bics <- length( unlist( motif.clusts ) )
      if ( 'motif' %in% names( filter ) ) {
        bys <- bys[ bys %chin% filter$motif ]
        bics1 <- lapply( bics1, function( b ) b[ b %chin% filter$motif ] )
        n.bics <- length( bics1 )
      }
    }
    pvals <- phyper( tab2, length( out ), n.bics, sapply( bics1, length )[ names( tab2 ) ], lower=F )
    ##print(cbind(tab2,length(out),n.bics,sapply( bics1, length )[ names( tab2 ) ], pvals ))
    ## if ( p.correct ) {
    ##   pvals[ tab2 <= 1 ] <- 1
    ##   pvals <- pvals * length( pvals )
    ##   pvals[ pvals > 1 ] <- 1
    ## }
    if ( ! is.na( q.val ) ) {
      if ( q.val == TRUE ) q.val <- "BH"
      ##require( qvalue )
      pv <- c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) )
      ##pv <- pv + 1e-200; pv[ pv > 1 ] <- 1
      ##qvals <- try( qvalue( pv, robust=T ) )
      ##if ( "try-error" %in% class( qvals ) || ( ! is.list( qvals ) && all( qvals == 0 ) ) )
      ##  qvals <- pvals * 0 ## I have found this means all p-values are very low
      ##else qvals <- qvals$qvalues
      qvals <- p.adjust( pv, q.val )
    }
  }
  if ( p.val ) {
    out2 <- cbind( out2, p.value=c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) ) )
    if ( ! is.na( q.val ) ) out2 <- cbind( out2, q.value=c( qvals, rep( 1, nrow( out2 ) - length( qvals ) ) ) )
  }
  if ( tf.weights ) {
    if ( targetType == 'tf' ) {
      if ( verbose ) cat( "Computing tf weight quantiles for", sum( out2[ ,'count' ] > 2 ), "tfs\n" )
      weight.tab <- t( sapply( rownames( out2 ), function( tf ) {
        if ( out2[ tf, 'count' ] <= 2 ) return( c( NA, NA, NA, NA, NA ) )
        wts <- unlist( lapply( out.orig, function( i ) i$coeffs[ tf ] ) )
        wts <- wts[ ! is.na( wts ) ]
        wts2 <- unlist( lapply( out.orig, function( i ) i$poss.reg[ tf ] ) )
        wts2 <- wts2[ ! is.na( wts2 ) ]
        c( length( wts ), length( wts2 ), quantile( wts, c(0.1,0.5,0.9) ) )
      } ) )
      colnames( weight.tab ) <- c( "n.inf", "n.poss", "10%", "50%", "90%" )
      out2 <- cbind( out2, weight.tab )
    } else if ( orig.srcType == 'tf' ) {
      if ( srcType == 'bicluster' ) {
        if ( verbose ) cat( "Computing tf weight quantiles for", sum( out2[ ,'count' ] > 2 ), targetType, "s\n" )
        coefs <- get.tfs( bic=bys )
        weight.tab <- t( sapply( rownames( out2 ), function( g ) {
          if ( out2[ g, 'count' ] <= 2 ) return( c( NA, NA, NA, NA, NA ) )
          bics <- get.biclusters( gene=g )[[ 1 ]]
          bics <- bics[ bics %chin% names( coefs ) ]
          if ( length( bics ) <= 2 ) return( c( NA, NA, NA, NA, NA ) )
          cc <- coefs[ bics ]
          wts <- unlist( lapply( cc, function( i ) i$coeffs[ orig.src ] ) )
          wts <- wts[ ! is.na( wts ) ]
          wts2 <- unlist( lapply( cc, function( i ) i$poss.reg[ orig.src ] ) )
          wts2 <- wts2[ ! is.na( wts2 ) ]
          c( length( wts ), length( wts2 ), quantile( wts, c(0.1,0.5,0.9) ) )
        } ) )
        colnames( weight.tab ) <- c( "n.inf", "n.poss", "10%", "50%", "90%" )
        out2 <- cbind( out2, weight.tab )
      }
    }
  }
  out2 <- as.data.frame( out2 )
  attr( out2, "total" ) <- length( bys )
  out2
}

## agglom <- function( src='MOTC_1', srcType='motif.cluster', targetType='gene', path='motif,bicluster', p.val=T, 
##                    q.val="BH", cond.filter=NULL, cond.filter.frac=0.5, verbose=F ) {
##   orig.src <- src; orig.srcType <- srcType
##   if ( ! is.null( path ) ) {
##   by <- unlist( strsplit( path, ',' ) )
##   for ( b in by ) {
##     text <- paste( 'bys <- unlist( get.', b,   's( ', srcType, "=src ) )", sep="" )
##     if ( verbose ) print( text )
##     eval( parse( text=text ) )
##     if ( ! is.null( cond.filter ) && b == 'bicluster' ) {
##       if ( verbose ) print( 'cc <- get.conditions( biclust=bys )' )
##       cc <- get.conditions( biclust=bys )
##       tmp <- sapply( cc, function( i ) mean( cond.filter %in% i ) >= cond.filter.frac ) ## 50% of cond.filter in the biclusters
##       if ( ! any( tmp ) ) warning( "No biclusters pass the cond.filter! Perhaps too many conditions?" )
##       else bys <- bys[ tmp ]
##     }
##     src <- bys
##     srcType <- b
##   }
##   } else {
##     bys <- src; by=srcType
##   }
##   text <- paste( 'out <- get.', targetType, 's( ', by[ length( by ) ],      '=bys )', sep="" )
##   if ( verbose ) print( text )
##   eval( parse( text=text ) )
##   if ( targetType == 'tf' ) {
##     out.orig <- out
##     out <- lapply( out, lapply, names )
##   }
##   tab <- sort( table( unlist( out ) ), decreasing=T )
##   out2 <- t( t( tab ) )
##   colnames( out2 ) <- "count"
##   tab2 <- tab[ tab > 1 ]
##   pvals <- qvals <- numeric()
##   if ( p.val && length( tab2 ) > 0 && by[ length( by ) ] == 'bicluster' ) { ## Right now this is only valid for path='..,..,bicluster'
##     bics1 <- as.list( names( tab ) ); names( bics1 ) <- bics1
##     if ( targetType != "bicluster" ) {
##       text <- paste( 'bics1 <- get.biclusters( ', targetType,      '=names( tab2 ) )', sep="" )
##       if ( verbose ) print( text )
##       eval( parse( text=text ) )
##     }
##     n.bics <- e$k.clust
##     if ( ! is.null( cond.filter ) ) {
##       bics <- table( unlist( get.biclusters( cond=cond.filter ) ) )
##       bics <- bics[ bics >= length( cond.filter ) * cond.filter.frac ] ##as.integer( gsub( "BIC_", "",
##       bics1 <- lapply( bics1, function( i ) i[ i %in% names( bics ) ] )
##       n.bics <- length( bics )
##     }
##     pvals <- phyper( tab2, length( out ), n.bics, sapply( bics1, length )[ names( tab2 ) ], lower=F )
##     if ( ! is.na( q.val ) ) {
##       if ( q.val == TRUE ) q.val <- "BH"
##       pv <- c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) )
##       qvals <- p.adjust( pv, q.val )
##     }
##   }
##   if ( p.val ) {
##     out2 <- cbind( out2, p.value=c( pvals, rep( 1, nrow( out2 ) - length( pvals ) ) ) )
##     if ( ! is.na( q.val ) ) out2 <- cbind( out2, q.value=c( qvals, rep( 1, nrow( out2 ) - length( qvals ) ) ) )
##   }
##   if ( tf.weights ) {
##     if ( targetType == 'tf' ) {
##       if ( verbose ) cat( "Computing tf weight quantiles for", sum( out2[ ,'count' ] > 2 ), "tfs\n" )
##       weight.tab <- t( sapply( rownames( out2 ), function( tf ) {
##         if ( out2[ tf, 'count' ] <= 2 ) return( c( NA, NA, NA, NA, NA ) )
##         wts <- unlist( lapply( out.orig, function( i ) i$coeffs[ tf ] ) )
##         wts <- wts[ ! is.na( wts ) ]
##         wts2 <- unlist( lapply( out.orig, function( i ) i$poss.reg[ tf ] ) )
##         wts2 <- wts2[ ! is.na( wts2 ) ]
##         c( length( wts ), length( wts2 ), quantile( wts, c(0.1,0.5,0.9) ) )
##       } ) )
##       colnames( weight.tab ) <- c( "n.inf", "n.poss", "10%", "50%", "90%" )
##       out2 <- cbind( out2, weight.tab )
##     } else if ( orig.srcType == 'tf' ) {
##       if ( srcType == 'bicluster' ) {
##         if ( verbose ) cat( "Computing tf weight quantiles for", sum( out2[ ,'count' ] > 2 ), targetType, "s\n" )
##         coefs <- get.tfs( bic=bys )
##         weight.tab <- t( sapply( rownames( out2 ), function( g ) {
##           if ( out2[ g, 'count' ] <= 2 ) return( c( NA, NA, NA, NA, NA ) )
##           bics <- get.biclusters( gene=g )[[ 1 ]]
##           bics <- bics[ bics %in% names( coefs ) ]
##           if ( length( bics ) <= 2 ) return( c( NA, NA, NA, NA, NA ) )
##           cc <- coefs[ bics ]
##           wts <- unlist( lapply( cc, function( i ) i$coeffs[ orig.src ] ) )
##           wts <- wts[ ! is.na( wts ) ]
##           wts2 <- unlist( lapply( cc, function( i ) i$poss.reg[ orig.src ] ) )
##           wts2 <- wts2[ ! is.na( wts2 ) ]
##           c( length( wts ), length( wts2 ), quantile( wts, c(0.1,0.5,0.9) ) )
##         } ) )
##         colnames( weight.tab ) <- c( "n.inf", "n.poss", "10%", "50%", "90%" )
##         out2 <- cbind( out2, weight.tab )
##       }
##     }
##   }
##   out2 <- as.data.frame( out2 )
##   attr( out2, "total" ) <- length( bys )
##   out2
## }

get.tfs <- function( biclusters, weight.cutoffs=c(0.01, 2), boot.results=F, ... ) {
  if ( ! missing( biclusters ) ) {
    inds <- as.integer( sapply( strsplit( biclusters, "_" ), "[", 2 ) )
    if ( ! boot.results || ! exists( 'e.coeffs.big' ) ) {
      out <- lapply( inds, function( i ) { tmp <- try( e.coeffs[[ i ]] )
                                           if ( is.null( tmp ) || class( tmp ) == 'try-error' ) return( NULL );
                                           tmp <- tmp[[ 1 ]]; if ( is.null( tmp ) ) return( NULL );
                                           list( coeffs=tmp$coeffs[ abs( tmp$coeffs ) %betw% weight.cutoffs ],
                                                poss.reg=tmp$possibly.regulates ) } )
    } else {
      out <- lapply( inds, function( i ) { tmp <- try( e.coeffs.big[[ i ]] )
                                           if ( is.null( tmp ) || class( tmp ) == 'try-error' ) return( NULL );
                                           tmp <- tmp[[ 1 ]]; if ( is.null( tmp ) ) return( NULL );
                                           qq <- tmp$coef.quantiles[ ,'50%' ]
                                           list( coeffs=qq[ abs( qq ) %betw% weight.cutoffs ],
                                                poss.reg=tmp$possibly.regulates ) } )
    }
    names( out ) <- biclusters
  }

  out
}

get.predicted.dynamics <- function( biclusters, ... ) {
  if ( ! missing( biclusters ) ) {
    inds <- as.integer( sapply( strsplit( biclusters, "_" ), "[", 2 ) )
    out <- lapply( inds, function( i ) e.coeffs.big[[ i ]][[ 1 ]]$pred.ts[ 1, ] )
    names( out ) <- biclusters
  }

  out
}

## Agglomerate from two "sources", to biclusters, and compute the p-value of the bicluster overlap, e.g.
##   agglom2( "MOTC_1", "motif,bicluster", "MOTC_170", "motif,bicluster" )
agglom2 <- function( src="MOTC_1", srcType=guess.type( src ), srcPath="motif,bicluster",
                    target="MOTC_170", targetType=guess.type( target ), targetPath="motif,bicluster",
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

if ( FALSE ) {
  for(cond in colnames(e$ratios$ratios)){
    bc <- as.integer(gsub('BIC_','',out$get.biclusters(cond=cond)[[1]]))
    cat(cond,length(bc),"\n")
    if(cond%in%names(nets))next
    nets[[cond]]<-Matrix(out$gene.gene.network(bc,T))
    print(range(nets[[cond]],na.rm=T))
  }
  nets2<-sapply(nets,as.vector)
  d=dist(t(nets2))
}

## Biclust.filter = a list of integers (k's) which are GOOD biclusters to use for the network (e.g. low residual)
## do this:
##    resids=unlist(lapply(e$clusterStack,'[[','resid'))
##    biclust.filter=which(resids<=0.5&!is.na(resids))
##    can use e.values as well:
##    evals=unlist(lapply(e$clusterStack,function(i)min(i$e.val,na.rm=T)))
## or filter by condition
## bc <- as.integer(gsub('BIC_','',out$get.biclusters(cond=x)[[1]]))
## motif.based - default is to use co-bicluster membership for edges; if motif.based=TRUE,
##    use co-motif membership for edges
gene.gene.network <- function( biclust.filter=NA, matrix.only=F, motif.based=F ) {
  all.rm <- lapply( rownames( e$ratios$ratios ), function( g ) integer() )
  names( all.rm ) <- rownames( e$ratios$ratios )
  if ( ! motif.based ) {
    for ( k in 1:e$k.clust ) {
      rows <- e$clusterStack[[ k ]]$rows; for ( r in rows ) all.rm[[ r ]] <- c( all.rm[[ r ]], k )
    }
    if ( ! is.na( biclust.filter ) ) all.rm <- mclapply( all.rm, function( i ) i[ i %in% biclust.filter ] )
  } else { ## motif.based -- use out$genes.to.motifs
    all.rm <- genes.to.motifs
    all.rm <- all.rm[ rownames(e$ratios$ratios) ]
  }
  row.ov <- do.call( rbind, mclapply( all.rm, function( i ) sapply( all.rm, function( j ) sum( i %in% j ) ) ) )
  ## Matrix of the total # of (unique) clusters that genes i OR j are in (union)
  all.len <- lapply( all.rm, length )
  n.clust <- do.call( rbind, mclapply( all.len, function( i ) sapply( all.len, function( j ) min( c( i, j ) ) ) ) )
  rm( all.rm, all.len )

  row.ov[ lower.tri( row.ov, diag=T ) ] <- 0
  n.clust[ lower.tri( n.clust, diag=T ) ] <- 1
  rownames( row.ov ) <- rownames( n.clust ) <- rownames( e$ratios$ratios )
  if ( matrix.only ) return( row.ov / n.clust )
  
  row.ov[ n.clust < 10 ] <- 0; n.clust[ n.clust < 10 ] <- Inf
  r.sif <- which( row.ov / n.clust > 0.2, arr=T )
  r.sif <- data.frame( g1=rownames( r.sif ), g2=rownames( row.ov )[ r.sif[ ,2 ] ], type=rep( "bic", nrow( r.sif ) ),
                      weight=row.ov[ r.sif ] / n.clust[ r.sif ] )
  r.sif <- r.sif[ order( r.sif$weight, decreasing=T ), ]; rownames( r.sif ) <- NULL

  gs <- unique( c( as.character( r.sif$g1 ), as.character( r.sif$g2 ) ) )
  r.na <- data.frame( g=gs, short=e$get.long.names( gs, short=T )[ gs ] )
  sh <- as.character( r.na$short ); sh[ which( sh == "" ) ] <- as.character( r.na$g[ which( sh == "" ) ] )
  r.na$short <- as.factor( sh ); rm( sh, gs )
  list( sif=r.sif, na=r.na )
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
                         filter=list(condition=conditions) )
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

  out2 <- apply.func( 1:mc.length, function( mc ) { ##length( motif.clusts ), function( mc ) {
    mc <- paste( "MOTC", mc, sep="_" )
    print( mc )
    tmp <- agglom( src=mc, srcType='motif.cluster', targetType='gene', path='motif',
                         q.val=q.val, filter=list(condition=conditions) )
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
  if ( ! exists( "p.scans" ) ) {
    p.scans <- pssm.scans
    require( data.table )
    if ( ! is.data.table( p.scans ) ) {
      p.scans$mots <- abs( p.scans$mots ); gc() ## Don't use strand or p-values for these funcs, so remove em
      p.scans <- p.scans[ ,c( "bic", "mots", "gene", "posns" ) ]; gc() ##, "fwd", "pvals" ) ]; gc()
      p.scans <- as.data.table( p.scans ); gc()
      setkey( p.scans, bic, mots, gene, posns )
    }
  }
  
  if ( motifs != "ALL" ) {
    mots <- strsplit( gsub( "MOT_", "", motifs ), "_" )
    bi <- as.integer( sapply( mots, "[", 1 ) )
    mo <- as.integer( sapply( mots, "[", 2 ) )
    scans <- p.scans[ J( c( bi, bi ), c( mo, -mo ) ), allow.cart=T ] ##p.scans[ J( bi, mo ) ] ##, names( seqs ) ) ]
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

get.pssm.scans <- function( motifs="ALL", seqs=e$genome.info$genome.seqs,
                           seq.type=names( e$meme.scores )[ 1 ], p.cutoff='0.0001', force=(motifs!='ALL') ) {
  out.p.cutoff <- as.numeric( p.cutoff ) / 10
  tmpf <- sprintf( 'filehash/pssm_scans_%s.tsv', as.character( out.p.cutoff ) )
  if ( force || ! file.exists( sprintf( "%s.bz2", tmpf ) ) ) {
    mast.cmd <- './progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.0001 -ev 9999 -comp' ## Different from mast cmd used for cMonkey; e.g. doesn't have -revcomp, also set the p-value cutoff to 1e-4 (-mt 0.0001)
    mast.cmd <- sprintf( './progs/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt %s -ev 9999 -comp', p.cutoff )

    ks <- 1:e$k.clust
    if ( motifs != "ALL" && ! is.null( motifs ) )
      ks <- as.integer( gsub( "BIC_", "", unlist( get.biclusters( motif=motifs ) ) ) )
    dir.create( 'filehash' )
    dir.create( sprintf( 'filehash/pssm_scans_%s', as.character( out.p.cutoff ) ) )

    qqq <- FALSE
    if ( is.null( e$genome.info$genome.seqs ) ) { ## might be the case if discard.genome==TRUE during cMonkey run
      genome.info <- e$get.genome.info()
      e$genome.info$genome.seqs <- genome.info$genome.seqs
      rm( genome.info ); gc()
      qqq <- TRUE
    }
    if ( is.list( e$genome.info$genome.seqs ) ) e$genome.info$genome.seqs <- unlist( e$genome.info$genome.seqs )
    if ( qqq ) e$genome.info$bg.list[[ 1 ]] <- e$mkBgFile( e$genome.info$all.upstream.seqs, order=e$bg.order[ 1 ],
                                                          bgfname=e$genome.info$bg.fname[ 1 ],
                                                          use.rev.comp=grepl( "-revcomp", e$meme.cmd[ 1 ] ) )
    
    tab <- mclapply( ks, function( k ) {
      tmpf <- sprintf( 'filehash/pssm_scans_%s/%08d.tsv', as.character( out.p.cutoff ), k )
      if ( file.exists( sprintf( "%s.bz2", tmpf ) ) && file.info( sprintf( "%s.bz2", tmpf ) )[ 'size' ] > 1 ) return( NULL )
      else if ( file.exists( sprintf( "%s.NULL", tmpf ) ) ) return( NULL )
      cat( "SCANNING MOTIFS across genome for bicluster:", k, "\n" )
      ##if ( is.null( e$meme.scores[[ seq.type ]][[ k ]] ) || e$meme.scores[[ seq.type ]][[ k ]] == "" ||
      ##    is.null( e$meme.scores[[ seq.type ]][[ k ]]$meme.out ) ) return( NULL )
      tmp1 <- try( e$cluster.meme.motif.lines( k, logodds=T ) )
      if ( class( tmp1 ) == 'try-error' || length( tmp1 ) < 10 ) { system( sprintf( "touch %s.NULL", tmpf ) ); return( NULL ) } ## No motifs
      tmp3 <- try( e$runMast( tmp1, mast.cmd, names( seqs ), seqs, verbose=F, seq.type=seq.type,
                             bg.list=e$genome.info$bg.list[[ 1 ]], bg.fname=e$genome.info$bg.fname[ 1 ], unlink=T ) ) ##NULL, unlink=T ) ##)
      if ( class( tmp3 ) == 'try-error' || length( tmp3 ) <= 0 ) { system( sprintf( "touch %s.NULL", tmpf ) ); return( NULL ) }
      tmp4 <- try( e$getMastPValuesAndEValues( tmp3, names( seqs ) )[[ 2 ]] )
      if ( class( tmp4 ) == 'try-error' || nrow( tmp4 ) <= 0 ) { system( sprintf( "touch %s.NULL", tmpf ) ); return( NULL ) }
      tmp4 <- subset( tmp4, pvals <= out.p.cutoff )
      if ( class( tmp4 ) == 'try-error' || nrow( tmp4 ) <= 0 ) { system( sprintf( "touch %s.NULL", tmpf ) ); return( NULL ) }
      tmp4 <- cbind( tmp4, bic=rep( k, nrow( tmp4 ) ) )
      if ( class( tmp4 ) == 'try-error' ) { system( sprintf( "touch %s.NULL", tmpf ) ); return( NULL ) }
      write.table( tmp4, col.names=F, quote=F, sep='\t', file=tmpf, append=F, row.names=F )
      system( sprintf( "bzip2 -fv %s", tmpf ) )
      ##tmp4
      NULL
    }, mc.preschedule=F )

    ##for ( i in 1:length( tab ) ) if ( ! is.null( tab[[ i ]] ) && ! class( tab[[ i ]] ) == 'try-error' ) {
    ##  write.table( tab[[ i ]], col.names=(i==1), quote=F, sep='\t', file=tmpf, append=(i!=1), row.names=F )
    ##}
    ##system( sprintf( "bzip2 -fv %s", tmpf ) )
    system( sprintf( "find ./filehash/pssm_scans_%s/ -name '*.tsv.bz2' -print | sort | xargs bunzip2 -c | bzip2 -c >%s", as.character( out.p.cutoff ),
                    sprintf( "%s.bz2", tmpf ) ) )
  }
  cat( "Reading", sprintf( "%s.bz2", tmpf ), "\n" )
  tab <- read.delim( bzfile( sprintf( "%s.bz2", tmpf ) ), sep='\t', head=F )
  #tab <- fread( bzfile( sprintf( "%s.bz2", tmpf ) ), head=F ) ## fread cant read from connections yet
  tab <- as.data.table( tab )
  setnames( tab, c('gene', 'pvals', 'posns', 'mots', 'bic') )
  setkey( tab, bic, mots, gene, posns )
  tab
}

get.motif.coding.fracs <- function( motifs='ALL', verbose=F, p.cutoff=1e-6 ) {
  coding.seqs <- lapply( names( e$genome.info$genome.seqs ), function( n )
                        rep( FALSE, nchar( e$genome.info$genome.seqs[ n ] ) ) )
  names( coding.seqs ) <- names( e$genome.info$genome.seqs )
  tmp <- subset( e$genome.info$feature.tab, type %in% c( 'CDS', 'rRNA', 'tRNA' ) )
  tmp$where <- as.character( tmp$contig )
  tmp$Start <- as.integer( as.character( tmp$start_pos ) )
  tmp$Stop <- as.integer( as.character( tmp$end_pos ) )
  for ( i in 1:nrow( tmp ) ) {
    ttmp <- tmp[ i, ]
    if ( ! is.na( ttmp$where ) && ! is.na( ttmp$Start ) && ! is.na( ttmp$Stop ) )
      coding.seqs[[ ttmp$where ]][ ttmp$Start:ttmp$Stop ] <- TRUE
    rm( ttmp )
  }
  
  ##mots <- rep( 1:e$k.clust, 2 ) ##strsplit( gsub( "MOT_", "", motifs ), "_" )
  scans <- as.data.table( subset( pssm.scans, pvals <= p.cutoff ) ) ##.2[ J( c( bi, bi ), c( mo, -mo ), chr ) ] ##subset( pssm.scans, bic %in% bi )
  setkey( scans, 'bic', 'mots', 'gene', 'posns' )
  if ( motifs == 'ALL' ) {
    ##motifs <- unlist( get.motifs( paste( 'BIC', 1:e$k.clust, sep='_' ) ) )
    all.motifs <- unlist( lapply( 1:nrow( out$motif.widths ),
                             function( i ) { ii <- which( out$motif.widths[ i, ] > 0 )
                                             if ( length( ii ) <= 0 ) return( NULL )
                                             paste( 'MOT', i, ii, sep='_' ) } ) ) ## This is faster than the above!
    ## Use only the motifs that are left in the pssm.scans after filtering the pssm.scans by p-value
    tmp <- unique( data.frame( mots=abs( scans$mots ), bic=scans$bic ) )
    bi <- tmp$bic
    mo <- tmp$mots
    motifs <- paste( 'MOT', bi, mo, sep='_' )
    rm( tmp )
  } else {
    mots <- strsplit( gsub( "MOT_", "", motifs ), "_" )
    bi <- as.integer( sapply( mots, "[", 1 ) )
    mo <- as.integer( sapply( mots, "[", 2 ) )
    scans <- scans[ J( c( bi, bi ), c( mo, -mo ) ), allow.cart=T ]  ## this speeds it up a lot
  }
  scans <- scans[ ! is.na( scans$posns ), ]
  setkey( scans, 'bic', 'mots', 'gene', 'posns' )
  in.coding <- rep( NA, nrow( scans ) )
  for ( cc in names( coding.seqs ) ) {
    print( cc ); 
    in.coding[ scans$gene == cc ] <- coding.seqs[[ cc ]][ scans[ scans$gene == cc, ]$posns ]
  }
  scans$in.coding <- in.coding
  scans <- as.data.table( scans )
  setkey( scans, 'bic', 'mots', 'gene', 'posns' )
  if ( verbose ) { cat( "Scanning for bad motifs..." ); print( length( bi ) ) }
  gc()

  frac.in.coding <- do.call( c, mclapply( 1:length( bi ), function( i ) {
    sc <- scans[ bic == bi[i] & abs( mots ) == mo[i] ] ## J( c( bi, bi ), c( mo, -mo ) ) ]
    if ( nrow( sc ) <= 0 || all( is.na( sc$posns ) ) ) return( NA ) ##next
    hits <- sc$in.coding
    if ( verbose && i %% 100 == 1 ) cat( i, ' ' )
    mean( hits, na.rm=T )
  }, mc.preschedule=F ) )

  if ( verbose ) cat( '\n' )
  ##motifs.orig <- motifs
  ##motifs <- motifs[ ! is.na( frac.in.coding ) & frac.in.coding < mean( unlist( coding.seqs ) ) ]
  motif.coding.fracs <- frac.in.coding
  names( motif.coding.fracs ) <- motifs
  if ( exists( 'all.motifs' ) ) {
    motif.coding.fracs <- motif.coding.fracs[ all.motifs ]
    names( motif.coding.fracs ) <- all.motifs
  }
  mean.coding.frac <- sum( sapply( coding.seqs, sum, na.rm=T ) ) / sum( sapply( coding.seqs, length ) ) ##mean( unlist( coding.seqs ), na.rm=T )
  ##rm( coding.seqs, scans, tmp, frac.in.coding, in.coding )
  out <- list( all.fracs=motif.coding.fracs, mean.fracs=mean.coding.frac, p.cutoff=p.cutoff )
  return( out )
} 

plot.genes.in.region <- function( center, where, window, new.plot=T, yoff=0, yscale=1, plot.axis=T, org=NA, gene.names=T,
                                 main="", gene.coords=NULL, op.shift=F, ... ) {
  
  if ( is.null( gene.coords ) ) gene.coords <- e$get.gene.coords( NULL, op.shift=op.shift )
  gene.coords <- subset( gene.coords, ! is.na( start_pos ) & ! is.na( end_pos ) )
  coords <- gene.coords

  window <- round( c( center - window / 2, center + window / 2 ) )
  w.expand <- window + c( -1000, 1000 ) ##round( c( window[ 1 ] - diff( window ) / 5, window[ 2 ] + diff( window ) / 5 ) )
  in.window <- function( x, w ) { x >= w[ 1 ] & x <= w[ 2 ] }

  starts <- as.integer( as.character( coords[ ,'start_pos' ] ) )
  stops <- as.integer( as.character( coords[ ,'end_pos' ] ) )
  genes.in <- coords[ tolower( as.character( coords[ ,'contig' ] ) ) == tolower( where ) &
                     ( in.window( starts, w.expand ) | ## WTF!!??
                      in.window( stops, w.expand ) ), ,drop=F ]

  if ( nrow( genes.in ) <= 0 ) {
    if ( new.plot ) {
      plot( 0, 0, xlim=window, ylim=c(-1,1), ann=F, xaxt="n", yaxt="n", bty="n" )
      if ( plot.axis ) axis( side=1, pos=0, mgp=c(0,0.5,0) )
    }
    return()
  }
  starts <- as.integer( as.vector( genes.in[ ,'start_pos' ] ) )
  stops <- as.integer( as.vector( genes.in[ ,'end_pos' ] ) )
  genes.x <- apply( cbind( starts, stops ), 1, mean )
  dirs <- ifelse( genes.in$strand == "D", 1, -1 )
  ##genes.x[ ( ! ( genes.x %betw% window ) ) & ( starts %betw% window ) ] <- 
    ##starts[ ( ! ( genes.x %betw% window ) ) & ( starts %betw% window ) ] - diff( window ) / 10 * dirs
  genes.x[ genes.x > window[ 2 ] - 10 & starts %betw% window ] <- window[ 2 ] - 10
  ##genes.x[ ( ! ( genes.x %betw% window ) ) & ( stops %betw% window ) ] <-
  ##  stops[ ( ! ( genes.x %betw% window ) ) & ( stops %betw% window ) ] + diff( window ) / 10 * dirs
  genes.x[ genes.x < window[ 1 ] + 10 & stops %betw% window ] <- window[ 1 ] + 10

  genes.is.fwd <- as.character( genes.in[ ,'strand' ] ) == "D"
  
  if ( new.plot ) {
    plot( 0, 0, xlim=window, ylim=c(-1,1), ann=F, xaxt="n", yaxt="n", bty="n" )
    if ( plot.axis ) axis( side=1, pos=0, mgp=c(0,0.5,0) )
  }
  for ( i in 1:nrow( genes.in ) ) {
    start <- as.integer( as.vector( genes.in[ i, "start_pos" ] ) )
    end <- as.integer( as.vector( genes.in[ i, "end_pos" ] ) )
    name <- genes.in[ i, "names" ]
    if ( genes.in[ i, "strand" ] == "D" ) {
      rect( start, yoff+yscale*0.2, end, yoff+yscale*1, col="yellow", border="black", lwd=1 )
    } else if ( genes.in[ i, "strand" ] == "R" ) {
      rect( start, yoff+yscale*-1, end, yoff+yscale*(-0.2), col="orange", border="black", lwd=1 )
    }
  }
  if ( gene.names ) {
    regex <- if ( ! is.null( e$genome.info$gene.regex ) ) e$genome.info$gene.regex else e$genome.info$gene.prefix
    for ( i in 1:nrow( genes.in ) ) {
      start <- as.integer( as.vector( genes.in[ i, "start_pos" ] ) )
      end <- as.integer( as.vector( genes.in[ i, "end_pos" ] ) )
      name <- as.character( genes.in[ i, "names" ] )
      name <- grep( paste( "^", regex, sep='' ), e$get.synonyms( name )[[ 1 ]], val=T, ignore=T )
      if ( length( name ) > 1 ) name <- name[ name %in% rownames( e$ratios$ratios ) ]
      if ( length( name ) <= 0 ) next
      if ( genes.in[ i, "strand" ] == "D" ) {
        text( genes.x[ i ], yoff+yscale*0.6, labels=name, adj=c(0.5,0.5), col="black", ... )
      } else if ( genes.in[ i, "strand" ] == "R" ) {
        text( genes.x[ i ], yoff+yscale*(-0.6), labels=name, adj=c(0.5,0.5), col="black", ... )
      }
    }
  }
}

## Expand motif hits to genes in same operons as those listed, since operon genes were not searched by meme (other than 1st gene)
get.meme.hits.table <- function( motifs="ALL", seqs=e$genome.info$genome.seqs, force=F,
                                seq.type=names( e$meme.scores )[ 1 ], p.cutoff=1e-4, expand.operons=T ) {
  ## o.list <- e$operon.list()  
  if ( motifs == 'ALL' ) {
    motifs <- unlist( lapply( 1:nrow( out$motif.widths ),
                             function( i ) { ii <- which( out$motif.widths[ i, ] > 0 )
                                             if ( length( ii ) <= 0 ) return( NULL )
                                             paste( 'MOT', i, ii, sep='_' ) } ) )
  }
  tmpf <- sprintf( 'filehash/meme_hits_%s.tsv', as.character( p.cutoff ) )
  if ( force || ! file.exists( sprintf( "%s.bz2", tmpf ) ) ) {
    tab <- lapply( motifs, function( m ) {
      print( m )
      m.info <- out$get.motif.info( m )[[ 1 ]]
      if ( is.null( m.info ) || is.null( m.info$posns ) ) return( NULL )
      m.info <- m.info$posns
      m.info <- subset( m.info, p.value <= p.cutoff )
      if ( nrow( m.info ) <= 0 ) return( NULL )
      ## if ( expand.operons ) {
      ##   g <- as.character( m.info$gene )
      ##   o.hits <- which( sapply( o.list, function( i ) any( i %in% g ) ) )
      ##   o.genes <- unlist( o.list[ o.hits ] ); names( o.genes ) <- NULL
      ##   o.genes <- o.genes[ ! o.genes %in% g ]
      ## }
      k.m <- as.integer( strsplit( m, '_' )[[ 1 ]][ 2:3 ] )
      m.info$bic <- k.m[ 1 ]
      m.info$mot <- k.m[ 2 ]
      m.info
    } )
    dir.create( 'filehash' )
    for ( i in 1:length( tab ) ) if ( ! is.null( tab[[ i ]] ) && ! class( tab[[ i ]] ) == 'try-error' ) {
      write.table( tab[[ i ]], col.names=(i==1), quote=F, sep='\t', file=tmpf, append=(i!=1), row.names=F )
    }
    system( sprintf( "bzip2 -fv %s", tmpf ) )
  }
  cat( "Reading", sprintf( "%s.bz2", tmpf ), "\n" )
  tab <- read.delim( bzfile( sprintf( "%s.bz2", tmpf ) ), sep='\t', head=T )
  as.data.table( tab )
}

## Compute equiv. of 'pssm.scans' table from 'meme.hits' table
get.meme.genome.positions <- function( motifs='ALL', genes=NULL, force=F ) {
  if ( motifs == 'ALL' ) motifs <- unlist( lapply( 1:nrow( out$motif.widths ),
         function( i ) { ii <- which( out$motif.widths[ i, ] > 0 )
                         if ( length( ii ) <= 0 ) return( NULL )
                         paste( 'MOT', i, ii, sep='_' ) } ) ) ## This is faster than the above!

  out.p.cutoff <- 0.001

  tmpf <- sprintf( "./filehash/meme_posns_%s.tsv", as.character( out.p.cutoff ) )
  if ( force || ! file.exists( sprintf( "%s.bz2", tmpf ) ) ) {
    ##    use the MEME hits (not MAST scan hits) to get positions across the genome
    if ( force || ! exists( 'meme.hits' ) ) { ## Meme hits table computed from get.meme.hits.table()
      tab <- get.meme.hits.table( motifs, p.cutoff=out.p.cutoff, force=force )
    } else {
      tab.mots <- paste( 'MOT', meme.hits$bic, meme.hits$mot, sep='_' )
      tab <- subset( meme.hits, paste( 'MOT', bic, mot, sep='_' ) %chin% motifs )
    }
    if ( ! is.null( genes ) ) tab <- subset( tab, gene %chin% genes )

    dir.create( 'filehash' )
    dir.create( sprintf( 'filehash/meme_posns_%s', as.character( out.p.cutoff ) ) )
    
    ##out.tab <- NULL
    all.genes <- sort( unique( as.character( tab$gene ) ) )
    ##for ( g in all.genes ) {
    tmp <- mclapply( all.genes, function( g ) {
      ##print( g )
      tmpf <- sprintf( 'filehash/meme_posns_%s/%s.tsv', as.character( out.p.cutoff ), g )
      if ( file.exists( sprintf( "%s.bz2", tmpf ) ) && file.info( sprintf( "%s.bz2", tmpf ) )[ 'size' ] > 1 ) next
      sites <- subset( tab, gene == g )
      coo <- e$get.gene.coords( g, op.shift=e$operon.shift )
      g.seq <- substr( e$genome.info$genome.seqs[ as.character( coo$contig ) ], coo$start - 1000, coo$end + 1000 )
      ss <- table( as.character( sites$site ) ) ## Do it only for unique copies (to speed it up)
      seqs <- names( ss )
      seqs <- gsub( 'N', '[GATC]', seqs ) ## replace N's with GATCs, then need to use perl=T
      seqs <- gsub( 'X', '[GATC]', seqs ) ## replace X's with GATCs, then need to use perl=T
      locs <- lapply( seqs, gregexpr, g.seq ) ##e$genome.info$genome.seqs[ chr ], perl=T )
      
      seqs <- e$rev.comp( names( ss ) )
      seqs <- gsub( 'N', '[GATC]', seqs ) ## replace N's with GATCs, then need to use perl=T
      seqs <- gsub( 'X', '[GATC]', seqs ) ## replace X's with GATCs, then need to use perl=T
      locs.rev <- lapply( seqs, gregexpr, g.seq ) ##e$genome.info$genome.seqs[ chr ], perl=T )
      genome.posns <- integer()

      coo.st <- coo$start - 1000 - 1
      
      for ( i in 1:length( locs ) ) {
        if ( i %% 100 == 0 ) cat( g, i, "\n" )
        l <- NULL
        if ( ! is.null( locs[[ i ]] ) && attr(locs[[i]][[1]],'match.length') != -1 ) l <- locs[[ i ]]
        else if ( ! is.null( locs.rev[[ i ]] ) && attr(locs.rev[[i]][[1]],'match.length') != -1 ) l <- locs.rev[[ i ]]
        tab.inds <- which( sites$site == names( ss )[ i ] )
        if ( is.null( l ) ) { genome.posns[ tab.inds ] <- NA; next }
        l <- l[[ 1 ]]
        if ( l[ 1 ] == -1 ) { genome.posns[ tab.inds ] <- NA; next }
        for ( j in 1:length( l ) ) genome.posns[ tab.inds ] <- l[ j ] + coo.st ## Store the global start posns of each site
      }
      if ( length( genome.posns ) < nrow( sites ) ) { cat( "WARNING:", length( genome.posns ), nrow( sites ), "\n" ); next }
      sites$genome.posns <- genome.posns
      sites$chr <- coo$contig
      sites$gene.strand <- coo$strand
      ##out.tab <- rbind( out.tab, sites )
      write.table( sites, col.names=(g==all.genes[1]), quote=F, sep='\t', file=tmpf, append=F, row.names=F )
      system( sprintf( "bzip2 -fv %s", tmpf ) )
      NULL
    }, mc.preschedule=F )
    ##as.data.table( out.tab )
    system( sprintf( "find ./filehash/meme_posns_%s/ -name '*.tsv.bz2' -print | sort | xargs bunzip2 -c | bzip2 -c >%s", as.character( out.p.cutoff ),
                    sprintf( "%s.bz2", tmpf ) ) )
  }
  cat( "Reading", sprintf( "%s.bz2", tmpf ), "\n" )
  tab <- read.delim( bzfile( sprintf( "%s.bz2", tmpf ) ), sep='\t', head=T )
  ##tab <- fread( bzfile( sprintf( "%s.bz2", tmpf ) ) )
  as.data.table( tab )
  ##tab
}  

plot.motif.clusters <- function( inds=1:mc.length, ... ) { ##length( motif.clusts ) ) {
  par( mfrow=c( 4, 4 ) )
  for ( i in inds ) {
    print( i )
    tmp <- get.motif.cluster.info( paste( "MOTC", i, sep="_" ) )[[ 1 ]]
    e$viewPssm( attr( tmp, 'combined.pssm' ), main=paste( i, length( attr( tmp, 'mot.names' ) ) ), ... )
  }
}  

plot.motif.instances <- function( motif.cluster, inds=NA, aligned=T, ... ) {
  tmp <- get.motif.cluster.info( motif.cluster )[[ 1 ]]
  if ( aligned ) {
    pssms <- attr( tmp, 'aligned.pssms' )
    if ( ! is.na( inds ) ) pssms <- pssms[ inds ]
  } else {
    mots <- get.motifs( motif.clust=motif.cluster )[[ 1 ]]
    if ( is.na( inds ) ) inds <- 1:length( mots )
    pssms <- lapply( mots[ inds ], function( m ) pssm <- get.motif.info( m )[[ 1 ]]$pssm )
  }
  for ( p in pssms[ inds ] ) e$viewPssm( p, ... )
}

## Use resampling to get significance of peak heights from plot.promoter.architecture()
## Note the height to compare is the top value on the right axis of the top plot of plot.promoter.architecture()
## or run tmp=plot.promoter.architecture() with count.all=T and get it from max(tmp$counts[,'ALL'])
## Default is resample over only noncoding regions.
## Another option (not yet implemented): try resampling over coding/noncoding with same frequency as in input motifs
get.motif.scan.peak.height.significances <- function( motifs='ALL', p.value.cutoff=1e-5, type=c('mast','meme','fimo'),
                                                     motif.filter=NULL, 
                                                     sampling.bg=c('noncoding_only','both_same_as_motifs','both')[2] ) {
  if ( is.null( motifs ) || motifs == 'ALL' ) { ## assume 'ALL' assumes remove.bad=TRUE
    if ( exists( 'coding.fracs' ) )
      motifs <- names( which( ! is.na( coding.fracs$all.fracs ) & coding.fracs$all.fracs < coding.fracs$mean.fracs - 0.01 ) )
    else motifs <- unlist( lapply( 1:nrow( out$motif.widths ), function( i ) {
      ii <- which( out$motif.widths[ i, ] > 0 )
      if ( length( ii ) <= 0 ) return( NULL )
      paste( 'MOT', i, ii, sep='_' ) } ) )
    print( length( motifs ) )
  }

  if ( ! is.null( motif.filter ) ) {
    motifs <- motifs[ motifs %chin% motif.filter ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (3)" )
  }
  
  mots <- strsplit( gsub( "MOT_", "", motifs ), "_" )
  bi <- as.integer( sapply( mots, "[", 1 ) )
  mo <- as.integer( sapply( mots, "[", 2 ) )
  if ( type == 'mast' ) {
    scans <- pssm.scans[ pvals <= p.value.cutoff ]
    scans <- scans[ J( c( bi, bi ), c( mo, -mo ) ), allow.cart=T ] ##subset( pssm.scans, bic %in% bi )
    scans <- scans[ ! is.na( scans$posns ), ]
    setnames( scans, 5, 'Start' )
    scans$mots <- abs( scans$mots )
    scans$Stop <- scans$Start + motif.widths[ cbind( scans$bic, abs( scans$mots ) ) ] - 1
  } else if ( type == 'meme' ) {
    scans <- subset( meme.hits, p.value <= p.value.cutoff )
    scans <- scans[ J( c( bi, bi ), c( mo, -mo ) ), allow.cart=T ]
    scans <- scans[ ! is.na( scans$genome.posns ), ]
    setnames( scans, 2, 'mots' ) ##colnames( scans )[ 2 ] <- 'mots'
    setnames( scans, 8, 'posns' ) ##colnames( scans )[ 8 ] <- 'posns'
    scans$posns <- scans$posns + 1 ## oops - a 1 bp offset needs to be fixed
    scans$mots[ scans$strand == '-' ] <- -scans$mots[ scans$strand == '-' ] ## flip if it's rev comped
    scans$mots[ scans$gene.strand == 'R' ] <- -scans$mots[ scans$gene.strand == 'R' ] ## flip for rev-comp genes too
  } else if ( type == 'fimo' ) {
    scans <- subset( fimo.out, `p-value` <= p.value.cutoff )
    print(dim(scans))
    scans <- scans[ J( bi, mo ), allow.cart=T ]
    print(dim(scans))
    setnames( scans, 2, 'mots' ) ##colnames( scans )[ 2 ] <- 'mots'
    scans <- scans[ ! is.na( Seq ) ]
    scans$mots[ scans$Strand == '-' ] <- -scans$mots[ scans$Strand == '-' ]
    scans$posns <- ifelse( scans$Strand == '+', scans$Start, scans$Stop ) + 1 ## necessary but I dont know why
    setnames( scans, 3, 'gene' ) ##colnames( scans )[ 3 ] <- 'gene'
    setnames( scans, 7, 'pvals' ) ##colnames( scans )[ 7 ] <- 'pvals'
  }
  setkey( scans, 'bic', 'mots' )
  scans <- scans[ ! is.na( gene ) ]

  non.coding.seqs <- lapply( names( e$genome.info$genome.seqs ), function( n )
                        rep( TRUE, nchar( e$genome.info$genome.seqs[ n ] ) ) )
  names( non.coding.seqs ) <- names( e$genome.info$genome.seqs )
  tmp <- subset( e$genome.info$feature.tab, type %in% c( 'CDS', 'rRNA', 'tRNA' ) )
  tmp$where <- as.character( tmp$contig )
  tmp$Start <- as.integer( as.character( tmp$start_pos ) )
  tmp$Stop <- as.integer( as.character( tmp$end_pos ) )
  for ( i in 1:nrow( tmp ) ) {
    ttmp <- tmp[ i, ]
    if ( ! is.na( ttmp$where ) && ! is.na( ttmp$Start ) && ! is.na( ttmp$Stop ) )
      non.coding.seqs[[ ttmp$where ]][ ttmp$Start:ttmp$Stop ] <- FALSE
    rm( ttmp )
  }
  coding.seqs <- lapply( non.coding.seqs, function( i ) which( i == FALSE ) )
  non.coding.seqs <- lapply( non.coding.seqs, function( i ) which( i == TRUE ) )
  names( non.coding.seqs ) <- names( coding.seqs ) <- names( e$genome.info$genome.seqs )

  ## sample positions for each motif from the non coding regions
  counts <- lapply( names( e$genome.info$genome.seqs ), function( n )
                        rep( 0, nchar( e$genome.info$genome.seqs[ n ] ) ) )
  names( counts ) <- names( e$genome.info$genome.seqs )

  for ( chr in names( counts ) ) {
    if ( sampling.bg == 'noncoding_only' ) { ## sample only from noncoding regions
      sc <- scans[ gene == chr ] ##& Start %in% non.coding.seqs[[ chr ]] & Stop %in% non.coding.seqs[[ chr ]] ]
      posns <- sample( non.coding.seqs[[ chr ]], nrow( sc ), replace=T )
    } else if ( sampling.bg == 'both_same_as_motifs' ) { ## sample from noncoding and coding with same frequency as seen in input motifs
      sc <- scans[ gene == chr ]
      n.in <- sum( sc$Start %in% non.coding.seqs[[ chr ]] & sc$Stop %in% non.coding.seqs[[ chr ]], na.rm=T )
      n.out <- nrow(sc) - n.in
      posns <- c( sample( non.coding.seqs[[ chr ]], n.in, replace=T ), sample( coding.seqs[[ chr ]], n.out, replace=T ) )
    } else if ( sampling.bg == 'both' ) { ## sample from noncoding and coding uniformly across genome
      sc <- scans[ gene == chr ]
      posns <- sample.int( max(c(max(coding.seqs[[chr]]),max(non.coding.seqs[[chr]]))), nrow( sc ), replace=T )
    }
    wids <- abs( sc$Start - sc$Stop ) ## motif widths
    cat( chr, dim(sc), length(posns), '\n' )
    for ( i in 1:length( posns ) ) {
      if ( i %% 1000000 == 0 ) cat( i, length(posns), "\n" )
      wid <- wids[ i ]
      pos <- posns[ i ]
      if ( is.na( wid ) || is.na( pos ) ) next
      inds <- pos-1+(1:wid)
      counts[[ chr ]][ inds ] <- counts[[ chr ]][ inds ] + 1
    }
  }

  qout <- unlist( lapply( names( counts ), function( chr ) counts[[ chr ]][ non.coding.seqs[[ chr ]] ] ) )
  print( quantile( qout, c( 0.9, 0.95, 0.99, 0.999, 0.99999 ) ) )
  cat( max( qout ), 1/length( qout ), "\n" )
  qout
}

## gene can be a gene name or coords e.g. c(Chr=34100,34300) to plot a fixed region (window/shift/op.shift are ignored)
##    or c(Chr=34100) in which case window is used to get second coord
## can use e.g. biclust.filter=paste("BIC",which(evals<=10&!is.na(evals)),sep='_')    ## (see gene.gene.network() comments for getting "evals" vector)
##  NOTE: doesn't seem that residual filtering or e.value filtering changes things too much! Huh?
## To only include motifs in motif.clusters, use motif.filter=unlist(motif.clusts)
## To get out the motif scan data as a big genome-size PSSM, do
##    tmp=out$plot.promoter.architecture(c(NC_002937.3=1,3570858),p.val=1e-6,type='fimo',dont.plot=T)
## Add "plot.sig.line=0.99" e.g. to plot line, above which p-value is < 0.1
## "top.only" is synonymous with "for the paper" -- default for paper is p.val=1e-5, motif.filt=unlist(out$motif.clusts[1:out$mc.length]), count.all=T
plot.promoter.architecture <- function( gene, window=125, shift=75, e.value.cutoff=Inf, p.value.cutoff=1e-5,
                                       op.shift=F, include.bad=F, verbose=T, biclust.filter=NULL, motif.filter=NULL,
                                       type=c('mast','meme','fimo'), dont.plot=F, top.only=F, plot.genome.seq=T,
                                       count.all=F, plot.sig.line=NA, ... ) {
  if ( is.character( gene ) ) {
    coo <- e$get.gene.coords( gene, op.shift=op.shift )
    if ( is.null( coo ) ) return( NULL )
    print( coo )
    st.st <- c( coo$start_pos, coo$end_pos )
    if ( coo$strand == "R" ) st.st <- rev( st.st )
    st.st <- st.st[ 1 ] + c( -window, +window ) + ( if ( coo$strand == 'R' ) +shift else -shift )
    chr <- as.character( coo$contig )
    names( st.st )[ 1 ] <- chr
  } else {
    if ( length( gene ) >= 2 ) st.st <- gene
    else if ( length( gene ) == 1 ) st.st <- gene + c( -1, 1 ) * window / 2
  }

  if ( type == 'mast' ) { ## use mast pssm scans
    if ( ! exists( "pssm.scans" ) ) pssm.scans <- get.pssm.scans()
    if ( p.value.cutoff < max( pssm.scans$pvals, na.rm=T ) )
      pssm.scans <- pssm.scans[ pvals <= p.value.cutoff ]
    if ( ! is.data.table( pssm.scans ) ) {
      pssm.scans <- as.data.table( pssm.scans ); gc()
      setkey( pssm.scans, bic, mots, gene, posns )
    }
  } else if ( type == 'meme' ) { ## use meme hits from meme.scores (direct output of cmonkey)
    if ( ! exists( "meme.hits" ) ) {
      meme.hits <- get.meme.genome.positions() ##get.meme.hits.table()
    }
    meme.hits <- subset( meme.hits, p.value <= p.value.cutoff )
    meme.hits <- as.data.table( meme.hits ); gc()
    setkey( meme.hits, bic, mot, gene, genome.posns )
  } else if ( type == 'fimo' ) {
    if ( ! exists( 'fimo.out' ) ) load( 'filehash/fimo_out_1e-05.RData' )
    fimo.out <- subset( fimo.out, `p-value` <= p.value.cutoff )
  }
  if ( st.st[ 1 ] < 1 ) st.st[ 1 ] <- 1
  chr <- names( st.st )[ 1 ]
  if ( is.null( chr ) ) {
    chr <- names( which( sapply( e$genome.info$genome.seqs, nchar ) ==
                        max( sapply( e$genome.info$genome.seqs, nchar ) ) ) )
    names( st.st )[ 1 ] <- chr
  }

  if ( type == 'mast' ) { ## use mast pssm scans
    scans <- pssm.scans[ gene == chr & posns %betw% ( st.st + c( -100, 100 ) ) ]
    scans <- unique( scans )
    motifs <- unique( paste( "MOT", scans$bic, abs( scans$mots ), sep="_" ) )
  } else if ( type == 'meme' ) {
    chr2 <- chr; rm( chr )
    scans <- meme.hits[ chr == chr2 & genome.posns %betw% ( st.st + c( -1000, 1000 ) ) ]
    chr <- chr2; rm( chr2 )
    scans <- unique( scans )
    motifs <- unique( paste( "MOT", scans$bic, abs( scans$mot ), sep="_" ) )
  } else if ( type == 'fimo' ) {
    scans <- fimo.out[ Seq == chr & Start %betw% ( st.st + c( -100, 100 ) ) ]
    scans <- unique( scans )
    motifs <- unique( paste( 'MOT', scans$bic, scans$mot, sep="_" ) )
  }

  if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria (1)!" )
  if ( verbose ) cat( length( motifs ), 'motifs.\n' )

  if ( ! is.null( biclust.filter ) ) {
    motifs <- motifs[ motifs %chin% unlist( get.motifs( biclust=biclust.filter ) ) ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (2)" )
  }

  if ( ! is.null( motif.filter ) ) {
    motifs <- motifs[ motifs %chin% motif.filter ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (3)" )
  }

  if ( ! include.bad && exists( "bad.clusts" ) ) {
    ##bad.mcs <- lapply( motif.clusts, function( i ) coding.fracs$all.fracs[ i ] )
    bad.ms <- unique( unlist( get.motifs( motif.clust=bad.clusts, expand=F ) ) )
    motifs <- motifs[ ! motifs %chin% bad.ms ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (4)" )
  }
  
  if ( ! is.infinite( e.value.cutoff ) && ! is.na( e.value.cutoff ) ) {
    minfo <- get.motif.info( motifs=motifs )
    e.vals <- do.call( c, lapply( minfo, function( tmp ) {
      if ( is.null( tmp ) ) return( NA )
      return( tmp$e.value )
    } ) )
    motifs <- motifs[ e.vals <= e.value.cutoff ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (5)" )
  }

  if ( ! include.bad ) {    
    if ( exists( "coding.fracs" ) ) { ## should have been pre-computed via get.motif.coding.fracs() for all motifs,
      frac.in.coding <- coding.fracs$all.fracs[ motifs ] ## in ensemble.analysis()
    } else {
      coding.fracs <- get.motif.coding.fracs( motifs, verbose=T )
      frac.in.coding <- coding.fracs$all.fracs
    }
    motifs.orig <- motifs
    motifs <- motifs[ ! is.na( frac.in.coding ) & frac.in.coding < coding.fracs$mean.fracs - 0.01 ]
    if ( length( motifs ) <= 0 ) stop( "No motifs pass criteria! (6)" )
    rm( coding.seqs, scans, in.coding )
  }
  
  if ( verbose ) cat( length( motifs ), 'motifs remain.\n' )
  mots <- strsplit( gsub( "MOT_", "", motifs ), "_" )
  bi <- as.integer( sapply( mots, "[", 1 ) )
  mo <- as.integer( sapply( mots, "[", 2 ) )

  if ( type == 'mast' ) {
    scans <- pssm.scans[ J( c( bi, bi ), c( mo, -mo ), chr ), allow.cart=T ] ##subset( pssm.scans, bic %in% bi )
    scans <- scans[ ! is.na( scans$posns ), ]
    scans <- scans[ scans$posns %betw% ( st.st + c( -500, 500 ) ), ]
  } else if ( type == 'meme' ) {
    scans <- meme.hits[ J( c( bi, bi ), c( mo, -mo ) ), allow.cart=T ]
    scans <- scans[ ! is.na( scans$genome.posns ), ]
    scans <- scans[ scans$genome.posns %betw% ( st.st + c( -500, 500 ) ), ]
    setnames( scans, 2, 'mots' ) ##colnames( scans )[ 2 ] <- 'mots'
    setnames( scans, 8, 'posns' ) ##colnames( scans )[ 8 ] <- 'posns'
    scans$posns <- scans$posns + 1 ## oops - a 1 bp offset needs to be fixed
    scans$mots[ scans$strand == '-' ] <- -scans$mots[ scans$strand == '-' ] ## flip if it's rev comped
    scans$mots[ scans$gene.strand == 'R' ] <- -scans$mots[ scans$gene.strand == 'R' ] ## flip for rev-comp genes too
  } else if ( type == 'fimo' ) {
    scans <- fimo.out[ J( bi, mo, chr ), allow.cart=T ]
    scans <- scans[ scans$Start %betw% ( st.st + c( -500, 500 ) ), ]
    setnames( scans, 2, 'mots' ) ##colnames( scans )[ 2 ] <- 'mots'
    scans$mots[ scans$Strand == '-' ] <- -scans$mots[ scans$Strand == '-' ]
    scans$posns <- ifelse( scans$Strand == '+', scans$Start, scans$Stop ) + 1 ## necessary but I dont know why
    setnames( scans, 3, 'gene' ) ##colnames( scans )[ 3 ] <- 'gene'
    setnames( scans, 7, 'pvals' ) ##colnames( scans )[ 7 ] <- 'pvals'
  }
  setkey( scans, 'bic', 'mots' )
  
  getEntropy <- function( pssm ) {
    pssm[ pssm == 0 ] <- 0.00001
    entropy <- apply( pssm, 1, function( i ) -sum( i * log2( i ) ) )
    return( entropy )
  }

  seq <- substr( e$genome.info$genome.seqs[ names( st.st )[ 1 ] ], st.st[ 1 ], st.st[ 2 ] )
  
  mat <- matrix( 0, nrow=diff( st.st ) + 1, ncol=4 )
  rownames( mat ) <- as.character( st.st[ 1 ]:st.st[ 2 ] )
  colnames( mat ) <- e$col.let
  if ( ! dont.plot ) mat2 <- mat * 0

  scans <- scans[ scans$posns %betw% ( st.st + c( -100, 100 ) ) ]

  if ( exists( "motif.clusts" ) ) {
    mcs <- get.motif.clusters( motif=motifs )
    mot.tab <- sort( table( unlist( mcs ) ) )##[ sort( table( unlist( mcs ) ) ) > 2 ]
    mot.tab <- rev( mot.tab[ mot.tab > 2 ] )
    print( mot.tab )
  } else {
    mot.tab <- character()
  }

  if ( length( mot.tab ) > 10 ) mot.tab <- mot.tab[ 1:9 ] ##-(0:9) + length( mot.tab ) ]
  if ( count.all ) mot.tab <- c( ALL=0, mot.tab )
  counts <- matrix( 0, nrow=nrow( mat ), ncol=length( mot.tab ) )
  colnames( counts ) <- names( mot.tab )
  rownames( counts ) <- rownames( mat )
  if ( ! dont.plot ) counts2 <- counts

  if ( nrow( scans ) > 0 ) {
    bics <- unique( scans$bic )
    for ( k in bics ) {
      if ( verbose ) {
        wh <- which( bics == k )
        if ( wh %% 100 == 1 ) cat( k, wh, length( bics ), "\n" )
      }
      sc <- scans[ bic == k ]
      mots <- unique( abs( sc$mots ) )
      for ( m in mots ) {
        width <- motif.widths[ k, m ]
        if ( width <= 0 ) next
        if ( exists( "mcs" ) ) {
          mc <- mcs[[ paste( 'MOT', k, m, sep='_' ) ]]
          mc <- mc[ mc %chin% colnames( counts ) ]
        }
        
        pssm.orig <- get.motif.info( paste( 'MOT', k, m, sep='_' ) )[[ 1 ]]$pssm
        pssm.rev <- pssm.orig[ nrow( pssm.orig ):1, 4:1 ]
        
        entr <- getEntropy( pssm.orig )
        scale.e.orig <- (2 - entr) / 2
        scale.e.rev <- rev( scale.e.orig )
        if ( ! dont.plot ) {
          pssm.orig2 <- pssm.orig * scale.e.orig
          pssm.rev2 <- pssm.rev * scale.e.rev
        }
        sc2 <- sc[ abs( sc$mots ) == m ]

        inds.pssm <- 1:nrow( pssm.orig )        
        for ( i in 1:nrow( sc2 ) ) {
          mot <- sc2$mots[ i ]
          if ( sign( mot ) == 1 ) {
            pssm <- pssm.orig
            if ( ! dont.plot ) {
              pssm2 <- pssm.orig2
              scale.e <- scale.e.orig
            }
          } else if ( sign( mot ) == -1 ) {
            pssm <- pssm.rev
            if ( ! dont.plot ) {
              pssm2 <- pssm.rev2
              scale.e <- scale.e.rev
            }
          }

          posn <- sc2$posns[ i ]          
          ##inds <- as.character( ( posn - 1 ):( posn - 2 + width ) )
          inds <- ( posn - 1 ):( posn - 2 + width )

          inds.1 <- inds.pssm[ inds %betw% st.st ] ##%in% rownames( mat ) ]
          if ( length( inds.1 ) <= 0 ) next
          inds <- inds[ inds %betw% st.st ] ##%in% rownames( mat ) ]
          #inds <- as.character( inds )
          inds <- inds - st.st[1] + 1 ## this instead of the above (commented out) line? is faster
          mat[ inds, ] <- mat[ inds, ] + pssm[ inds.1, ]
          if ( ! dont.plot ) mat2[ inds, ] <- mat2[ inds, ] + pssm2[ inds.1, ]
          if ( count.all ) counts[ inds, 'ALL' ] <- counts[ inds, 'ALL' ] + 1 ## count up all!!
          if ( exists( "mc" ) && length( mc ) > 0 && ! is.null( mc ) ) {
            counts[ inds, mc ] <- counts[ inds, mc ] + 1
            if ( ! dont.plot ) counts2[ inds, mc ] <- counts2[ inds, mc ] + scale.e[ inds.1 ]
          }
        }
      }
    }
  }

  if ( dont.plot ) return( invisible( list( motifs=motifs, scans=scans, mat=mat, counts=counts,
                                           st.st=st.st, mot.tab=mot.tab ) ) )

  qout <- NULL
  if ( ! is.na(plot.sig.line) ) qout <- out$get.motif.scan.peak.height.significances( type=type, p.val=p.value.cutoff, motif.filt=motif.filter, ... )
  if ( ! is.infinite( plot.sig.line ) ) sig.line <- quantile( qout, plot.sig.line )
  else sig.line <- max( qout ) ## max is about p-value = 1e-6
  
  ## Note if plotting to x11 window, these options help with speed. Antialiasing doesnt seem to slow it down:
  ## X11.options(type="dbcairo"); options(X11updates=0.25)
  if ( ! top.only ) par( mfrow=c( 4, 1 ) ) ##, 
  else layout(matrix(c(1,2),2,1),heights=c(.7,.3))
      ##xpd=NA, mai=rep( 0, 4 ), mar=rep( 0, 4 ), mgp=rep( 0, 3 ), oma=rep( 0, 4 ), omd=rep( 0, 4 ), omi=rep( 0, 4 ) )

  yr <- e$viewPssm( mat, scale.e=apply( mat, 1, sum ), no.axis.labels=T, no.render=top.only, ... ) ##, border=NA ) ## can add 'border=NA' arg. to the
  axis( 1, labels=rownames( mat )[ seq( 1, nrow( mat ), by=10 ) ], at=seq( 1, nrow( mat ), by=10 ), line=-0.3 )
  cmax <- max( apply( mat, 1, sum ), na.rm=T ) ##max( counts, na.rm=T ) #* 1.2
  axis( 4, at=c(0,max(yr)/2,max(yr)), labels=as.character(c('0',round(cmax/2),round(cmax))), pos=nrow( mat ) + 1 )

  if ( ! is.na(plot.sig.line) ) {
    sig.line <- max(yr)*sig.line/round(cmax)
    lines( c( 1, nrow(mat) ), rep( sig.line, 2 ), col='red', lty=3, lwd=3 )
  }

  if ( ncol( counts ) > 0 ) {       ## plot.promoter.architecture.2() call.
    counts[ counts <= 0 ] <- NA
    counts2[ counts2 <= 0 ] <- NA
    for ( i in 2:( nrow( counts ) - 1 ) ) {
      counts[ i, is.na( counts[ i, ] ) & ! is.na( counts[ i-1, ] ) & counts[ i-1, ] > 0 ] <- 0
      counts[ i, is.na( counts[ i, ] ) & ! is.na( counts[ i+1, ] ) & counts[ i+1, ] > 0 ] <- 0
      counts2[ i, is.na( counts2[ i, ] ) & ! is.na( counts2[ i-1, ] ) & counts2[ i-1, ] > 0 ] <- 0
      counts2[ i, is.na( counts2[ i, ] ) & ! is.na( counts2[ i+1, ] ) & counts2[ i+1, ] > 0 ] <- 0
    }

    leg.count <- apply( counts, 2, max, na.rm=T )
    if ( exists( "motif.clusts" ) ) leg.count <- mot.tab[ colnames( counts ) ]
    counts.orig <- counts
    if ( top.only ) {
      counts <- counts[ ,leg.count==0 | leg.count >= max(leg.count)/10 ]
      leg.count <- leg.count[ leg.count==0 | leg.count >= max(leg.count)/10 ]
    }

    ##counts1a <- counts * max( yr ) / ( max( counts, na.rm=T ) * 1.2 )
    counts1a <- counts # * max( yr ) / cmax #( max( counts, na.rm=T ) * 1.2 )
    if ( ncol( counts1a ) > 0 ) matlines( counts1a, typ='l', lwd=3, col=1:ncol( counts ), lty=((1:ncol( counts )) %/% 9) + 1 )

    legend( 'topleft', legend=paste( colnames( counts ), leg.count ), lwd=3,
           col=1:ncol( counts ), lty=((1:ncol( counts )) %/% 9) + 1, cex=0.6, horiz=F, trace=F, seg.len=5 )
  }
  ##plot( st.st, c( -2, 2 ), typ='n' ) ##, mar=c(0,0,0,0), xaxs='i', yaxs='i' )
  try( plot.genes.in.region( mean( st.st[1:2] ), chr, diff( range( st.st[1:2] ) ), new=T, yscale=0.5, plot.axis=!top.only ) )

  if ( plot.genome.seq ) {
    mat3 <- mat * 0
    mat3[ cbind( rownames( mat3 ), strsplit( seq, '' )[[ 1 ]] ) ] <- 1
    ##par( xpd=NA )
    e$viewPssm( mat3, scale.e=rep( 0.2, nrow( mat3 ) ), new=F, 
               xoff=min( as.integer( rownames( mat ) ) ) - 1, yoff=-0.1, no.axis.labels=T, ... )
    par( xpd=FALSE )
  }
  
  if ( top.only ) return( invisible( list( motifs=motifs, scans=scans, mat=mat, counts=counts,
                                          st.st=st.st, mot.tab=mot.tab ) ) )

  ##e$viewPssm( mat2, scale.e=apply( mat2, 1, sum ), no.axis.labels=T, ... ) ## OLD
  ## NEW, MAKING IT MORE LIKE A MOTIF LOGO
  yr <- e$viewPssm( mat2 + max(mat2,na.rm=T)/50, scale.e=NA, no.axis.labels=T, min.height=0, ... )
  pssm <- mat2 + 1e-10; pssm <- t( apply( pssm, 1, function( i ) i / sum( i ) ) )
  entropy <- apply( pssm, 1, function( i ) -sum( i * log2( i ) ) )
  bits <- 2 - entropy
  br <- range( bits, na.rm=T ) ## draw labels at 0, 1, 2
  print(max(bits))
  axis( 2, at=c(0,max(yr)/(max(bits)/2)/2, max(yr)/(max(bits)/2)), labels=c('0','1','2'), pos=0, xpd=NA )

  cmax <- max( counts, na.rm=T ) #* 1.2
  counts1a <- counts * max( yr ) / cmax
  if ( ncol( counts1a ) > 0 ) matlines( counts1a, typ='l', lwd=3, col=1:ncol( counts ), lty=(1:ncol( counts ) %/% 8) + 1 )
  axis( 4, at=c(0,max(yr)/2,max(yr)), labels=as.character(c('0',round(cmax/2),round(cmax))), pos=nrow( mat ) + 1 )
  ## END NEW STUFF
  
  yr <- e$viewPssm( mat2 * mat3, scale.e=apply( mat2 * mat3, 1, sum ), no.axis.labels=T, ... ) ##, border=NA )
  if ( ncol( counts2 ) > 0 ) {
    counts2 <- counts * max( yr ) / cmax
    matlines( counts2, typ='l', lwd=3, col=1:ncol( counts ), lty=(1:ncol( counts ) %/% 8) + 1 )
  }
  axis( 1, labels=rownames( mat )[ seq( 1, nrow( mat ), by=10 ) ], at=seq( 1, nrow( mat ), by=10 ), line=-0.3 )
  
  invisible( list( motifs=motifs, scans=scans, mat=mat, mat2=mat2, mat3=mat3, counts=counts, st.st=st.st,
                  mot.tab=mot.tab ) )
}

## Note, this can be done on all motifs, or a specific set. E.g. to do it on a set of motifs that hit a
##   specific region of the genome:
## qqq=out$plot.promoter.architecture('MPN332',verbose=T)
## locator(1) ## get coord - returns 140 (local coord)
## rownames(qqq$mat)[140]   ## convert to genome coord - 390267
## m=out$get.motifs(position=390267)
## qqq=out$cluster.motifs('mcl',m,mcl.I=6.0,min.gene=NA)
## OR:
## qqq=out$cluster.motifs('hclust',m,min.gene=NA,k.cut=0.99)
## qqq=out$cluster.motifs('hclust',m,min.gene=1,expand=T,k.cut=0.99)
## e$viewPssm(attr(qqq$tt.out2[[1]],'combined.pssm'))
## par(mfrow=c(4,4));for(i in 1:16)e$viewPssm(attr(qqq$tt.out2[[1]],'aligned.pssms')[[i]])

cluster.motifs <- function( cluster.option, motifs='ALL', min.gene.overlap=1, ## NA for all-vs-all motif comparison
                           e.value.cutoff=Inf, p.value.cutoff=0.001, resid.cutoff=Inf, n.cutoff=10,
                           expand=F, include.bad=F, find.bad=NA, #get.motif.ccs=NA,
                           ##parallelize.on.compute.cluster=T,
                           in.tt.out=NULL, mcl.I=3.0,
                           mcl.cmd='./progs/mcl-10-201/local/bin/mcl mcltemp2 -o %s --abc -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi 4.0',
                           ... ) {
  if ( file.exists( "cmonkey-motif-other.R" ) ) {
    try( sys.source( "cmonkey-motif-other.R", envir=e ) )
    try( sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e ) )
  }
  out2 <- list()
  in.args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments
               sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## nifty trick, eh?
  out2$args <- in.args
  out2$e <- e
  if ( length( e$meme.scores ) > 1 ) { ## HACK - combined different motif types into one (assumes all are promoters)
    e$meme.scores.orig <- e$meme.scores
    e$meme.scores <- list(); e$meme.scores[[ names( e$meme.scores.orig )[ 1 ] ]] <- list()
    for ( k in 1:e$k.clust ) {
      for ( i in names( e$meme.scores.orig ) ) {
        if ( length( e$meme.scores.orig[[ i ]] ) >= k && ! is.null( e$meme.scores.orig[[ i ]][[ k ]] ) ) {
          e$meme.scores[[ 1 ]][[ k ]] <- e$meme.scores.orig[[ i ]][[ k ]]
          next
        }
      }
    }
  }

  save.touts <- TRUE
  motifs.all <- FALSE
  if ( is.null( motifs ) || motifs == 'ALL' ) {
    motifs.all <- TRUE
    if ( is.na( find.bad ) ) find.bad <- TRUE
    #if ( is.na( get.motif.ccs ) ) get.motif.ccs <- TRUE
    if ( include.bad == FALSE && exists( 'coding.fracs' ) )
      motifs <- names( which( ! is.na( coding.fracs$all.fracs ) & coding.fracs$all.fracs < coding.fracs$mean.fracs - 0.01 ) )
    else motifs <- unlist( lapply( 1:nrow( out$motif.widths ), function( i ) {
      ii <- which( out$motif.widths[ i, ] > 0 )
      if ( length( ii ) <= 0 ) return( NULL )
      paste( 'MOT', i, ii, sep='_' ) } ) )
  } else {
    if ( is.na( find.bad ) ) find.bad <- FALSE
    #if ( is.na( get.motif.ccs ) ) get.motif.ccs <- FALSE
    save.touts <- FALSE
  }

  tt.out <- NULL
  cat( "Getting alignments for", length( motifs ), "motifs...\n" )
  ks <- as.integer( sapply( strsplit( motifs, '_', ), '[', 2 ) )
  mots <- as.integer( sapply( strsplit( motifs, '_', ), '[', 3 ) )
  if ( ! is.null( in.tt.out ) ) {
    tt.out <- in.tt.out ## Input tomtom data frame
  } else if ( file.exists( 'filehash/tt.out.RData' ) ) {
    print( 'Loading pre-computed tt.out data frame.' )
    if ( ! exists( 'tt.out' ) ) load( 'filehash/tt.out.RData' )
    if ( ! motifs.all ) {
      m1 <- paste( 'MOT', tt.out$biclust1, tt.out$motif1, sep='_' )
      m2 <- paste( 'MOT', tt.out$biclust2, tt.out$motif2, sep='_' )
      if ( ! expand ) tt.out <- subset( tt.out, m1 %in% motifs & m2 %in% motifs )
      else tt.out <- subset( tt.out, m1 %in% motifs | m2 %in% motifs )
    }
  } else {
    save.touts.file <- FALSE
    if ( save.touts ) save.touts.file <- 'filehash/touts'
    ##sprintf( 'filehash/touts_%08d',ks.i )##%s',paste(paste(range(ks.i),collapse='_'),paste(range(ks.j),collapse='_'),sep='_') )
    print( save.touts.file )
    if ( save.touts ) {
      ##if ( parallelize.on.compute.cluster && file.exists( save.touts.file ) ) next
      ##else
      dir.create( save.touts.file )
    }
    if ( motifs.all ) { ## For each bicluster k, compare vs. ks less than that one; that way we can add more biclusters
      ##  to the ensemble and only compare the new ones against the previously existing ones.
      ##seqs <- round( seq( 1, length( ks ), length=10 ) ) ## make sure this makes the # per tomtom run about <16,000
      ##for ( i in seqs[ -10 ] ) {
        ##ks.i <- i:min( length( ks ), seqs[ which( seqs == i ) + 1 ] )
        ##for ( j in seqs[ which( seqs == i ):( length( seqs ) - 1 ) ] ) {
          ##ks.j <- j:min( length( ks ), seqs[ which( seqs == j ) + 1 ] )
      ##for ( ks.i in ks ) {
      tmp <- mclapply( unique( ks ), function( ks.i ) {
        ks.j <- unique( ks )[ unique( ks ) < ks.i ]
        ## want at lst 100 compari motifs
        if ( length( ks.j ) < 100 ) ks.j <- unique( c( ks.j, sample( unique( ks )[ unique( ks ) > ks.i ], 100-length(ks.j) ) ) ) 
        cat("TOMTOM-ing:", ks.i, 'VS.', ks.j, "\n")
        ##tmp <-
        e$motif.similarities.tomtom( ks[ ks%in%ks.i ], ks[ ks%in%ks.j ], mots[ ks%in%ks.i ], mots[ ks%in%ks.j ],
                                    min.gene.overlap=min.gene.overlap,
                                    save.touts=save.touts.file, desymm=F, e.value.cutoff=e.value.cutoff, ... )
        ##if ( save.touts && ! parallelize.on.compute.cluster ) tt.out <- rbind( tt.out, tmp )
        ##if ( ! save.touts ) tt.out <- rbind( tt.out, tmp )
        ##rm( tmp ); gc()
        NULL
      }, mc.preschedule=F )
      ##  if ( save.touts && ! parallelize.on.compute.cluster ) save( tt.out, file='filehash/tt.out.RData' )
      ##}
      ##if ( save.touts && ! parallelize.on.compute.cluster ) save( tt.out, file='filehash/tt.out.RData' )
    } else {
      ks <- as.integer( sapply( strsplit( sort( motifs ), '_' ), '[', 2 ) )
      mots <- as.integer( sapply( strsplit( sort( motifs ), '_' ), '[', 3 ) )
      cat(length(ks),length(mots),"\n")
      tt.out <- e$motif.similarities.tomtom( ks, ks, mots, mots, min.gene.overlap=min.gene.overlap,
                                            save.touts=F, desymm=F, e.value.cutoff=e.value.cutoff, ... )
    }
  }

  if ( is.null( tt.out ) ) {
    ##system( 'bunzip2 -c filehash/touts_*_*_*_*/*_tout.tsv.bz2 | bzip2 -c >filehash/tt_outs.tsv.bz2' )
    ## Skip the first (header) lines of each file.
    ##system( "find ./filehash/touts_*_*_*_*/ -name '*_tout.tsv.bz2' -print | sort | parallel '(bunzip2 -vc {} | tail -n +2)' |bzip2 -c > filehash/tt_outs.tsv.bz2" )
    #system( sprintf( "find %s -name '*_tout.tsv.bz2' -print | sort | parallel '(bunzip2 -vc {} | tail -n +2)' |bzip2 -c > filehash/tt_outs.tsv.bz2",
    #                save.touts.file ) )
    system( sprintf( "bunzip2 -vc %s/*_tout.tsv.bz2 | grep -v '#Query ID' | bzip2 -c >filehash/tt_outs.tsv.bz2", save.touts.file ) )
    tt.out <- read.delim( bzfile( 'filehash/tt_outs.tsv.bz2' ), head=F )
    tmp <- system( sprintf( "find %s -name '*_tout.tsv.bz2' -print | head -1", save.touts.file ), intern=T )
                  ##tmp <- system(          "find ./filehash/touts_*_*_*_*/ -name '*_tout.tsv.bz2' -print | head -1", intern=T )
    tmp <- read.delim( bzfile( tmp ), sep='\t', head=T, check.names=F )
    colnames( tt.out ) <- colnames( tmp ); rm( tmp )
    tout <- tt.out; rm( tt.out )
    colnames( tout ) <- gsub( '#', '', colnames( tout ) )
    q.id <- t( sapply( strsplit( as.character( tout[ ,1 ] ), "_", fixed=T ), '[', 2:5 ) ) ##function( i ) i[ 2:5 ] )
    ##q.res.ev <- do.call( rbind, strsplit( as.character( tout[ ,1 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) )
    t.id <- t( sapply( strsplit( as.character( tout[ ,2 ], fixed=T ), "_" ), '[', 2:5 ) )
    ##t.res.ev <- do.call( rbind, strsplit( as.character( tout[ ,2 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) )
    tout2 <- data.frame( biclust1=as.integer( q.id[ ,1 ] ), motif1=as.integer( q.id[ ,2 ] ),
                        resid1=as.numeric( q.id[ ,3 ] ), e.value1=as.numeric( q.id[ ,4 ] ),
                        biclust2=as.integer( t.id[ ,1 ] ), motif2=as.integer( t.id[ ,2 ] ),
                        resid2=as.numeric( t.id[ ,3 ] ), e.value2=as.numeric( t.id[ ,4 ] ),
                        offset=as.integer( as.character( tout[[ 'Optimal offset' ]] ) ),
                        p.value=as.numeric( as.character( tout[[ 'p-value' ]] ) ),
                        q.value=as.numeric( as.character( tout[[ 'q-value' ]] ) ),
                        overlap=as.integer( as.character( tout[[ 'Overlap' ]] ) ),
                        consensus1=as.factor( as.character( tout[[ 'Query consensus' ]] ) ),
                        consensus2=as.factor( ifelse( tout[[ 'Orientation' ]] == "-",
                          e$rev.comp( as.character( tout[[ 'Target consensus' ]] ) ),
                          as.character( tout[[ 'Target consensus' ]] ) ) ),
                        orientation=tout[[ 'Orientation' ]] ) ##,
    tout2 <- subset( tout2, ! ( biclust1 == biclust2 & motif1 == motif2 ) )
    tout2 <- tout2[ order( tout2$q.value, tout2$p.value ), ]
    tt.out <- tout2; rm( tout, tout2 )
  }
  cat( nrow( tt.out ), 'total alignments.\n' )
##  tt.out <- subset( tt.out, bic == 
  ##     if ( ! is.na( min.gene.overlap ) ) {
  ##       cat( "Removing hits between motifs in clusters with <", min.gene.overlap, "overlapping genes\n" )
  ##       g1 <- sapply( as.integer( as.character( tt.out$biclust1 ) ), e$get.rows )
  ##       g2 <- sapply( as.integer( as.character( tt.out$biclust2 ) ), e$get.rows )
  ##       n.ov <- sapply( 1:length( g1 ), function( i ) sum( g1[[ i ]] %in% g2[[ i ]] ) )
  ##       tt.out <- subset( tt.out, n.ov >= min.gene.overlap )
  ##     }
  out2$tt.out <- tt.out
  if ( cluster.option != 'mcl' ) {
    e$parallel.cores <- 1
    tt.out2 <- e$cluster.tomtom.results( tt.out, e.value.cutoff=e.value.cutoff, p.value.cutoff=p.value.cutoff,
                                        resid.cutoff=resid.cutoff, n.cutoff=n.cutoff, ... )
    tmp <- tt.out2[ sapply( tt.out2, function( i ) length( attr( i, "mot.names" ) ) ) > 0 ]
    attributes( tmp ) <- attributes( tt.out2 )
    out2$tt.out2 <- tmp; rm( tmp, tt.out2 )
  } else { ## TODO: use mcl to do the clustering...? see cmPostProc2_mot_metaclustering2.R
    ##tmp <- subset( tt.out, p.value <= 0.01 ) ##p.value <= 0.01 )
    out2$mcl.I <- mcl.I
    if ( is.na( e.value.cutoff ) ) e.value.cutoff <- Inf
    if ( is.na( resid.cutoff ) ) resid.cutoff <- Inf
    if ( is.na( p.value.cutoff ) ) p.value.cutoff <- Inf
    if ( ! is.infinite( e.value.cutoff ) || ! is.infinite( resid.cutoff ) || ! is.infinite( p.value.cutoff ) )
      tmp <- subset( tt.out, p.value <= p.value.cutoff & e.value1 <= e.value.cutoff & e.value2 <= e.value.cutoff &
                    resid1 <= resid.cutoff & resid2 <= resid.cutoff )
    cat( nrow( tmp ), "alignments pass the cutoff filters.\n" )
    n1 <- paste( as.integer( as.character( tmp$biclust1 ) ), tmp$motif1, sep="_" )
    n2 <- paste( as.integer( as.character( tmp$biclust2 ) ), tmp$motif2, sep="_" )
    pv <- log10( tmp$p.value )
    rm( tmp ); gc()
    require( Matrix )
    rnames <- 1:length( unique( c( n1, n2 ) ) )
    names( rnames ) <- unique( c( n1, n2 ) )
    xx <- Matrix( data=0, nrow=length( rnames ), ncol=length( rnames ), sparse=T )
    xx[ cbind( rnames[ n1 ], rnames[ n2 ] ) ] <- pv
    rm( n1, n2, pv )

    tmp <- which( xx != 0, arr=T )
    ##w <- 1 - 10^x[ tmp ] ##; rm( x ); gc()
    require( igraph0 )
    gr <- graph.edgelist( tmp, directed=F ); rm( tmp ); gc()
    gr <- simplify( gr ); gc()
    ## Look into MCL http://www.micans.org/mcl/ for clustering network!!!
    ## I=4.5 worked best for Halo egrin so we'll use it for default.
    ## Actually, I=4.0 works better for the 10-cmonkey-run ecoli egrin
    try( unlink( "mcltemp*" ) ); write.graph( gr, file="mcltemp", format="ncol" ); rm( gr ); gc()
    ## NOTE: look in to hierarchical mcl:
    ##   ./progs/mcl-10-201/local/bin/mclcm zzz --subcluster --mplex y -a "-I 1.2" -- "--abc"
    ##   ./progs/mcl-10-201/local/bin/mcxdump -imx-tree mcl.cone --newick -o NEWICKFILE
    ## or (?)   ./progs/mcl-10-201/local/bin/clm order mcl.cone
    system( "cp mcltemp mcltemp2; awk '{print $2,$1}' mcltemp >>mcltemp2" ) ## Need to symmetrize it!
    outfile <- sprintf( "out.mcltemp2.I%s", gsub('.','',sprintf("%.1f",mcl.I),fixed=T) )
    ##cmd <- sprintf( "./progs/mcl-10-201/local/bin/mcl mcltemp2 -o %s --abc -I %.1f -v all -te 3 -S 200000 -scheme 6 --analyze=y -pi 2.0", mcl.I )
    cmd <- sprintf( mcl.cmd, outfile, mcl.I )
    out2$mcl.cmd <- cmd
    print( cmd )
    system( cmd )
    system( sprintf( 'gzip -fv %s', outfile ) )
    if ( save.touts ) {
      system( sprintf( 'mv -v %s.gz filehash/', outfile ) )
      outfile <- sprintf( 'filehash/%s', outfile )
    }
    clusts <- lapply( strsplit( readLines( gzfile( sprintf( "%s.gz", outfile ) ) ), "\t" ), as.integer )
    cat( "GOT", length( clusts ), "motif clusters with", length( unlist( clusts ) ), "motifs.\n" )
    clusts <- clusts[ sapply( clusts, length ) >= 3 ] ## eliminate tiny clusters
    clusts <- lapply( clusts, function( i ) paste( "MOT", names( rnames )[ which( rnames %in% i ) ], sep="_" ) )
    cat( "GOT", sum( sapply( clusts, length ) >= 3 ), "motif clusters (length > 3) with",
        length( unlist( clusts[ sapply( clusts, length ) >= 3 ] ) ), "motifs.\n" )
    cat( "GOT", sum( sapply( clusts, length ) >= 10 ), "motif clusters (length > 10) with",
        length( unlist( clusts[ sapply( clusts, length ) >= 10 ] ) ), "motifs.\n" )
    out2$motif.clusts <- clusts
    ##out2$mc.length <- length( clusts )
    mc.length <- out2$mc.length <- max( which( sapply( clusts, length ) >= n.cutoff ) ) ## ignore small clusters
    if ( is.infinite( mc.length ) ) mc.length <- out2$mc.length <- max( which( sapply( clusts, length ) >= 3 ) )

    tmp <- subset( tt.out, p.value <= p.value.cutoff & e.value1 <= e.value.cutoff & e.value2 <= e.value.cutoff &
                  resid1 <= resid.cutoff & resid2 <= resid.cutoff )
    m1 <- paste( as.integer( as.character( tmp$biclust1 ) ), tmp$motif1, sep="_" )
    m2 <- paste( as.integer( as.character( tmp$biclust2 ) ), tmp$motif2, sep="_" )
    tt.out2 <- list()
    for ( i in 1:mc.length ) { ##( clusts ) ) {
      print( i )
      cl <- clusts[[ i ]]
      ccl <- gsub( "MOT_", "", cl )
      tttt <- subset( tmp, m1 %in% ccl & m2 %in% ccl )
      if ( nrow( tttt ) <= 1 ) next
      ## This actually doesnt cluster them since k.cut=1, but computes the aligned motifs pretty well enough...
      tttt.out <- e$cluster.tomtom.results( tttt, p.cutoff=Inf, k.cut=1, return.aligned.pssms=TRUE )
      tt.out2[[ i ]] <- tttt.out[[ 1 ]]
    }
    out2$tt.out2 <- tt.out2; rm( tttt.out, tmp, tttt, tt.out2 ); gc()
  }


  if ( find.bad ) { ## ID motif.clusters whose most of their motifs are in coding regions ==> BAD
    ## motif.clusts <- out2$motif.clusts ##get.motifs( motif.clust=paste( "MOTC", 1:length( out2$tt.out2 ), sep="_" ) )
    ## chr <- names( e$genome.info$genome.seqs )[ 1 ]
    ## coding.seqs <- lapply( names( e$genome.info$genome.seqs ), function( n )
    ##                       rep( FALSE, nchar( e$genome.info$genome.seqs[ n ] ) ) )
    ## names( coding.seqs ) <- names( e$genome.info$genome.seqs )
    ## tmp <- subset( e$genome.info$feature.tab, type %in% c( 'CDS', 'tRNA', 'rRNA' ) )
    ## tmp$contig <- as.character( tmp$contig )
    ## tmp$start_pos <- as.integer( as.character( tmp$start_pos ) )
    ## tmp$end_pos <- as.integer( as.character( tmp$end_pos ) )
    ## for ( i in 1:nrow( tmp ) ) {
    ##   ttmp <- tmp[ i, ]
    ##   if ( ! is.na( ttmp$contig ) && ! is.na( ttmp$start_pos ) && ! is.na( ttmp$end_pos ) )
    ##     coding.seqs[[ ttmp$contig ]][ ttmp$start_pos:ttmp$end_pos ] <- TRUE
    ## }

    ##fracs <- list()

    bad.clusts <- character()
    all.mn <- out$coding.fracs$mean.frac
    for ( i in 1:out2$mc.length ) { ##length( motif.clusts ) ) {
      mn <- mean( out$coding.fracs$all.fracs[ out2$motif.clusts[[ i ]] ], na.rm=T )
      sd <- sd( out$coding.fracs$all.fracs[ out2$motif.clusts[[ i ]] ], na.rm=T )
      if ( is.na( mn ) || is.na( sd ) || mn + sd * 2 > all.mn ) bad.clusts <- c( bad.clusts, paste( 'MOTC', i, sep='_' ) )
      cat( i, paste( 'MOTC', i, sep='_' ) %in% bad.clusts, "\n" )
    ##   motc <- paste( "MOTC", i, sep="_" )
    ##   if ( length( fracs ) >= i && ! is.null( fracs[[ i ]] ) &&
    ##       ( is.null( fracs[[ i ]]$pval ) || ! is.na( fracs[[ i ]]$pval ) ) ) next
    ##   print( motc )
    ##   motifs <- unlist( out2$motif.clusts[ i ] ) ##get.motifs( motif.clust=motc ) )
    ##   if ( length( motifs ) <= 3 ) next
    ##   posns <- posns.orig <- get.motif.positions( motifs, seqs=e$genome.info$genome.seqs )
    ##   if ( length( posns ) <= 1 ) next
    ##   locs <- which( posns > 0 ) 
    ##   posns <- posns[ locs ]
    ##   tmp <- unlist( sapply( 1:length( posns ), function( i ) rep( locs[ i ], posns[ i ] ) ) )
    ##   widths <- motif.widths[ t( sapply( strsplit( motifs, "_" ), function( i ) as.integer( i[ 2:3 ] ) ) ) ]
      
    ##   dens <- density( c( 0, tmp, max( tmp ) + 10 ), bw=mean( widths ) / 8, n=max( locs ) + 10 )
    ##   lower <- c( diff( dens$y ), NA )
    ##   upper <- c( NA, diff( dens$y ) )
    ##   deriv.zero <- abs( lower + upper ) < sd( lower + upper, na.rm=T )
    ##   deriv.2nd.neg <- lower < 0 & upper > 0
    ##   good <- deriv.zero & deriv.2nd.neg & dens$y * length( tmp ) > 10 & ! is.na( lower ) & ! is.na( upper )
    ##   if ( sum( good ) <= 1 )
    ##     good <- deriv.zero & deriv.2nd.neg & dens$y * length( tmp ) > 4.9 & ! is.na( lower ) & ! is.na( upper )
      
    ##   pks <- dens$x[ good ]
    ##   pk.ht <- dens$y[ good ] * length( tmp )
    ##   if ( length( pks ) <= 0 ) next
    ##   ttmp <- coding.seqs[[ 1 ]][ round( pks ) ]
    ##   fracs[[ i ]] <- list( pks=pks, means=mean( ttmp ) )
    ##   cat( i, length( pks ), fracs[[ i ]]$means, "\n" )
    ## }
    
    ## cs <- unlist( coding.seqs )
    ## for ( i in 1:length( fracs ) ) {
    ##   tmp <- mclapply( 1:1000, function( k ) mean( sample( cs, length( fracs[[ i ]]$pks ) ) ) )
    ##   fracs[[ i ]]$pval <- mean( tmp <= fracs[[ i ]]$mean )
    ##   cat( i, length( fracs[[ i ]]$pks ), fracs[[ i ]]$mean, fracs[[ i ]]$pval, "\n" )
    }
    ##bad.clusts <- which( sapply( fracs, function( i ) is.null( i ) || i[ 'Observed In', ] < 0.1 ) )
    ##bad.clusts <- which( sapply( fracs, function( i ) is.na( i$pval ) || i$pval > 0.1 ) )
    ##bad.clusts <- paste( 'MOTC', bad.clusts, sep='_' )
    ##names( bad.clusts ) <- NULL
    print( bad.clusts )
    out2$bad.clusts <- bad.clusts
  }      
    
  out2
}  

get.motif.cluster.clusters <- function( cutoff=0.3, e.value.cutoff=Inf, resid.cutoff=Inf, p.value.cutoff=Inf ) {
  ## re-cluster the motif clusters into motif.cluster.clusters to find 
  ## "stragglers" and smaller mot.clusters that werent added to any bigger mot.cluster
  ## but probably should be, (same as in original Halo ensemble)
  if ( is.na( e.value.cutoff ) ) e.value.cutoff <- Inf
  if ( is.na( resid.cutoff ) ) resid.cutoff <- Inf
  if ( is.na( p.value.cutoff ) ) p.value.cutoff <- Inf
  tmp <- tt.out
  if ( ! is.infinite( e.value.cutoff ) || ! is.infinite( resid.cutoff ) || ! is.infinite( p.value.cutoff ) )
    tmp <- subset( tt.out, e.value1 <= e.value.cutoff & e.value2 <= e.value.cutoff &
                  resid1 <= resid.cutoff & resid2 <= resid.cutoff & p.value <= p.value.cutoff )
  n1 <- paste( as.integer( as.character( tmp$biclust1 ) ), tmp$motif1, sep="_" )
  n2 <- paste( as.integer( as.character( tmp$biclust2 ) ), tmp$motif2, sep="_" )
  pv <- log10( tmp$p.value )
  rm( tmp ); gc()
  require( Matrix )
  rnames <- 1:length( unique( c( n1, n2 ) ) )
  names( rnames ) <- unique( c( n1, n2 ) )
  xx <- Matrix( data=0, nrow=length( rnames ), ncol=length( rnames ), sparse=T )
  xx[ cbind( rnames[ n1 ], rnames[ n2 ] ) ] <- pv
  rm( n1, n2, pv )

  ttmp <- motif.clusts ##get.motifs( motif.clust=paste( "MOTC", 1:length( out$tt.out2 ), sep="_" ) )
  ov.means <- Matrix( data=0, nrow=length( ttmp ), ncol=length( ttmp ), sparse=T )
  ##cutoff <- 0.30
  for ( i in 1:( mc.length - 1 ) ) { ##length( ttmp ) - 1 ) ) {
    m.i <- rnames[ gsub( "MOT_", "", ttmp[[ i ]] ) ]
    m.i <- m.i[ ! is.na( m.i ) ]
    if ( length( m.i ) <= 0 ) next
    x.tmp <- xx[ rnames[ m.i ], ]
    for ( j in (i+1):length( ttmp ) ) {
      m.j <- rnames[ gsub( "MOT_", "", ttmp[[ j ]] ) ]
      m.j <- m.j[ ! is.na( m.j ) ]
      if ( length( m.j ) <= 0 ) next
      xx.tmp <- x.tmp[ ,rnames[ m.j ] ] ##xx[ rnames[ m.i ], rnames[ m.j ] ]
      ov.means[ i, j ] <- mean( xx.tmp != 0 )
    }
    cat( "Similar:", i, max( ov.means[ i, ] ), which( ov.means[ i, ] > cutoff ), "\n" )
  }

  ov.means <- as.matrix( ov.means )
  ov.means[ lower.tri( ov.means ) ] <- t( ov.means )[ lower.tri( t( ov.means ) ) ]
  motif.cluster.clusters <- rep( 0, ncol( ov.means ) )
  if ( ! exists( 'mc.length' ) ) mc.length <- max( which( sapply( ttmp, length ) > 3 ) ) ## should be stored in "out" already
  if ( ! any( ov.means > cutoff ) ) cutoff <- cutoff - 0.05 ##0.25
  ind <- 1
  for ( i in 1:mc.length ) {
    tt <- which( ov.means[ i, ] > cutoff )
    ttmp2 <- tt; ttmp2 <- ttmp2[ ttmp2 != i ]
    if ( length( ttmp2 ) <= 0 ) next
    best.match <- ttmp2[ which.max( ov.means[ i, ttmp2 ] ) ]
    if ( motif.cluster.clusters[ i ] != 0 ) {
      motif.cluster.clusters[ tt ] <- motif.cluster.clusters[ i ]
    } else if ( best.match <= mc.length && motif.cluster.clusters[ best.match ] != 0 ) {
      motif.cluster.clusters[ tt ] <- motif.cluster.clusters[ best.match ]
    } else {
      motif.cluster.clusters[ tt ] <- ind
      ind <- ind + 1
    }
    cat( i, tt, motif.cluster.clusters[ i ], "\n" )
  }
  out2 <- new.env()
  out2$ov.means <- ov.means
  out2$motif.cluster.clusters <- motif.cluster.clusters
  m=c(which(motif.cluster.clusters[1:mc.length]==0),which(motif.cluster.clusters!=0))
  cat( "Total unique motif clusters:", length(table(motif.cluster.clusters[m]))+table(motif.cluster.clusters[m])['0']-1, "\n" )
  out2
}

if ( FALSE ) {
  total.cond.count=table(unlist(lapply(e$clusterStack,'[[','cols')))
  cc=lapply(1:out$mc.length,function(mc){
    mc=paste("MOTC",mc,sep='_')
    print(mc)
    conds=table(unlist(out$get.conditions(biclust=unlist(out$get.biclusters(motif=unlist(out$get.motifs(motif.clust=mc)))))))
    conds=conds[names(total.cond.count)]
    conds[is.na(conds)]<-0
    names(conds)<-names(total.cond.count)
    conds=conds/total.cond.count
    conds
  } )
  cors=cor(t(do.call(rbind,cc)));diag(cors)<-NA
}

try( source( "cmonkey-ensemble-funcs2.R" ) )
