## Updated motif clustering (faster) that uses 'motif shadows'

new.cluster.motifs <- function( motifs='ALL', p.cutoff=1e-6, include.bad=F, n.cutoff=10,
                               mcl.I=1.2, mcl.pi=2.0, distance.weight.cutoff=0.999999, 
                               mcl.cmd=paste( './progs/mcl-10-201/local/bin/mcl new.mcltemp -o %s --abc',
                                            ' -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1' ),
                               meme.the.clusters=T, improve.the.clusters=F, plot.them=F, in.clusts=NULL ) {
  
  in.args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments
               sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## nifty trick, eh?

  nc <- sapply( e$genome.info$genome.seqs, nchar )
  if ( ! is.null( in.clusts ) ) {
  clusts <- in.clusts
  } else {
  if ( exists( 'fimo.out' ) ) p.scans <- subset( fimo.out, `p-value` <= p.cutoff )
  else if ( exists( 'pssm.scans' ) ) p.scans <- subset( pssm.scans, pvals <= p.cutoff )
  else {
    load( sprintf( "filehash/fimo_out_%s.RData", as.character( p.cutoff ) ) )
    p.scans <- subset( fimo.out, `p-value` <= p.cutoff )
  }
  if ( motifs == 'ALL' ) {
    motifs <- unlist( lapply( 1:nrow( motif.widths ), function( i ) {
      ii <- which( motif.widths[ i, ] > 0 )
      if ( length( ii ) <= 0 ) return( NULL )
      paste( 'MOT', i, ii, sep='_' ) } ) )
  }
  cat( length( motifs ), "motifs\n" )

  ## Load unfortunately-named "m" list of motif shadow overlaps
  if ( file.exists( sprintf( "filehash/new_motif_shadows_%s.RData", as.character( p.cutoff ) ) ) ) {
    cat( "Loading motif distances pre-computed...\n" )
    load( sprintf( "filehash/new_motif_shadows_%s.RData", as.character( p.cutoff ) ) ) ## loads 'm'
  } else { ## Idea is to do it for ALL motifs, save those, and subset those results.
    m.tmp <- lapply( nc, function( i ) rep( FALSE, length=i ) ); names( m.tmp ) <- names( nc )
    m <- mclapply( 1:length( motifs ), function( i ) {
      bi <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 2 ])
      mo <- as.integer( strsplit( motifs[i], '_' )[[ 1 ]][ 3 ])
      ##wi <- motif.widths[ bi, mo ]
      scans <- p.scans[ J( c( bi, bi ), c( mo, -mo ) ) ]
      if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
      scans <- scans[ ! is.na( scans$Start ), ] ##posns ), ]
      if ( nrow(scans) <= 0 ) return( lapply( m.tmp, function( i ) integer(0) ) ) ##Matrix( m.tmp ) )
      ##ii <- 0:(wi-1)
      for ( iii in 1:nrow( scans ) ) { ##posn in scans$posns ) { ##1:nrow( scans ) ) {
        ##posn <- scans$posns[ iii ]
        inds <- scans$Start[ iii ]:scans$Stop[ iii ] ##posn + ii
        chr <- levels( scans$Seq )[ scans$Seq[ iii ] ] ## gene
        ##if ( posn < 20 || posn > nc-wi-10 ) inds <- inds[ inds > 0 & inds <= nc[ chr ] ]
        m.tmp[[ chr ]][ inds ] <- TRUE ##scans$posns[j]:(scans$posns[j]+wi-1),1] = 1
      }
      cat( i, length( motifs ), motifs[ i ], nrow( scans ), "\n" ) ##sum( m.tmp ), "\n" )
      lapply( m.tmp, which ) ##which( m.tmp ) ##Matrix( m.tmp )
    }, mc.preschedule=F )
    rm( m.tmp, p.scans ); gc()
    
    names( m ) <- motifs
    save( m, file=sprintf( "filehash/new_motif_shadows_%s.RData", as.character( p.cutoff ) ) )
    ##posn.hit <- rep( 0, nc ) ## Idea to remove posns from distance calcs where there are 0 (or few) hits; won't work.
    ##for ( mm in m ) { tmp <- which( mm ); posn.hit[ tmp ] <- posn.hit[ tmp ] + 1 }
    ##for ( i in 1:length( m ) ) m[[ i ]] <- which( m[[ i ]][ ,1 ] )
    ##stop()
    ##if ( ! get.dists ) return( m )
  }

  cat( length( motifs ), "motifs\n" )
  if ( ! include.bad && exists( "coding.fracs" ) ) {
    frac.in.coding <- coding.fracs$all.fracs[ motifs ] ## in ensemble.analysis()
    motifs <- motifs[ ! is.na( frac.in.coding ) & frac.in.coding < coding.fracs$mean.fracs - 0.01 ]
    cat( length( motifs ), "motifs\n" )
    m <- m[ motifs ]
    rm( frac.in.coding )
  }
  
  ##dists <- Matrix( 0, nrow=length( motifs ), ncol=length( motifs ) )
  nc.cumsum <- c( 0, cumsum( nc ) )[ 1:length( nc ) ]; names( nc.cumsum ) <- names( nc )
  for ( i in 1:length( m ) ) {
    x <- m[[ i ]]
    if ( length( unlist( x ) ) <= 0 ) next
    for ( ii in 1:length( x ) ) x[[ ii ]] <- x[[ ii ]] + nc.cumsum[ ii ]
    m[[ i ]] <- unlist( x ); names( m[[ i ]] ) <- NULL
  }

  ##require( bigmemory ); ##require( Matrix ) ## ; require( bit )

  dirname <- sprintf( "filehash/new_motif_shadows_%s_%d", as.character( p.cutoff ), length( motifs ) )
  dir.create( dirname )
  if ( ! file.exists( sprintf( "%s.tsv.bz2", dirname ) ) ) {
    tmp <- mclapply( 1:( length( m ) - 1 ), function( i ) {
      ## out <- Matrix( 0, nrow=1, ncol=length(m) ) ## Matrix( 0, nrow=1, ncol=length(m) ) ##
      fname <- sprintf( "%s/%08d.tsv.bz2", dirname, i )
      if ( file.exists( fname ) ) return()
      print( i )
      x <- m[[ i ]]
      if ( length( x ) <= 0 ) return()
      out <- rep( 1, length( m ) )
      for ( j in (i+1):length( m ) ) {
        y <- m[[ j ]]
        if ( length( y ) <= 0 ) next
        tmp <- x %in% y
        if ( ! any( tmp ) ) next
        tmp <- ( sum( ! tmp ) + sum( ! ( y %in% x ) ) ) / length( unique( c( x, y ) ) ) ##sum( x != y ) / sum( x | y )
        out[ j ] <- tmp
      }
      tmp <- which( out < 1 )
      if ( length( tmp ) <= 0 ) return()
      write.table( cbind( i, tmp, out[tmp] ), quote=F, sep='\t', row.names=F, col.names=F, file=bzfile( fname ) )
      return()
    }, mc.preschedule=F )

    system( sprintf( "find ./%s/ -name '*.tsv.bz2' -print | sort | xargs bunzip2 -c | bzip2 -c >%s.tsv.bz2",
                    dirname, dirname ) )
  }

  cat( "Going to run mcl now...\n" )
  system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
                  dirname, distance.weight.cutoff ) )
  system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
                  dirname, distance.weight.cutoff ) )

  outfile <- sprintf( "new.mcltemp.I%s.pi%s", gsub('.','',sprintf("%.1f",mcl.I),fixed=T),
                     gsub('.','',sprintf("%.1f",mcl.pi),fixed=T) )
  mcl.cmd <- sprintf( mcl.cmd, outfile, mcl.I, mcl.pi )
  print( mcl.cmd )
  mcl.out <- system( mcl.cmd, intern=T, ignore.stderr=F )
  if ( ! improve.the.clusters ) unlink( 'new.mcltemp' )
  print( mcl.out )

  system( sprintf( 'gzip -fv %s', outfile ) )
  clusts <- lapply( strsplit( readLines( gzfile( sprintf( "%s.gz", outfile ) ) ), "\t" ), as.integer )
  cat( "GOT", length( clusts ), "motif clusters with", length( unlist( clusts ) ), "motifs.\n" )
  clusts <- clusts[ sapply( clusts, length ) >= 3 ] ## eliminate tiny clusters
  clusts <- lapply( clusts, function( i ) motifs[ i ] )
  cat( "GOT", sum( sapply( clusts, length ) >= 3 ), "motif clusters (length > 3) with",
      length( unlist( clusts[ sapply( clusts, length ) >= 3 ] ) ), "motifs.\n" )
  cat( "GOT", sum( sapply( clusts, length ) >= 10 ), "motif clusters (length > 10) with",
      length( unlist( clusts[ sapply( clusts, length ) >= 10 ] ) ), "motifs.\n" )
  mc.length <- max( which( sapply( clusts, length ) >= n.cutoff ) ) ## ignore small clusters
  if ( is.infinite( mc.length ) ) mc.length <- max( which( sapply( clusts, length ) >= 3 ) )
  attr( clusts, 'mcl.cmd' ) <- mcl.cmd
  attr( clusts, 'mcl.out' ) <- mcl.out
  attr( clusts, 'mc.length' ) <- mc.length
  attr( clusts, 'in.args' ) <- in.args
  attr( clusts, 'motifs' ) <- motifs
  save( clusts, file=sprintf( "filehash/new_motif_shadows_%s_clusts.RData", as.character( p.cutoff ) ) )
  }

  if ( improve.the.clusters ) {
    out2 <- cluster.the.new.motif.clusters( clusts )
    out$motif.cluster.clusters <- out2$motif.cluster.clusters
  }

  if ( meme.the.clusters ) { ## requires "m" to be in memory - saved as filehash/new_motif_shadows_P_CUTOFF.RData
    full.genome <- paste( e$genome.info$genome.seqs, collapse="", sep="" )
    seq.type <- 1
    bg.list <- e$genome.info$bg.list[[ seq.type ]]
    bgo <- e$bg.order[ seq.type ]
    bg.fname <- e$my.tempfile( "meme.tmp", suf=".bg" )
    tmp <- e$mkBgFile( e$genome.info$genome.seqs, order=bgo, bgfname=bg.fname, input.list=bg.list, use.rev.comp=T )

    require( Matrix )
    if ( ! exists('attrs') ) attrs <- list() ##mclapply( 1:mc.length, function( i ) {
    for ( i in (length(attrs)+1):mc.length ) {
      print(i)
      if ( length( unlist( m[ clusts[[ i ]] ] ) ) <= 0 ) next ##return( NULL ) ##next
      thresh.frac = 10
      tmp <- table( unlist( m[ clusts[[ i ]] ] ) )
      if ( length( tmp ) <= 1 ) return( NULL ) ##next
      tmp2 <- tmp[ tmp >= max( tmp ) / thresh.frac ]
      if ( length( tmp2 ) <= 1 ) return( NULL ) ##next
      m.tmp <- matrix( 0, nrow=sum( nc ), ncol=1 )
      m.tmp[ as.integer( names( tmp2 ) ), 1 ] <- as.integer( tmp2 )
      starts <- which( m.tmp > 0 & c( 0, m.tmp[ -length(m.tmp) ] ) == 0 )
      cat(i,length(starts),"\n")
      while( length( starts ) > 300 || length( starts ) < 5 ) {
        was.toobig <- length( starts ) > 300
        if ( length( starts ) > 300 ) thresh.frac <- thresh.frac / 1.1
        else if ( length( starts ) < 5 ) thresh.frac <- thresh.frac * 1.1
        tmp2 <- tmp[ tmp >= max( tmp ) / thresh.frac ]
        m.tmp <- matrix( 0, nrow=sum( nc ), ncol=1 )
        m.tmp[ as.integer( names( tmp2 ) ), 1 ] <- as.integer( tmp2 )
        starts <- which( m.tmp > 0 & c( 0, m.tmp[ -length(m.tmp) ] ) == 0 )
        if ( was.toobig && length( starts ) < 30 ) {
          thresh.frac <- thresh.frac * 1.1 ## go back up if we get too small
          tmp2 <- tmp[ tmp >= max( tmp ) / thresh.frac ]
          m.tmp <- matrix( 0, nrow=sum( nc ), ncol=1 )
          m.tmp[ as.integer( names( tmp2 ) ), 1 ] <- as.integer( tmp2 )
          starts <- which( m.tmp > 0 & c( 0, m.tmp[ -length(m.tmp) ] ) == 0 )
          break
          cat(i,thresh.frac,length(starts),"\n")
        }
        gc()
      }
      ends <- which( m.tmp > 0 & c( m.tmp[ -1 ], 0 ) == 0 )
      if ( length( ends ) != length( starts ) ) next ##return( NULL ) ##stop( "ERROR 1" )
      tmp <- which( ends - starts < 10 )
      starts[ tmp ] <- starts[ tmp ] - 5
      ends[ tmp ] <- ends[ tmp ] + 5
      sseqs <- substring( full.genome, starts, ends )
      tmp <- table( nchar( sseqs ) )
      tmp <- as.integer( names( tmp )[ c( min( which( tmp > 1 ) ), max( which( tmp > 3 ) ) ) ] )
      if ( any( is.na( tmp ) ) ) tmp <- range( nchar( sseqs ) )
      if ( tmp[2] > 40 ) tmp[2] <- 40
      starts <- starts - 10; ends <- ends + 10 ## expand out the actual seqs that are included
      sseqs <- substring( full.genome, starts, ends )
      meme.out <- 'Segmentation fault'
      cat( i, length(sseqs), "\n" )
      while( any( grepl( 'Segmentation fault', meme.out ) ) ) {
        meme.out <- e$runMeme( as.character( 1:length( sseqs ) ), sseqs, verbose=T,
                              cmd=sprintf( paste( "./progs/meme $fname -bfile %s -time 600 -dna -revcomp -maxsize",
                                "9999999 -nmotifs 2 -evt 1e9 -minw %d -maxw %d -mod zoops -nostatus -text" ),
                                bg.fname, tmp[1], tmp[2] ), filter=F )
        tmp[2] <- tmp[2] - 2
      }
      meme.out <- e$getMemeMotifInfo( meme.out )
      cat(i,length(sseqs),meme.out[[1]]$e.value,"\n")
      meme.out[[ 1 ]]$posns$abs.start <-
        starts[ as.integer( as.character( meme.out[[ 1 ]]$posns$gene ) ) ] + meme.out[[ 1 ]]$posns$start
      ##attr( clusts[[ i ]], 'meme.out' ) <- meme.out
      ##attr( clusts[[ i ]], 'sseqs' ) <- sseqs
      attrs[[i]] <- list( meme.out=meme.out, sseqs=sseqs )
    }#, mc.preschedule=F )
    
    for ( i in 1:mc.length ) if ( ! is.null( attrs[[ i ]] ) && class( attrs[[ i ]] ) != 'try-error' ) {
      attr( clusts[[ i ]], 'meme.out' ) <- attrs[[ i ]]$meme.out
      attr( clusts[[ i ]], 'sseqs' ) <- attrs[[ i ]]$sseqs
    }
    ##save( clusts, file=sprintf( "filehash/new_motif_shadows_%s_clusts.RData", as.character( p.cutoff ) ) )

    if ( get.combined.pssms ) { ## run tomtom between the motifs in each cluster, construct combined pssm for each
      for ( i in 1:length( clusts ) ) {
        tmp <- out$cluster.motifs( 'hclust', motifs=clusts[[i]], min.gene.overlap=1,
                              e.value.cutoff=Inf, p.value.cutoff=0.0001, resid.cutoff=Inf, n.cutoff=1,
                              expand=F, include.bad=T, find.bad=NA, in.tt.out=NULL, k.cut=1 )
        attr( clusts[[i]], 'tt.out' ) <- tmp$tt.out
        attr( clusts[[i]], 'tt.out2' ) <- tmp$tt.out2[[1]]
      }
    }
  }

  
  if ( plot.them && ( meme.the.clusters || get.combined.pssms ) ) {
    pdf( "Rplots.pdf" ); par( mfrow=c(4,4) )
    for ( i in 1:mc.length ) {
      print(i)
      if ( ! get.combined.pssms ) {
        mo <- attr(clusts[[i]],'meme.out')[[1]]
        e$viewPssm(mo$pssm,main=paste(i, length(clusts[[i]]), length(attr(clusts[[i]],'sseqs')), mo$e.value))
      } else {
        pssm <- attr( attr(clusts[[i]], 'tt.out2'), 'combined.pssm' )
        e$viewPssm( pssm, main=paste(i, length(clusts[[i]])) )
      }
    }
    dev.off()
  }
  ## for ( i in 1:mc.length ) {
  ##   tmp <- cluster.motifs( 'hclust', k.cut=0.99, motifs=clusts[[i]], min.gene.overlap=NA,
  ##                         e.value.cutoff=Inf, p.value.cutoff=0.001, resid.cutoff=Inf, n.cutoff=10,
  ##                         expand=F, include.bad=F, find.bad=NA, get.motif.ccs=NA )
  ## }  
  clusts
}

cluster.the.new.motif.clusters <- function( c, cutoff=0.15 ) {
  args <- attr( c, 'in.args' )
  motifs <- args$motifs
  if ( motifs == 'ALL' ) {
    motifs <- unlist( lapply( 1:nrow( motif.widths ), function( i ) {
      ii <- which( motif.widths[ i, ] > 0 )
      if ( length( ii ) <= 0 ) return( NULL )
      paste( 'MOT', i, ii, sep='_' ) } ) )
  }
  cat( length( motifs ), "motifs\n" )
  frac.in.coding <- coding.fracs$all.fracs[ motifs ] ## in ensemble.analysis()
  motifs <- motifs[ ! is.na( frac.in.coding ) & frac.in.coding < coding.fracs$mean.fracs - 0.01 ]
  cat( length( motifs ), "motifs\n" )

  dirname <- sprintf( "filehash/new_motif_shadows_%s_%d", as.character( args$p.cutoff ), length( motifs ) )
  system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
                  dirname, args$distance.weight.cutoff ) )
  system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
                  dirname, args$distance.weight.cutoff ) )

  tmp <- read.delim( 'new.mcltemp', sep=' ', head=F, row.names=NULL )
  require( Matrix )
  m <- Matrix( 0, nrow=length( motifs ), ncol=length( motifs ) )
  m[ as.matrix( tmp[ ,1:2 ] ) ] <- tmp[ ,3 ]
  m[ as.matrix( tmp[ ,2:1 ] ) ] <- tmp[ ,3 ]
  rm( tmp ); gc()

  ov.means <- Matrix( data=0, nrow=length( c ), ncol=length( c ), sparse=T )
  rnames <- 1:length( motifs ); names( rnames ) <- motifs
  ##cutoff <- 0.30
  for ( i in 1:( attr( c, 'mc.length' ) - 1 ) ) { ##length( ttmp ) - 1 ) ) {
    m.i <- rnames[ c[[ i ]] ]
    m.i <- m.i[ ! is.na( m.i ) ]
    if ( length( m.i ) <= 0 ) next
    x.tmp <- m[ m.i, ]
    for ( j in (i+1):length( c ) ) {
      m.j <- rnames[ c[[ j ]] ]
      m.j <- m.j[ ! is.na( m.j ) ]
      if ( length( m.j ) <= 0 ) next
      xx.tmp <- x.tmp[ ,m.j ] ##xx[ rnames[ m.i ], rnames[ m.j ] ]
      ov.means[ i, j ] <- mean( xx.tmp != 0 )
    }
    cat( "Similar:", i, max( ov.means[ i, ] ), which( ov.means[ i, ] > cutoff ), "\n" )
  }
  
  ov.means <- as.matrix( ov.means )
  ov.means[ lower.tri( ov.means ) ] <- t( ov.means )[ lower.tri( t( ov.means ) ) ]
  mc.length <- attr( c, 'mc.length' )
  if ( ! any( ov.means > cutoff ) ) cutoff <- cutoff - 0.05 ##0.25
  ind <- 1
  motif.cluster.clusters <- rep( 0, ncol( ov.means ) )
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
      tt <- c( i, tt ); tt <- tt[ motif.cluster.clusters[ tt ] == 0 ]
      motif.cluster.clusters[ tt ] <- ind
      ind <- ind + 1
    }
    cat( i, tt, motif.cluster.clusters[ i ], "\n" )
  }
  m=c(which(motif.cluster.clusters[1:mc.length]==0),which(motif.cluster.clusters!=0))
  cat( "Total unique motif clusters:", length(table(motif.cluster.clusters[m]))+table(motif.cluster.clusters[m])['0']-1, "\n" )
  out2 <- new.env()
  out2$ov.means <- ov.means
  out2$motif.cluster.clusters <- motif.cluster.clusters

  unlink( 'new.mcltemp' )
  out2
}

fimo.all.motifs <- function( p.cutoff=1e-6 ) { ## Note this creates files w/ no cutoff, then filters it at the end
  in.args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments
               sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## nifty trick, eh?

  dir.create( 'filehash/fimo_out' )
  seqs.file <- e$my.tempfile( 'fimo_seqs', tmpdir='./filehash/fimo_out' )
  writeLines( paste( paste( ">", names( e$genome.info$genome.seqs ), sep="" ), e$genome.info$genome.seqs, sep="\n" ),
             con=seqs.file )
  inds <- c( seq( 1, e$k.clust, by=100 ), e$k.clust )
  files <- mclapply( 1:( length( inds ) - 1 ), function( i ) {
    mots.file <- e$all.motifs.to.mast.file( ks=inds[i]:inds[i+1], seq.type=names(e$mot.weights)[1],
                                           e.value.cutoff=Inf, resid.cutoff=Inf )
    out.file <- e$my.tempfile( 'fimo_out', tmpdir='./filehash/fimo_out')

    ## Assume using customized version of meme_4.3.0 with MAX_MOTIFS in motif.h changed to 24000
    cmd <- sprintf( './progs/fimo --max-stored-scores 9999999999 --text --verbosity 2 %s %s |bzip2 -c >%s.bz2',
                   mots.file, seqs.file, out.file )
    print( cmd )
    out <- system( cmd, intern=T )
    out.file ## save it to tempfiles because storing it all up on each processor takes too much RAM
  }, mc.preschedule=F )

  ##tmp <- lapply( files, function( f ) read.delim( bzfile( sprintf( "%s.bz2", f ) ) ) )
  ##lines <- do.call( c, lines )

  fimo.out <- system( sprintf( 'bunzip2 -c filehash/fimo_out/fimo_out_*.bz2 | awk \'($6<=%s){print}\'',
                         as.character( p.cutoff ) ), intern=T )
  fimo.out <- as.data.frame( do.call( rbind, lapply( fimo.out, function( i ) strsplit( i, '\t' )[[ 1 ]] ) ) )
  colnames( fimo.out ) <- strsplit( system( sprintf( 'bunzip2 -c %s.bz2 | head -1', files[1] ), intern=T ), '\t' )[[ 1 ]]
  fimo.out$Start <- as.integer( as.character( fimo.out$Start ) )
  fimo.out$Stop <- as.integer( as.character( fimo.out$Stop ) )
  fimo.out$`Log-odds` <- as.numeric( as.character( fimo.out$`Log-odds` ) )
  fimo.out$`p-value` <- as.numeric( as.character( fimo.out$`p-value` ) )
  tmp <- do.call( rbind, strsplit( as.character( fimo.out$Motif ), '_' ) )
  fimo.out$bic <- as.integer( tmp[ ,2 ] )
  fimo.out$mot <- as.integer( tmp[ ,3 ] )
  fimo.out$Strand <- substr( tmp[ ,1 ], 1, 1 )
  rm( tmp )
  fimo.out$Motif <- NULL
  fimo.out <- as.data.table( fimo.out )
  setkey( fimo.out, bic, mot, Seq, Start )
  attr( fimo.out, 'in.args' ) <- in.args
  save( fimo.out, file=sprintf( 'filehash/fimo_out_%s.RData', as.character( p.cutoff ) ) )
  fimo.out
}

evaluate.all.clusterings <- function( ... ) {
  load( 'filehash/fimo_out_1e-06.RData' )
  load( 'filehash/new_motif_shadows_1e-06.RData' )

  retry.new.clustering <- function( mcl.I=1.2, mcl.pi=2.0, distance.weight.cutoff=0.999999, 
                                   mcl.cmd=paste( './progs/mcl-10-201/local/bin/mcl new.mcltemp -o %s --abc',
                                     ' -I %.1f -v all -te 3 -S 200000 -scheme 7 --analyze=y -pi %.1f 2>&1' ),
                                   meme.the.clusters=F,
                                   dirname='filehash/new_motif_shadows_1e-06_52977' ) {

    outfile <- sprintf( "filehash/new.mcltemp.I%s.pi%s.dcut%s", gsub('.','',sprintf("%.1f",mcl.I),fixed=T),
                       gsub('.','',sprintf("%.1f",mcl.pi),fixed=T),
                       sprintf("%.3f",distance.weight.cutoff) )

    system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $1,$2,1-$3}' > new.mcltemp", 
                    dirname, distance.weight.cutoff ) )
    system( sprintf( "bunzip2 -c %s.tsv.bz2 | awk '($3<=%.2f){print $2,$1,1-$3}' >> new.mcltemp", 
                    dirname, distance.weight.cutoff ) )
    mcl.cmd <- sprintf( "%s  2>&1", mcl.cmd )
    mcl.cmd <- sprintf( mcl.cmd, outfile, mcl.I, mcl.pi )
    print( mcl.cmd )
    mcl.out <- system( mcl.cmd, intern=T, ignore.stderr=F )
    
    system( sprintf( 'gzip -fv %s', outfile ) )
    clusts <- lapply( strsplit( readLines( gzfile( sprintf( "%s.gz", outfile ) ) ), "\t" ), as.integer )
    cat( "GOT", length( clusts ), "motif clusters with", length( unlist( clusts ) ), "motifs.\n" )
    clusts <- clusts[ sapply( clusts, length ) >= 3 ] ## eliminate tiny clusters
    ##clusts <- lapply( clusts, function( i ) motifs[ i ] )
    cat( "GOT", sum( sapply( clusts, length ) >= 3 ), "motif clusters (length > 3) with",
        length( unlist( clusts[ sapply( clusts, length ) >= 3 ] ) ), "motifs.\n" )
    cat( "GOT", sum( sapply( clusts, length ) >= 10 ), "motif clusters (length > 10) with",
        length( unlist( clusts[ sapply( clusts, length ) >= 10 ] ) ), "motifs.\n" )
    attr( clusts, 'mcl.cmd' ) <- mcl.cmd
    attr( clusts, 'mcl.out' ) <- mcl.out
    clusts
  }
  out$retry.new.clustering <- retry.new.clustering
  environment( out$retry.new.clustering ) <- out
  
  tmpl <- list()
  try( load( 'tmpl.RData' ) )
  sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out, chdir=T )
  for ( mcl.I in rev( c( 1.2, 1.5, 2.0, 2.5, 3.0, 4.5, 6.0 ) ) ) {
    for ( mcl.pi in c( 1, 2, 5, 10, 20 ) ) {
      for ( dc in c( 0.8, 0.9, 0.95, 0.98, 0.99, 0.999, 0.999999 ) ) {
        if ( paste( mcl.I, mcl.pi, dc ) %in% names( tmpl ) ) next
        print( paste( mcl.I, mcl.pi, dc ) )
        tmp2 <- capture.output( tmp <- out$retry.new.clustering( mcl.I=mcl.I, mcl.pi=mcl.pi, distance.weight.cutoff=dc, ... ) )
        attr( tmp, 'captured.output' ) <- tmp2
        tmpl[[ paste( mcl.I, mcl.pi, dc ) ]] <- tmp
        print( paste( mcl.I, mcl.pi, dc ) ); print( attr( tmp, 'mcl.cmd' ) )
        print( attr( tmp, 'mcl.out' ) ); print( tmp2 ); print( date() )
        save( tmpl, file='tmpl.RData' )
      }
    }
  }
  save(tmpl,file='tmpl.RData')

  ## Calculate in vs. out densities of each clustering
  motifs <- unlist( lapply( 1:nrow( out$motif.widths ), function( i ) {
    ii <- which( out$motif.widths[ i, ] > 0 )
    if ( length( ii ) <= 0 ) return( NULL )
    paste( 'MOT', i, ii, sep='_' ) } ) )
  frac.in.coding <- out$coding.fracs$all.fracs[ motifs ] ## in ensemble.analysis()
  motifs <- motifs[ ! is.na( frac.in.coding ) & frac.in.coding < out$coding.fracs$mean.fracs - 0.01 ]
  
  require(data.table)
  df=data.table(read.delim(bzfile(sprintf('filehash/new_motif_shadows_1e-06_%d.tsv.bz2',length(motifs))),head=F))
  setkey(df,V1,V2)
  bad.inds<-integer()
  for(i in 1:length(motifs)){ ## some motifs w/ no hits got files anyway, but that have all 0's -- we need to remove these.
    if(!file.exists(sprintf('filehash/new_motif_shadows_1e-06_%d/%08d.tsv.bz2',length(motifs),i))) {bad.inds=c(bad.inds,i);next}
    tmp=read.delim(bzfile(sprintf('filehash/new_motif_shadows_1e-06_%d/%08d.tsv.bz2',length(motifs),i)),head=F)
    if(all(tmp$V3==0)){bad.inds=c(bad.inds,i);print(i)}
  }
  cat(length(bad.inds),'motifs with no hits to genome...?\n')
  df=df[!(V1%in%bad.inds)&!(V2%in%bad.inds)]
  rm(frac.in.coding);gc()

  require(Matrix)
  m=Matrix(0,nrow=length(motifs),ncol=length(motifs)); m=forceSymmetric(m)
  m[cbind(df$V1,df$V2)] <- 1-df$V3 ## Note m is a SIMILARITY matrix (i.e., 1-distance) so 1 is similar!
  m[cbind(df$V2,df$V1)] <- 1-df$V3 ## this is to keep it sparse
  save(m,df,file='m_and_df.RData')

  max.len <- 400 ## Make sure a clustering doesn't cheat w/ too many small clusters. TODO: how to set this organism-specifically.
  m.inds=1:length(motifs)
  dens.out=list()
  for( i in names(tmpl) ) {
    print(i)
    if ( ! is.null(dens.out[[i]]) ) next
    ii<-i; i=tmpl[[i]]
    mc.len=sum(sapply(i,length)>=10)
    if(mc.len>max.len)mc.len=max.len
    cluster.densities=mclapply(i[1:mc.len],function(j) {
      mm=m[j,j] ##as.matrix(m[j,j])
      d=c(sum(mm),nrow(mm))##/nrow(mm)/(nrow(mm)-1)
    } )
    between.densities=mclapply(i[1:mc.len],function(j){
      not.j=m.inds[m.inds%in%unlist(i[1:mc.len]) & !m.inds%in%j]
      mm=m[j,not.j]
      d=c(sum(mm),nrow(mm),ncol(mm))##/nrow(mm)/ncol(mm)
    })
    out.density={
      non.hits=m.inds[!m.inds%in%unlist(i[1:mc.len])]
      mm=m[non.hits,non.hits]
      d=c(sum(mm),nrow(mm))##/nrow(mm)/(nrow(mm)-1)
    }
    dens.out[[ii]]=list(cluster.dens=cluster.densities,between.dens=between.densities,out.dens=out.density)
    gc()
  }  

  save(dens.out,file='dens.out.RData')

  tmp2 = t( sapply( names(dens.out), function(i) { ##for ( i in names(dens.out) ) {
    cl.size = sapply(dens.out[[i]]$cluster.dens,function(j)j[2])
    max.ind = if ( length(cl.size) > 400 ) 400 else length(cl.size)
    in.dens = sapply(dens.out[[i]]$cluster.dens[1:max.ind],function(j)j[1]/j[2]/j[2])
    betw.dens = sapply(dens.out[[i]]$between.dens[1:max.ind],function(j)j[1]/j[2]/j[3])
    in.dens[in.dens==0&betw.dens==0] = betw.dens[in.dens==0&betw.dens==0] = NA ## bad cluster!
    out.dens = dens.out[[i]]$out.dens[1]/dens.out[[i]]$out.dens[1]/dens.out[[i]]$out.dens[1]
    ##in.dens.vs.out.dens = in.dens/sum(sapply(dens.out[[i]]$cluster.dens[1:max.ind],function(j)j[2]))^2 /
    c(sum(cl.size[1:max.ind]),max.ind,weighted.mean(in.dens,cl.size[1:max.ind],na.rm=T),
      weighted.mean(betw.dens,cl.size[1:max.ind],na.rm=T),out.dens)
  } ) )

  tmp <- t(sapply(strsplit(names(tmpl)," "),as.numeric))
  dc <- table( tmp[,3] )
  dc[] <- 1:length(dc)
  pi <- table( tmp[,2] )
  pi[] <- 0:(length(pi)-1)
  mcl.I <- table( tmp[,1] )
  mcl.I[] <- 1:length(mcl.I)
  mc.len <- sapply(tmpl, function(i)sum(sapply(i,length)>=10) )
  
  par(mfrow=c(3,3))
  plot(tmp2[,1],log10(tmp2[,3]),xlab='cl.size',ylab='log10(in.dens)',cex=0.7,col=dc[as.character(tmp[,3])],
       pch=pi[as.character(tmp[,2])])
  legend( 'bottomleft', legend=names(dc), pch=20, cex=0.7, col=1:length(dc), title='distance.cut' )  
  plot(tmp2[,1],log10(tmp2[,4]),xlab='cl.size',ylab='log10(betw.dens)',cex=0.7,col=dc[as.character(tmp[,3])],
       pch=pi[as.character(tmp[,2])])       
  legend( 'bottomleft', legend=names(pi), cex=0.7, pch=0:(length(pi)-1), title='mcl.PI' )  
  plot(tmp2[,1],log10(tmp2[,5]),xlab='cl.size',ylab='log10(out.dens)',cex=0.7,col=dc[as.character(tmp[,3])],
       pch=pi[as.character(tmp[,2])])
  plot(tmp2[,2],log10(tmp2[,3]/tmp2[,4]),xlab='mc.len',ylab='log10(in.dens/betw.dens)',cex=0.7,col=dc[as.character(tmp[,3])],
              pch=pi[as.character(tmp[,2])])
  plot(tmp2[,1],log10(tmp2[,3]/tmp2[,4]),xlab='cl.size',ylab='log10(in.dens/betw.dens)',cex=0.7,col=dc[as.character(tmp[,3])],
       pch=pi[as.character(tmp[,2])])       
  plot(tmp2[,1],log10(tmp2[,3]/tmp2[,4]),xlab='cl.size',ylab='log10(in.dens/betw.dens)',cex=0.7,col=mcl.I[as.character(tmp[,1])],
       pch=pi[as.character(tmp[,2])])       
  legend( 'bottomleft', legend=names(mcl.I), pch=20, cex=0.7, col=1:length(mcl.I), title='mcl.I' )  
  plot(mc.len,log10(tmp2[,3]/tmp2[,4]),xlab='mc.len',ylab='log10(in.dens/betw.dens)',cex=0.7,col=mcl.I[as.character(tmp[,1])],
       pch=pi[as.character(tmp[,2])])
  graphics.to.pdf( "evaluate_all_clusterings.pdf" )

  clusts <- tmpl[[which.max(log10(tmp2[,3]/tmp2[,4]))]] ## find the best
  for(i in 1:length(clusts))clusts[[i]]<-motifs[clusts[[i]]]

  mc.length <- sum( sapply(clusts, length) >=10 )
  for ( i in 1:mc.length ) {
    tmp <- out$cluster.motifs( 'hclust', motifs=clusts[[i]], min.gene.overlap=1,
                              e.value.cutoff=Inf, p.value.cutoff=0.0001, resid.cutoff=Inf, n.cutoff=1,
                              expand=F, include.bad=T, find.bad=NA, in.tt.out=NULL, k.cut=1 )
    attr( clusts[[i]], 'tt.out' ) <- tmp$tt.out
    attr( clusts[[i]], 'tt.out2' ) <- tmp$tt.out2[[1]]
  }

  pdf( "evaluate_all_clusterings_clusters.pdf" )
  par(mfrow=c(4,4))
  for ( i in 1:mc.length ) {
    pssm <- attr( attr(clusts[[i]], 'tt.out2'), 'combined.pssm' )
    e$viewPssm( pssm, main=paste(i, length(clusts[[i]])) )
  }
  dev.off()

  attr( clusts, 'mc.length' ) <- mc.length
  attr( clusts, 'motifs' ) <- motifs
  
  ##tmpl
}


ensemble.to.database <- function(out) {
  dir.create( 'filehash/DATABASES/' )
  ## DATA - ratios, string, operons, annotations
  tmp <- data.table( round(out$e$ratios$ratios*1000)/1000 )
  tmp <- data.table( gene=rownames(out$e$ratios$ratios), tmp )
  write.table( tmp, file='filehash/DATABASES/ratios.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/ratios.tsv', wait=F )
  tmp <- data.table( out$e$genome.info$feature.tab )
  write.table( tmp, file='filehash/DATABASES/feature_tab.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/feature_tab.tsv', wait=F )
  tmp <- data.table( out$e$genome.info$feature.names )
  write.table( tmp, file='filehash/DATABASES/feature_names.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/feature_names.tsv', wait=F )
  tmp <- data.table( out$e$genome.info$operons )
  write.table( tmp, file='filehash/DATABASES/operons.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/operons.tsv', wait=F )
  tmp <- data.table( gene=names(out$e$genome.info$all.upstream.seqs[[1]]), sequence=out$e$genome.info$all.upstream.seqs[[1]] )
  write.table( tmp, file='filehash/DATABASES/upstream_seqs.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/upstream_seqs.tsv', wait=F )
  if ( ! is.null( out$e$networks$operons ) ) {
    tmp <- data.table( out$e$networks$operons )
    write.table( tmp, file='filehash/DATABASES/operon_network.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
    rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/operon_network.tsv', wait=F )
  }
  if ( ! is.null( out$e$networks$string ) ) {
    tmp <- data.table( out$e$networks$string )
    write.table( tmp, file='filehash/DATABASES/string_network.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
    rm( tmp ); system( 'bzip2 -fv filehash/DATABASES/string_network.tsv', wait=F )
  }
  
  ## BICLUSTERS - info, genes, conditions, motifs
  bzcon1 <- 'filehash/DATABASES/bicluster.tsv'
  bzcon2 <- 'filehash/DATABASES/bicluster_genes.tsv'
  bzcon3 <- 'filehash/DATABASES/bicluster_conds.tsv'
  bzcon4 <- 'filehash/DATABASES/bicluster_motifs.tsv'
  wrote <- FALSE
  for ( k in 1:out$e$k.clust ) {
    if ( k %% 100 == 0 ) print(k)
    tab1 <- tab2 <- tab3 <- tab4 <- NULL
    bb <- paste('BIC',k,sep='_')
    clust <- out$get.bicluster.info(bb)[[1]]
    if ( ! is.null(clust) ) {
      fname <- names(out$e$fnames.to.cluster[which(out$e$fnames.to.cluster==k)])
      k.orig <- clust$k
      tab1 <- data.table( bic=k, nrow=clust$nrow, ncol=clust$ncol, resid=clust$resid, pclust=clust$p.clust,
                         eval=min(clust$e.val,na.rm=T), k_orig=k.orig, fname=fname )
      if ( ! is.null( tab1 ) ) write.table( tab1, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    }
    genes <- clust$rows
    if ( ! is.null( genes ) ) tab2 <- data.table( bic=k, gene=genes )
    if ( ! is.null( tab2 ) ) write.table( tab2, bzcon2, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    conds <- clust$cols
    if ( ! is.null( conds ) ) tab3 <- data.table( bic=k, cond=conds )
    if ( ! is.null( tab3 ) ) write.table( tab3, bzcon3, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    mots <- out$get.motifs(bicluster=bb)[[1]]
    if ( ! is.null( mots ) ) tab4 <- data.table( bic=k, mot=gsub('MOT_','',mots) )
    if ( ! is.null( tab4 ) ) write.table( tab4, bzcon4, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    wrote <- TRUE
  }
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon2 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon3 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon4 ), wait=F )

  ## MOTIFS - info; meme positions; mast positions; pssms
  bzcon1 <- 'filehash/DATABASES/motif.tsv'
  bzcon2 <- 'filehash/DATABASES/motif_meme_posn.tsv'
  bzcon3 <- 'filehash/DATABASES/motif_mast_posn.tsv'
  bzcon4 <- 'filehash/DATABASES/motif_pssm.tsv'
  wrote <- FALSE
  for ( k in 1:out$e$k.clust ) {
    bic <- paste('BIC',k,sep='_')
    mots <- out$get.motifs(bicluster=bic)[[1]]
    for ( m in mots ) {
      tab1 <- NULL
      if ( k %% 100 == 0 ) print(m)
      minf <- out$get.motif.info(m)[[1]]
      mm <- gsub( 'MOT_', '', m )
      if ( ! is.null( minf ) ) {
        cf <- out$coding.fracs$all.fracs[m]
        tab <- data.table( mot=mm, width=minf$width, llr=minf$llr, eval=minf$e.value, sites=minf$sites,
                          coding=cf, good=(cf < out$coding.fracs$mean.fracs - 0.01) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        tab <- NULL; if ( nrow(minf$posns) > 0 ) tab <- as.data.table( cbind( mot=mm, minf$posns ) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon2, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        tab <- NULL; if ( nrow(minf$mast) > 0 ) tab <- as.data.table( cbind( mot=mm, minf$mast ) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon3, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        pssm <- minf$pssm; colnames( pssm ) <- out$e$col.let; pssm <- as.data.table( pssm )
        tab <- data.table( cbind( mot=mm, ind=1:nrow(pssm), round(pssm*1000)/1000 ) )
        if ( ! is.null( tab ) ) write.table( tab, bzcon4, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
        wrote <- TRUE
      }
    }
  }
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon2 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon3 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon4 ), wait=F )

  ## MOTIF CLUSTS -- motif members and combined pssms
  bzcon1 <- 'filehash/DATABASES/motif_clust.tsv'
  bzcon2 <- 'filehash/DATABASES/motif_clust_combined_pssm.tsv'
  bzcon3 <- 'filehash/DATABASES/motif_clust_aligned_pssms.tsv'
  wrote <- FALSE
  for ( k in 1:out$mc.length ) {
    if ( k %% 100 == 0 ) print(k)
    mc <- paste('MOTC', k, sep='_')
    mots <- out$get.motifs(motif.clust=mc)[[1]]
    tab <- data.table( motc=k, mot=gsub('MOT_','',mots) )
    if ( ! is.null( tab ) ) write.table( tab, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    mc.info <- out$get.motif.cluster.info(mc)[[1]]
    pssm <- attr(mc.info, 'combined.pssm' )
    colnames( pssm ) <- out$e$col.let; pssm <- as.data.table( pssm )
    tab <- data.table( cbind( motc=k, ind=1:nrow(pssm), round(pssm*1000)/1000 ) )
    if ( ! is.null( tab ) ) write.table( tab, bzcon2, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
    aligned.pssms <- attr(mc.info, 'aligned.pssms')
    wrote2 <- wrote
    for ( m in names( aligned.pssms ) ) {
      pssm <- aligned.pssms[[m]]
      colnames( pssm ) <- out$e$col.let; pssm <- as.data.table( pssm )
      tab <- data.table( cbind( motc=k, motif=gsub('MOT_','',m), ind=1:nrow(pssm), round(pssm*1000)/1000 ) )
      if ( ! is.null( tab ) ) write.table( tab, bzcon3, quote=F, sep='\t', row.names=F, col.names=!wrote2, append=wrote2 )
      wrote2 <- TRUE
    }
    wrote <- TRUE
  }
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon2 ), wait=F )
  system( sprintf( 'bzip2 -fv -9 %s', bzcon3 ), wait=F )

  ## MOTIF SCANS
  tmp <- out$pssm.scans; setnames( tmp, 'gene', 'p.value', 'posn', 'mot', 'bic' )
  tmp$strand <- ifelse( tmp$mot > 0, '+', '-' ); tmp$mot <- abs( tmp$mot )
  write.table( tmp, file='filehash/DATABASES/motif_scans.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  system( sprintf( 'bzip2 -fv -9 %s', 'filehash/DATABASES/motif_scans.tsv' ), wait=F )
  rm( tmp ); gc()

  load( "filehash/fimo_out_1e-06.RData" )
  write.table( fimo.out, file='filehash/DATABASES/motif_scans_fimo.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  system( sprintf( 'bzip2 -fv -9 %s', 'filehash/DATABASES/motif_scans_fimo.tsv' ), wait=F )
  rm( fimo.out ); gc()

  ## MEME HITS
  write.table( out$meme.hits, file='filehash/DATABASES/motif_meme_hits.tsv', quote=F, sep='\t', row.names=F, col.names=T, append=F )
  system( sprintf( 'bzip2 -fv -9 %s', 'filehash/DATABASES/motif_meme_hits.tsv' ), wait=F )

  ## MOTIF GENOME-WIDE LOGO INFO
  genomes <- lapply( names( out$e$genome.info$genome.seqs ), function(i) {
    nc <- nchar( out$e$genome.info$genome.seqs[i] )
    output <- c( 1, nc ); names( output ) <- c( i, '' ); output } )
  options( cores=4 )
  bzcon1 <- 'filehash/DATABASES/motif_promoter_counts.tsv'
  wrote <- FALSE
  ## THIS GETS THE COUNTS/PSSMs FOR MOTIFS IN MOTIF CLUSTERS ONLY:
  for ( mc in 1:out$mc.length ) {
    for ( i in genomes ) {
      mots <- out$motif.clusts[[mc]]
      tmp <- out$plot.promoter.architecture( i, type='fimo', dont.plot=T, motif.filter=mots, include.bad=T )
      ## tmp$mat is regular pssm; tmp$mat2 is pssm scaled by info-content.
      not.zeroes <- apply( cbind( tmp$counts[,1], tmp$mat ), 1, function(i) any(i!=0) )
      tab <- data.table( motc=mc, loc=as.integer(names(tmp$counts[not.zeroes,])), count=tmp$counts[not.zeroes,1],
                        round(tmp$mat[not.zeroes,]*1000)/1000 )
      cat( mc, nrow(tab), "\n" )
      if ( ! is.null( tab ) ) write.table( tab, bzcon1, quote=F, sep='\t', row.names=F, col.names=!wrote, append=wrote )
      wrote <- TRUE
    }
  }

  ## NOW DO IT FOR ALL MOTIFS (WHETHER IN MOTIF CLUSTERS OR NOT)
  ## NOTE we use '0' as the value for motc here.
  tmp <- out$plot.promoter.architecture( i, type='fimo', dont.plot=T, include.bad=F, count.all=T )
  not.zeroes <- apply( cbind( tmp$counts[,1], tmp$mat ), 1, function(i) any(i!=0) )
  tab <- data.table( motc=0, loc=as.integer(names(tmp$counts[not.zeroes,1])), count=tmp$counts[not.zeroes,1],
                    round(tmp$mat[not.zeroes,]*1000)/1000 )
  write.table( tab, bzcon1, quote=F, sep='\t', row.names=F, col.names=F, append=T )
  
  system( sprintf( 'bzip2 -fv -9 %s', bzcon1 ), wait=F )

  ## INFERELATOR COEFFS (TBD)


  ## TRICK TO READ THESE IN FAST (using data.table package):
  ## fread(paste(readLines(bzfile('bicluster.tsv.bz2')),collapse='\n')) ## never mind - this crashes. Need to bunzip2 it first, then use fread()
  ## or use sqldf
  ## read.csv.sql(bzfile('bicluster.tsv.bz2'),sep='\t')
  ## so this:
  ## tf <- e$my.tempfile(); system(sprintf('bunzip2 -c filehash/DATABASES/motif_promoter_counts.tsv.bz2 >%s', tf))
  ## tab <- fread(tf); unlink(tf)
  ## or, e.g.:
  ## tf <- file(tf); tmp <- sqldf("select ind,count from tf where motc=0",file.format=list(sep='\t',head=T),verbose=T)
}

