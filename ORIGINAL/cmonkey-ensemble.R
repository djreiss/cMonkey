source( "cmonkey-motif-other.R" )

## Create an ensemble from a single run by setting, e.g.
##   cm.func.each.iter <- function() if ( iter %% 25 == 0 ) save.image( sprintf( "zzz_hpy_%04d.RData", iter ) )
## Then run this with glob="zzz_hpy_*.RData"
## This can also be used for RData output from multiple runs.
## Then, to get bicluster #24 in RData file from iter=1900, for example:
## which(names(e$fnames.to.cluster)=="./zzz_hpyu_1900.RData"&e$fnames.to.cluster==24)
## TESTING: can add new RData files to an env that was already built using this function.
cmonkey.ensemble <- function( glob, env=NULL, filter=NULL, filehash=T, plot=F ) {
  if ( length( glob ) == 1 ) files <- list.files( patt=glob2rx( glob ), full=T )
  else files <- glob
  if ( length( files ) <= 0 ) return( NULL )
  all.ratios <- if ( is.null( env ) ) NULL else env$ratios$ratios
  fnames.to.cluster <- if ( is.null( env ) ) numeric() else env$fnames.to.cluster
  files <- files[ ! files %in% names( fnames.to.cluster ) ] ## if adding to previously built env, only add new files
  if ( filehash && file.exists( "filehash" ) ) cat( "WARNING: filehash/ exists -- press CTRL-C and delete it via unlink('filehash',recurs=T)\n" )
  for ( f in files ) {
    print( f )
    load( f )
    if ( ! exists( "e" ) ) next
    if ( ! is.null( filter ) ) {
      tmp <- filter( e )
      cat( "FILTER:", tmp, "\n" )
      if ( tmp == FALSE ) next
    }
    if ( plot ) {
      print( tail( e$stats ) )
      e$plotStats()
    }
    e$parallel.cores <- 1
    e$clusterStack <- e$get.clusterStack( force=T )

    if ( TRUE ) { ## Fix for old versions of cmonkey
      sys.source("~/scratch//biclust/cmonkey-motif-other.R",envir=e,chdir=T)
      print( names( e$meme.scores ) )
      if ( length( e$meme.scores ) == 1 && names( e$meme.scores ) == 'upstream' ) { 
        names( e$meme.scores ) <- names( e$mot.weights ) <- names( e$n.motifs ) <- names( e$motif.upstream.scan ) <-
          names( e$motif.upstream.search ) <- names( e$meme.cmd ) <- names( e$mast.cmd ) <- names( e$bg.order ) <-
            names( e$operon.shift ) <- names( e$genome.info$bg.list ) <- 'upstream meme'
        sys.source("~/scratch//biclust/cmonkey.R",envir=e,chdir=T)
        e$genome.info$bg.fname<-e$my.tempfile( "meme.tmp", suf=".bg" )
        names(e$genome.info$bg.fname)<-'upstream meme'
      } else if ( length( e$meme.scores ) == 2 ) { ## wtf! some weird halo runs there
        e$meme.scores[[ 1 ]] <- e$meme.scores[[ 'upstream meme' ]]
        e$meme.scores[[ 2 ]] <- NULL
      }
    }

    good.ks <- 1:e$k.clust

    print( names( e$meme.scores ) )
    if ( is.null( env ) ) {
      env <- e
      all.ratios <- e$get.cluster.matrix()
      tmp <- e$stats ##get.stats() ##stats[ nrow( e$stats ), ]
      env$stats <- list()
      env$stats[[ f ]] <- tmp
      tmp.f <- good.ks; names( tmp.f ) <- rep( f, length( tmp.f ) )
      fnames.to.cluster <- tmp.f
      env$from.files <- f
      if ( filehash ) {
        require( filehash ) ## Must use DJR customized version that allows integer indices!
        dir.create( 'filehash' )
        if ( dbCreate( "filehash/clusterStack.dump", type='RDS' ) )
          tmpa <- dbInit( "filehash/clusterStack.dump", type='RDS' ) ##dumpList( env$clusterStack, "filehash/clusterStack.dump", type='RDS' )
        for ( k in 1:env$k.clust ) dbInsert( tmpa, sprintf( "%08d", k ), env$clusterStack[[ k ]] )
        env$clusterStack <- tmpa; rm( tmpa )
        for ( i in 1:length( env$meme.scores ) ) { ##names( env$meme.scores ) ) {
          ind <- i ##which( names( env$meme.scores ) == i )
          tmpa <- env$meme.scores[[ i ]]
          tmpa$all.pv <- NULL
          if ( dbCreate( sprintf( "filehash/meme.scores.%d.dump", ind ), type='RDS' ) )
            tmpa <- dbInit( sprintf( "filehash/meme.scores.%d.dump", ind ), type='RDS' ) ##dumpList( tmpa, sprintf( "filehash/meme.scores.%d.dump", ind, type='RDS' )
          for ( k in 1:env$k.clust ) dbInsert( tmpa, sprintf( "%08d", k ), env$meme.scores[[ i ]][[ k ]] )
          env$meme.scores[[ i ]] <- tmpa; rm( tmpa ); gc()
        }
      }
      ##env<<-env
      cat(length(e$clusterStack),length(e$meme.scores[[1]]),"\n")
      cat(length(env$clusterStack),length(env$meme.scores[[1]]),"\n")
      cat(dim(all.ratios),"\n")
      rm( e )
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

    tmp.f <- good.ks + max( c( 0, fnames.to.cluster ) ) ##);
    names( tmp.f ) <- rep( f, length( tmp.f ) )
    fnames.to.cluster <- c( fnames.to.cluster, tmp.f )
    
    ## No row.membership anymore (version >=4.9.0)!
    if ( as.integer( gsub( '.', '', env$cmonkey.version, fixed=T ) ) < 490 ) {
      env$get.rows <- function( k ) clusterStack[[ k ]]$rows; environment( env$get.rows ) <- env
      if ( FALSE && exists( 'row.membership', envir=e ) ) {
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

      env$get.cols <- function( k ) clusterStack[[ k ]]$cols; environment( env$get.cols ) <- env
      if ( FALSE && exists( 'col.membership', envir=e ) ) {
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
    }
    
    if ( length( e$clusterStack ) < e$k.clust ) e$clusterStack[[ e$k.clust ]] <- ""
    if ( ! filehash ) {
      env$clusterStack <- c( env$clusterStack, e$clusterStack[ good.ks ] )
    } else {
      ll <- env$k.clust - e$k.clust ##length( env$clusterStack )
      for ( k in good.ks ) { ##env$clusterStack[[ length( env$clusterStack ) + 1 ]] <- e$clusterStack[[ k ]]
        dbInsert( env$clusterStack, sprintf( "%08d", ll + k ), e$clusterStack[[ k ]] )
      }
    }

    ## tmp <- e$get.stats() ##stats[ nrow( e$stats ), ]
    ## if ( nrow( e$stats ) > 0 && ! all( colnames( e$stats ) %in% colnames( tmp ) ) ) {
    ##   tmp2 <- e$stats[ 1, ] * NA
    ##   tmp2[ , colnames( tmp ) ] <- tmp; tmp <- tmp2
    ## } else if ( nrow( e$stats ) <= 0 || ! all( colnames( tmp ) %in% colnames( e$stats ) ) ) {
    ##   tmp <- tmp[ ,colnames( e$stats ) ]
    ## }
    ## if ( ! all( colnames( env$stats ) %in% colnames( tmp ) ) ) {
    ##   tmp2 <- rbind( env$stats, env$stats[ 1, ] )
    ##   tmp2[ nrow( tmp2 ), ] <- NA
    ##   tmp2[ nrow( tmp2 ), colnames( tmp ) ] <- tmp
    ##   tmp <- tmp2[ nrow( tmp2 ), ]
    ## } else if ( ! all( colnames( tmp ) %in% colnames( env$stats ) ) ) {
    ##   tmp <- tmp[ ,colnames( env$stats ) ]
    ## }
    ## env$stats <- rbind( env$stats, tmp )
    ## rm( tmp, tmp2 )

    if ( ! is.null( e$meme.scores ) && ! is.null( e$meme.scores[[ 1 ]] ) ) {
      ##e$meme.scores[[ 1 ]] <- 
      for ( i in 1:length( e$meme.scores ) ) { ##names( e$meme.scores ) ) {
        ##if ( ! i %in% names( env$meme.scores ) ) env$meme.scores[[ i ]][ tmp.f ] <- as.list( rep( "", env$k.clust ) )
        if ( length( e$meme.scores[[ i ]] ) < e$k.clust ) e$meme.scores[[ i ]][[ e$k.clust ]] <- ""
        if ( i > length( env$meme.scores ) ) {
          if ( dbCreate( sprintf( "filehash/meme.scores.%d.dump", i ), type='RDS' ) ) {
            tmpa <- dbInit( sprintf( "filehash/meme.scores.%d.dump", i ), type='RDS' )
            env$meme.scores[[ i ]] <- tmpa; rm( tmpa ); gc()
          }
        }
        if ( ! filehash ) {
          env$meme.scores[[ i ]][ tmp.f ] <- e$meme.scores[[ i ]][ good.ks ]
          names( env$meme.scores[[ i ]] ) <- NULL
        } else {
          ll <- env$k.clust - e$k.clust ##length( env$meme.scores[[ i ]] )
          tmp <- try( env$meme.scores[[i]][[ll]] )
          while( class( tmp ) == 'try-error' || ( class( tmp ) == 'character' && tmp == '' ) ) {
            ll <- ll - 1
            tmp <- try( env$meme.scores[[i]][[ll]] )
          }
          for ( k in good.ks ) {
            dbInsert( env$meme.scores[[ i ]], sprintf( "%08d", ll + k ), e$meme.scores[[ i ]][[ k ]] )
          }
        }
      }
    }
    env$stats[[ f ]] <- e$stats ## get.stats() ## stats[ nrow( e$stats ) ]
    cat(length(e$clusterStack),length(e$meme.scores[[1]]),"\n")
    cat(length(env$clusterStack),length(env$meme.scores[[1]]),"\n")
    cat(dim(all.ratios),"\n")
    cat(env$clusterStack[[length(env$clusterStack)]]$k,
        env$meme.scores[[1]][[length(env$meme.scores[[1]])-1]]$k,"\n")
    rm( e )
    rm( n.scores, net.scores, row.scores, r.scores, rr.scores, mot.scores, m.scores, row.memb, col.memb, envir=env )
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
  for ( n in names( env$meme.scores ) ) env$meme.scores[[ n ]]$all.pv <- ''
  for ( n in names( env$meme.scores ) ) env$meme.scores[[ n ]]$all.pv <- NULL
  
  ##rm( n.scores, net.scores, row.scores, r.scores, rr.scores, mot.scores, m.scores, row.memb, col.memb, envir=env )
  cat( "Got", env$k.clust, "clusters for the ensemble.\n" )
  rm( row.membership, col.membership, row.scores, net.scores, mot.scores,
     m.scores, n.scores, r.scores, rr.scores, row.memb, col.scores,
     c.scores, cc.scores, col.memb, cmonkey.init, write.project, old.row.membership,
     pre.adjusted.row.membership, old.col.membership, envir=env )
  ##rm( net.scores, row.scores, mot.scores, envir=env )
  env
}

## TESTING: using it on an "ensemble" from a single run (as described above) -- should use 'remove.exact.dupes=T' to remove
##   exact duplicate biclusters, but this currently does not work)
## Set e$dupe.ks to the output of this func; then those clusters will be ignored in subsequent analyses (TODO)
##    THIS IS NOT TO BE USED:::
get.all.dupes <- function( e, compare.cols=F, max.no.dupe.count=1000 ) {
  dupe.inds <- integer()
  for ( k1 in e$k.clust:2 ) {
    if ( k1 %in% dupe.inds ) next
    print(k1)
    cc <- e$clusterStack[[ k1 ]]
    r1 <- cc$rows
    c1 <- cc$cols
    res1 <- cc$resid
    if ( length( r1 ) <= 1 || length( c1 ) <= 1 || is.na( res1 ) ) next
    no.dupe.count <- 0
    for ( k2 in (k1-1):1 ) {
      if ( k2 %in% dupe.inds ) next
      if ( no.dupe.count > max.no.dupe.count ) { cat( k2, no.dupe.count, "\n" ); break } ## Don't need to count all the way down to 1 if we have already counted down this many and NOT seen a dupe
      cc <- e$clusterStack[[ k2 ]]
      res2 <- cc$resid
      if ( is.na( res2 ) ) next
      r2 <- cc$rows
      if ( length( r1 ) != length( r2 ) ) { no.dupe.count <- no.dupe.count + 1; next }
      if ( ! all( r1 %in% r2 ) ) { no.dupe.count <- no.dupe.count + 1; next }
      if ( compare.cols ) {
        c2 <- cc$cols
        if ( length( c1 ) != length( c2 ) ) { no.dupe.count <- no.dupe.count + 1; next }
        if ( ! all( c1 %in% c2 ) ) { no.dupe.count <- no.dupe.count + 1; next }
      }
      dupe.inds <- unique( c( dupe.inds, k2 ) )
      cat( "DUPE:", k2, k1, length( dupe.inds ), no.dupe.count, "\n" )
      no.dupe.count <- 0
    }
  }
  dupe.inds
}  

## cluster.motifs can be 'mcl' to use mcl instead of hclust
cmonkey.ensemble.analysis <- function( e, make.sif=T, cluster.motifs=F, min.gene.overlap=1,
                                      e.value.cutoff=Inf, p.value.cutoff=0.001, resid.cutoff=Inf,
                                      find.bad=T, get.motif.ccs=T, verbose=T, filehash=F, n.cutoff=10,
                                      ... ) { ## e from cmonkey.ensemble()
  ##if ( force.single ) mclapply <- lapply ## even 2 cores is too much for my mac
  ##else
  require( parallel )
  out <- new.env() ##list( sif=r.sif, na=r.na )
  out$args <- c( mget( names( formals() ), env=as.environment( -1 ) ), ## Store the function call's arguments, incl
                sapply( as.list( substitute( { ... } )[ -1 ] ), deparse ) ) ## expanding the '...' nifty trick, eh?

  out$e <- e
  try( sys.source( "cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )    
  try( sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )    
  
  if ( make.sif ) {
    ## Get a matrix of the # of times (total) that each pair of genes is in the same cluster
    net <- out$gene.gene.network()
    out$sif <- net$sif
    out$na <- net$na
    rm( net ); gc()
  }

  got.new.clusts <- integer()
  
  ## Get table of all motif widths
  motif.widths <- NULL
  if ( file.exists( "filehash/motif.widths.RData" ) ) load( "filehash/motif.widths.RData" )
  if ( is.null( motif.widths ) || nrow( motif.widths ) < e$k.clust ) {  ## Append if necessary
    seq.type <- 1 ##'upstream meme'
    ms <- e$meme.scores[[ seq.type ]]
    nr <- if ( ! is.null( motif.widths ) ) nrow( motif.widths ) else 0
    got.new.clusts <- (nr+1):e$k.clust
    widths1 <- do.call( rbind, mclapply( (nr+1):e$k.clust, function( i ) {
      if ( verbose && i %% 100 == 0 ) print( i )
      out <- c( 0, 0, 0 )
      if ( is.null( ms[[ i ]] ) || ms[[ i ]] == '' ) return( out )
      mm <- ms[[ i ]]$meme.out
      for ( j in 1:length( mm ) ) if ( ! is.null( mm[[ j ]]$width ) ) out[ j ] <- mm[[ j ]]$width
      out
    } ) )
    ##if ( ncol( widths1 ) > 2 ) widths1[ widths1[ ,2 ] == 0, 3:ncol( widths1 ) ] <- 0
    if ( all( widths1[ ,3 ] == 0 ) ) widths1 <- widths1[ ,1:2 ]
    motif.widths <- rbind( motif.widths, widths1 ); rm( widths1 )
    save( motif.widths, file="filehash/motif.widths.RData" )
  }
  out$motif.widths <- motif.widths

  ##try( unload( "filehash" ) )
  ##require( filehashRO )

  if ( find.bad ) { ## Find 'bad.motifs' -- motifs predominantly in coding regions 
    if ( ! exists( "pssm.scans" ) ) {
      if ( ! exists( "get.pssm.scans", envir=out ) ) {
        try( sys.source( "cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )
        try( sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )
      }

      ## Get mast pssm scans over entire genome - takes a while!
      require( data.table )
      pssm.scans <- out$get.pssm.scans( force=length( got.new.clusts ) > 0 ) ## Only for new clusters
      out$pssm.scans <- pssm.scans; rm( pssm.scans ); gc()
      
      ## p.scans <- pssm.scans
      ## require( data.table )
      ## ##p.scans$mots <- abs( p.scans$mots ); gc() ## Don't use strand or p-values for these funcs, so remove em
      ## ##p.scans <- p.scans[ ,c( "bic", "mots", "gene", "posns" ) ]; gc() ##, "fwd", "pvals" ) ]; gc()
      ## p.scans <- as.data.table( p.scans ); gc()
      ## setkey( p.scans, bic, mots, gene, posns )
    }
    
    if ( file.exists( "filehash/coding.fracs.RData" ) ) {
      load( "filehash/coding.fracs.RData", envir=out )
    }

    if ( ! file.exists( "filehash/coding.fracs.RData" ) || length( got.new.clusts ) > 0 ) { ## Append the fracs for the new clusters
      if ( length( got.new.clusts ) == 0 ) { got.new.clusts <- 1:e$k.clust; rm( coding.fracs, envir=out ) }
      mots <- unlist( lapply( got.new.clusts, function( i ) { ii <- which( out$motif.widths[ i, ] > 0 )
                                                              if ( length( ii ) <= 0 ) return( NULL )
                                                              paste( 'MOT', i, ii, sep='_' ) } ) )
      p.cutoff <- 1e-6 ## default for get.motif.coding.fracs()
      if ( exists( 'coding.fracs', envir=out ) ) p.cutoff <- out$coding.fracs$p.cutoff
      coding.fracs <- out$get.motif.coding.fracs( motifs=mots, verbose=verbose, p.cutoff=p.cutoff )
      if ( exists( 'coding.fracs', envir=out ) ) {
        out$coding.fracs$all.fracs <- c( out$coding.fracs$all.fracs, coding.fracs$all.fracs )
        out$coding.fracs$mean.fracs <- coding.fracs$mean.fracs
        coding.fracs <- out$coding.fracs
      }
      save( coding.fracs, file="filehash/coding.fracs.RData" )
      out$coding.fracs <- coding.fracs; rm( coding.fracs )
    }
    
    if ( length( got.new.clusts ) > 0 )
      system( "rm -r filehash/meme_hits_0.001.tsv.bz2 filehash/meme_posns_0.001*" ) ## pretty fast so just re-do completely
    meme.hits <- out$get.meme.genome.positions() ##get.meme.hits.table()
    setkey( meme.hits, gene )
    out$meme.hits <- meme.hits; rm( meme.hits ); gc()
  }

  ## Good params: p.value.cutoff=LOW (e.g. 1e-9) and mcl.I=4.5
  if ( cluster.motifs != FALSE ) {
    if ( n.cutoff > length( unique( names( out$e$fnames.to.cluster ) ) ) )
      n.cutoff <- length( unique( names( out$e$fnames.to.cluster ) ) )

    if ( FALSE ) { ## new "motif shadows"-based clustering
      try( sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs2.R", envir=out, chdir=T ) ) 
      fimo.out <- out$fimo.all.motifs( p.value.cutoff )
      fimo.out$Site <- as.character( fimo.out$Site ) ## seems to really slow it down if it is a factor
      out2 <- out$new.cluster.motifs( p.cutoff=p.value.cutoff, n.cutoff=n.cutoff, meme.the.clusters=F )
      out$motif.clusts <- out2; rm( out2 )
      out2 <- out$cluster.the.new.motif.clusters( out$motif.clusts )
      out$motif.cluster.clusters <- out2$motif.cluster.clusters; rm( out2 )
    } else {
      ## Orig tomtom-based clustering
      out2 <- out$cluster.motifs( cluster.option=cluster.motifs, e.value.cutoff=e.value.cutoff,
                                 p.value.cutoff=p.value.cutoff, resid.cutoff=resid.cutoff,
                                 min.gene.overlap=min.gene.overlap,
                                 get.motif.ccs=get.motif.ccs, find.bad=find.bad, n.cutoff=n.cutoff, ... )
      for ( i in names( out2 ) ) out[[ i ]] <- out2[[ i ]]
      rm( out2 ); gc()
    
      if ( get.motif.ccs ) { ## re-cluster the motif clusters into motif.cluster.clusters to find
        cat( 'Getting motif cluster clusters!\n' )
        out2a <- out$get.motif.cluster.clusters( cutoff=0.3 ) ## Note this cutoff worked for Halo but I used 0.2 for eco and yeast.
        out$ov.means <- Matrix( out2a$ov.means )
        out$motif.cluster.clusters <- out2a$motif.cluster.clusters
        rm( out2a ); gc()
      }
    
      if ( filehash ) {
        mclI <- gsub( '.', '', sprintf( "I%.1f", out$mcl.I ), fixed=T )
        print( mclI )
        try( dbCreate( sprintf( "filehash/tt.out2.%s", mclI ), type='RDS' ) )
        tt.out2 <- dbInit( sprintf( "filehash/tt.out2.%s", mclI ), type='RDS' )
        for ( i in 1:length( out$tt.out2 ) ) { print( i ); tt.out2[[ i ]] <- out$tt.out2[[ i ]] }
        out$tt.out2 <- tt.out2
        rm( tt.out2 ); gc()
      }
    }
  }
  
  out
}

## out is output from cmonkey.ensemble.analysis()
cluster.the.biclusters <- function( out, weight=c( g=0.6, m=0.1, c=0.3 ), mcl.I=4.5 ) { 
  e <- out$e

  if ( ! file.exists( 'clusterTheBiclusters.txt.gz' ) ) {
    ag <- lapply( 1:e$k.clust, function( k ) e$get.rows( k ) )
    tmp <- unique( unlist( ag ) ); gene.lookup <- 1:length( tmp ); names( gene.lookup ) <- tmp; rm( tmp )
    ag <- lapply( ag, function( i ) { out <- gene.lookup[ i ]; names( out ) <- NULL; out } )
    nrow <- sapply( ag, length )
    ac <- lapply( 1:e$k.clust, function( k ) e$get.cols( k ) )
    tmp <- unique( unlist( ac ) ); cond.lookup <- 1:length( tmp ); names( cond.lookup ) <- tmp; rm( tmp )
    ac <- lapply( ac, function( i ) { out <- cond.lookup[ i ]; names( out ) <- NULL; out } )
    ncol <- sapply( ac, length )
    if ( exists( "tt.out", envir=out ) ) {
      am <- mclapply( paste( "BIC", 1:e$k.clust, sep="_" ), function( k ) unlist( out$get.motifs( k ), use.names=F ) )
      nmot <- sapply( am, length )
    }
    
    ## require( Matrix )
    ## biclust.nw.g1 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )
    ## biclust.nw.g2 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )
    ## biclust.nw.c1 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )
    ## biclust.nw.c2 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )
    ## biclust.nw.m1 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )
    ## biclust.nw.m2 <- Matrix( data=0, nrow=e$k.clust, ncol=e$k.clust )

    ##fname <- 'clusterTheBiclusters_temp.txt'
    ##dir.create( "clusterTheBiclusters_temp" )
    tmp <- mclapply( 1:e$k.clust, function( k1 ) { ##for ( k1 in 1:e$k.clust ) {
      print( k1 )
      ##fname <- sprintf( "clusterTheBiclusters_temp/%08d.txt.gz", k1 )
      ##if ( file.exists( fname ) ) return( NULL )
      
      gg <- ag[[ k1 ]]
      ##nr <- ifelse( nrow < nrow[ k1 ], nrow, nrow[ k1 ] )
      nr <- unlist( lapply( ag, function( gs ) length( unique( c( gs, gg ) ) ) ) )
      ##biclust.nw.g1[ ,k1 ] <- nr
      nn1 <- unlist( lapply( ag, function( gs ) sum( gs %in% gg ) ), use.names=F ) ## / nr
      ##biclust.nw.g2[ ,k1 ] <- nn1
      ##nn1[ nn1 < quantile( nn1, 0.999 ) ] <- 0

      cc <- ac[[ k1 ]]
      ##nc <- ifelse( ncol < ncol[ k1 ], ncol, ncol[ k1 ] )
      nc <- unlist( lapply( ac, function( cs ) length( unique( c( cs, cc ) ) ) ) )
      ##biclust.nw.c1[ ,k1 ] <- nc
      nn2 <- unlist( lapply( ac, function( cs ) sum( cs %in% cc ) ), use.names=F ) ## / nc
      ##biclust.nw.c2[ ,k1 ] <- nn2
      ##nn2[ nn2 < quantile( nn2, 0.999 ) ] <- 0

      if ( exists( "am" ) ) {
        mm <- am[[ k1 ]]
        nn3 <- 0
        if ( ! is.null( mm ) ) {
          ##nm <- ifelse( nmot < nmot[ k1 ], nmot, nmot[ k1 ] )
          nm <- unlist( lapply( am, function( ms ) length( unique( c( ms, mm ) ) ) ) )
          ##biclust.nw.m1[ ,k1 ] <- nm
          nn3 <- unlist( lapply( am, function( ms ) sum( ms %in% mm ) ), use.names=F ) ##/ nm
          ##biclust.nw.m2[ ,k1 ] <- nn3
        } ##else {
        ##  nn3 <- rep( 0, length( am ) )
        ##}
      }
      
      ##nn <- nn1 + nn2 + nn3 ##nn1 * weight[ 'g' ] + nn2 * weight[ 'c' ] + nn3 * weight[ 'm' ]
      ##lock.file( fname )
      nn0 <- which( nn1 > 0 )
      out <- NULL
      if ( length( nn0 ) > 0 ) {
        if ( exists( "am" ) ) {
          ##write.table(
          out <- data.frame( k1, ##which( nn >= 0.05 ), nn1[ nn >= 0.05 ], nn2[ nn >= 0.05 ], nn3[ nn >= 0.05 ],
                            nn0, nn1[ nn0 ], nr[ nn0 ], nn2[ nn0 ], nc[ nn0 ], nn3[ nn0 ], nm[ nn0 ] ) ##, file=gzfile( fname ), quote=F, sep=" ",
          ##row.names=F, col.names=F, append=TRUE )
        } else {
          out <- data.frame( k1, nn0, nn1[ nn0 ], nr[ nn0 ], nn2[ nn0 ], nc[ nn0 ] )
        }
      }
      ## ##unlock.file( fname )
      ## ##print( file.info( "clusterTheBiclusters_temp.txt" )$size )
      ## ##system( sprintf( "gzip -v %s", fname ) )
      out ##NULL
    }, mc.preschedule=F )

    lapply( 1:length( tmp ), function( i ) write.table( tmp[[ i ]], file="filehash/clusterTheBiclusters.txt",
                                                       quote=F, sep=" ", row.names=F, col.names=F, append=(i!=1) ) )
    system( "gzip -v filehash/clusterTheBiclusters.txt" )
  }
  
  ## if ( ! exists( "weight" ) ) {
  ##   weight <- c( g=0.6, m=0.1, c=0.3 ) 
  ##   weight <- c( g=0.4, m=0.2, c=0.1 ) ## GOOD!
  ##   mcl.I <- 4.5
  print( weight ); print( mcl.I )
  ##}
  
  ## how to combine the columns into a single weight (OLD):
  formula <- sprintf( '$3*%.2f+$4*%.2f+$5*%.2f', weight['g'], weight['c'], weight['m'] )
  formula <- sprintf( '($3/$4-($4-$3)/$4)*%.2f+($5/$6-($6-$5)/$6)*%.2f+($7/$8-($8-$7)/$8)*%.2f',
                     weight['g'], weight['c'], weight['m'] )
  formula <- sprintf( '($3/$4*%.2f)+($5/$6*%.2f)+($7/$8*%.2f)', weight['g'], weight['c'], weight['m'] )

  if ( ! exists( "am" ) ) {
    formula <- sprintf( '($3/$4*%.2f)+($5/$6*%.2f)', weight['g'], weight['c'] )
  }

  out.fname <- sprintf( "clusterTheBiclusters.mcl.g.%.1f.m.%.1f.c.%.1f.new.I%s.out",
                       weight['g'], weight['m'], weight['c'],
                       gsub('.','',sprintf("%.1f",mcl.I),fixed=T) )
  print( out.fname )
  
  ## mcl liks symmetric graphs, so symmetrize it using the -tf option
  system( sprintf( "gunzip -c filehash/clusterTheBiclusters.txt.gz | awk '{print $1,$2,%s}' | awk '($3>0.1){print}' | ./progs/mcl-10-201/local/bin/mcl - --abc -o %s -I %.1f -v all -te 3 -S 200000 -tf '#max()'", formula, out.fname, mcl.I ) )

  system( sprintf( "gzip -v %s", out.fname ) )
  biclust.clusts <- lapply( strsplit( readLines( gzfile( sprintf( "%s.gz", out.fname ) ) ), "\t" ), as.integer )
  biclust.clusts <- biclust.clusts[ sapply( biclust.clusts, length ) >= 5 ]
  biclust.clusts <- lapply( biclust.clusts, function( i ) paste( "BIC", i, sep="_" ) )
  attr( biclust.clusts, 'mcl.I' ) <- mcl.I

  if ( FALSE ) {
    bc.genes <- mclapply( biclust.clusts, function( i ) out$agglom( i, srcType='bicluster', 'gene', path=NULL ) )
    bc.genes <- lapply( bc.genes, function( i ) subset( i, p.value <= 0.01 ) )
    bc.conds <- mclapply( biclust.clusts, function( i ) out$agglom( i, srcType='bicluster', 'condition', path=NULL ) )
    bc.conds <- lapply( bc.conds, function( i ) subset( i, p.value <= 0.01 ) )
    bc.motifs <- mclapply( biclust.clusts, function( i ) out$agglom( i, srcType='bicluster', 'motif.cluster', path='motif' ) )
    bc.motifs <- lapply( bc.motifs, function( i ) subset( i, p.value <= 0.01 ) )
    save( biclust.clusts, bc.genes, bc.conds, bc.motifs, file=sprintf( "%s.RData", out.fname ) )

    tt.outs <- list()
    for ( ind in ( length( tt.outs ) + 1 ):length( biclust.clusts ) ) {
      ks <- as.integer( gsub( "BIC_", "", biclust.clusts[[ ind ]] ) )
      cat( ind, length( ks ), "\n" )
      e$parallel.cores <- 4
      tt.out <- e$motif.similarities.tomtom( ks, ks, desymm=F, e.value.cutoff=e.value.cutoff ) ##, ... )
      e$parallel.cores <- 1
      tt.out2 <- try( e$cluster.tomtom.results( tt.out, k.cut=0.99,
                                          e.value.cutoff=e.value.cutoff, p.value.cutoff=p.value.cutoff,
                                          resid.cutoff=resid.cutoff, n.cutoff=n.cutoff ) ) ##, ... )
      if ( class( tt.out2 ) != 'try-error' ) {
        tmp <- tt.out2[ sapply( tt.out2, function( i ) length( attr( i, "mot.names" ) ) ) > 0 ]
        attributes( tmp ) <- attributes( tt.out2 )
        tt.out2 <- tmp; rm( tmp )
      }
      tt.outs[[ ind ]] <- list( tt.out=tt.out, tt.out2=tt.out2 )
    }
    
    for ( ind in 1:length( tt.outs ) ) {
      print( ind )
      par( mfrow=c( 4, 4 ) )
      tt.out2 <- tt.outs[[ ind ]]$tt.out2
      if ( length( tt.out2 ) <= 0 ) next
      for ( i in 1:length( tt.out2 ) ) {
        print( i )
        tmp <- tt.out2[[ i ]]
        out$e$viewPssm( attr( tmp, 'combined.pssm' ), main=paste( i, length( attr( tmp, 'mot.names' ) ) ) )
      }
    }
  }
  
  biclust.clusts
}

## NOTE: used min+2se up until Apr. 28 !!!
do.nwInf <- function( ..., all.tfs, cv.choose="min+4se", n.boot=100, make.shortcuts=F ) {
  source( "~/scratch/biclust/nwInf/runnit.R", chdir=T )
  dir.create( "filehash/nwInf_coeffs", recursive=T )
  DEBUG <- TRUE
  debug.on()

  for ( k in 1:e$k.clust ) { ## re-do for ones that failed
    fname <- sprintf( "filehash/nwInf_coeffs/%08d.RData", k )
    print( fname )
    if ( file.exists( fname ) ) { ##&& file.info( fname )$size > 1 ) next
      tmp <- try( load( fname ) )
      if ( class( tmp ) != "try-error" ) next
    }
    DEBUG <- FALSE
    coeffs <- try( runnit.wrapper.halo( e, ks=k, predictors=all.tfs, cv.choose=cv.choose, tf.groups=Inf, alpha=0.8,
                                       n.boot=1, tau=0, r.cutoff=Inf, r.filter=Inf, weighted=T, funcs=NA,
                                       aic.filter=NA, plot=F, ... ) )
    if ( class( coeffs ) != 'try-error' ) save( coeffs, file=fname )
    rm( coeffs )
  }

  require( filehash )
  try( dbCreate( "filehash/all.coeffs.dump", type="RDS" ) ) ## indiv. file for each entry
  try( dbCreate( "filehash/all.coeffs.sm.dump" ) ) ## store small-ified coeffs in all.coeffs.sm, but not plottable
  all.coeffs <- dbInit( "filehash/all.coeffs.dump", type="RDS" ) ## indiv. file for each entry
  all.coeffs.sm <- dbInit( "filehash/all.coeffs.sm.dump" )
  for ( i in 1:e$k.clust ) {
    fname <- sprintf( "filehash/nwInf_coeffs/%08d.RData", i )
    rm( coeffs )
    if ( file.exists( fname ) ) {
      print( fname )
      load( fname )
      if ( class( coeffs ) == "try-error" ) message( "BAD! ", i )
      all.coeffs[[ i ]] <- coeffs
      for ( n in c( "cluster.conds", "coeffs.boot", "coef.quantiles", "all.inputs", "plot.info",
                   "pred.ss", "pred.ts", "observed", "n.boot", "boot.opt" ) ) coeffs[[ 1 ]][[ n ]] <- NULL
      all.coeffs.sm[[ i ]] <- coeffs
    }
  }

  coeffs.to.biclust <- NULL
  if ( make.shortcuts ) {
    if ( ! file.exists( "filehash/coeffs.to.biclust.RData" ) && exists( "all.coeffs.sm" ) ) {
      coeffs.to.biclust <- list()
      for ( k in 1:e$k.clust ) {
        print( k )
        tmp <- all.coeffs.sm[[ k ]][[ 1 ]]
        for ( coe in names( tmp$coeffs ) ) coeffs.to.biclust[[ coe ]] <- c( coeffs.to.biclust[[ coe ]], k )
        if ( ! is.null( tmp$possibly.regulates ) )
          for ( coe in names( tmp$possibly.regulates ) ) coeffs.to.biclust[[ coe ]] <- c( coeffs.to.biclust[[ coe ]], k )
      }
      coeffs.to.biclust <- lapply( coeffs.to.biclust, unique )
      save( coeffs.to.biclust, file="filehash/coeffs.to.biclust.RData" )
    } else {
      try( load( "filehash/coeffs.to.biclust.RData" ) )
    }
  }
  
  invisible( list( all.coeffs.sm, all.coeffs, coeffs.to.biclust ) )
}

re.init <- function( make.shortcuts=T, filehash="filehash" ) {
  require( filehashRO )
  require( data.table )
  require( parallel )
  X11.options(type="dbcairo"); options(X11updates=0.25)
  try( source( "~/scratch/biclust/cmonkey-ensemble.R", chdir=T ) )
  e$meme.scores[[ 1 ]] <- dbInit( sprintf('%s/meme.scores.1.dump',filehash),
                                 type=if ( file.info( sprintf('%s/meme.scores.1.dump',filehash) )[ ,'isdir' ] ) 'RDS' else NULL )
  e$clusterStack <- dbInit( sprintf('%s/clusterStack.dump',filehash),
                           type=if ( file.info( sprintf('%s/clusterStack.dump',filehash) )[ ,'isdir' ] ) 'RDS' else NULL )
  if ( exists( 'out' ) ) {
    out$e <- e
    if ( ! is.null( out$mcl.I ) ) {
      mclI <- gsub( '.', '', sprintf( "I%.1f", out$mcl.I ), fixed=T )
      if ( file.exists( sprintf( "%s/tt.out2.%s",filehash, mclI ) ) ) {
        print( mclI )
        out$tt.out2 <- dbInit( sprintf( "%s/tt.out2.%s",filehash, mclI ), type='RDS' )
      }
    }
    ##if ( file.exists( sprintf("%s/motif.cluster.info" ) )
    ##  out$motif.cluster.info <- dbInit( sprintf("%s/motif.cluster.info", type="RDS" )
    if ( file.exists( sprintf("%s/all.coeffs.sm.dump",filehash) ) ) out$e.coeffs <- dbInit( sprintf("%s/all.coeffs.sm.dump",filehash) )
    if ( file.exists( sprintf("%s/all.coeffs.dump.RDS",filehash) ) ) {
      out$e.coeffs.big <- dbInit( sprintf("%s/all.coeffs.dump.RDS",filehash), type="RDS" ) ## indiv. file for each entry
      ##try( source( "~/scratch/biclust/nwInf/runnit.R", chdir=T ) ) ## Why was this here (6/4/12) ???
    }
    if ( file.exists( sprintf("%s/coeffs.to.biclust.RData",filehash) ) ) {
      load( sprintf("%s/coeffs.to.biclust.RData",filehash) )
      out$coeffs.to.biclust <- coeffs.to.biclust; rm( coeffs.to.biclust )
    }
  }
  if ( exists( 'out' ) ) try( sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )
  try( sys.source( "~/scratch/biclust/cmonkey-funcs.R", envir=e ) )
  if ( exists( 'out' ) ) try( sys.source( "~/scratch/biclust/cmonkey-funcs.R", envir=out$e, chdir=T ) )
  try( sys.source( "~/scratch/biclust/cmonkey-plotting.R", envir=e ) )
  if ( exists( 'out' ) ) try( sys.source( "~/scratch/biclust/cmonkey-plotting.R", envir=out$e, chdir=T ) )
  try( sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e ) )
  if ( exists( 'out' ) ) try( sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=out$e, chdir=T ) )

  if ( make.shortcuts && exists( 'out' ) ) {
    if ( ! exists( 'genes.to.biclust', envir=out ) ) {
      print( "genes.to.biclust" ) ## map from genes to biclusters to speed up get.biclusters(genes)
      if ( ! file.exists( sprintf("%s/genes.to.biclust.RData",filehash) ) ) {
        genes.to.biclust <- list()
        for ( k in 1:e$k.clust ) {
          if ( k %% 1000 == 0 ) print( k )
          tmp <- out$get.genes( bic=paste( "BIC", k, sep="_" ) )[[ 1 ]]
          for ( g in tmp ) genes.to.biclust[[ g ]] <- c( genes.to.biclust[[ g ]], k )
        }
        save( genes.to.biclust, file=sprintf("%s/genes.to.biclust.RData",filehash) )
      } else {
        load( sprintf("%s/genes.to.biclust.RData",filehash) )
      }
      out$genes.to.biclust <- genes.to.biclust
    }

    ## Map genes to motifs through meme.hits; but need to expand operons! (NEW as of May, 2013!)
    if ( ! exists( 'genes.to.motifs', envir=out ) ) {
      print( 'genes.to.motifs' )
      if ( ! file.exists( sprintf("%s/genes.to.motif.RData",filehash) ) ) {
        setkey( out$meme.hits, gene )
        o.list <- out$e$operon.list(); o.genes <- unlist( o.list )
        genes.to.motifs <- lapply( names( out$genes.to.biclust ), function( g ) {
          mh <- out$meme.hits[ J( g ) ]
          mots <- paste( 'MOT', mh$bic, mh$mot, sep='_' )
          cat( g, length(mots), "\n" )
          if ( g %in% o.genes ) { ## it's an operon gene; include other genes in the operon (if they're in the same bicluster)
            new.mots <- character()
            bics <- unlist( out$get.biclusters( gene=g ) )
            gg <- o.list[[ which( sapply( o.list, function( i ) any( i == g ) ) ) ]]
            gg <- gg[ gg != g ]
            for ( ggg in gg ) {
              mh2 <- out$meme.hits[ J( ggg ) ]
              mots2 <- paste( 'MOT', mh2$bic, mh2$mot, sep='_' )
              bics2 <- unlist( out$get.biclusters( motif=mots2 ) )
              mots2 <- mots2[ bics2 %in% bics ] ## get motifs for biclusters that also contain 'g'
              new.mots <- unique( c( new.mots, mots2 ) )
            }
            mots <- unique( c( mots, new.mots ) )
            cat( '\t', g, length(mots), "\n" )
          }
          return( mots )
        } )
        names( genes.to.motifs ) <- names( out$genes.to.biclust )
        save( genes.to.motifs, file=sprintf("%s/genes.to.motif.RData",filehash) )
      } else {
        load( sprintf("%s/genes.to.motif.RData",filehash) )
        out$genes.to.motifs <- genes.to.motifs
      }
    }

    if ( ! exists( 'motifs.to.motif.clusts', envir=out ) && exists( 'mc.length', envir=out ) ) {
      print( "motifs.to.motif.clusts" ) ## map from motifs to motif.clusts to speed up get.motif.clusters(motifs)
      if ( ! file.exists( sprintf("%s/motifs.to.motif.clusts.RData",filehash) ) ) {
        motifs.to.motif.clusts <- list()
        for ( k in 1:out$mc.length ) {
          tmp <- out$motif.clusts[[ k ]]
          for ( g in tmp ) motifs.to.motif.clusts[[ g ]] <-
            c( motifs.to.motif.clusts[[ g ]], paste( "MOTC", k, sep='_' ) )
        }
        save( motifs.to.motif.clusts, file=sprintf("%s/motifs.to.motif.clusts.RData",filehash) )
      } else {
        load( sprintf("%s/motifs.to.motif.clusts.RData",filehash) )
      }
      out$motifs.to.motif.clusts <- motifs.to.motif.clusts
    }

    if ( ! exists( 'motifs.to.genes', envir=out ) ) {
      print( "motifs.to.genes" ) ## map from motifs to genes to speed up get.genes(motif)
      if ( ! file.exists( sprintf("%s/motifs.to.genes.RData",filehash) ) ) {
        motifs.to.genes <- list()
        for ( g in rownames( e$ratios$ratios ) ) {
          print( g )
          tmp <- out$get.motifs( gene=g )[[ 1 ]]
          for ( m in tmp ) motifs.to.genes[[ m ]] <-
            c( motifs.to.genes[[ m ]], paste( g, sep='_' ) )
        }
        save( motifs.to.genes, file=sprintf("%s/motifs.to.genes.RData",filehash) )
      } else {
        load( sprintf("%s/motifs.to.genes.RData",filehash) )
      }
      out$motifs.to.genes <- motifs.to.genes
    }
        
    if ( ! exists( 'conds.to.biclust', envir=out ) ) {
      print( "conds.to.biclust" ) ## map from conds to biclusters to speed up get.biclusters(conditions)
      if ( ! file.exists( sprintf("%s/conds.to.biclust.RData",filehash) ) ) {
        conds.to.biclust <- list()
        for ( k in 1:e$k.clust ) {
          if ( k %% 1000 == 0 ) print( k )
          tmp <- out$get.conditions( bic=paste( "BIC", k, sep="_" ) )[[ 1 ]]
          for ( g in tmp ) conds.to.biclust[[ g ]] <- c( conds.to.biclust[[ g ]], k )
        }
        conds.to.biclust <- lapply( conds.to.biclust, unique )
        save( conds.to.biclust, file=sprintf("%s/conds.to.biclust.RData",filehash) )
      } else {
        load( sprintf("%s/conds.to.biclust.RData",filehash) )
      }
      out$conds.to.biclust <- conds.to.biclust
    }

    if ( ! exists( 'total.cond.count', envir=out ) ) {
      if ( ! file.exists( sprintf( '%s/total.cond.count.RData', filehash ) ) ) {
        total.cond.count <- table( unlist( lapply( e$clusterStack, '[[', 'cols' ) ) )
        save( total.cond.count, file=sprintf( '%s/total.cond.count.RData', filehash ) )
      }
      load( sprintf( '%s/total.cond.count.RData', filehash ), envir=out )
    }
  }
}

## Code for running multiple cmonkey runs on subsets of data, varying different parameters, saving each run
## assumes 'ratios' is already loaded
## code to run mpn ensemble:
## rm(list=ls());data(mpn,package='cMonkey.data');ratios=mpn$ratios;rm(mpn)
## source("cmonkey-ensemble.R");run.ensemble('mpn')
## To test code in current R session, set the 'iter' parameter; note all variables set get copied to .GlobalEnv !!
## Create a 'cleanup.function(fname)' to be run at end of cmonkey run - e.g. upload file to s3.
run.ensemble <- function( org, iter ) {
  if ( ! missing( iter ) ) {
    source( "cmonkey.R" ) ##require( cMonkey )
    dir.create( sprintf( 'ENSEMBLE_%s', org ) )
    if ( ! exists( 'ratios' ) ) load( sprintf( 'ENSEMBLE_%s/zzz_%s_ratios.RData', org, org ) )
    nc <- ncol( ratios )
    ratios <- ratios[ ,sample( 1:ncol( ratios ), sample( round(nc*0.3):round(nc*0.6), 1 ) ) ]
    print( dim( ratios ) )
    n.clust.per.row <- sample( 1:3, 1 )
    k.clust <- round( nrow( ratios ) * n.clust.per.row / 20 )
    k.clust <- sample( round(k.clust*0.6):round(k.clust*1.4), 1 )
    if ( sample( 1:2, 1 ) == 1 ) {
      motif.upstream.scan <- c( 0, 250 )
      motif.upstream.search <- c( 0, 150 )
      mot.weights <- c( `upstream.noncod meme`=1 )
    } else {
      motif.upstream.scan <- c( sample( (-50):0, 1 ), sample( 150:250, 1 ) )
      motif.upstream.search <- c( sample( (-20):0, 1 ), sample( 100:200, 1 ) )
    }
    net.weights <- c( string=runif( 1 ) * 0.5 + 0.2, operons=runif( 1 ) * 0.5 + 0.2 )
    operon.shift <- sample( 1:2, 1 ) == 1
    if ( sample( 1:2, 1 ) == 1 ) net.weights <- net.weights[ -which( names( net.weights ) == 'string' ) ]
    if ( sample( 1:2, 1 ) == 1 ) net.weights <- net.weights[ -which( names( net.weights ) == 'operons' ) ]
    if ( sample( 1:3, 1 ) == 1 && length( net.weights ) > 0 )
      seed.method <- c( rows=paste( 'net=', sample( names( net.weights ), 1 ), ':5', sep='' ), cols='rnd' )
    else if ( sample( 1:3, 1 ) == 1 ) seed.method <- c( rows='rnd', cols='rnd' ) ## otherwise default (kmeans)
    bg.order <- sample( 0:3, 1 )
    maxw <- sample( 12:30, 1 )
    meme.cmd <- "./progs//meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $none"
    meme.cmd <- gsub( "maxw 24", sprintf( "maxw %d", maxw ), meme.cmd )
    n.motifs <- sample( 1:3, 1 )
    n.iter <- 2000
    row.scaling <- sample( 4:8, 1 ) ## default 6
    net.scaling <- seq( 0, runif( 1 ) * 0.5 + 0.1, length=n.iter*3/4 ) ## default 0.5
    for ( i in ls() ) assign( i, get( i ), envir=.GlobalEnv )
    e <- cmonkey.init( organism=org, plot.iters=0 )
    if ( runif( 1 ) <= 0.3 ) { ## don't uniquify the sequences! (not legit, but often works better!)
      e$filter.sequences.orig <- e$filter.seqs
      e$filter.seqs <- function( seqs, start.stops=NULL, 
                                seq.type=paste( c("upstream","upstream.noncod","upstream.noncod.same.strand",
                                  "downstream","gene")[ 1 ], "meme" ), distance=motif.upstream.search[[ seq.type ]],
                                uniquify=T, remove.repeats=T, remove.atgs=T, mask.overlapping.rgns=F,
                                blast.overlapping.rgns=F, verbose=F, ... ) {
        filter.seqs.orig( seqs, start.stops, seq.type, distance, uniquify=F, remove.repeats=T, remove.atgs=T,
                         mask.overlapping.rgns=F, blast.overlapping.rgns=F, verbose=F, ... )
      }
      environment( e$filter.seqs ) <- e
    }
    cmonkey( e, dont.init=T )
    e <<- e ## Doesn't get saved unless I do this!
    save.image( file=sprintf( "ENSEMBLE_%s/zzz_%s_%03d.RData", org, org, as.integer( iter ) ) )
    rm( list=ls() )
    load( sprintf( "ENSEMBLE_%s/zzz_%s_%03d.RData", org, org, as.integer( iter ) ) )
    print( exists( 'e' ) )
    if ( exists( 'cleanup.function' ) )
      cleanup.function( file=sprintf( "ENSEMBLE_%s/zzz_%s_%03d.RData", org, org, as.integer( iter ) ) )
    stop()
  } else {
    if ( ! exists( 'ratios' ) ) stop( 'missing ratios!' )
    dir.create( sprintf( 'ENSEMBLE_%s', org ) )
    save( ratios, file=sprintf( 'ENSEMBLE_%s/zzz_%s_ratios.RData', org, org ) )
    for ( iter in 1:100 ) {
      cat( iter, "\n" )
      if ( file.exists( sprintf( "ENSEMBLE_%s/zzz_%s_%03d.RData", org, org, iter ) ) ) next
      cmd <- sprintf( "R CMD BATCH --no-save --no-restore --iter=%d --org=%s cmonkey-ensemble.R ENSEMBLE_%s/zzz_%s_%03d.Rout",
                     iter, org, org, org, iter )
      print( cmd )
      print( date() )
      system( cmd )
    }
  }
}

if ( length( grep( "--iter=", commandArgs(), fixed=T ) > 0 ) ) {
  iter <- commandArgs(); iter <- grep( "--iter=", iter, fixed=T, val=T )
  iter <- strsplit( iter, "=", fixed=T )[[ 1 ]][ 2 ]
  org <- commandArgs(); org <- grep( "--org=", org, fixed=T, val=T )
  org <- strsplit( org, "=", fixed=T )[[ 1 ]][ 2 ]
  run.ensemble( org, iter )
}

full.cmonkey.ensemble.analysis <- function( org='mpn', glob=sprintf( 'zzz_%s_???.RData', org ), cluster.motifs='mcl',
                                           p.value.cutoff=1e-5, ... ) {
  if ( ! exists( 'e', envir=.GlobalEnv ) ) e <- cmonkey.ensemble( glob, filehash=T )
  else e <- cmonkey.ensemble( glob, env=e, filehash=T )

  ## example for Halo:
  ## e=cmonkey.ensemble(list.files('jobs',recursive=T,full=T),filehash=T)
  rm( row.membership, col.membership, m.scores, n.scores, r.scores, rr.scores, row.memb, col.scores,
     c.scores, cc.scores, col.memb, cmonkey.init, write.project, old.row.membership,
     pre.adjusted.row.membership, old.col.membership, envir=e )

  if ( FALSE ) { ## TODO -- all subsequent analyses should ignore all ks in dupe.ks
    e$dupe.ks <- get.all.dupes( e )
  }
  
  save( e, file=sprintf( 'zzz_%s_ensemble_%d.RData', org, length(unique(names(e$fnames.to.cluster))) ) )

  sys.source( "~/scratch/biclust/cmonkey-motif-other.R", envir=e )
  e$progs.dir <- './progs/'
  e$parallel.cores=4
  ##e$genome.info$bg.fname <- NULL
  if ( is.null( e$genome.info$genome.seqs ) ) { ## might be the case if discard.genome==TRUE
    genome.info <- e$get.genome.info()
    e$genome.info$genome.seqs <- genome.info$genome.seqs
    rm( genome.info ); gc()
  }
  if ( is.list( e$genome.info$genome.seqs ) ) e$genome.info$genome.seqs <- unlist( e$genome.info$genome.seqs )
  if ( is.null( e$genome.info$bg.list ) ) e$genome.info$bg.list[[ 1 ]] <- e$mkBgFile( e$genome.info$all.upstream.seqs, order=e$bg.order[ 1 ],
                                                                                      bgfname=e$genome.info$bg.fname[ 1 ],
                                                                                      use.rev.comp=grepl( "-revcomp", e$meme.cmd[ 1 ] ) )
  save( e, file=sprintf( 'zzz_%s_ensemble_%d.RData', org, length(unique(names(e$fnames.to.cluster))) ) )
  
  save( e, file=sprintf( 'zzz_%s_ensemble_%d.RData', org, length(unique(names(e$fnames.to.cluster))) ) )

  if ( e$parallel.cores >= 4 ) e$parallel.cores <- e$parallel.cores.meme <- ceiling( e$parallel.cores / 2 )
  options( cores=e$parallel.cores )

  if ( FALSE ) {
    out1 <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs='mcl', min.gene.overlap=1 )
    out2 <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs='mcl', min.gene.overlap=1, mcl.I=4.5 )
    out3 <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs='mcl', min.gene.overlap=1, mcl.I=4.5,
                                      p.value.cutoff=0.001 )
    out <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs='mcl', min.gene.overlap=1, p.value.cutoff=1e-6, mcl.I=1.5 ) ## Final used for Ecoli
    out.hc <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs=T, min.gene.overlap=1 ) ## hclust requires too much memory
    out <- cmonkey.ensemble.analysis( e, make.sif=F )
    out <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs='mcl', min.gene.overlap=1, p.value.cutoff=1e-5 )
  }

  e$parallel.cores <- e$parallel.cores.motif <- parallel:::detectCores()
  options( cores=parallel:::detectCores() )

  if ( ! 'n.cutoff' %in% names(list(...)) ) { ## default is 10 but make it number.of.cmonkey.runs/10 (min 3)
    n.cutoff <- max( 3, round( length( unique( names( e$fnames.to.cluster ) ) ) / 10 ) )
  } else {
    n.cutoff <- list(...)$n.cutoff
  }
  cat( "Using n.cutoff =", n.cutoff, "\n" )
  
  out <- cmonkey.ensemble.analysis( e, make.sif=F, cluster.motifs=cluster.motifs, min.gene.overlap=1,
                                   p.value.cutoff=p.value.cutoff, n.cutoff=n.cutoff, ... )
  try( sys.source( "cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )
  try( sys.source( "~/scratch/biclust/cmonkey-ensemble-funcs.R", envir=out, chdir=T ) )
  save( e, out, file=sprintf( 'zzz_%s_ensemble_%d_w_pssm.scans.RData', org, length(unique(names(e$fnames.to.cluster))) ) )
  
  pdf( sprintf( 'zzz_%s_ensemble_%d_w_pssm.scans.pdf', org, length(unique(names(e$fnames.to.cluster))) ) )
  par( mfrow=c( 4, 4 ) )
  out$plot.motif.clusters( 1:out$mc.length )
  dev.off()

  ##if ( FALSE ) { ## to re-init the offline databases
  e <<- e
  out <<- out
  try( source( "cmonkey-ensemble.R" ) )
  try( source( "~/scratch/biclust/cmonkey-ensemble.R", chdir=T ) )
  re.init() ## see that func for more
  ##}

  out <- get( 'out', envir=.GlobalEnv )
  out$e <- e
  save( e, out, file=sprintf( 'zzz_%s_ensemble_%d_w_pssm.scans.RData', org,
                  length(unique(names(e$fnames.to.cluster))) ) )
  out
}

if ( FALSE ) {
  out = full.cmonkey.ensemble.analysis( 'mpn' )
}

## PLot a 'bicluster.cluster'' as output from 'cluster.the.biclusters()'
plot.regulon <- function( out, biclust.clusts, ind, p.cutoff.r=0.001, p.cutoff.c=0.1 ) {
  ##source( "~/scratch/biclust/heatmap3.R" )
## Modified version of heatmap2 to plot 'regulons' output by 'cluster.the.biclusters()'
  heatmap3 <-
    function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
              distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                               "row", "column", "none"), symm = FALSE, scale = c("none", 
                                              "row", "column"), na.rm = TRUE, revC = identical(Colv, 
                                                                                "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
              scale != "none", col = "heat.colors", colsep, rowsep, 
              sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
              notecol = "cyan", na.color = par("bg"), trace = c("column", 
                                                        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
              vline = median(breaks), linecol = tracecol, margins = c(5, 
                                                            5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
              cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
              key = TRUE, keysize = 1.5, density.info = c("histogram", 
                                           "density", "none"), denscol = tracecol, symkey = min(x < 
                                                                                     0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
              xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, 
              ...) 
      {
        scale01 <- function(x, low = min(x), high = max(x)) {
          x <- (x - low)/(high - low)
          x
        }
        retval <- list()
        scale <- if (symm && missing(scale)) 
          "none"
        else match.arg(scale)
        dendrogram <- match.arg(dendrogram)
        trace <- match.arg(trace)
        density.info <- match.arg(density.info)
        if (length(col) == 1 && is.character(col)) 
          col <- get(col, mode = "function")
        if (!missing(breaks) && (scale != "none")) 
          warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
                  "specified can produce unpredictable results.", "Please consider using only one or the other.")
        if (is.null(Rowv) || is.na(Rowv)) 
          Rowv <- FALSE
        if (is.null(Colv) || is.na(Colv)) 
          Colv <- FALSE
        else if (Colv == "Rowv" && !isTRUE(Rowv)) 
          Colv <- FALSE
        if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
          stop("`x' must be a numeric matrix")
        nr <- di[1]
        nc <- di[2]
        if (nr <= 1 || nc <= 1) 
          stop("`x' must have at least 2 rows and 2 columns")
        if (!is.numeric(margins) || length(margins) != 2) 
          stop("`margins' must be a numeric vector of length 2")
        if (missing(cellnote)) 
          cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
        if (!inherits(Rowv, "dendrogram")) {
          if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
                                                       c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
              dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                    dendrogram, "'. Omitting row dendogram.")
          }
        }
        if (!inherits(Colv, "dendrogram")) {
          if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
                                                       c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
              dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                    dendrogram, "'. Omitting column dendogram.")
          }
        }
        if (inherits(Rowv, "dendrogram")) {
          ddr <- Rowv
          rowInd <- order.dendrogram(ddr)
        }
        else if (is.integer(Rowv)) {
          hcr <- hclustfun(distfun(x))
          ddr <- as.dendrogram(hcr)
          ddr <- reorder(ddr, Rowv)
          rowInd <- order.dendrogram(ddr)
          if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Rowv)) {
          Rowv <- rowMeans(x, na.rm = na.rm)
          hcr <- hclustfun(distfun(x))
          ddr <- as.dendrogram(hcr)
          ddr <- reorder(ddr, Rowv)
          rowInd <- order.dendrogram(ddr)
          if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
        }
        else {
          rowInd <- nr:1
        }
        if (inherits(Colv, "dendrogram")) {
          ddc <- Colv
          colInd <- order.dendrogram(ddc)
        }
        else if (identical(Colv, "Rowv")) {
          if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
          if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
          }
          else colInd <- rowInd
        }
        else if (is.integer(Colv)) {
          hcc <- hclustfun(distfun(if (symm) 
                                   x
          else t(x)))
          ddc <- as.dendrogram(hcc)
          ddc <- reorder(ddc, Colv)
          colInd <- order.dendrogram(ddc)
          if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
        }
        else if (isTRUE(Colv)) {
          Colv <- colMeans(x, na.rm = na.rm)
          hcc <- hclustfun(distfun(if (symm) 
                                   x
          else t(x)))
          ddc <- as.dendrogram(hcc)
          ddc <- reorder(ddc, Colv)
          colInd <- order.dendrogram(ddc)
          if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
        }
        else {
          colInd <- 1:nc
        }
        retval$rowInd <- rowInd
        retval$colInd <- colInd
        retval$call <- match.call()
        x <- x[rowInd, colInd]
        x.unscaled <- x
        cellnote <- cellnote[rowInd, colInd]
        if (is.null(labRow)) 
          labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
          else rownames(x)
        else labRow <- labRow[rowInd]
        if (is.null(labCol)) 
          labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
          else colnames(x)
        else labCol <- labCol[colInd]
        if (scale == "row") {
          retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
          x <- sweep(x, 1, rm)
          retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
          x <- sweep(x, 1, sx, "/")
        }
        else if (scale == "column") {
          retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
          x <- sweep(x, 2, rm)
          retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
          x <- sweep(x, 2, sx, "/")
        }
        if (missing(breaks) || is.null(breaks) || length(breaks) < 
            1) {
          if (missing(col) || is.function(col)) 
            breaks <- 16
          else breaks <- length(col) + 1
        }
        if (length(breaks) == 1) {
          if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                          length = breaks)
          else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
          }
        }
        nbr <- length(breaks)
        ncol <- length(breaks) - 1
        if (class(col) == "function") 
          col <- col(ncol)
        min.breaks <- min(breaks)
        max.breaks <- max(breaks)
        x[x < min.breaks] <- min.breaks
        x[x > max.breaks] <- max.breaks
        if (missing(lhei) || is.null(lhei)) 
          lhei <- c(keysize, 4)
        if (missing(lwid) || is.null(lwid)) 
          lwid <- c(keysize, 4)
        if (missing(lmat) || is.null(lmat)) {
          lmat <- rbind(4:3, 2:1)
          if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
              stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                          1)
            lhei <- c(lhei[1], 0.2, lhei[2])
          }
          if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
              stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                               1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
          }
          lmat[is.na(lmat)] <- 0
        }
        if (length(lhei) != nrow(lmat)) 
          stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
        if (length(lwid) != ncol(lmat)) 
          stop("lwid must have length = ncol(lmat) =", ncol(lmat))
        op <- par(no.readonly = TRUE)
                                        #on.exit(par(op))
                                        #print(lmat)
                                        #print(lhei)
                                        #print(lwid)
                                        #print(margins)
        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
        if (!missing(RowSideColors)) {
          par(mar = c(margins[1], 0, 0, 0.5))
          image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        }
        if (!missing(ColSideColors)) {
          par(mar = c(0.5, 0, 0, margins[2]))
          image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
        }
        par(mar = c(margins[1], 0, 0, margins[2]))
        x <- t(x)
        cellnote <- t(cellnote)
        if (revC) {
          iy <- nr:1
          if (exists("ddr")) 
            ddr <- rev(ddr)
          x <- x[, iy]
          cellnote <- cellnote[, iy]
        }
        else iy <- 1:nr
        image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
              c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
              breaks = breaks, ...)
        retval$carpet <- x
        if (exists("ddr")) 
          retval$rowDendrogram <- ddr
        if (exists("ddc")) 
          retval$colDendrogram <- ddc
        retval$breaks <- breaks
        retval$col <- col
        if (!invalid(na.color) & any(is.na(x))) {
          mmat <- ifelse(is.na(x), 1, NA)
          image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
                col = na.color, add = TRUE)
        }
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
             cex.axis = cexCol)
        if (!is.null(xlab)) 
          mtext(xlab, side = 1, line = margins[1] - 1.25)
        axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
             cex.axis = cexRow)
        if (!is.null(ylab)) 
          mtext(ylab, side = 4, line = margins[2] - 1.25)
        if (!missing(add.expr)) 
          eval(substitute(add.expr))
        if (!missing(colsep)) 
          for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
                                                          length(csep)), xright = csep + 0.5 + sepwidth[1], 
                                    ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
                                    col = sepcolor, border = sepcolor)
        if (!missing(rowsep)) 
          for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                      1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                         1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                                    col = sepcolor, border = sepcolor)
        min.scale <- min(breaks)
        max.scale <- max(breaks)
        x.scaled <- scale01(t(x), min.scale, max.scale)
        if (trace %in% c("both", "column")) {
          retval$vline <- vline
          vline.vals <- scale01(vline, min.scale, max.scale)
          for (i in colInd) {
            if (!is.null(vline)) {
              abline(v = i - 0.5 + vline.vals, col = linecol, 
                     lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
          }
        }
        if (trace %in% c("both", "row")) {
          retval$hline <- hline
          hline.vals <- scale01(hline, min.scale, max.scale)
          for (i in rowInd) {
            if (!is.null(hline)) {
              abline(h = i + hline, col = linecol, lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
          }
        }
        if (!missing(cellnote)) 
          text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
               col = notecol, cex = notecex)
        par(mar = c(margins[1], 0, 0, 0))
        if (dendrogram %in% c("both", "row")) {
          plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
        }
        else plot.new()
        par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
        if (dendrogram %in% c("both", "column")) {
          plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
        }
        else plot.new()
        if (!is.null(main)) 
          title(main, cex.main = 1.5 * op[["cex.main"]])
        if (key) {
          par(mar = c(5, 4, 2, 1), cex = 0.75)
          tmpbreaks <- breaks
          if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
          }
          else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
          }
          z <- seq(min.raw, max.raw, length = length(col))
          image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
                xaxt = "n", yaxt = "n")
          par(usr = c(0, 1, 0, 1))
          lv <- pretty(breaks)
          xv <- scale01(as.numeric(lv), min.raw, max.raw)
          axis(1, at = xv, labels = lv)
          if (scale == "row") 
            mtext(side = 1, "Row Z-Score", line = 2)
          else if (scale == "column") 
            mtext(side = 1, "Column Z-Score", line = 2)
          else mtext(side = 1, "Value", line = 2)
          if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                  lwd = 1)
            axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
            title("Color Key\nand Density Plot")
            par(cex = 0.5)
            mtext(side = 2, "Density", line = 2)
          }
          else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                  col = denscol)
            axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
            title("Color Key\nand Histogram")
            par(cex = 0.5)
            mtext(side = 2, "Count", line = 2)
          }
          else title("Color Key")
        }
        else plot.new()
        retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                        high = retval$breaks[-1], color = retval$col)
        invisible(retval)
      }

  
  ## bonferroni cutoff w/ p-val=0.01
  ## expects biclust.clusts, bc.genes, bc.conds, bc.motifs to exist (see below)
  genes <- out$agglom( biclust.clusts[[ ind ]], srcType='bicluster', 'gene', NULL )
  genes <- subset( genes, q.value <= p.cutoff.r / ( length( biclust.clusts ) * nrow( e$ratios$ratios ) ) )
  genes <- genes[ order( genes$p.value ), ]
  conds <- out$agglom( biclust.clusts[[ ind ]], srcType='bicluster', 'condition', NULL )
  conds <- subset( conds, p.value <= p.cutoff.r ) ## .c/ ( length( biclust.clusts ) * ncol( e$ratios$ratios ) ) )
  conds <- conds[ order( conds$p.value ), ]
  rats <- e$ratios$ratios[ rownames( genes ), rownames( conds ) ]
  print( dim( rats ) )
  print( range( cor( t( rats ), use='pairwise' ), na.rm=T ) )
  mots <- out$agglom( biclust.clusts[[ ind ]], srcType='bicluster', 'motif.cluster', 'motif' )
  mots <- subset( mots, q.value <= p.cutoff.c / ( length( biclust.clusts ) * out$mc.length ) )
  ##mots <- mots[ order( mots$p.value ), ,drop=F ]
  print( mots )
  require( gplots )
  r.col <- colorRampPalette( c( "blue", "blue", "green", "yellow", "yellow" ) )
  tmp.r <- range( -round( log10(genes$p.value+1e-299) ) )
  g.col <- colorRampPalette( c( "white", "blue" ) )( tmp.r[ 2 ] )[ -round( log10(genes$p.value+1e-299) ) ]
  tmp.c <- range( -round( log10(conds$p.value+1e-299) ) )
  c.col <- colorRampPalette( c( "white", "blue" ) )( tmp.r[ 2 ] )[ -round( log10(conds$p.value+1e-299) ) ]
  lmat <- matrix( c( 6, 0, 0, 5,
                    0, 0, 0, 2,
                    4, 7, 1, 3 ), nrow=3, byrow=T )
  heatmap3( rats, trace="none", Rowv=F, Colv=T, cexRow=0.5, cexCol=0.5,
           col=r.col( 100 ), density.info="none", RowSideColors=g.col, ColSideColors=c.col,
           scale="none", dendrogram="none", lmat=lmat, lwid=c( 1.5, 0.5, 0.2, 4.0 ), lhei=c( 0.8, 0.2, 4.0 ) )
  if ( nrow( mots ) > 0 ) {
    m.mat <- matrix( 1, nrow=nrow( genes ), ncol=nrow( mots ) )
    rownames( m.mat ) <- rownames( genes ); colnames( m.mat ) <- rownames( mots )
    for ( g in rownames( genes ) ) {
      if ( FALSE ) { ## Only use motif.clusters linked from gene->motif->motif.cluster
        tmp <- out$get.motifs( gene=g )[[ 1 ]]
        tmp2 <- out$get.biclusters( motif=tmp )
        is.in <- sapply( tmp2, function( i ) i %in% biclust.clusts[[ ind ]] )
        m2 <- names( is.in[ is.in == TRUE ] )
        if ( length( m2 ) > 0 ) m <- out$agglom( m2, srcType='motif', 'motif.cluster', NULL )
        else m <- NULL
      } else {     ## This uses motif.clusters linked from gene->bicluster->motif->motif.cluster (less specific)
        b <- out$get.biclusters( gene=g )[[ 1 ]]
        b <- b[ b %in% biclust.clusts[[ ind ]] ]
        if ( length( b ) > 0 ) m <- out$agglom( b, srcType='bicluster', 'motif.cluster', 'motif' )
        else m <- NULL
      }
      if ( ! is.null( m ) ) {
        m <- m[ rownames( mots )[ rownames( mots ) %in% rownames( m ) ], ]
        m.mat[ g, rownames( m ) ] <- m$p.value
      }
    }
    m.col <- colorRampPalette( c( "white", "red" ) )( 299 )#[ -round( log10(genes$p.value+1e-299) ) ]
    par( mar=c( 4.5, 0, 0, 0.5 ) ) ## why 4.5 and not 5.0? I dont know but it seems to work.
    m.mat <- m.mat[ nrow( m.mat ):1, ,drop=F ]
    image( 1:ncol( m.mat ), 1:nrow( m.mat ), t( -log10( m.mat + 1e-299 ) ), col=m.col, axes=F, xlab="", ylab="" )
    axis( 1, ncol( m.mat ):1, labels=rev( colnames( m.mat ) ), las=2, line=-0.5, tick=0, cex.axis=0.5 )
  }
}

