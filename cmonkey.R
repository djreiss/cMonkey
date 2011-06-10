##################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

<<<<<<< HEAD
cm.version <- "4.8.2"
=======
cm.version <- "4.8.6"
>>>>>>> 95ebf69154f46fab7b09fd3687841cb03a88c68c


cmonkey <- function( env=NULL, ... ) {
  if ( ( ( is.null( list( ... )$dont.init ) || ! list( ... )$dont.init ) &&
      ( is.null( env$dont.init ) || ! env$dont.init ) && ( ! exists( "dont.init" ) || ! dont.init ) ) ||
      is.null( env ) || is.null( env$genome.info ) ) {
    env <- cmonkey.init( env, ... ) ##cog.org, rsat.species,
  } else {
    if ( sink.number() > 0 ) for ( i in 1:sink.number() ) try( sink(), silent=T ) ## Sometimes R gets into a funky sunky state
    if ( env$save.logfile != FALSE ) sink( env$save.logfile, split=T, append=T )
  }  
  ##for ( i in ls( envir=env ) ) assign( i, get( i, envir=env ) )
  cat( "\33[31mTIME STARTED:", env$time.started, "\33[0m\n" )

  ## SEEDING -- seed new clusters if no row.membership exists.
  ##    ALSO re-set the random seed if env$cmonkey.params$rnd.seed doesn't exist (e.g. for multiple runs on the
  ##       same pre-initialized env.) -- need to this via:
  ##       rm(list="rnd.seed",envir=env$cmonkey.params);rm(list="rnd.seed",envir=env)
  if ( ( ! exists( "row.membership", envir=env ) || all( env$row.membership == 0 ) ||
        nrow( env$row.membership ) != attr( env$ratios, "nrow" ) ||
        nrow( env$col.membership ) != attr( env$ratios, "ncol" ) ) && exists( 'ratios', envir=env ) )
    env$cmonkey.re.seed( env )
  
  ## iter <- env$iter
  ## while( iter <= env$n.iter ) {
  ##   env$iter <- iter
  ##   env$cmonkey.one.iter( env )
  ##   iter <- iter + 1
  ##   ##if ( attr( env$ratios, "nrow" ) > env$big.run ) gc() ## Clean up if we require lots of memory
  ## }

  while( env$iter <= env$n.iter ) {
    iter <- env$iter
    env$cmonkey.one.iter( env )
    ##env$iter <- env$iter + 1
    ##if ( attr( env$ratios, "nrow" ) > env$big.run ) gc() ## Clean up if we require lots of memory
  }

  ##for ( i in ls( envir=env ) ) if ( i != "env" ) assign( i, get( i, envir=env ) )
  
  if ( ! is.na( env$plot.iters ) && ( iter %in% env$plot.iters || ( iter - 1 ) %in% env$plot.iters ) )
    try( env$plotStats( iter, plot.clust=env$favorite.cluster() ) ) ## Plot final results ... can be set for your given organism

  env$iter <- iter <- env$iter - 1
  print( env$cluster.summary() ) ## Print out summary stats of each final cluster

  if ( env$post.adjust == TRUE ) { ## post-adjustment of all clusters???
    env$pre.adjusted.row.membership <- env$row.membership
    ##env2 <-
    env$adjust.all.clusters( env, expand.only=F )
    ##env <- env2; rm( env2 )
    gc()
##     row.membership.orig <- row.membership
##     adjust.all.clusters( expand=3 )
  }

  ## Get rid of extraneous links to old environments that are lying around (not really necessary anymore but doesnt hurt)
  parent.env( env ) <- globalenv(); parent.env( env$cmonkey.params ) <- env
  
  env$clusterStack <- env$get.clusterStack( ks=1:env$k.clust, force=T ) 
  print( env$cluster.summary() ) ## Print out summary stats of each final cluster
  env$set.param( "time.ended", date(), env$cmonkey.params )
  env$time.ended <- date()
  cat( "\33[31mTIME ENDED:", env$time.ended, "\33[0m\n" )

  ##if ( sink.number() > 0 ) for ( i in 1:sink.number() ) try( sink(), silent=T ) ## Sometimes R gets into a funky sunky state

##   for ( i in ls() ) if ( ! i %in% c( "i", "env" ) ) {
##     tmp <- get( i )
##     if ( class( tmp ) != "function" ) assign( i, tmp, envir=env )
##   }
  invisible( env )
}

DEBUG <- function( ... ) {
}

install.binaries <- function( meme.version="4.3.0",
                       url=sprintf( "http://meme.nbcr.net/downloads/old_versions/meme_%s.tar.gz", meme.version ) ) {
  cwd <- setwd( system.file( package="cMonkey" ) ); on.exit( setwd( cwd ) )
  if ( ! exists( "progs" ) ) dir.create( "progs" )
  setwd( "progs/" )
  cMonkey:::dlf( sprintf( "meme_%s.tar.gz", meme.version ), url )
  system( sprintf( "tar -xzf meme_%s.tar.gz", meme.version ) ); unlink( sprintf( "meme_%s.tar.gz", meme.version ) )
  setwd( sprintf( "meme_%s", meme.version ) ); dir.create( "local" )
  system( sprintf( "./configure --prefix=%s/local/ --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial --enable-build-libxml2 --enable-build-libxslt --disable-shared --enable-static --with-gnu-ld", getwd() ) )
  system( "make" ); system( "make install" )
  setwd( ".." )
  system( sprintf( "ln -s meme_%s/local/bin/meme", meme.version ) )
  system( sprintf( "ln -s meme_%s/local/bin/mast", meme.version ) )
  system( sprintf( "ln -s meme_%s/local/bin/dust", meme.version ) )
  ## Download and install blast executables? From:
  ## ftp.ncbi.nih.gov/blast/executables/LATEST
  setwd( cwd )
}
