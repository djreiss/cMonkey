###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

motif.one.cluster <- function( k, seq.type=names( mot.weights )[ 1 ], verbose=F, ##mask=T, blast=F, ##addl.args="", 
                             ##pseudocount=1/length(get.rows(k)), ms=meme.scores[[ seq.type ]][[ k ]], 
                             ##min.seqs=cluster.rows.allowed[ 1 ], max.seqs=cluster.rows.allowed[ 2 ],
                             ##pal.opt="non",
                             ... ) {
  st <- strsplit( seq.type, " " )[[ 1 ]]
  out <- NULL
  if ( st[ 2 ] == "meme" ) out <- meme.one.cluster( k, seq.type=seq.type, verbose, ... )
  invisible( out )
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
      else cat( k, if ( attr( out$meme.out, "is.pal" ) ) "pal" else "non", 
               sapply( out$meme.out[ 1:min( 3, length( out$meme.out ) ) ], "[[", "e.value" ),
               if ( ! is.null( out$pv.ev ) )
               mean( log10( out$pv.ev[[ 1 ]][ rownames( out$pv.ev[[ 1 ]] ) %in% get.rows( k ), "p.value" ] ), na.rm=T )
               else 'Inf', '\t', pssm.to.string( out$meme.out[[ 1 ]]$pssm ), "\n" )
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
  for ( k in 1:k.clust ) {
    m <- out.ms[[ k ]]
    if ( ! is.null( m ) && ! is.null( m$pv.ev ) ) out.ms[[ k ]]$pv.ev[[ 1 ]] <- NULL
  }
  
  ##if ( require( ff ) && object.size( out.ev ) / 1048600 > big.memory ) { ## Matrices > 50MB get filebacked
  ##  dir.create( cmonkey.filename, recursive=T, show=F )
  ##  out.ms$all.ev <- as.ff( out.ev, filename=paste( cmonkey.filename, "/all.ev.", seq.type, sep="" ), overwrite=T )
  ##} else out.ms$all.ev <- out.ev
  
  attr( out.ms, "seq.type" ) <- seq.type
  invisible( out.ms )
}  

meme.one.cluster <- function( k, seq.type=names( mot.weights )[ 1 ], verbose=F, ##mask=T, blast=F, ##addl.args="", 
                             ##pseudocount=1/length(get.rows(k)), ##ms=meme.scores[[ seq.type ]][[ k ]], 
                             ##min.seqs=cluster.rows.allowed[ 1 ], max.seqs=cluster.rows.allowed[ 2 ],
                             ##pal.opt="non",
                             ... ) { 
  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, verbose=verbose, ... ) ##mask=mask, blast=blast,
  min.seqs <- cluster.rows.allowed[ 1 ]; max.seqs <- cluster.rows.allowed[ 2 ]
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) ##meme.seqs.allowed[[ seq.type ]][ 1 ] )
    return( list( k=k ) ) 
  ##uniq <- uniquify.seqs[ seq.type ]
  ##if ( uniq ) seqs <- seqs[ ! get.dup.seqs( seqs ) ] 
  if ( length( seqs ) < min.seqs || ##meme.seqs.allowed[[ seq.type ]][ 1 ] ||
      length( seqs ) > max.seqs ) ##meme.seqs.allowed[[ seq.type ]][ 2 ] )
    return( list( k=k ) ) 
  meme.out <- mast.out <- NULL
  all.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  if ( is.null( all.seqs ) ) all.seqs <- get.sequences( "all", seq.type=seq.type,
                                                       distance=motif.upstream.scan[[ seq.type ]], filter=F, ... )
  bg.list <- genome.info$bg.list[[ seq.type ]]
  if ( is.null( bg.list ) && ! is.na( bg.order[ seq.type ] ) ) {
    tmp.seqs <- all.seqs[ ! names( all.seqs ) %in% rows ]
    ##if ( uniq ) tmp.seqs <- tmp.seqs[ ! get.dup.seqs( seqs ) ] 
    capture.output( bg.list <- mkBgFile( tmp.seqs, order=bg.order[ seq.type ],
                                        use.rev.comp=grepl( "-revcomp", meme.cmd[ seq.type ] ) ) ) ##addl.args[ i ] ) ) )
    rm( tmp.seqs )
  }

##  pal.opt <- motif.palindrome.option[ seq.type ]
  ##addl.args <- paste( addl.args, meme.addl.args[ seq.type ] ) 
  ##addl.args <- sprintf( addl.args, n.motifs[[ seq.type ]][ iter ], min.motif.width[[ seq.type ]][ iter ], max.motif.width[[ seq.type ]][ iter ] )
  cmd <- sprintf( meme.cmd[ seq.type ], n.motifs[[ seq.type ]][ iter ] ) ##, min.motif.width[[ seq.type ]][ iter ], max.motif.width[[ seq.type ]][ iter ] )

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )

  ms <- meme.scores[[ seq.type ]][[ k ]]
  get.meme.consensus <- function( cmd, min.iter=500, max.e.value=0.1, ... ) {
    ##consensus <- meme.consensus[ seq.type ]
    if ( grepl( "-cons $compute", cmd, fixed=T ) ) {
      ##if ( is.na( consensus ) ) consensus <- if ( iter > 500 ) "compute" else ""
      ##if ( consensus == "compute"
      if ( iter > min.iter && ! is.null( ms ) && ! is.null( ms$meme.out ) ) {
        e.val <- sapply( 1:length( ms$meme.out ), function( i ) ms$meme.out[[ i ]]$e.value )
        if ( min( e.val, na.rm=T ) < max.e.value ) {
          best <- which.min( e.val )
          consensus <- toupper( pssm.to.string( ms$meme.out[[ best ]]$pssm ) )
          cmd <- gsub( "$compute", consensus, cmd, fixed=T )
        }
      }
    }
    ##if ( ! is.null( consensus ) && ! is.na( consensus ) && consensus != "compute" && consensus != "" )
    ##  addl.args <- paste( addl.args, "-cons", consensus )
    if ( grepl( "-cons $compute", cmd, fixed=T ) ) cmd <- gsub( "-cons $compute", "", cmd, fixed=T )
    else if ( grepl( "-cons $none", cmd, fixed=T ) ) cmd <- gsub( "-cons $none", "", cmd, fixed=T )
    cmd
  }

  cmd <- get.meme.consensus( cmd, ... )

  ##if ( ( ! "psps" %in% names( list( ... ) ) || is.null( psps ) ) && exists( "get.sequence.psps" ) )
    psps <- NULL

  ## Options: normal use: have either "-pal=non" or no pal option in meme.cmd
  ##          force pal: have either "-pal=pal" or "-pal" option in meme.cmd
  ##   Use best (e-value) of pal/non-pal: have "-pal=both" in meme.cmd
  pal.opt <- "non"
  if ( grepl( "-pal=non", cmd ) ) { ## Can have "-pal=non" or nothing for this case
    cmd <- gsub( "-pal=non", "", cmd )
  }
  
  run.meme <- function( sgenes, seqs, cmd, seq.type, ... ) { ## "non", "pal", or "both" (try both and use one w/ best E-value)
    if ( pal.opt == "non" ) {
      out <- runMeme( sgenes, seqs, cmd, seq.type=seq.type, ... )
    }
    out
  }

  ##run.meme <- runMeme
  ## if ( pal.opt == "both" ) run.meme <- runMemePalNonPal
  ## else if ( pal.opt == "pal" ) run.meme <- runMemePal

  if ( verbose ) {
    meme.out <- try( run.meme( names( seqs ), seqs, cmd, ##nmotif=n.motifs[[ seq.type ]][ iter ],
                              verbose=verbose, bg.list=bg.list, psps=psps, seq.type=seq.type, ... ) )
  } else {
    capture.output( meme.out <- try( run.meme( names( seqs ), seqs, cmd, ##nmotif=n.motifs[[ seq.type ]][ iter ], 
                                              verbose=verbose, bg.list=bg.list, psps=psps, seq.type=seq.type, ... ) ) )
  }

    meme.out2 <- getMemeMotifInfo( meme.out )
    attr( meme.out2, "meme.command.line" ) <- attr( meme.out, "meme.command.line" )
    ##meme.out2$is.pal <- FALSE
    attr( meme.out2, "is.pal" ) <- pal.opt == "pal" ##FALSE
    ##attr( meme.out2, "pal.nonpal.evals" ) <- NA
  
  if ( length( meme.out2 ) <= 0 ) return( list( k=k ) ) ## No significant motif was found
    
  if ( verbose ) mast.out <- try( runMast( meme.out, names( all.seqs ), all.seqs, verbose=verbose,
                                          ##addl.args=mast.addl.args[ seq.type ],
                                          seq.type=seq.type, bg.list=bg.list, ... ) ) 
  else capture.output( mast.out <- try( runMast( meme.out, names( all.seqs ), all.seqs,
                                                verbose=verbose, ##addl.args=mast.addl.args[ seq.type ],
                                                seq.type=seq.type, bg.list=bg.list, ... ) ) ) 

  pv.ev <- get.pv.ev.single( mast.out, rows )

  ## pv.ev <- NULL
  ## if ( length( grep( "Error reading log-odds matrix file", mast.out ) ) <= 0 && class( meme.out ) != "try-error" &&
  ##     class( mast.out ) != "try-error" && length( meme.out2 ) > 0 && length( mast.out ) > 0 ) {
  ##   pv.ev <- getMastPValuesAndEValues( mast.out, get.p.values=rows )
  ##   attr( pv.ev, "mast.command.line" ) <- attr( mast.out, "mast.command.line" )
  ##   if ( length( pv.ev ) > 0 && nrow( pv.ev[[ 1 ]] ) == 0 && nrow( pv.ev[[ 2 ]] ) == 0 ) {
  ##     pv.ev <- NULL
  ##   } else { ## New code allowing pv.ev data frames to be stored as named (possibly big-memory) matrix
  ##     for ( i in 1 ) { ##:length( pv.ev ) ) { ## For now don't bother w/ 2nd data.frame - it's small.
  ##       tmp <- as.matrix( pv.ev[[ i ]][ ,2:ncol( pv.ev[[ i ]] ) ] )
  ##       rownames( tmp ) <- pv.ev[[ i ]][ ,1 ]
  ##       pv.ev[[ i ]] <- tmp
  ##     }
  ##   }
  ## }
  invisible( list( k=k, meme.out=meme.out2, pv.ev=pv.ev ) )
}

get.pv.ev.single <- function( mast.out, rows ) {
  pv.ev <- NULL
  if ( length( grep( "Error reading log-odds matrix file", mast.out ) ) <= 0 &&
      class( mast.out ) != "try-error" && length( mast.out ) > 0 ) {
    pv.ev <- getMastPValuesAndEValues( mast.out, get.p.values=rows )
    attr( pv.ev, "mast.command.line" ) <- attr( mast.out, "mast.command.line" )
    if ( length( pv.ev ) > 0 && nrow( pv.ev[[ 1 ]] ) == 0 && nrow( pv.ev[[ 2 ]] ) == 0 ) {
      pv.ev <- NULL
    } else { ## New code allowing pv.ev data frames to be stored as named (possibly big-memory) matrix
      for ( i in 1 ) { ##:length( pv.ev ) ) { ## For now don't bother w/ 2nd data.frame - it's small.
        tmp <- as.matrix( pv.ev[[ i ]][ ,2:ncol( pv.ev[[ i ]] ) ] )
        rownames( tmp ) <- pv.ev[[ i ]][ ,1 ]
        pv.ev[[ i ]] <- tmp
      }
    }
  }
  pv.ev
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


## "seq.weights" is a named numeric vector that allows weighting of each *entire* input sequence.
##    They MUST be 0 < p <= 1.
## "psps" is a named list that gives "position-specific prior probabilities" for each residue in each sequence.
##    These must also be in the range [0,1] and sum to <= 1. Note that they say: "The last w-1 numbers for each
##    entry should be 0 (shown in blue in the example), since a motif of that width cannot start in those positions."
## See: http://meme.sdsc.edu/meme4_3_0/doc/meme.html
mkTempMemeFiles <- function( sgenes, seqs, fname="meme.tmp.fst", bgseqs=NULL, bgfname=NULL, filter.seqs=T,
                            bg.list=NULL, force.overwrite=F, seq.type=names( mot.weights )[ 1 ], seq.weights=NULL, psps=NULL, ... ) {

  if ( ! file.exists( fname ) || file.info( fname )$size == 0 || force.overwrite ) {
    sgenes <- sgenes[ ! ( is.na( seqs ) | is.null( seqs ) | seqs == "" ) ]
    seqs <- seqs[ ! ( is.na( seqs ) | is.null( seqs ) | seqs == "" ) ]
    max.width <- as.integer( strsplit( meme.cmd[ seq.type ], " " )[[ 1 ]][ which( strsplit( meme.cmd[ seq.type ],
                                                                                " " )[[ 1 ]] == "-maxw" ) + 1 ] )
    if ( filter.seqs ) {
      sgenes <- sgenes[ nchar( seqs ) >= max.width ] ##motif.width.range[[ seq.type ]][ 2 ] ]
      seqs <- seqs[ nchar( seqs ) >= max.width ] ##motif.width.range[[ seq.type ]][ 2 ] ]
    }
    lengths <- sum( nchar( seqs ) ) + length( seqs ) * 3
    if ( ! is.null( seq.weights ) ) {
      seq.weights <- seq.weights[ sgenes ]
      seq.weights[ is.na( seq.weights ) ] <- 0
      cat( paste( ">WEIGHTS", paste( seq.weights, collapse=" " ) ), paste( ">", sgenes, "\n", seqs, sep="" ),
          file=fname, sep="\n" )
    } else {
      cat( paste( ">", sgenes, "\n", seqs, sep="" ), file=fname, sep="\n" )
    }

  }
  
  if ( force.overwrite || ( ! is.null( bgfname ) && ( ! file.exists( bgfname ) || file.info( bgfname )$size <= 0 ) ) ) {
    if ( ! is.null( bg.list ) ) mkBgFile( input.list=bg.list, order=bg.list$order, bgfname=bgfname )
    else if ( ! is.null( bgseqs ) ) mkBgFile( bgseqs, order=0, bgfname=bgfname )
  }

  length( seqs )
}

## Thread-safe tempfile (so parallelized calls dont use same filename)
## Make sure to set .options.multicore=list( set.seed=T ) to make sure this is true!
my.tempfile <- function( pattern="file", tmpdir=tempdir(), suffix="", n.rnd.char=20 ) {
  ## f <- paste( tempfile( pattern, tmpdir ), "_", k, "_", iter, suffix, sep="" )
  ## while( file.exists( f ) ) f <- paste( tempfile( pattern, tmpdir ), "_", k, "_", iter, suffix, sep="" )
  ## cat( "", file=f, append=F )
  ## f
  file.path( paste( tmpdir, "/", pattern, "_", paste( sample( c( LETTERS, letters, 0:9, 0:9, 0:9, 0:9 ),
                                                             n.rnd.char ), collapse="" ), suffix, sep="" ) )
}

runMeme <- function( sgenes, seqs, cmd=meme.cmd[ names( mot.weights )[ 1 ] ], bgseqs=NULL, bgfname=NULL, bg.list=NULL, 
                    nmotif=1, unlink=T, verbose=T, seq.weights=NULL, psps=NULL, ... ) { 
  fname <- my.tempfile( "meme.tmp", suf=".fst" ) 
  bgfname <- my.tempfile( "meme.tmp", suf=".bg" ) 
  tmp <- mkTempMemeFiles( sgenes, seqs, fname=fname, bgseqs=bgseqs, bg.list=bg.list,
                         bgfname=bgfname, seq.weights=seq.weights, psps=psps, ... )
  if ( tmp <= 0 ) return( NULL )

  if ( is.null( bgfname ) || ! file.exists( bgfname ) ) cmd <- gsub( '-bfile $bgFname', '', cmd, fixed=T )
  else cmd <- gsub( '$bgFname', bgfname, cmd, fixed=T )
  
  if ( is.null( psps ) ) cmd <- gsub( '-psp $pspFname', '', cmd, fixed=T )
  else cmd <- gsub( '$pspFname', sprintf( "%s.psp", fname ), cmd, fixed=T )

  cmd <- gsub( "$fname", fname, cmd, fixed=T )
  
  if ( verbose ) cat( cmd, "\n" )
  output <- system.time.limit( cmd )
  attr( output, "meme.command.line" ) <- cmd

  if ( unlink ) unlink( c( fname, bgfname, sprintf( "%s.psp", fname ) ) )
  return( output )
}

## Meme (very) occasionally sits indefinitely without running Even when -time option is set
## (I dont know when or why - is it when there are no sequences? or too many? or messed up input?).
## This func is a generic "system" call that will force it to be killed after "tlimit" seconds.
## NOTE: probably only runs on UNIX-y systems. NOTE: doesnt work, so we now use pipe() which seems a bit more stable than system()
system.time.limit <- function( cmd, tlimit=600 ) { ## default 10 minutes
  out <- readLines( pipe( cmd, "rt" ) )
  ##closeAllConnections()
  out
}

## memeOut can be direct text output of meme or parsed output containing pssms (which can be modified)
runMast <- function( memeOut, genes, seqs, bgseqs=NULL, bg.list=NULL, ##e.value.cutoff=99999, p.value.cutoff=0.5, 
                    ##motif.e.value.cutoff=99999, seq.weights=NULL, psps=NULL, 
                    unlink=T, verbose=F, ... ) {
  fname <- my.tempfile( "mast.tmp", suf=".fst" ) 
  bgfname <- my.tempfile( "mast.tmp", suf=".bg" ) 
  memeOutFname <- my.tempfile( "meme.tmp", suf=".out" ) 

  cat( memeOut, sep="\n", file=memeOutFname )
  tmp <- mkTempMemeFiles( genes, seqs, fname=fname, bgseqs=bgseqs, bg.list=bg.list,
                         bgfname=bgfname, seq.weights=NULL, psps=NULL, ... )
  if ( tmp <= 0 ) return( NULL )

  cmd <- mast.cmd
  if ( is.null( bgfname ) || ! file.exists( bgfname ) ) cmd <- gsub( '-bfile $bgFname', '', cmd, fixed=T )
  else cmd <- gsub( '$bgFname', bgfname, cmd, fixed=T )
  
  ##bgtmp <- paste( "-bfile", bgfname )
  ##if ( is.null( bgfname ) || ! file.exists( bgfname ) ) bgtmp <- ""
  ##if ( ! exists( "mast.cmd" ) ) mast.cmd <- "./progs/mast"
  
  ##cmd <- paste( mast.cmd, memeOutFname, "-d", fname, bgtmp, "-nostatus -stdout -text", "-brief", addl.args )
  cmd <- gsub( "$memeOutFname", memeOutFname, cmd, fixed=T )
  cmd <- gsub( "$fname", fname, cmd, fixed=T )
  
  ##cmd <- paste( cmd, addl.args )
  if ( verbose ) cat( cmd, "\n" )
  output <- system.time.limit( cmd )  
  attr( output, "mast.command.line" ) <- cmd

  if ( unlink ) unlink( c( memeOutFname, fname, bgfname ) )
  output
}

getMemeMotifInfo <- function( memeOutput ) {
  out <- list()
  lines <- grep( "^MOTIF\\s+\\d", memeOutput, perl=T )
  if ( length( lines ) <= 0 ) lines <- grep( "^MOTIF\\s+", memeOutput, perl=T )
  if ( length( lines ) > 0 ) {
    pssms <- getMemeMotifPssm( memeOutput, n.motif=length( lines ) )
    splitted <- strsplit( memeOutput[ lines ], "[\\t\\s]+", perl=T )
    for ( i in 1:length( lines ) ) {
### MOTIF  1	width=  17   sites= 15   llr=163   E-value=5.5e-002
      splt <- splitted[[ i ]]
      motif <- as.integer( splt[ 2 ] )
      width <- as.integer( splt[ 5 ] )
      sites <- as.integer( splt[ 8 ] )
      llr <- as.integer( splt[ 11 ] )
      e.value <- as.numeric( sub( "\\+", "", splt[ 14 ] ) )
      pssm <- pssms[[ motif ]]$pssm

      l2 <- grep( paste( "Motif", motif, "sites sorted by position p-value" ), memeOutput ) + 4
      l3 <- grep( "--------------------------------------------------------------------------------",
                 memeOutput[ (l2+1):length( memeOutput ) ] )[ 1 ] + l2 - 1
      posns <- do.call( rbind, strsplit( memeOutput[ l2:l3 ], "[\\t\\s]+", perl=T ) )[ ,c( 1:4, 6 ) ]
      colnames( posns ) <- c( "gene", "strand", "start", "p.value", "site" )
      posns <- data.frame( gene=posns[ ,"gene" ], strand=posns[ ,"strand" ], start=as.integer( posns[ ,"start" ] ),
                          p.value=as.numeric( posns[ ,"p.value" ] ), site=posns[ ,"site" ] )
      
      out[[ motif ]] <- list( width=width, sites=sites, llr=llr, e.value=e.value,
                             pssm=pssm, posns=posns )
    }
  }
  out
}

getMastPValuesAndEValues <- function( mastOutput, get.p.values=NULL ) {
  lines <- grep( "COMBINED P-VALUE", mastOutput )
  if ( length( lines ) > 0 ) {
    splitted <- strsplit( mastOutput[ lines ], "[\\t\\s]+", perl=T )
    out <- t( sapply( 1:length( lines ), function( i ) {
      gene <- mastOutput[ lines[ i ] - 2 ]
      splt <- splitted[[ i ]]
      p.val <- splt[ 8 ] 
      e.val <- splt[ 11 ] 
      c( gene=gene, p.value=p.val, e.value=e.val ) 
    } ) )
    out <- data.frame( gene=out[ ,"gene" ], p.value=as.numeric( out[ ,"p.value" ] ),
                      e.value=as.numeric( out[ ,"e.value" ] ) )
  } 

  out2 <- data.frame()
  if ( ! is.null( get.p.values ) && ! is.na( get.p.values ) ) {
    tmp <- get.mast.pvals( mastOutput, in.genes=get.p.values )
    for ( g in names( tmp ) ) {
      pv <- as.numeric( tmp[[ g ]]$pvals )
      pos <- as.integer( tmp[[ g ]]$posns )
      mots <- as.integer( tmp[[ g ]]$mots )
      if ( ! all( c( length( pv ), length( pos ) ) == length( mots ) ) ) ## BAD HACK In case the parsing messed up
        pv <- c( pv, rep( pv[1], length( pos ) - length( pv ) ) )
      out2 <- rbind( out2, data.frame( gene=g, pvals=pv, posns=pos, mots=mots ) )
    }
  }
  
  return( list( out, out2 ) )
}

getMemeMotifPssm <- function( memeOut, n.motif=1 ) {
  pssms <- list()
  for ( i in 1:n.motif ) {
    m.line1 <- grep( ##paste( "Motif ", i, " position-specific probability matrix", sep="" ), memeOut )
                    sprintf( "Motif %d position-specific probability matrix", i ), memeOut )
    if ( length( m.line1 ) > 0 ) {
      m.desc <- strsplit( memeOut[ m.line1 + 2 ], " " )[[ 1 ]]
      winLen <- as.numeric( m.desc[ 6 ] )
      e.val  <- as.numeric( m.desc[ 10 ] )
      pssm <- do.call( rbind, strsplit( memeOut[ m.line1 + 2 + 1:winLen ], "\\s+", perl=T ) )[ ,2:5 ]
      pssm <- matrix( as.numeric( pssm ), nrow=winLen, ncol=4, byrow=F )
      pssms[[ i ]] <- list( pssm=pssm, e.val=e.val )
    } else {
      pssms[[ i ]] <- list( pssm=NULL, e.val=99999 )
    }
  }
  return ( pssms )
}

get.mast.pvals <- function( mast.output, in.genes=NULL ) {
  space.pad <- function( lines, length ) {
    nc <- nchar( lines )
    nc[ nc >= length ] <- 0
    spaces <- sapply( 1:length( lines ), function( i ) paste( rep( " ", length - nc[ i ] ), sep="", collapse="" ) )
    paste( lines, spaces )
  }
  
  out <- list()
  start <- grep( "SECTION III: ANNOTATED SEQUENCES", mast.output )
  if ( length( start ) == 0 || is.na( start ) ) return( out )
  
  end <- grep( "\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*",
              mast.output[ (start+3):length(mast.output) ] ) + start + 3
  line.starts <- grep( "LENGTH = ", mast.output[ (start+2):(start+1+end) ] ) + start + 1
  if ( is.null( line.starts ) || length( line.starts ) == 0 ) return( out )

  for ( i in 1:length( line.starts ) ) {
    l <- line.starts[ i ]
    gene <- mast.output[ l - 2 ]
    if ( is.null( gene ) || is.na( gene ) || ( ! is.null( in.genes ) && ! ( gene %in% in.genes ) ) ) next

    l.next <- line.starts[ i + 1 ] - 2
    if ( i >= length( line.starts ) ) l.next <- end
    
    if ( l.next - l <= 5 ) next

    submast <- mast.output[ l:( l.next - 1 ) ]
    l.start <- which( submast == "" )[ 1 ] + 1
    if ( submast[ l.start ] == "" ) l.start <- l.start + 1
    q <- list()
    for ( i in 1:6 ) q[[ i ]] <- space.pad( submast[ seq( (l.start+i-1), length( submast ), by=6 ) ], 80 )
  
    seq.starts <- as.integer( sapply( strsplit( q[[ 5 ]], " " ), "[", 1 ) )
    
    char.skip <- which( strsplit( q[[ 5 ]][ 1 ], "" )[[ 1 ]] %in% c( 'G', 'A', 'T', 'C', 'N', 'X' ) )[ 1 ]
    mots <- unlist( strsplit( gsub( "[\\[\\]\\<\\>]", "", paste( substr( q[[ 1 ]], char.skip, 80 ), collapse="" ),
                                   perl=T ), "\\s+", perl=T ) )
    mots <- as.integer( mots[ ! is.na( as.integer( mots ) ) ] )
    mots <- mots[ ! is.na( mots ) ]
    p.vals <- strsplit( paste( substr( q[[ 2 ]], char.skip, 80 ), collapse="" ), "\\s+" )[[ 1 ]]
    p.vals <- as.numeric( p.vals[ ! is.na( as.numeric( p.vals ) ) ] )

    posns <- integer()
    for ( i in 1:length( q[[ 1 ]] ) ) {
      posns <- c( posns, which( strsplit( substr( q[[ 1 ]][ i ], char.skip, 80 ), "" )[[ 1 ]] %in%
                               c( '[', '<' ) ) + seq.starts[ i ] )
    }
    
    out[[ gene ]] <- list( pvals=p.vals, mots=mots, posns=posns )
  }

  return( out )
}

mkBgFile <- function( bgseqs=NULL, order=0, bgfname=NULL, input.list=NULL, use.rev.comp=T, verbose=T ) {
  if ( ! is.null( input.list ) && ! is.null( bgfname ) ) { ## Already ran this function, so lets just output the list to a file
    tmp <- unlist( input.list[ 2:length( input.list ) ] )
    tmp2 <- sprintf( "%.8f", tmp ) ##sapply( tmp, function( i ) sprintf( "%.8f", i ) )
    names( tmp2 ) <- names( tmp )
    write.table( tmp2, row.names=names( tmp2 ),
                col.names=paste( "#", order, "th order Markov background model" ), quote=F, file=bgfname )
    return( input.list )
  }

  ## Some seqs (e.g. hpy genome) have degenerate codes; remove those by sampling
  repl <- list( R=c("G","A"), Y=c("T","C"), K=c("G","T"), M=c("A","C"), S=c("G","C"), W=c("A","T"),
               N=c("G","A","T","C") )
  bad.seqs <- grep( "[^GATCX]", bgseqs, perl=T )
  if ( length( bad.seqs ) > 0 ) {
    if ( verbose ) message( length( bad.seqs ), " sequences with degenerate residues...fixing." )
    for ( i in bad.seqs ) {
      tmp <- strsplit( bgseqs[ i ], character(0) )[[ 1 ]]
      inds <- grep( "[^GATCX]", tmp, perl=T )
      for ( ind in inds ) tmp[ ind ] <- sample( repl[[ tmp[ ind ] ]], 1 )
      bgseqs[ i ] <- paste( tmp, collapse="" )
    }
  }
  
  if ( verbose ) cat( "Calculating", order, "th order background Markov model from", length( bgseqs ), "sequences\n" )
  if ( use.rev.comp && verbose ) cat( "Using reverse-complement too.\n" )
  
  bgseqs <- unique( bgseqs )
  if ( use.rev.comp ) bgseqs <- unique( c( bgseqs, rev.comp( bgseqs ) ) )
  
  mc <- get.parallel( order + 1 )
  apply.func <- lapply
  
  tmp <- mc$apply( 0:order, function( ord, mc.cores ) {
    out <- list()
    if ( verbose ) cat( "Calculating", ord, "th order part of background Markov model from", length( bgseqs ), "sequences\n" )
    
    if ( ord == 0 ) {
      all.substrings <- unlist( strsplit( bgseqs, character( 0 ) ), use.names=F )
    } else {
      all.substrings <- sapply( 1:( max( nchar( bgseqs ) ) - ord ), function( i ) substr( bgseqs, i, i + ord ) ) 
      all.substrings <- as.vector( all.substrings )
    }
    all.substrings <- all.substrings[ ! is.na( all.substrings ) & all.substrings != "" &
                                     nchar( all.substrings ) == ord+1 ]
    
    counts <- table( as.factor( all.substrings ) )
    counts <- sort( counts )
    counts <- counts / length( all.substrings )
    counts <- counts[ grep( "N", names( counts ), val=T, invert=T ) ] ## New R-2.9.0 grep param "invert" - faster?
    out <- as.list( counts )
    for ( i in names( out ) ) {
      names( out[[ i ]] ) <- NULL
      if ( verbose && ord <= 3 ) cat( "FREQ:", i, "=", counts[ i ], "\n")
    }
    out
  }, mc.cores=min( order + 1, mc$par ) )

  out <- list()
  out$order <- order
  for ( i in 1:length( tmp ) ) for ( j in 1:length( tmp[[ i ]] ) )
    out[[ names( tmp[[ i ]] )[ j ] ]] <- tmp[[ i ]][[ j ]]

  if ( ! is.null( bgfname ) && ! file.exists( bgfname ) ) {
    cat( "Writing to file:", bgfname, "\n" )
    tmp <- unlist( out )
    tmp <- tmp[ 2:length( tmp ) ]
    tmp2 <- sprintf( "%.8f", tmp ) ##sapply( tmp, function( i ) sprintf( "%.8f", i ) ) 
    names( tmp2 ) <- names( out )[ 2:length( out ) ]
    write.table( tmp2, row.names=names( tmp2 ),
                col.names=paste( "#", order, "th order Markov background model" ), quote=F, file=bgfname )
  }
  invisible( out )
}

col.let <- c( "A", "C", "G", "T" )

pssm.to.string <- function( pssm, cutoff.1=0.7, cutoff.2=0.4 ) {
  maxes <- max.col( pssm )
  letters <- col.let[ maxes ]
  values <- pssm[ cbind( 1:nrow( pssm ), maxes ) ]

  letters[ letters == "A" & values < cutoff.1 ] <- "a"
  letters[ letters == "C" & values < cutoff.1 ] <- "c"
  letters[ letters == "G" & values < cutoff.1 ] <- "g"
  letters[ letters == "T" & values < cutoff.1 ] <- "t"
  letters[ values < cutoff.2 ] <- "n"
  return( paste( letters, collapse="" ) )
}

pssm.to.consensus <- function( pssm, cutoff.1=0.8, cutoff.2=0.6, ##cutoff.3=0.4,
                              regex=T ) {
  c1 <- apply( pssm, 1, function( i ) which( i > cutoff.1 ) )
  c2 <- apply( pssm, 1, function( i ) which( i > cutoff.2 ) )
##  c3 <- apply( pssm, 1, function( i ) which( i > cutoff.3 ) )
  ca <- rbind( sapply( c1, length ), sapply( c2, length ) ) ##, sapply( c3, length ) )
  cond <- function( l ) paste( "[", paste( l, collapse="" ), "]", sep="" )
  deg.codes <- c( AG='R', GT='K', CG='S', CT='Y', AC='M', AT='W', CGT='B', ACT='H', AGT='D', ACG='V', ACGT='N' )
  out <- character()
  for ( i in 1:length( c1 ) ) {
    if ( ca[ 1, i ] == 1 ) {
      out[ i ] <- col.let[ c1[[ i ]] ]
    } else if ( ca[ 2, i ] >= 2 ) {
      if ( sum( pssm[ i, c2[[ i ]] ] ) > cutoff.1 ) out[ i ] <- cond( col.let[ c2[[ i ]] ] )
      else if ( sum( pssm[ i, c2[[ i ]] ] ) > cutoff.2 ) out[ i ] <- cond( tolower( col.let[ c2[[ i ]] ] ) )
    } else if ( ca[ 2, i ] >= 1 ) {
      if ( pssm[ i, c2[[ i ]] ] > cutoff.2 ) out[ i ] <- tolower( col.let[[ c2[[ i ]] ]] )
##    } else if ( ca[ 3, i ] >= 1 ) {
##      if ( sum( pssm[ i, c3[[ i ]] ] ) > cutoff.1 ) out[ i ] <- cond( col.let[ c3[[ i ]] ] )
##      else if ( sum( pssm[ i, c3[[ i ]] ] ) > cutoff.2 ) out[ i ] <- cond( tolower( col.let[ c3[[ i ]] ] ) )
    }
    if ( is.na( out[ i ] ) ) out[ i ] <- "[ACGT]"
    if ( ! regex && nchar( out[ i ] ) > 1 ) { ## Use IAUC degenerate codes
      tmp <- substr( out[ i ], 2, nchar( out[ i ] ) - 1 )
      if ( tmp %in% names( deg.codes ) ) out[ i ] <- deg.codes[ tmp ]
      else if ( toupper( tmp ) %in% names( deg.codes ) ) out[ i ] <- tolower( deg.codes[ toupper( tmp ) ] )
    }
  }
  paste( out, collapse="" )
}

read.fasta <- function( fname, lines=NULL ) {
  if ( is.null( lines ) ) lines <- readLines( fname )
  lines <- lines[ lines != "" ]
  starts <- grep( "^>", lines, perl=T )
  if ( length( starts ) > 1 ) stops <- c( starts[ 2:length( starts ) ], length( lines ) + 1 )
  else stops <- length( lines ) + 1
  seqs <- sapply( 1:length( starts ), function( i ) paste( lines[ ( starts[ i ] + 1 ):( stops[ i ] - 1 ) ],
                                                          collapse="", sep="" ) )
  names( seqs ) <- gsub( "^>", "", lines[ starts ], perl=T )
  closeAllConnections()
  seqs
}

remove.low.complexity <- function( seqs, length=8, entropy.cutoff=0.6, repl="N", use.dust=T,
                                  seq.type=names( mot.weights )[ 1 ] ) {
  write.fasta <- function( seqs, fname )
    writeLines( paste( paste( ">", names( seqs ), sep="" ), seqs, sep="\n" ), con=fname )

  if ( use.dust ) { ## Use "dust" on seqs - output to fasta file, run dust,
    ##if ( ! exists( "dust.cmd" ) ) dust.cmd <- "./progs/dust"
    ##if ( ! file.exists( dust.cmd ) ) warning( paste( "For best results, install", dust.cmd ), call.=F )
    ##!else {
    seqs <- seqs[ ! is.null( seqs ) & ! is.na( seqs ) ]
    max.width <- as.integer( strsplit( meme.cmd[ seq.type ], " " )[[ 1 ]][ which( strsplit( meme.cmd[ seq.type ],
                                                                                " " )[[ 1 ]] == "-maxw" ) + 1 ] )
    seqs <- seqs[ nchar( seqs ) >= max.width ] ##motif.width.range[[ seq.type ]][ 2 ] ]
    if ( length( seqs ) > 0 ) {
      fname <- my.tempfile( "dust", suf=".fst" ) ## parse input (fasta file format)
      write.fasta( seqs, fname )
      cmd <- gsub( "$fname", fname, dust.cmd, fixed=T )
      fst <- system.time.limit( paste( cmd, ##fname,
                                      "2>/dev/null" ), tlimit=60 ) ## prevent stderr output
      unlink( fname )
      if ( length( fst ) <= 1 ) cat( "WARNING: you probably don't have 'dust' installed.\n" )
      else seqs <- read.fasta( NULL, fst )
      return( seqs )
    }
    ##}
  }

}

rev.comp <- function( seqs ) { ## Fast reverse-complement
  sapply( seqs, function( seq ) paste( rev( strsplit( toupper( chartr( "ATCG", "tagc", seq ) ), "" )[[ 1 ]] ),
                                      collapse="" ) )
}

## seq.type can be any of those listed or e.g. 'file=asdfg.fst' then get seqs from fasta file asdfg.fst
get.sequences <- function( k, seq.type=paste( c("upstream","upstream.noncod","upstream.noncod.same.strand",
                                "downstream","gene")[ 1 ], "meme" ), verbose=F, filter=T,
                          distance=motif.upstream.search[[ seq.type ]], ... ) {
  if ( length( k ) <= 0 ) return( NULL )
  if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
  else rows <- k
  if ( is.null( rows ) ) return( NULL )
  start.stops <- NULL
  if ( is.na( seq.type ) || strsplit( seq.type, " " )[[ 1 ]][ 1 ] == "gene" ) op.shift <- FALSE
  seq.type <- seq.type; n.seq.type <- strsplit( seq.type, " " )[[ 1 ]][ 1 ]

  if ( substr( seq.type, 1, 5 ) == "file=" ) { ## if seq.type is e.g. 'file=asdfg.fst' then get seqs from fasta file
    seqs <- read.fasta( fname=strsplit( n.seq.type, "=" )[[ 1 ]][ 2 ] )
    seqs <- seqs[ rows ]; names( seqs ) <- rows
  } else {
    if ( is.null( genome.info$feature.tab ) ) stop( "Motif searching is on but no feature.tab!" )
    op.shift <- operon.shift[ seq.type ]
    
    coos <- get.gene.coords( rows, op.shift )
    if ( is.null( coos ) || nrow( coos ) <= 0 ) return( NULL )
    coos <- subset( coos, ! is.na( start_pos ) & ! is.na( end_pos ) )
    if ( is.null( coos ) || nrow( coos ) <= 0 ) return( NULL )
    seqs <- character()

    if ( n.seq.type %in% c( "upstream.noncod", "upstream.noncod.same.strand" ) ) {
      all.coos <- genome.info$feature.tab[ ,c( "id", "name", "contig", "strand", "start_pos", "end_pos" ) ]
      all.coos <- subset( all.coos, name %in% unlist( genome.info$synonyms ) )
    }

    ##len <- distance ##motif.upstream.search[[ n.seq.type ]]
    ##for ( i in 1:nrow( coos ) ) {
    mc <- get.parallel( nrow( coos ) )
    tmp <- mc$apply( 1:nrow( coos ), function( i ) {
      if ( n.seq.type == "gene" ) { ## Get the gene's seq
        st.st <- coos[ i, c( "start_pos", "end_pos" ), drop=F ]
      } else if ( n.seq.type == "upstream" ) { ## Get upstream
        st.st <- if ( coos$strand[ i ] == "D" ) ## "D" is forward, "R" is reverse (what does D stand for?)
          c( coos$start_pos[ i ] - 1 - distance[ 2 ], coos$start_pos[ i ] - 1 - distance[ 1 ] )
        else c( coos$end_pos[ i ] + 1 + distance[ 1 ], coos$end_pos[ i ] + 1 + distance[ 2 ] )
      } else if ( n.seq.type == "downstream" ) {
        st.st <- if ( coos$strand[ i ] == "D" ) ## "D" is forward, "R" is reverse (what does D stand for?)
          c( coos$end_pos[ i ] + 1 + distance[ 1 ], coos$end_pos[ i ] + 1 + distance[ 2 ] )
        else c( coos$start_pos[ i ] - 1 - distance[ 2 ], coos$start_pos[ i ] - 1 - distance[ 1 ] )
      } else if ( n.seq.type %in% c( "upstream.noncod", "upstream.noncod.same.strand" ) ) { ## Get upstream seq. but only up to previous gene
        cc <- all.coos[ as.character( all.coos$contig ) == as.character( coos$contig[ i ] ) &
                       abs( all.coos$start_pos - coos$start_pos[ i ] ) <= 100000, ]
        if ( n.seq.type == "upstream.noncod.same.strand" )
          cc <- all.coos[ as.character( all.coos$strand ) == as.character( coos$strand[ i ] ), ]
        if ( coos$strand[ i ] == "D" ) {
          nearest <- max( cc$end_pos[ cc$end_pos < coos$start_pos[ i ] ] )
          st.st <- c( nearest, coos$start_pos[ i ] - distance[ 1 ] - 1 )
        } else if ( coos$strand[ i ] == "R" ) {
          nearest <- min( cc$start_pos[ cc$start_pos > coos$end_pos[ i ] ] )
          st.st <- c( coos$end_pos[ i ] + distance[ 1 ] + 1, nearest )
        }
      }
      seq <- substr( genome.info$genome.seqs[[ as.character( coos$contig[ i ] ) ]], st.st[ 1 ], st.st[ 2 ] )
      if ( coos$strand[ i ] == "R" ) seq <- rev.comp( seq )
      if ( nchar( seq ) > abs( diff( distance ) ) ) {
        if ( coos$strand[ i ] == "D" ) seq <- substr( seq, 1, abs( diff( distance ) ) )
        else seq <- rev.comp( substr( rev.comp( seq ), 1, abs( diff( distance ) ) ) )
      }
      ##seqs[ as.character( coos$names[ i ] ) ] <- seq
      ##start.stops <- rbind( start.stops, data.frame( start=st.st[ 1 ], end=st.st[ 2 ],
      ##                                              strand=as.character( coos$strand[ i ] ),
      ##                                              contig=as.character( coos$contig[ i ] ) ) )
      ##rownames( start.stops ) <- ##[ nrow( start.stops ) ] <- as.character( coos$names[ i ] )
      ##  make.unique( c( rownames( start.stops )[ -nrow( start.stops ) ], as.character( coos$names[ i ] ) ) )
      out <- list( seq=seq, name=as.character( coos$names[ i ] ),
                  start.stops=data.frame( start=st.st[ 1 ], end=st.st[ 2 ],
                    strand=as.character( coos$strand[ i ] ),
                    contig=as.character( coos$contig[ i ] ) ) )
      out
    } )

    for ( i in tmp ) {
      seqs[ i$name ] <- i$seq
      start.stops <- rbind( start.stops, i$start.stops )
      ##rownames( start.stops ) <- make.unique( c( rownames( start.stops )[ -nrow( start.stops ) ], i$name ) )
      rownames( start.stops )[ nrow( start.stops ) ] <- i$name
    }
    rownames( start.stops ) <- names( seqs ) <- make.unique( rownames( start.stops ) )
    
    rows <- rows[ rows %in% names( seqs ) ]
    start.stops <- start.stops[ rows, ,drop=F ]
    seqs <- seqs[ rows ]
    names( seqs ) <- rownames( start.stops ) <- rows
  }
  
  if ( any( is.na( seqs ) ) ) {
    warning( "Warning: could not find '", n.seq.type, "' sequences for all input genes", call.=F )
    if ( ! is.null( start.stops ) ) start.stops <- start.stops[ ! is.na( seqs ), ]
    seqs <- seqs[ ! is.na( seqs ) ]
  }

  if ( filter ) seqs <- filter.sequences( seqs, start.stops, seq.type, distance, verbose=verbose, ... )
  
  attr( seqs, "start.stops" ) <- start.stops
  invisible( seqs )
}

get.dup.seqs <- function( seqs ) {
  out <- duplicated( seqs ) 
  names( out ) <- names( seqs )
  out
}

## TODO: need to provide options to filter/N out OTHER (upstream gene) ATGs and coding regions
filter.sequences <- function( seqs, start.stops=NULL, 
                             seq.type=paste( c("upstream","upstream.noncod","upstream.noncod.same.strand",
                               "downstream","gene")[ 1 ], "meme" ), distance=motif.upstream.search[[ seq.type ]],
                             uniquify=T, remove.repeats=T, remove.atgs=T, mask.overlapping.rgns=F,
                             blast.overlapping.rgns=F, verbose=F ) {

  if ( uniquify ) seqs <- seqs[ ! get.dup.seqs( seqs ) ]

  ##remove.repeats <- remove.low.complexity.subseqs[ seq.type ]
  if ( remove.repeats && ##! no.remove.repeats &&
      length( grep( "NNNNNN", seqs ) ) <= 1 ) {
    if ( verbose ) cat( "Removing low-complexity regions from sequences.\n" )
    seqs.new <- remove.low.complexity( seqs, seq.type=seq.type ) ## uses "dust" by default, now.
    if ( length( seqs.new ) == length( seqs ) ) seqs <- seqs.new
    else warning( "Remove low complexity failed - skipping!" )
    rm( seqs.new )
  }

  if ( remove.atgs && any( distance < 0 ) ) {
    tmp <- names( seqs )
    substr( seqs, distance[ 2 ] + 1, distance[ 2 ] + 4 ) <- "NNNN" ## Mask out ATGs
    names( seqs ) <- tmp
  }

  if ( mask.overlapping.rgns ) {
    if ( is.null( start.stops ) ) start.stops <- attr( seqs, "start.stops" )
    if ( ! is.null( start.stops ) ) { ## Mask out one copy of regions that are overlapping with other seqs (e.g. from divergent promoters) using coordinates (i.e. must ACTUALLY be overlapping on genome!)
      overlaps <- apply( start.stops, 1, function( i ) subset( start.stops, i[ 4 ] == contig &
                                                              ( i[ 1 ] >= start & i[ 1 ] <= end ) | ( i[ 2 ] >= start & i[ 2 ] <= end ) ) )
      overlaps <- lapply( names( overlaps ), function( g ) overlaps[[ g ]][ rownames( overlaps[[ g ]] ) != g, ] )
      names( overlaps ) <- rownames( start.stops )
      is.overlapping <- sapply( overlaps, nrow ); overlaps <- overlaps[ is.overlapping > 0 ]
      for ( i in names( overlaps ) ) {
        if ( nrow( overlaps[[ i ]] ) <= 0 ) next
        seq1 <- seqs[ i ]
        if ( start.stops[ i, 3 ] == "R" ) seq1 <- rev.comp( seq1 )
        ss1 <- sapply( 20:nchar( seq1 ), function( i ) substr( seq1, 1, i ) )
        ss2 <- sapply( 1:( nchar( seq1 ) - 20 ), function( i ) substr( seq1, i, nchar( seq1 ) ) )
        for ( j in 1:nrow( overlaps[[ i ]] ) ) {
          seq2 <- seqs[ rownames( overlaps[[ i ]] )[ j ] ]
          if ( overlaps[[ i ]][ j, 3 ] == "R" ) seq2 <- rev.comp( seq2 )

          g1 <- sapply( sapply( ss1, grep, seq2 ), length )
          rgn <- c( 1, nchar( seq1 ) )
          if ( all( g1 > 0 ) ) {
          } else if ( any( g1 > 0 ) ) {
            ind <- which( diff( g1 ) != 0 )
            rgn <- c( 1, ind - 1 )
          } else {
            g2 <- sapply( sapply( ss2, grep, seq2 ), length )
            if ( any( g2 > 0 ) ) {
              ind <- which( diff( g2 ) != 0 )
              rgn <- c( ind + 1, nchar( seq1 ) )
            }
          }

          if ( verbose ) cat( sprintf( "Masking region %d-%d of sequence %s (%s)\n", rgn[ 1 ], rgn[ 2 ], i,
                                      rownames( overlaps[[ i ]] )[ j ] ) )
          substr( seq1, rgn[ 1 ], rgn[ 2 ] ) <- paste( rep( "N", rgn[ 2 ] - rgn[ 1 ] + 1 ), collapse="" )
          seq <- seq1
          ##seq <- strsplit( seq1, "" )[[ 1 ]]
          ##seq[ rgn[ 1 ]:rgn[ 2 ] ] <- "N"
          ##seq <- paste( seq, collapse="" )
          if ( start.stops[ i, 3 ] == "R" ) seq <- rev.comp( seq )
          seqs[ i ] <- seq
          other.ov <- rownames( overlaps[[ i ]] )[ j ]
          overlaps[[ other.ov ]] <- overlaps[[ other.ov ]][ rownames( overlaps[[ other.ov ]] ) != i, ,drop=F ]
        }
      }
    }
  }

  
  if ( ! is.null( start.stops ) ) attr( seqs, "start.stops" ) <- start.stops[ names( seqs ), ,drop=F ]
  seqs
}
