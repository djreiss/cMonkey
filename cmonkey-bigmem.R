## Slurp an env's ff's and filehashes back into memory, save out image,
##    and optionally restore them back to the original ff's
save.cmonkey.env <- function( env=NULL, file=NULL, verbose=T ) { ##restore=T, 
  ## If 'env' is null, search for all env's in globalenv() and try to save those
  if ( is.null( env ) ) {
    for ( i in ls( env=globalenv() ) ) 
      if ( is.environment( get( i, globalenv() ) ) && "cmonkey" %in% class( get( i, globalenv() ) ) )
        save.cmonkey.env( get( i, globalenv() ), file, verbose ) ##restore, 
    return( invisible() )
  }

  if ( is.null( file ) ) file <- paste( env$cmonkey.filename, ".RData", sep="" )
  if ( verbose ) message( "Saving environment to ", file )
  save( env, file=file ) ## WARNING: assumes env exists in global environment!
  invisible( env )
}

matrix.reference <- function( m, ... ) {
    return( m ) ## This is for official package - don't deal with ref's yet (too slow)
}

list.reference <- function( l, file, ... ) {
  if ( ! big.memory || ! require( filehash ) ) return( l )
  if ( ! is.null( l ) && ( "filehashDB1" %in% class( l ) || length( l ) <= 0 ) ) return( l ) ##l <- dbInit( file )
  if ( big.memory && ! file.exists( cmonkey.filename ) ) dir.create( cmonkey.filename, recursive=T, show=F )
  if ( big.memory.verbose ) try( message( "Filehashing: ", file ), silent=T )
  if ( file.exists( file ) ) unlink( file )
  if ( is.null( l ) ) { dbCreate( file, ... ); l <- dbInit( file, ... ) } ## No list given, just create a database and initialize it.
  else l <- dumpList( l, file, ... )
  l
}

ffify.env <- function( env ) { ## Make all internal big matrices and lists disk-based
  for ( i in c( "row.scores", "mot.scores", "net.scores", "r.scores", "rr.scores", "col.scores",
               "c.scores", "cc.scores", "net.scores", "row.memb", "col.memb" ) ) {
    if ( exists( i, envir=env ) ) {
      tmp <- matrix.reference( env[[ i ]], backingfile=i, backingpath=env$cmonkey.filename )
      assign( i, tmp, envir=env )
    }
  }
  for ( i in names( env$ratios ) ) {
    env$ratios[[ i ]] <- matrix.reference( env$ratios[[ i ]], backingfile=paste( "ratios.", i, sep="" ),
                                          backingpath=env$cmonkey.filename )
  }
  for ( i in names( env$meme.scores ) ) {
    file <- paste( env$cmonkey.filename, "/meme.scores.", i, sep="" )
    if ( exists( "meme.scores", envir=env ) ) { 
      env$meme.scores[[ i ]] <- list.reference( env$meme.scores[[ i ]], file ) ##dumpList
    }
  }
  for ( i in c( "clusterStack", "genome.info", "networks" ) ) {
    file <- paste( env$cmonkey.filename, "/", i, sep="" )
    if ( exists( i, envir=env ) ) env[[ i ]] <- list.reference( env[[ i ]], file ) ## dumpList
  }
  invisible( env )
}

## NOTE: this should work whether the objects are on disk or not!
un.ffify.env <- function( env ) {
  for ( i in c( "row.scores", "mot.scores", "net.scores", "r.scores", "rr.scores", "col.scores", "net.scores",
               "cc.scores", "row.memb", "col.memb" ) ) 
    if ( exists( i, envir=env ) ) env[[ i ]] <- env[[ i ]][,]
  for ( i in names( env$ratios ) ) if ( ! is.null( env$ratios[[ i ]] ) ) ##&& is.big.matrix( env$ratios[[ i ]] ) )
    env$ratios[[ i ]] <- env$ratios[[ i ]][,]
  for ( i in names( env$meme.scores ) ) if ( ! is.null( env$meme.scores[[ i ]] ) )
    env$meme.scores[[ i ]] <- as.list( env$meme.scores[[ i ]] )
  for ( i in c( "clusterStack", "genome.info", "networks" ) ) if ( exists( i, envir=env ) )
    env[[ i ]] <- as.list( env[[ i ]] )
  invisible( env )
}
