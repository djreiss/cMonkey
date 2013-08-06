##################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

## Get genome data from RSAT (see http://rsat.ulb.ac.be/rsat/supported-organisms.cgi for all organisms)
## Use "http://rsat.ccb.sickkids.ca/" instead of http://rsat.ulb.ac.be/rsat/ (faster/more reliable?)
## Questions: use feature.tab instead of cds.tab? (includes trnas, etc) YES!!! ##Also get trnas via trna.tab?
##   NOTE: some org's (e.g. Human) don't have a feature.tab - need to use cds.tab then.

#' Load a ratios matrix either from a file or a variable
#'  
#' @param ratios  The input ratios matrix or the file name for a tab delinated file containging one
#' @param preprocess  Set to true to preprocess the ratios matrix (DEFAULT: T)
#' @param verbose  Set to FALSE to supress the output (DEFAULT: T)
#' @export
#' @usage ratios <- load.ratios( ratios, preprocess=T, verbose=T )
load.ratios <- function( ratios ) {
  if ( is.null( ratios ) ) return( NULL )

  if ( is.character( ratios ) && file.exists( ratios ) ) {
    cat( "Loading ratios file", ratios, "\n" )
    ratios <- read.delim( file=gzfile( ratios ), sep="\t", as.is=T, header=T, check.names=F )
  }

  if ( is.matrix( ratios ) || is.data.frame( ratios ) ) {
    if ( class( ratios[ ,1 ] ) == "character" ) {
      ratios <- ratios[ ! duplicated( ratios[ ,1 ] ), ] ## Remove second probe obs. for a few genes (?)
      rownames( ratios ) <- attr( ratios, "rnames" ) <- ratios[ ,1 ]; ratios <- ratios[ ,-1 ]
    }
    if ( class( ratios[ ,1 ] ) == "character" ) ratios <- ratios[ ,-1 ]
  }

  cat( "Original ratios matrix is", paste( dim( ratios ), collapse="x" ), "\n" )
  ## First filter out rows/cols with too many (>= 3/4 of total) NAs or zeroes (e.g. for sce data)
  if ( ! is.matrix( ratios ) ) ratios <- as.matrix( ratios )
  if ( is.null( attr( ratios, "isPreProcessed" ) ) || attr( ratios, "isPreProcessed" ) == FALSE ) {
    ratios <- preprocess.ratios( ratios )
    attr( ratios, "isPreProcessed" ) <- TRUE
  }
  cat( "Processed ratios matrix is", paste( dim( ratios ), collapse="x" ), "\n" )

  closeAllConnections()
  ratios
}

#' Scale the ratios matrix so that the profiles will be easily comparable
#'  Will apply "scale" to the matrix and may perform other actions.
#'  
#' @param ratios  The input ratios matrix
#' @param filter  Set to TRUE to filter any rows or columns that doesn't change (DEFAULT: T)
#' @param normalize  Set to TRUE to row center by median AND scale by sd (DEFAULT: T)
#' @param col.groups  A vector with a group number for each column (DEFAULT: NULL, i.e. set all column to group 1)
#' @param frac.cutoff  Used with "filter" to determine what constitutes "doesn't change" (DEFAULT: 0.98)
#' @param verbose  Set to FALSE to supress the output (DEFAULT: T)
#' @export
#' @usage ratios <- preprocess.ratios( ratios, filter=T, normalize=T, col.groups=NULL, frac.cutoff=0.98, verbose=T )
preprocess.ratios <- function( ratios, filter=T, normalize=T, col.groups=NULL, frac.cutoff=0.98 ) {
  if ( is.null( col.groups ) ) col.groups <- rep( 1, ncol( ratios ) )
  if ( is.null( names( col.groups ) ) ) names( col.groups ) <- colnames( ratios )
  ## for ( cg in unique( col.groups ) ) {
  ##   cols <- names( which( col.groups == cg ) )
  ##   ratios[ ,cols ] <- t( scale( t( ratios[ ,cols ] ), center=apply( ratios[ ,cols ,drop=F ], 1, median, na.rm=T ),
  ##                               scale=F ) )
  ## }
  if ( filter ) {
    cat( "Filtering out nochange rows/cols from ratios matrix...\n" )
    tmp1 <- apply( ratios, 1, function( i ) mean( is.na( i ) | abs( i ) <= 0.17 ) ) < frac.cutoff
    tmp2 <- apply( ratios, 2, function( i ) mean( is.na( i ) | abs( i ) <= 0.1 ) ) < frac.cutoff
    ratios <- ratios[ tmp1, ,drop=F ]
    ratios <- ratios[ ,tmp2 ,drop=F ]
    cat( "Filtered ratios matrix is", paste( dim( ratios ), collapse="x" ), "\n" )
    col.groups <- col.groups[ tmp2 ]
  }
  if ( normalize ) {
    for ( cg in unique( col.groups ) ) {
      cols <- names( which( col.groups == cg ) )
      cat( "Converting ratios matrix", cg, "to z-scores...\n" )
      ratios[ ,cols ] <- t( scale( t( ratios[ ,cols ,drop=F ] ),
                                  center=apply( ratios[ ,cols ,drop=F ], 1, median, na.rm=T ),
                 scale=apply( ratios[ ,cols ,drop=F ], 1, sd, na.rm=T ) ) ) ## Row center by median AND scale by sd!
    }
  }
  ratios
}


#' Download a file from the internet.  A wrapper for download.file()
#'  
#' @param f  The file name
#' @param url  The URL to download the file from
#' @param msg  A message to display on the console (DEFAULT: NULL )
#' @param mode  The mode with which to write the file. Useful values are "w", "wb" (binary), "a" (append) and "ab". (DEFAULT: 'wb' )
#' @param quiet  Set to FALSE to supress the output (DEFAULT: FALSE )
#' @param ...  Additional parameters for download.file()
#' @export
#' @usage err <- dlf( f, url, msg=NULL, mode="wb", quiet=F)
dlf <- function( f, url, msg=NULL, mode="wb", quiet=F, ... ) {
  err <- 0
  if ( mode == "ab" || ! file.exists( f ) || file.info( f )$size == 0 ) {
    if ( ! file.exists( dirname( f ) ) ) try( dir.create( dirname( f ), recursive=T ) )
    if ( ! is.null( msg ) ) cat( msg, "\n" )
    err <- try( download.file( url, destfile=f, mode=mode, quiet=quiet, ... ) )
  }
  closeAllConnections()
  err
}

#' Get the genome information.
#' Returns a list: ( species=rsat.species, genome.seqs=genome.seqs, feature.tab=feature.tab, 
#'       feature.names=feature.names,org.id=org.id, taxon.id=taxon.id, synonyms=synonyms )
#'  
#' @param fetch.upstream  (DEFAULT: F)
#' @param fetch.predicted.operons  (DEFAULT: "rsat")
#' @export
#' @usage err <- get.genome.info( fetch.upstream=F, fetch.predicted.operons="rsat" )
get.genome.info <- function( fetch.upstream=F ) { ##, fetch.predicted.operons="rsat" ) { 
  rsat.url <- rsat.urls[ 1 ]
  feature.tab <- feature.names <- genome.seqs <- operons <- org.id <- synonyms <- NULL 
  genome.loc <- paste( rsat.url, "/data/genomes/", rsat.species, "/genome/", sep="" )

  fname <- paste( "data/", rsat.species, "/organism_names.tab", sep="" )
  err <- dlf( fname, paste( genome.loc, "/organism_names.tab", sep="" ) )
  if ( class( err ) == "try-error" ) {
    tmp.url <- paste( rsat.url, "/data/genomes/", rsat.species, "_EnsEMBL/genome/organism_names.tab", sep="" )
    err <- dlf( fname, tmp.url ) 
    if ( class( err ) != "try-error" ) genome.loc <- paste( rsat.url, "/data/genomes/", rsat.species, "_EnsEMBL/genome/", sep="" )
  }
  if ( ! file.exists( fname ) || file.info( fname )$size <= 0 )
    stop( paste( "Genome info for", rsat.species, "does not exist. Please check", genome.loc,
                "and let me know if I am wrong" ) )
  nskip <- sum( substr( readLines( gzfile( fname ), n=20 ), 1, 2 ) == "--" |
               readLines( gzfile( fname ), n=20 ) == "" )
  org.id <- read.delim( gzfile( fname ), head=F, as.is=T, skip=nskip ) 
  if ( ! exists( "taxon.id" ) || is.na( taxon.id ) || is.null( taxon.id ) ) taxon.id <- org.id$V1[ 1 ]
  cat( "Organism taxon id:", taxon.id, "\n" ) 
  closeAllConnections()
  
  ## If we are getting all sequences from a fasta or csv file, skip downloading of genome seq info
  if ( ! no.genome.info && ! all( grepl( "file=", names( mot.weights ) ) ) ) {
    fname <- paste( "data/", rsat.species, "/feature.tab", sep="" )
    use.cds <- FALSE
    err <- dlf( fname, paste( genome.loc, "feature.tab", sep="" ),
               paste( "Fetching genome annotation data from RSAT", rsat.url, "..." ) )
    if ( class( err ) == "try-error" ) {
      err <- dlf( fname, paste( genome.loc, "cds.tab", sep="" ) )
      use.cds <- TRUE
    }
    cat( "Loading genome annotation data...\n" )
    head <- readLines( gzfile( fname ), n=30 ); nskip <- length( grep( '^--', head ) )
    feature.tab <- read.delim( gzfile( fname ), skip=nskip, head=F, comment='', as.is=F ) ##skip=16, 
    closeAllConnections()
    head <- strsplit( gsub( '^-- ', '', head[ grep( '^-- id', head, perl=T ) ], perl=T ), "\t" )[[ 1 ]]
    colnames( feature.tab ) <- head[ 1:ncol( feature.tab ) ]
    fname <- paste( "data/", rsat.species, "/feature_names.tab", sep="" )
    err <- dlf( fname, paste( genome.loc, if ( ! use.cds ) "feature_names.tab" else "cds_names.tab", sep="" ) )
    nskip <- sum( substr( readLines( gzfile( fname ), n=20 ), 1, 2 ) == "--" )
    closeAllConnections()
    feature.names <- read.delim( gzfile( fname ), head=F, as.is=T, skip=nskip, row.names=NULL, comment='' ) 
    closeAllConnections()
    colnames( feature.names ) <- c( "id", "names", "type" ) 
    feature.names <- unique( feature.names )
    chroms <- unique( as.character( feature.tab$contig ) )
    chroms <- chroms[ ! is.na( chroms ) & chroms != "" ]

    
    if ( ! is.na( mot.iters[ 1 ] ) ) {
      genome.seqs <- list()
      ##genome.seqs <- toupper( lapply(
      for ( i in chroms ) { ##, function( i ) {
        ## NOTE this is still bad as it loads entire genome into memory FIRST -- even when big.memory==TRUE
        cat( "Loading genome sequence, chromosome", i, "\n" )
        fname <- paste( "data/", rsat.species, "/", i, ".raw", sep="" )
        err <- dlf( fname, paste( genome.loc, i, ".raw", sep="" ) )
        if ( class( err ) == "try-error" ) {
          ii <- gsub( ":", "_", i, fixed=T )
          err <- dlf( fname, paste( genome.loc, ii, ".raw", sep="" ) )
          if ( class( err ) == "try-error" ) {
            err <- dlf( fname, paste( genome.loc, gsub( ".[0-9]$", "", i ), ".raw", sep="" ) )
            if ( class( err ) == "try-error" ) cat( "ERROR reading genome sequence", i, "\n" )
            else fname <- paste( "data/", rsat.species, "/", gsub( ".[0-9]$", "", i ), ".raw", sep="" )
          } else fname <- paste( "data/", rsat.species, "/", ii, ".raw", sep="" )
        }
        out <- try( readLines( gzfile( fname ) ), silent=T )
        ##closeAllConnections()
        if ( class( out ) == "try-error" || length( out ) == 0 || is.na( out ) || out == "" ||
            out == "CHARACTER(0)" ) out <- try( readLines( fname ), silent=T )
        if ( class( out ) == "try-error" || length( out ) == 0 || is.na( out ) || out == "" ||
            out == "CHARACTER(0)" ) {
          cat( "ERROR reading genome sequence", i, "\n" )
          next
        }
        if ( length( out ) > 1 ) out <- paste( out, collapse='', sep='' ) ## probably multiple lines
        out <- toupper( out )
        print( nchar( out ) )
        ##out } ) )
        genome.seqs[[ i ]] <- out
      }
      ##names( genome.seqs ) <- chroms
      ##genome.seqs <- genome.seqs[ ! ( is.na( gs ) | gs == "" | gs == "CHARACTER(0)" ) ]
      if ( length( genome.seqs ) != length( chroms ) ) {
        cat( "WARNING: Could not read sequence for chromosomes", chroms[ ! chroms %in% names( genome.seqs ) ], "\n" )
        feature.tab <- subset( feature.tab, contig %in% names( genome.seqs ) )
      }
      if ( length( genome.seqs ) <= 0 ) genome.seqs <- NULL
    }

  }

  synonyms <- NULL
  if ( exists( "synonym.thesaurus" ) ) {
    ## Load thesaurus file and turn it into a list of character arrays which are keyed with
    ## a probe id, and then have the synonyms as the elements of the character array.
    ## If an element doesn't exist for a gene then leave it out of the character array.
    cat( "Loading synonym information from custom thesaurus...\n" )
    tmp <- read.csv( gzfile( synonym.thesaurus ), header=F )
    synonyms <- strsplit( toupper( as.character( tmp[ ,2 ] ) ), ';' )
    names( synonyms ) <- toupper( as.character( tmp[ ,1 ] ) )
    rm( tmp )
  } else if ( exists( "ratios" ) && ! is.null( ratios ) ) {
    cat( "Gathering all \"standard\" orf names and other synonyms for all probe names...\n" )
    synonyms <- get.synonyms( attr( ratios, "rnames" ), feature.names, verbose=T ) 
  }

  if ( ! is.null( synonyms ) ) {
    cat( "Mean number of synonyms per probe:", mean( sapply( synonyms, length ), na.rm=T ), "\n" )
    ##synonyms <- tmp
    is.bad <- sapply( names( synonyms ), function( i ) length( synonyms[[ i ]] ) == 0 ||
                     substr( synonyms[[ i ]][ 1 ], 1, 5 ) == "Error" )
    if ( sum( is.bad ) > 0 ) {
      cat( "These", sum( is.bad ), "probe names have no matching ORF annotation:\n" )
      print( names( which( is.bad ) ) )
    }
    rm( is.bad ) ##tmp
  }
  
  closeAllConnections()
  
  invisible( list( species=rsat.species, genome.seqs=genome.seqs, feature.tab=feature.tab, feature.names=feature.names,
                  org.id=org.id, taxon.id=taxon.id, synonyms=synonyms ) ) 
}

## Can get operon predictions from RSAT (http://rsat.ccb.sickkids.ca/infer-operons.cgi).
## See http://rsat.ccb.sickkids.ca/help.infer-operons.html for how these are computed.
## or MicrobesOnline (list is at http://www.microbesonline.org/operons/OperonList.html).
get.operon.predictions <- function( fetch.predicted.operons="microbes.online", org.id=genome.info$org.id$V1[ 1 ] ) {
  operons <- NULL
  if ( fetch.predicted.operons == "rsat" ) {
    rsat.url <- rsat.urls[ 1 ]
    cat( "Using operon predictions from RSAT...\n" )
    fname <- paste( "data/", rsat.species, "/rsat_operon_predictions.html", sep="" )
    err <- dlf( fname, paste( rsat.url, "/infer-operons.cgi?organism=",
                             rsat.species, "&genes=all&return_leader=on&return_operon=on&return_query=on&",
                             "output=display&dist_thr=55", sep="" ) ) ##, "Fetching operon predictions from RSAT..." )
    operons <- readLines( gzfile( fname ) )
    start <- which( operons == "<INPUT type=\"hidden\" NAME=\"gene_selection\" VALUE=\"#lead\toperon\tquery" ) + 1
    end <- which( operons == "<INPUT type=\"hidden\" NAME=\"feattype\" VALUE=\"\">" ) - 2
    operons <- do.call( rbind, strsplit( operons[ start:end ], "\t+", perl=T ) )
    colnames( operons ) <- c( "lead", "operon", "query" ) 
    operons <- as.data.frame( operons )
  } else if ( fetch.predicted.operons == "microbes.online" ) {
    cat( "Using operon predictions from MicrobesOnline...\n" )
    fname <- paste( "data/", rsat.species, "/microbesonline_operons_gnc", org.id, ".named", sep="" )
    err <- dlf( fname, paste( "http://www.microbesonline.org/operons/gnc", org.id, ".named", sep="" ) ) ##,
               ##"Fetching operon predictions from MicrobesOnline..." )
    if ( org.id != taxon.id && ( ! file.exists( fname ) || file.info( fname )$size == 0 ) ) {
      fname <- paste( "data/", rsat.species, "/microbesonline_operons_gnc", taxon.id, ".named", sep="" )
      err <- dlf( fname, paste( "http://www.microbesonline.org/operons/gnc", taxon.id, ".named", sep="" ) ) ##, 
               ##"Fetching operon predictions from MicrobesOnline (2)..." )
    }
    if ( file.exists( fname ) ) cat( "Succesfully fetched operon predictions. Parsing...\n" )
    ops <- read.delim( gzfile( fname ) )
    ops2 <- ops ##subset( ops, bOp == "TRUE" & SysName1 != "" & SysName2 != "" )
    on <- as.character( ops2$SysName1 )[ -1 ]
    bOp <- ops2$bOp[ -1 ]
    on[ which( on == '' ) ] <- on[ which( on == '' ) - 1 ]
    ops2$SysName1 <- c( as.character( ops2$SysName1 )[ 1 ], on )
    on <- as.character( ops2$SysName2 )[ -1 ]
    on[ which( on == '' ) ] <- on[ which( on == '' ) - 1 ]
    ops2$SysName2 <- c( as.character( ops2$SysName2 )[ 1 ], on )
    ops2 <- subset( ops2, bOp == "TRUE" ) ##& SysName1 != "" & SysName2 != "" )    
    ##gns <- sort( unique( c( as.character( ops2$SysName1 ), as.character( ops2$SysName2 ) ) ) )
    ##gns <- gns[ gns != "" ]

    sn1 <- as.character( ops2$SysName1 )
    sn1[ sn1 == "" | is.na( sn1 ) ] <- as.character( ops2$Name1 )[ sn1 == "" | is.na( sn1 ) ]
    sn2 <- as.character( ops2$SysName2 )
    sn2[ sn2 == "" | is.na( sn2 ) ] <- as.character( ops2$Name2 )[ sn2 == "" | is.na( sn2 ) ]
    operons <- list( 0 )
    for ( i in 1:length( sn1 ) ) {
      sn1i <- sn1[ i ]
      if ( sn1i == '' ) next
      found <- which( sapply( operons, function( j ) sn1i %in% j ) )
      if ( length( found ) > 0 ) operons[[ found[ 1 ] ]] <- c( operons[[ found[ 1 ] ]], sn2[ i ] )
      else operons[[ length( operons ) + 1 ]] <- c( sn1i, sn2[ i ] )
    }
    operons <- operons[ -1 ]
    operons <- lapply( operons, unique )
    operons <- lapply( operons, function( i ) i[ i != '' ] )
    operons <- operons[ sapply( operons, length ) > 1 ]

    gns <- unlist( operons )
    search.names <- c( gns, as.character( genome.info$feature.names$id ) )
    if ( exists( "ratios" ) ) search.names <- c( attr( ratios, "rnames" ), search.names )
    mc <- get.parallel( length( operons ) )
    nms <- mc$apply( 1:length( operons ), function( i ) {
      s <- get.synonyms( operons[[ i ]] )
      s <- lapply( s, function( i ) i[ i %in% search.names ] )
      ids <- unlist( lapply( s, function( i ) i[ i %in% genome.info$feature.names$id ][ 1 ] ) )
      if ( length( ids ) <= 0 ) { 
        warning( paste( "No genome annotation for any genes in operon #", i, " -- don't know what to do!", call.=F ) )
        return( "" ) }
      ids[ is.na( ids ) ] <- names( ids )[ is.na( ids ) ]
      vngs <- unlist( lapply( s, function( i ) {
        out <- i[ ! i %in% genome.info$feature.names$id ]
        if ( length( out ) <= 0 ) ## && exists( "ratios" ) )
          out <- i[ i %in% search.names ] ##attr( ratios, "rnames" ) ]
        if ( length( out ) <= 0 ) out <- i[ genome.info$feature.names$id == i & genome.info$feature.names$id == "primary" ]
        if ( length( out ) <= 0 ) out <- i
        out
      } ) )
      coos <- get.gene.coords( ids, op.shift=F )
      vngs <- vngs[ ids %in% coos$names ]
      if ( is.null( coos ) || nrow( coos ) <= 0 ) { 
        warning( paste( "No genome annotation for any genes in operon #", i, " -- don't know what to do!", call.=F ) )
        return( "" ) }
      if ( mean( as.character( coos$strand ) == "D" ) > 0.6 ) head <- vngs[ which.min( coos$start_pos ) ]
      else if ( mean( as.character( coos$strand ) == "R" ) > 0.6 ) head <- vngs[ which.max( coos$end_pos ) ]
      else { head <- ""
        warning( paste( "About 50% of operon #", i, "are on opposite strands -- don't know what to do!", call.=F ) ) }
      head
    } )
    names( operons ) <- unlist( nms )
    operons <- operons[ names( operons ) != "" ]
    operons <- do.call( rbind, lapply( names( operons ), function( h ) data.frame( head=h, gene=operons[[ h ]] ) ) )
    operons <- subset( operons, head != "" )
  }
  if ( ! is.null( operons ) ) attr( operons, "source" ) <- fetch.predicted.operons
  closeAllConnections()
  operons
}

get.STRING.links <- function( org.id=genome.info$org.id$V1[ 1 ], detailed=T ) {
  ## Original file fetched from http://string-db.org/newstring_download/protein.links.v9.0.txt.gz
  ##   and parsed into individual organism and saved on bragi in ~cmonkey/data/STRING/
  ## line format is something like this (85962 is Hpy Org.id):
  ## 85962.HP0001 85962.HP0617 180
  ## My parser removed the "85962." from each line.
  ## Note old versions of string can be fetched from here: http://string.embl.de/server_versions.html
  ##   e.g. v8.0 from http://string80.embl.de/newstring_download/protein.links.v8.0.txt.gz
  
  ##if ( exists( "get.STRING.links.NEW" ) ) return( invisible( get.STRING.links.NEW( org.id ) ) )
  
  url <- sprintf( "http://baliga.systemsbiology.net/cmonkey/data/STRING/%d_STRING.tsv.gz", org.id )
  fname <- sprintf( "data/%s/string_links.tab.gz", rsat.species )
  if ( ( ! file.exists( fname ) || file.info( fname )$size <= 10 ) ) {
    err <- dlf( fname, url, paste( "Fetching STRING protein links file", url, "\n" ) )
    ##if ( class( err ) == "try-error" || ! file.exists( fname ) || file.info( fname )$size < 10 ) 
    ##  stop( paste( "Whoops, the file was not completely downloaded. Please try to download it yourself from",
    ##              url, "and place it in", fname, "then restart cMonkey.\n" ) )
  }
  
  if ( file.exists( fname ) && file.info( fname )$size > 10 ) {
    cat( "Loading EMBL STRING interaction links from local file", fname, "\n" )
    string.links <- read.delim( gzfile( fname ), sep="\t", head=F )
    colnames( string.links ) <- c( "protein1", "protein2", "combined_score" )
  } else {
    string.links <- get.string.links.NEW( org.id )
  }
  closeAllConnections()
  invisible( string.links )
}

## NEW String links function: use organism ID from uniprot, can download full list from here:
## http://www.uniprot.org/taxonomy/?query=*&compress=yes&format=tab
## Can query it (this is an idea to try in future!), e.g. via:
## http://www.uniprot.org/taxonomy/?query=pylori&format=tab
## In theory we could get the links for a single organism from the STRING API via:
## dlf( fname, paste("http://string-db.org/api/psi-mi-tab/interactionsList?identifiers=",
##                    paste(paste(org.id,rownames(ratios),sep="."),collapse="%0D"),sep="") )
## but apparently this is not allowed (but can do it in increments of 200 ...
## See STRING API info page: http://string.embl.de/help/index.jsp?topic=/org.string-db.docs/api.html
## NOTE for mammalian (at least mouse and human(?)) you do orgid.geneid ONLY if it's e.g.
##      9544.ENSMMUP00000032511  but if it's e.g. Apo3, then leave out the orgid prefix.
##      (need to implement this!)
## TODO: download parsed organism-specific files from bragi instead. See proc_string.pl ...
get.STRING.links.NEW <- function( org.id=genome.info$org.id$V1[ 1 ], all.genes=attr( ratios, "rnames" ),
                             score="score", min.score=02, string.url="http://string-db.org/" ) {
  if ( file.exists( paste( "data/", rsat.species, "/string_links_FALSE_", org.id, ".tab", sep="" ) ) ) {
    ## Has already-parsed (from original BIG BIG STRING flat file) 3-column space-delimited file
    string.links <- read.delim( paste( "data/", rsat.species, "/string_links_FALSE_", org.id, ".tab", sep="" ),
                               head=T, sep=" " )
    string.links$protein1 <- gsub( paste( org.id, ".", sep="" ), "", string.links$protein1 ) 
    string.links$protein2 <- gsub( paste( org.id, ".", sep="" ), "", string.links$protein2 )
    return( string.links )
  }

  ## if ( ! is.na( genome.info$gene.regex ) )
  ##   all.genes <- unique( grep( paste( "^", genome.info$gene.regex, sep="" ), genome.info$feature.names$names,
  ##                             val=T ) )
  ## else
  
  file <- sprintf( "data/%s/string_links_%s.tab", rsat.species, org.id )

  proc.string.df <- function( file ) {
    err <- try( tmp <- unique( read.delim( file, head=F, sep="" ) ) ) ## "" Allows tabs OR spaces
    if ( "try-catch" %in% class( err ) ) return( NULL )
    if ( ! exists( "tmp" ) || nrow( tmp ) <= 0 ) return( NULL )
    tmp2 <- strsplit( as.character( tmp$V15 ), "[:|]", perl=T )
    tmp2a <- sapply( tmp2, function( i ) which( i == score ) )
    tmp2b <- sapply( 1:length( tmp2 ), function( i ) if ( length( tmp2a[[ i ]] ) == 0 ) NA else
                    as.numeric( tmp2[[ i ]][ tmp2a[[ i ]] + 1 ] ) )
    string.links <- data.frame( protein1=gsub( paste( "string:", org.id, ".", sep="" ), "", tmp$V1 ),
                               protein2=gsub( paste( "string:", org.id, ".", sep="" ), "", tmp$V2 ),
                               combined_score=tmp2b )
    string.links <- unique( subset( string.links, ! is.na( combined_score ) ) )
    string.links
  }

  string.links <- NULL; tried <- character()
  if ( file.exists( file ) ) string.links <- proc.string.df( file )
  if ( file.exists( sprintf( "%s.tried", file ) ) ) tried <- readLines( sprintf( "%s.tried", file ) )

  ##tried <- unique( c( as.character( string.links$protein1 ), as.character( string.links$protein2 ) ) )
  tmp2 <- all.genes %in% tried ##sapply( all.genes, function( g ) length( grep( g, tmp, perl=T ) ) ) > 0
  if ( ! file.exists( file ) || ( ! is.null( string.links ) && any( ! tmp2 ) ) ) {
    if ( ! is.null( string.links ) ) all.genes <- all.genes[ ! tmp2 ]
    id.file <- tempfile()
    options( timeout=300 ) ## default is 60secs which may be too short
    for ( i in seq( 1, length( all.genes ), by=100 ) ) {
      ##ids <- paste( org.id, all.genes[ i:min( i + 99, length( all.genes ) ) ], sep="." )
      ids <- all.genes[ i:min( i + 99, length( all.genes ) ) ]
      ids <- c( ids, paste( org.id, all.genes[ i:min( i + 99, length( all.genes ) ) ], sep="." ) )
      cat( i, "of", length( all.genes ), "\n" )
      if ( org.id == 3702 ) { ## Note can resolve gene names (e.g. needed for Ath) via:
        ##ids <- all.genes[ i:min( i + 99, length( all.genes ) ) ]
        url <- paste( string.url, "api/tsv/resolveList?caller_identity=cMonkey&identifiers=",
                     URLencode( paste( ids, collapse="\r" ), reserved=T ), sep="" )
        dlf( id.file, url, mode="wb", quiet=T )
        ids <- unique( as.character( read.delim( id.file )$stringId ) )
        ##cat("HERE1:",length(ids),"\n")
        unlink( id.file )
      }
      url <- paste( string.url, "api/psi-mi-tab/interactionsList?required_score=", min.score,
                   "&caller_identity=cMonkey&network_graph=2&limit=99999&identifiers=",
               URLencode( paste( ids, collapse="\r" ), reserved=T ), "&species=", org.id, sep="" )
      if ( ! file.exists( file ) )
        dlf( file, url, mode="wb", msg="Fetching STRING protein links (piecewise)... this may take a while...",
            quiet=F )
      else dlf( file, url, mode="ab", quiet=F ) ## Append result to end of file
      ##cat("HERE2:",system(sprintf("wc -l %s",file),intern=T),"\n")
    }

    string.links <- proc.string.df( file )
    writeLines( unique( c( all.genes, tried ) ), sprintf( "%s.tried", file ) )
    options( timeout=60 ) ## default is 60secs which may be too short
  }
  invisible( string.links )
}


get.COG.code <- function( org, rows=attr( ratios, "rnames" ) ) {
  ## Can we use the updated mappings from STRING? Not sure...
  ## fname <- "data/COG.mappings.txt.gz"
  ## err <- dlf( fname, 'http://string.embl.de:8080/newstring_download/COG.mappings.v8.3.txt.gz',
  ##            "Fetching COG codes from EMBL STRING..." )
  
  up.rows <- toupper( rows )
  out <- rep( "-", length( rows ) ); names( out ) <- up.rows ## if none, return all "-"'s (named character vector)
  ## Note can see what each COG letter code corresponds to in
  ## ftp://ftp.ncbi.nih.gov/pub/COG/COG/fun.txt
  ## And full organism names from
  ## ftp://ftp.ncbi.nih.gov/pub/COG/COG/org.txt

  fname <- "data/COG_whog.txt"
  err <- dlf( fname, 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog', "Fetching COG codes from NCBI..." )
  lines <- readLines( gzfile( fname ) )
  closeAllConnections()
  hits <- grep( paste( org, '\\:|COG', sep="" ), lines )
  hpy.hits <- grep( paste( org, '\\:', sep="" ), lines[ hits ] )
  if ( length( hpy.hits ) <= 0 ) return( NULL )
  genes <- gsub( paste( '\\s+', org, '\\:\\s+', sep='' ), '', lines[ hits ][ hpy.hits ], perl=T )
  cogs <- lines[ hits ][ hpy.hits - 1 ]
  cog.codes <- sapply( strsplit( cogs, '[\\s+\\[\\]]', perl=T ), "[", 2 )
  cog.codes <- substr( cog.codes, 1, 1 )
  genes <- toupper( genes )

  mc <- get.parallel( length( genes ) ) ##attr( ratios, "nrow" ) )
  tmp <- mc$apply( 1:length( genes ), function( i ) { 
    gn <- strsplit( genes[ i ], " " )[[ 1 ]]
    if ( length( gn ) <= 0 ) next
    gn <- gn[ ! is.na( gn ) ]
    if ( ! all( gn %in% up.rows ) ) gn <- toupper( unlist( get.synonyms( gn, ignore=T ) ) ) 
    if ( sum( gn %in% up.rows ) <= 0 ) return( character() ) 
    gn <- gn[ gn %in% up.rows ]
    out[ up.rows %in% gn ] <- cog.codes[ i ] 
    out
  } ) 
  for ( t in tmp ) if ( length( t ) > 0 ) out[ t != "-" ] <- t[ t != "-" ]
  out[ out == "-" ] <- NA
  names( out ) <- rows 
  closeAllConnections()
  out
}
