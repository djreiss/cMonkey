blast.align <- function( seqs1, seqs2, addl.params="", full.out=F, verbose=F, unlink=T,
                        bl2seq.cmd=sprintf( "%s/bl2seq", progs.dir ),
                        formatdb.cmd=sprintf( "%s/formatdb", progs.dir ),
                        blast.cmd=sprintf( "%s/blastall", progs.dir ) ) {
  file1 <- my.tempfile( "blast." )
  file2 <- my.tempfile( "blast." )
  tmp.log.file <- my.tempfile( "formatdb.log" )
  if ( length( seqs1 ) == 1 && length( seqs2 ) == 1 ) {
    cat( ">seq1\n", seqs1, "\n", sep="", file=file1 )
    cat( ">seq2\n", seqs2, "\n", sep="", file=file2 )
    ##cmd <- paste( bl2seq.cmd, "-i", file1, "-j", file2, "-p blastn", addl.params )
    cmd <- sprintf( "%s -i %s -j %s -p blastn %s", bl2seq.cmd, file1, file2, addl.params )
    if ( ! full.out ) cmd <- paste( cmd, "-D 1" )
    if ( verbose ) cat( cmd, "\n" )
    output <- system( cmd, intern=TRUE, ignore.stderr=TRUE )
  } else {

    if ( length( seqs1 ) == 1 ) cat( ">seq1\n", seqs1, "\n", sep="", file=file1 )
    else cat( paste( ">", names( seqs1 ), "\n", seqs1, sep="" ), file=file1, sep="\n" )
    if ( length( seqs2 ) == 1 ) cat( ">seq2\n", seqs2, "\n", sep="", file=file2 )
    else cat( paste( ">", names( seqs2 ), "\n", seqs2, sep="" ), file=file2, sep="\n" )

    if ( length( seqs1 ) >= length( seqs2 ) ) {
      ##cmd <- paste( formatdb.cmd, "-l", tmp.log.file, "-i", file1, "-p F -o T" )
      cmd <- sprintf( "%s -l %s -i %s -p F -o T", formatdb.cmd, tmp.log.file, file1 )
      if ( verbose ) cat( cmd, "\n" )
      system( cmd )
      ##cmd <- paste( blast.cmd, "-p blastn -d", file1, "-i", file2, addl.params )
      cmd <- sprintf( "%s -p blastn -d %s -i %s %s", blast.cmd, file1, file2, addl.params ) 
    } else if ( length( seqs2 ) > length( seqs1 ) ) {
      ##cmd <- paste( formatdb.cmd, "-l", tmp.log.file, "-i", file2, "-p F -o T" )
      cmd <- sprintf( "%s -l %s -i %s -p F -o T", formatdb.cmd, tmp.log.file, file2 )
      if ( verbose ) cat( cmd, "\n" )
      system( cmd )
      ##cmd <- paste( blast.cmd, "-p blastn -d", file2, "-i", file1, addl.params )
      cmd <- sprintf( "%s -p blastn -d %s -i %s %s", blast.cmd, file2, file1, addl.params ) 
    }
    if ( ! full.out ) cmd <- paste( cmd, "-m 8" )
    if ( verbose ) cat( cmd, "\n" )
    output <- system( cmd, intern=TRUE, ignore.stderr=TRUE )
  }
  ## Fields are: Query id, Subject id, % identity, alignment length, mismatches, gap openings,
  ## q. start, q. end, s. start, s. end, e-value, bit score"
  ##if ( ! full.out ) output <- strsplit( output, "\t" )
  if ( unlink ) try( unlink( c( paste( file1, "*", sep="" ), paste( file2, "*", sep="" ), tmp.log.file ) ) )
  return( output )
}

## Assumes blast.align was run with full.out=F                   
parse.blast.out <- function( blast.out ) {
  if ( substr( blast.out[ 1 ], 1, 12 ) == "# BLASTN 2.2" ) blast.out <- blast.out[ -(1:3) ]
  out <- t( sapply( strsplit( blast.out, "\t" ), cbind ) )
  out <- data.frame( `Query id`=out[ ,1 ], `Subject id`=out[ ,2 ], `% identity`=as.numeric( out[ ,3 ] ),
                    `alignment length`=as.integer( out[ ,4 ] ), mismatches=as.integer( out[ ,5 ] ),
                    `gap openings`=as.integer( out[ ,6 ] ), `q. start`=as.integer( out[ ,7 ] ),
                    `q. end`=as.integer( out[ ,8 ] ), `s. start`=as.integer( out[ ,9 ] ),
                    `s. end`=as.integer( out[ ,10 ] ), `e-value`=as.numeric( out[ ,11 ] ),
                    `bit score`=as.numeric( out[ ,12 ] ) )
  out
}

blast.match.seqs <- function( seqs, match=NULL, e.cutoff=1 ) {
  if ( is.null( match ) ) match <- seqs
  out <- parse.blast.out( blast.align( seqs, match, paste( "-e", e.cutoff ) ) )
  ##good <- which( as.numeric( out[ ,"e.value" ] ) <= e.cutoff )
  ##out <- out[ good, c( "Query.id", "Subject.id" ) ]
  out <- subset( out, as.character( Query.id ) != as.character( Subject.id ) )
  out <- subset( out, as.character( Query.id ) %in% names( seqs ) )
  out <- subset( out, as.character( Subject.id ) %in% names( match ) )
  return( out[ order( out$e.value ), ] )
}

all.dna.seqs <- function( l, lett=c( "G", "A", "T", "C" ), as.matrix=F ) {
  n.lett <- length( lett )
  out <- sapply( 1:l, function( ll ) rep( as.vector( sapply( lett, function( i ) rep( i, n.lett^( ll - 1 ) ) ) ),
                                         n.lett^( l - ll ) ) )
  if ( as.matrix ) return( out )
  apply( out, 1, paste, collapse="" )
}

weeder.one.cluster <- function( k, seq.type="upstream weeder", n.motifs=4,
                               verbose=F, unlink=T, weeder.size="medium", ... ) {
  ntides <- c( "T", "G", "A", "C" )
  for ( w in c( 6, 8 ) ) {
    if ( ! file.exists( sprintf( "%s/FreqFiles", progs.dir ) ) ) dir.create( sprintf( "%s/FreqFiles", progs.dir ) )
    if ( ! file.exists( paste( sprintf( "%s/FreqFiles/", progs.dir ),
                              toupper( organism ), ".", w, ".", "freq", sep="" ) ) ) {
      seqs <- unique( genome.info$all.upstream.seqs[[ seq.type ]] )
      all.substrings <- as.vector( sapply( 1:( max( nchar( seqs ) ) - w + 1 ),
                                          function( i ) substr( seqs, i, i + w - 1 ) ) )  
      all.substrings <- all.substrings[ ! is.na( all.substrings ) & all.substrings != "" &
                                       nchar( all.substrings ) == w ]
      all.substrings <- all.substrings[ ! grepl( "[^GATC]", all.substrings ) ]
      hist.substrings <- table( as.factor( all.substrings ) )
      
      ##all.combos <- unique( combn( rep( ntides, w ), w, FUN=paste, sep="", collapse="" ) )
      all.combos <- all.dna.seqs( w, ntides )
      ##ss.counts <- mkBgFile( genome.info$all.upstream.seqs[[ seq.type ]], order=w, verbose=F )
      all.combos <- all.combos[ ! all.combos %in% names( hist.substrings ) ]
      tmp <- rep( 0, length( all.combos ) ); names( tmp ) <- all.combos
      hist.substrings <- c( hist.substrings, tmp ) + 1 ## Add pseudocount
      hist.substrings <- hist.substrings[ sort( names( hist.substrings ) ) ] ## must be in alphabetical order!
      write.table( hist.substrings, quote=F, sep=" ", col.names=F,
             file=paste( sprintf( "%s/FreqFiles/", progs.dir ), toupper( organism ), ".", w, ".", "freq", sep="" ) )
    }
  }

  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, ... ) ##mask=mask, blast=blast, 
  min.seqs <- cluster.rows.allowed[ 1 ]; max.seqs <- cluster.rows.allowed[ 2 ]
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )
  
  fst.file <- my.tempfile( sprintf( "weeder_%d_%d_", k, iter ), suf=".fst" )
  cat( paste( ">", names( seqs ), "\n", seqs, sep="" ), file=fst.file, sep="\n" )
  file.remove( paste( fst.file, c( "fst", "html", "mix", "wee" ), sep="." ) )
  cwd <- setwd( progs.dir ); on.exit( setwd( cwd ) ) ##"./progs" )
  ## Options: S: search both strands. Other options: A: assume motifs in ALL seqs; M: motif can occur > once per seq;
  ## Tn: save n highest scoring motifs (default 10)
  ## weeder.size can be "small" (w = 6 or 8); "medium" (6, 8 or 10), "large" (default; 6, 8, 10, 12; takes too long!)
  ## Documentation: http://159.149.109.9/modtools/downloads/weedermanual.pdf
  cmd <- sprintf( weeder.cmd, fst.file, toupper( organism ), weeder.size, n.motifs*5 )
  if ( verbose ) print( cmd )
  out <- system( cmd, intern=T, ignore.stderr=!verbose )
  setwd( cwd )
  if ( ! file.exists( paste( fst.file, "mix", sep="." ) ) && ! file.exists( paste( fst.file, "wee", sep="." ) ) )
    return( list( k=k ) ) 
  out <- c( readLines( paste( fst.file, "mix", sep="." ) ), readLines( paste( fst.file, "wee", sep="." ) ) )
  if ( unlink ) file.remove( c( fst.file, paste( fst.file, c( "fst", "html", "mix", "wee" ), sep="." ) ) )

  mot.scores <- as.data.frame( do.call( rbind, strsplit( grep( "^\\d+\\) ", out, perl=T, val=T )," " ) )[ ,2:4 ] )
  mot.scores[ ,2 ] <- as.numeric( as.character( mot.scores[ ,2 ] ) )
  mot.scores[ ,3 ] <- as.numeric( as.character( mot.scores[ ,3 ] ) )
  mot.scores[ ,3 ][ is.na( mot.scores[ ,3 ] ) ] <- 0
  
  is.highest.ranking.motif <- grep( "Interesting motifs (highest-ranking) seem to be", out, fixed=T )
  is.not.highest.ranking.motif <- grep( "Interesting motifs (not highest-ranking) can also be", out, fixed=T )
  starts <- grep( "Best occurrences", out, fixed=T )
  weeder.out <- list()
  for ( i in 1:length( starts ) ) {
    start <- starts[ i ]
    weeder.out[[ i ]] <- list()
    weeder.out[[ i ]]$motifs.redund <- strsplit( out[ start - 2 ], "\\s\\-\\s", perl=T )[[ 1 ]]
    weeder.out[[ i ]]$motifs <- out[ c( start - 5, start - 6 ) ]
    weeder.out[[ i ]]$is.highest.ranking <- start > is.highest.ranking.motif && start < is.not.highest.ranking.motif
    weeder.out[[ i ]]$score <- unique( subset( mot.scores, V1 %in% weeder.out[[ i ]]$motifs )$V2 )
    end <- start - 1 + min( which( out[ ( start + 1 ):length( out ) ] == "" ) )
    lines <- strsplit( out[ (start+1):end ], "\\s+", perl=T )
    lines <- do.call( rbind, lines )
    colnames( lines )[ 2:ncol( lines ) ] <- lines[ 1, 1:( ncol( lines ) - 1 ) ]
    lines <- as.data.frame( lines[ -1, -1 ] )
    lines$match <- gsub( "(", "", gsub( ")", "", lines$match, fixed=T ), fixed=T )
    weeder.out[[ i ]]$matches <- lines
    i2 <- end + min( grep( "Frequency Matrix", out[ end:length( out ) ] ) ) + 1
    lines <- strsplit( out[ i2:( i2 - 1 + min( which( out[ ( i2 + 1 ):length( out ) ] == "" ) ) ) ], "\t+", perl=T )
    lines <- do.call( rbind, lines )[ ,-1 ]
    counts1 <- do.call( rbind, strsplit( lines[ ,1 ], "\\s+", perl=T ) )[ ,-1 ]
    colnames( counts1 ) <- counts1[ 1, ]; counts1 <- counts1[ -1, ]
    counts1 <- apply( counts1, 2, as.integer )
    counts2 <- do.call( rbind, strsplit( lines[ ,2 ], "\\s+", perl=T ) )[ ,-1 ]
    colnames( counts2 ) <- counts2[ 1, ]; counts2 <- counts2[ -1, ]
    counts2 <- apply( counts2, 2, as.integer )
    weeder.out[[ i ]]$counts.all <- counts1
    weeder.out[[ i ]]$counts.best <- counts2
  }

  if ( length( weeder.out ) <= 0 ) {
    attr( weeder.out, "weeder.out" ) <- out
    return( list( k=k, weeder.out=weeder.out ) )
  }
  
  ##TODO: include motif length (longer better) in weeder.out motif sorting, then use n.motifs best motifs (can use 2 for iters 100-1000, then 4 for 1001-2000)
  ## Order motifs by increasing length
  weeder.out <- weeder.out[ order( sapply( weeder.out, function( i ) nchar( i$motifs[ 1 ] ) ) ) ]  
  
  ## Now, see if any smaller motifs are in the "motifs.redund" of bigger motifs
  n.mot <- length( weeder.out )
  if ( n.mot > 1 ) {
    m.redund <- matrix( 0, nrow=n.mot, ncol=n.mot )
    for ( i in 1:( n.mot - 1 ) ) for ( j in (i + 1):n.mot )
      m.redund[ i, j ] <- sum( weeder.out[[ i ]]$motifs %in% weeder.out[[ j ]]$motifs.redund )
    weeder.out <- weeder.out[ order( apply( m.redund, 1, sum, na.rm=T ), decreasing=F ) ]
    m.redund <- m.redund * 0 
    n.mot <- length( weeder.out )
    for ( i in 1:( n.mot - 1 ) ) for ( j in (i + 1):n.mot )
      m.redund[ i, j ] <- sum( weeder.out[[ i ]]$motifs %in% weeder.out[[ j ]]$motifs.redund )
    ## Reorder to have highest scoring, longer, non-redundant motifs first
    sum.redund <- apply( m.redund, 1, sum, na.rm=T )
    ##sum.redund <- max( sum.redund, na.rm=T ) - sum.redund
    m.length <- sapply( weeder.out, function( i ) nchar( i$motifs[ 1 ] ) )
    weeder.out <- weeder.out[ order( sapply( weeder.out, "[[", "is.highest.ranking" ), m.length,
                                    sum.redund, decreasing=T ) ]
    m.redund <- m.redund * 0 
    n.mot <- length( weeder.out )
    for ( i in 1:( n.mot - 1 ) ) for ( j in (i + 1):n.mot )
      m.redund[ i, j ] <- sum( weeder.out[[ i ]]$motifs %in% weeder.out[[ j ]]$motifs.redund )
    attr( weeder.out, "is.redund" ) <- m.redund
    weeder.out <- weeder.out[ 1:n.motifs ]
    weeder.out <- weeder.out[ ! sapply( weeder.out, is.null ) ]
  }
  attr( weeder.out, "weeder.out" ) <- out

  m.in <- character()
  for ( i in 1:length( weeder.out ) )
    m.in <- c( m.in, pssm.motif.lines( weeder.out[[ i ]]$counts.all, id=sprintf( "weeder_%d", i ), header=(i==1) ) )
  all.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  mast.out <- runMast( m.in, mast.cmd[ seq.type ], names( all.seqs ), all.seqs,
                      bg.list=genome.info$bg.list[[ seq.type ]], unlink=T, verbose=verbose )

  pv.ev <- get.pv.ev.single( mast.out, rows )

  meme.out <- list() ## Make a "meme.out" structure similar to getMemeMotifInfo() for MEME output
  for ( ii in 1:length( weeder.out ) ) {
    wo <- weeder.out[[ ii ]]
    pssm <- wo$counts.all
    pssm <- pssm + max( pssm, na.rm=T ) / 100 ## Pseudocount?
    for ( i in 1:nrow( pssm ) ) pssm[ i, ] <- pssm[ i, ] / sum( pssm[ i, ], na.rm=T )
    posns <- data.frame( gene=names( seqs )[ wo$matches$Seq ], strand=wo$matches$St,
                        start=as.integer( wo$matches$pos ),
                        p.value=(100-as.numeric( wo$matches$match )+0.001)/100,
                        site=gsub( '[\\[\\]]', '', wo$matches$oligo, perl=T ) )
    meme.out[[ ii ]] <- list( width=nrow( wo$counts.all ), sites=nrow( wo$matches ),
                            llr=wo$score, e.value=wo$score, pssm=pssm, posns=posns )
  }
  attr( meme.out, "is.pal" ) <- FALSE
  
  invisible( list( k=k, weeder.out=weeder.out, meme.out=meme.out, pv.ev=pv.ev ) )
  ##invisible( weeder.out )
}

pssm.motif.lines <- function( pssm, id, e.value=1, header=T, seq.type="upstream weeder" ) {
  meme.let <- c( "A", "C", "G", "T" )
  if ( missing( id ) ) id <- paste( pssm, collapse="" )
  lines <- character()
  if ( header ) {
    lines <- ##c( "MEME version 3.0", "",
               "ALPHABET= ACGT"##, "", "strands: + -", "",
               ##"Background letter frequencies (from dataset with add-one prior applied):" )
               ##)
    lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                        sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )
  }
    
  if ( is.null( colnames( pssm ) ) ) colnames( pssm ) <- col.let
  pssm <- pssm[ ,meme.let ] + max( pssm, na.rm=T ) / 100 ## Pseudocount?
  for ( i in 1:nrow( pssm ) ) pssm[ i, ] <- pssm[ i, ] / sum( pssm[ i, ], na.rm=T )
  idd <- gsub( "[_/]", ".", id )
  ##lines <- c( lines, "", sprintf( "MOTIF %s", idd ), sprintf( "BL   MOTIF %s width=0 seqs=0", idd ) )
  ##lines <- c( lines, sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
  ##                           20, e.value ) )
  lines <- c( lines, sprintf( "log-odds matrix: alength= 4 w= %d", nrow( pssm ) ) )
  for ( j in 1:nrow( pssm ) ) lines <- c( lines, paste( sprintf( "%5.3f", log2( pssm[ j, ] ) ),
                                                       collapse=" ", sep=" " ) )
  lines
}

## Note, with seq.type == "... spacer", need to set filter=F because SPACER doesn't like N's in sequences !
## seq.type can contain "spacer" (gapped) or "prism" (ungapped) or neither (assumes "spacer")
spacer.one.cluster <- function( k, seq.type="upstream spacer", hits.to.all=F, verbose=F, unlink=T, score.cutoff=-3,
                               ... ) {
  min.seqs <- cluster.rows.allowed[ 1 ]
  max.seqs <- cluster.rows.allowed[ 2 ]
  bg.fname <- paste( progs.dir, "/SPACER/data/genomes/", organism, ".upstream.raw", sep="" )
  if ( ! file.exists( bg.fname ) ) dir.create( sprintf( "%s/SPACER/data/genomes", progs.dir ) )
  writeLines( genome.info$all.upstream.seqs[[ seq.type ]], con=bg.fname )

  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  filter <- FALSE
  if ( "filter" %in% names( list( ... ) ) && list( ... )$filter && ! grepl( "spacer", seq.type )[ 1 ] ) filter <- TRUE
  seqs <- get.sequences( rows, seq.type=seq.type, filter=filter, ... ) ##mask=mask, blast=blast, 
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )
  
  bg.file <- paste( organism, ".upstream.raw", sep="" )
  tempfile <- my.tempfile( sprintf( "spacer_%d_%d_", k, iter ), suf=".fasta" )
  temp.out <- my.tempfile( sprintf( "spacer_out_%d_%d_", k, iter ), suf=".txt" )
  cat( paste( ">", names( seqs ), "\n", seqs, sep="" ), file=tempfile, sep="\n" )
  cwd <- setwd( sprintf( "%s/SPACER", progs.dir ) ); on.exit( setwd( cwd ) )
  ##cat( "java -mx500M -jar SPACER.jar -b", bg.file, "-o", temp.out, tempfile, "\n" ); stop()
  if ( ! grepl( "prism", seq.type )[ 1 ] ) { ##is.null( opt ) ) {
    ##cmd <- paste( "java -Xmx1000M -Xshare:off -jar SPACER.jar -b", bg.file, "-o", temp.out, tempfile )
    ##cmd <- sprintf( "java -Xmx1000M -Xshare:off -jar SPACER.jar -b %s -o %s %s", bg.file, temp.out, tempfile )
    cmd <- sprintf( spacer.cmd[ 1 ], bg.file, temp.out, tempfile )
  } else { ##if ( grepl( "prism", seq.type )[ 1 ] ) { ##opt == "PRISM" ) {
    ##cmd <- paste( "java -Xmx1000M -Xshare:off -jar PRISM.jar -b", bg.file, "-o", temp.out, tempfile )
    ##cmd <- sprintf( "java -Xmx1000M -Xshare:off -jar PRISM.jar -b %s -o %s %s", bg.file, temp.out, tempfile )
    cmd <- sprintf( gsub( "SPACER.jar", "PRISM.jar", spacer.cmd[ 1 ] ), bg.file, temp.out, tempfile )
  }
  if ( verbose ) print( cmd )
  out <- system( cmd, intern=T, ignore.stderr=!verbose )
  if ( ! file.exists( temp.out ) ) return( list( k=k ) )
  spacer.out <- readLines( temp.out )
  if ( unlink ) unlink( temp.out )
  out <- strsplit( spacer.out[ spacer.out != "" ], "[\\t\\,]", perl=T )

  all.fasta <- my.tempfile( sprintf( "spacer_all_%d_%d_", k, iter ), suf=".fasta" )
  bg.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  if ( hits.to.all && ! is.na( bg.seqs ) && length( out ) > 1 )
    cat( paste( ">", names( bg.seqs ), "\n", bg.seqs, sep="" ), file=all.fasta, sep="\n" )
  
  out2 <- out3 <- list()
  for ( i in 1:length( out ) ) {
    motif <- out[[ i ]][ 2 ]
    motif.score <- as.numeric( out[[ i ]][ 1 ] )
    if ( motif.score < score.cutoff ) break
    prob <- 2^( -motif.score )
    temp.out <- my.tempfile( sprintf( "spacer_out_%d_%d_", k, iter ), suf=".txt" )
    ##if ( is.null( opt ) ) {
    if ( ! grepl( "prism", seq.type )[ 1 ] ) { ##is.null( opt ) ) {
      ##cmd <- paste( "java -Xmx1000M -Xshare:off -jar SPACER.jar -l", motif, "-o", temp.out, tempfile )
      ##cmd <- sprintf( "java -Xmx1000M -Xshare:off -jar SPACER.jar -l %s -o %s %s", motif, temp.out, tempfile )
      cmd <- sprintf( spacer.cmd[ 2 ], motif, temp.out, tempfile )
    } else { ##if ( opt == "PRISM" ) {
      ##cmd <- paste( "java -Xmx1000M -Xshare:off -jar PRISM.jar -l", motif, "-o", temp.out, tempfile )
      ##cmd <- sprintf( "java -Xmx1000M -Xshare:off -jar PRISM.jar -l %s -o %s %s", motif, temp.out, tempfile )
      cmd <- sprintf( gsub( "SPACER.jar", "PRISM.jar", spacer.cmd[ 2 ] ), motif, temp.out, tempfile )
    }
    if ( verbose ) print( cmd )
    tmp <- system( cmd, intern=T, ignore.stderr=!verbose )
    tmp <- readLines( temp.out )
    spacer.out <- c( spacer.out, tmp )
    if ( unlink && file.exists( temp.out ) ) unlink( temp.out )

    tmp <- do.call( rbind, strsplit( tmp[ tmp != "" ], "," ) )
    out2[[ length( out2 ) + 1 ]] <- tmp
    if ( hits.to.all ) {
      tmp.out2 <- my.tempfile( sprintf( "spacer_out_%d_%d_", k, iter ) )
      if ( ! file.exists( all.fasta ) ) { ##"./zzz.all.fasta" ) ) {
        tmp <- bg.seqs[ ! bg.seqs %in% seqs ]
        cat( paste( ">", names( bg.seqs ), "\n", bg.seqs, sep="" ), file=all.fasta, sep="\n" )
      }
      if ( ! grepl( "prism", seq.type )[ 1 ] ) { ##is.null( opt ) ) {
      ##if ( is.null( opt ) ) {
        ##cmd <- paste( "java -mx500M -jar SPACER.jar -l", motif, "-o", tmp.out2, all.fasta )
        ##cmd <- sprintf( "java -mx1000M -Xshare:off -jar SPACER.jar -l %s -o %s %s", motif, tmp.out2, all.fasta )
        cmd <- sprintf( spacer.cmd[ 2 ], motif, tmp.out2, all.fasta )
      } else { ##if ( opt == "PRISM" ) {
        ##cmd <- paste( "java -mx500M -jar PRISM.jar -l", motif, "-o", tmp.out2, all.fasta )
        ##cmd <- sprintf( "java -mx1000M -Xshare:off -jar PRISM.jar -l %s -o %s %s", motif, tmp.out2, all.fasta )
        cmd <- sprintf( gsub( "SPACER.jar", "PRISM.jar", spacer.cmd[ 2 ] ), motif, temp.out2, all.fasta )
      }
      if ( verbose ) print( cmd )
      tmp <- system( cmd, intern=T, ignore.stderr=!verbose )
      tmp <- readLines( tmp.out2 )
      tmp <- do.call( rbind, strsplit( tmp[ tmp != "" ], "," ) )
      tmp[ ,1 ] <- names( bg.seqs )[ as.integer( as.character( tmp[ ,1 ] ) ) ]
      out3[[ length( out3 ) + 1 ]] <- tmp
      if ( unlink && file.exists( tmp.out2 ) ) unlink( tmp.out2 )
    }
  }

  if ( unlink && file.exists( tempfile ) ) file.remove( tempfile )
  if ( unlink && file.exists( all.fasta ) ) file.remove( all.fasta )
  
  setwd( cwd )
  out <- list( motifs=do.call( rbind, out ), hits=out2 )
  attr( out, "spacer.out" ) <- spacer.out
  if ( hits.to.all ) out$all.hits=out3
  ## Info on 3rd column of out$motifs:
  ## [r] scores higher when it is counted on only the top strand in the given fasta file AND the background sequences
  ## [R] scorer higher when it is counted on both strands (palindromes are counted only once).
  ## Also note for out$hits, position is a negative number counting from the right end of the geneID (starting with -)
  ## Also note that "score" is like a -log2(e-value).
  meme.let <- c( "A", "C", "G", "T" )
  out$motifs <- data.frame( score=as.numeric( out$motifs[ ,1 ] ), motif=out$motifs[ ,2 ], flag=out$motifs[ ,3 ] )

  out$motifs <- subset( out$motifs, score >= score.cutoff ) ## remove motifs with "e-value" > 100
  if ( length( out$hits ) <= 0 || nrow( out$motifs ) <= 0 ) return( list( k=k, spacer.out=out ) )
  
  for ( i in 1:length( out$hits ) ) {
    hits <- data.frame( seq=as.integer( out$hits[[ i ]][ ,1 ] ), posn=as.integer( out$hits[[ i ]][ ,2 ] ),
                       strand=out$hits[[ i ]][ ,3 ], site=out$hits[[ i ]][ ,4 ] )
    sites <- toupper( do.call( rbind, strsplit( as.character( hits$site ), "" ) ) )
    if ( any( ! sites %in% meme.let ) )
      sites[ ! sites %in% meme.let ] <- sample( meme.let, sum( ! sites %in% meme.let ) )
    ## Make a pssm from the enumerated sites:
    pssm <- matrix( 0, nrow=ncol( sites ), ncol=4 )
    rownames( pssm ) <- as.character( 1:nrow( pssm ) ); colnames( pssm ) <- meme.let
    for ( j in 1:nrow( sites ) ) pssm[ cbind( 1:nrow( pssm ), sites[ j, ] ) ] <-
      pssm[ cbind( 1:nrow( pssm ), sites[ j, ] ) ] + 1
    out$hits[[ i ]] <- list( sites=hits, pssm=pssm )
  }
  
  ## Now do same thing as for weeder: scan pssms on sequences to determine the p-values of their matches.
  ## TODO: should use SPACER itself to get gene scores since pssms don't capture variable spacing.
  m.in <- character()
  for ( i in 1:length( out$hits ) ) {
    if ( out$motifs$score[ i ] < score.cutoff ) next ## Dont include bad motifs
    m.in <- c( m.in, pssm.motif.lines( out$hits[[ i ]]$pssm, id=sprintf( "spacer_%d", i ), header=(i==1) ) )
  }
  all.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  mast.out <- runMast( m.in, mast.cmd[ seq.type ], names( all.seqs ), all.seqs,
                      bg.list=genome.info$bg.list[[ seq.type ]], unlink=T, verbose=verbose )
  
  pv.ev <- get.pv.ev.single( mast.out, rows )

  meme.out <- list() ## Make a "meme.out" structure similar to getMemeMotifInfo() for MEME output
  for ( ii in 1:length( out$hits ) ) {
    wo <- out$hits[[ ii ]]
    pssm <- wo$pssm
    pssm <- pssm + max( pssm, na.rm=T ) / 100 ## Pseudocount?
    for ( i in 1:nrow( pssm ) ) pssm[ i, ] <- pssm[ i, ] / sum( pssm[ i, ], na.rm=T )
    posns <- data.frame( gene=names( seqs )[ wo$sites$seq ], strand=wo$sites$strand,
                        start=as.integer( wo$sites$posn ),
                        p.value=NA, ##(100-as.numeric( wo$matches$match )+0.001)/100,
                        site=toupper( wo$sites$site ) )
    meme.out[[ ii ]] <- list( width=nrow( wo$pssm ), sites=nrow( wo$sites ),
                             llr=out$motifs$score[ ii ], e.value=2^(-out$motifs$score[ ii ] ),
                             pssm=pssm, posns=posns )
  }
  attr( meme.out, "is.pal" ) <- FALSE
  invisible( list( k=k, spacer.out=out, meme.out=meme.out, pv.ev=pv.ev ) )
  ##invisible( out )
}

## Note: other motif finding algo's in R are rGADEM and MotIV -- see
## http://manuals.bioinformatics.ucr.edu/home/ht-seq#TOC-Additional-Motif-Discovery-Packages
cosmo.one.cluster <- function( k, seq.type="upstream cosmo", n.motifs=2, minW=6, maxW=24, model="ZOOPS", 
                               verbose=F, unlink=T, ... ) { ##plot=F, 
  ## TODO: get minW, maxW, n.motifs, ZOOPS from an initialized param
  ## TODO: mask out sequences w/ first motif and run again if n.motifs > 1  DONE BUT cosmo doesn't seem to ignore 'N's
  require( cosmo )
  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, ... ) ##mask=mask, blast=blast, 
  min.seqs <- cluster.rows.allowed[ 1 ]; max.seqs <- cluster.rows.allowed[ 2 ]
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )

  bgs <- genome.info$all.upstream.seqs[[ seq.type ]]
  ss <- list(); for ( s in names( seqs ) ) ss[[ s ]] <- list( desc=s, seq=seqs[ s ] )
  bgsl <- list(); for ( s in names( bgs ) ) bgsl[[ s ]] <- list( desc=s, seq=bgs[ s ] )
  bgm <- bgModel( bgsl ) ##, order=0 ) ## Note can store this for all calls so dont need to recompute

  res <- list()
  for ( i in 1:n.motifs ) {
    res[[ i ]] <- cosmo( seqs=ss, minW=minW, maxW=maxW, models=model, transMat=bgm$transMat, minSites=length(seqs)/2,
                        maxSites=2*length(seqs), wCrit="eval", intCrit="eval", silent=!verbose, ... )
    pssm <- t( attr( attr( res[[ i ]], "pwm" ), "pwm" )[ col.let, ] )
    ##if ( plot ) viewPssm( pssm )
    if ( i < n.motifs ) { ## cosmo doesnt seem to ignore 'N's!!!
      sites <- attr( res[[ i ]], "motifs" )
      gene <- attr( sites, "seq" ); ## strand=ifelse( attr( sites, "orient" ) < 0, '-', '+' ),
      start <- attr( sites, "pos" )
      
      for ( i in 1:length( start ) )
        ss[[ gene[ i ] ]]$seq <- substr( ss[[ gene[ i ] ]]$seq, start[ i ], start[ i ] + nrow( pssm ) - 1 ) <-
          paste( rep( 'N', nrow( pssm ) ), collapse='' )
    }
  }

  ## Now do same thing as for weeder: scan pssm on sequences to determine the p-values of its matches.
  m.in <- character()
  for ( i in 1:length( res ) ) {
    pssm <- t( attr( attr( res[[ i ]], "pwm" ), "pwm" )[ col.let, ] )
    m.in <- c( m.in, pssm.motif.lines( pssm, id=sprintf( "cosmo_%d", i ), header=(i==1) ) )
  }
  all.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  mast.out <- runMast( m.in, names( all.seqs ), all.seqs, bg.list=genome.info$bg.list[[ seq.type ]],
                      unlink=T, verbose=verbose )
  
  pv.ev <- get.pv.ev.single( mast.out, rows )

  meme.out <- list() ## Make a "meme.out" structure similar to getMemeMotifInfo() for MEME output
  for ( ii in 1:length( res ) ) {
    pssm <- t( attr( attr( res[[ ii ]], "pwm" ), "pwm" )[ col.let, ] )
    pssm <- pssm + max( pssm, na.rm=T ) / 100 ## Pseudocount?
    for ( i in 1:nrow( pssm ) ) pssm[ i, ] <- pssm[ i, ] / sum( pssm[ i, ], na.rm=T )
    sites <- attr( res[[ ii ]], "motifs" )
    posns <- data.frame( gene=attr( sites, "seq" ), strand=ifelse( attr( sites, "orient" ) < 0, '-', '+' ),
                        start=attr( sites, "pos" ), p.value=1-attr( sites, "prob" )+1e-5,
                        site=attr( sites, "motif" ) )
    meme.out[[ ii ]] <- list( width=nrow( pssm ), sites=nrow( sites ),
                             llr=attr(res,"sel")["Width","critVal"], e.value=attr(res,"sel")["Width","critVal"],
                             pssm=pssm, posns=posns )
  }
  attr( meme.out, "is.pal" ) <- FALSE
  invisible( list( k=k, cosmo.out=res, meme.out=meme.out, pv.ev=pv.ev ) )
}

gibbs.one.cluster <- function( k, seq.type=names(mot.weights)[1], n.motifs=2, width=8,
                              verbose=F, unlink=T, ... ) {
  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, ... ) ##mask=mask, blast=blast, 
  min.seqs <- cluster.rows.allowed[ 1 ]; max.seqs <- cluster.rows.allowed[ 2 ]
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )

  bgs <- genome.info$all.upstream.seqs[[ seq.type ]]
  ##ss <- list(); for ( s in names( seqs ) ) ss[[ s ]] <- list( desc=s, seq=seqs[ s ] )
  ##bgsl <- list(); for ( s in names( bgs ) ) bgsl[[ s ]] <- list( desc=s, seq=bgs[ s ] )
  ##bgm <- bgModel( bgsl ) ##, order=0 ) ## Note can store this for all calls so dont need to recompute
  ##res <- cosmo( seqs=ss, minW=12, maxW=20, models=model, transMat=bgm$transMat, minSites=length(seqs)/2,
  ##             maxSites=2*length(seqs), wCrit="eval", intCrit="eval", ... )
  ##if ( plot ) plot( res ) ##viewPssm( t( attr( attr( res, "pwm" ), "pwm" )[ col.let, ] ) )
  ##res
  fname <- my.tempfile( "gibbs_" )
  cat( paste( ">", names( seqs ), "\n", seqs, sep="" ), file=fname, sep="\n" )
  bgtab <- data.frame( names( genome.info$bg.list[[ seq.type ]] ),
                      sapply( genome.info$bg.list[[ seq.type ]], "[", 1 ) )
  rownames( bgtab ) <- NULL; colnames( bgtab ) <- c( "V1", "V2" )
  bgtab <- bgtab[ -1, ]
  bgtab$V1 <- tolower( bgtab$V1 )
  bgfname <- my.tempfile( "gibbs_bg_" )
  write.table( bgtab, bgfname, sep="\t", quote=F, col.names=F, row.names=F )
  outfname <- my.tempfile( "gibbs_out_" )

  ## Run it in 'motif sampler' mode: Gibbs file lengths expect {flags}
  ## Note 'lengths' is possible motif widths separated by a comma, or see -d flag ???
  ## Not sure how to set the background model for this algo
  cmd <- sprintf( "%s/gibbs/Gibbs.i686-apple-darwin %s %d %d -n", progs.dir, fname, width,
                 floor( sqrt( length( seqs ) ) ) )
  if ( verbose ) print( cmd )
  out <- system( cmd, intern=T, ignore.stderr=!verbose )  
  out <- readLines( outfname )
  if ( unlink ) unlink( c( fname, bgfname, outfname ) )
  out
}

gadem.one.cluster <- function( k, seq.type=names(mot.weights)[1], n.motifs=2,
                               verbose=F, unlink=T, ... ) {
  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, ... ) ##mask=mask, blast=blast, 
  min.seqs <- cluster.rows.allowed[ 1 ]; max.seqs <- cluster.rows.allowed[ 2 ]
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  cat( k, "\t", Sys.getpid(), date(), "\t\t", seq.type, "\tSEQUENCES:", length( seqs ), "\n" )

  bgs <- genome.info$all.upstream.seqs[[ seq.type ]]
  ##ss <- list(); for ( s in names( seqs ) ) ss[[ s ]] <- list( desc=s, seq=seqs[ s ] )
  ##bgsl <- list(); for ( s in names( bgs ) ) bgsl[[ s ]] <- list( desc=s, seq=bgs[ s ] )
  ##bgm <- bgModel( bgsl ) ##, order=0 ) ## Note can store this for all calls so dont need to recompute
  ##res <- cosmo( seqs=ss, minW=12, maxW=20, models=model, transMat=bgm$transMat, minSites=length(seqs)/2,
  ##             maxSites=2*length(seqs), wCrit="eval", intCrit="eval", ... )
  ##if ( plot ) plot( res ) ##viewPssm( t( attr( attr( res, "pwm" ), "pwm" )[ col.let, ] ) )
  ##res
  fname <- my.tempfile( "gadem_" )
  cat( paste( ">", names( seqs ), "\n", seqs, sep="" ), file=fname, sep="\n" )
  bgtab <- data.frame( names( genome.info$bg.list[[ seq.type ]] ),
                      sapply( genome.info$bg.list[[ seq.type ]], "[", 1 ) )
  rownames( bgtab ) <- NULL; colnames( bgtab ) <- c( "V1", "V2" )
  bgtab <- bgtab[ -1, ]
  bgtab$V1 <- tolower( bgtab$V1 )
  bgfname <- my.tempfile( "gadem_bg_" )
  write.table( bgtab, bgfname, sep="\t", quote=F, col.names=F, row.names=F )
  outfname <- my.tempfile( "gadem_out_" )

  cmd <- sprintf( "%s/gadem -fseq %s -minN %d -bOrder %d -fbm %s -fout %s -pgf 0 -verbose %d", progs.dir, fname,
                 floor( sqrt( length( seqs ) ) ), bg.order[ seq.type ], bgfname, outfname, if ( verbose ) 1 else 0 )
  paste( cmd, "-pgf 0 -gen 20 -pop 100 -em 200 -fEM 0.9 -maxgap 18 -useScore 1" ) ## other options ?
  if ( verbose ) print( cmd )
  out <- system( cmd, intern=T, ignore.stderr=!verbose )  
  out <- readLines( outfname )
  if ( unlink ) unlink( c( fname, bgfname, outfname ) )
  out
}

## Can try using glam2 (gapped motif finder), e.g. command to find HPY flag. motif is:
##  system( "./progs/glam2 -2 -a 6 -b 25 -q 1 -o zzzseqs_glam.out n zzzseqs.fna" )
## Then scan the motif across all upstream seqs (get output, easily parsable):
##  qqq <- system( "./progs/glam2scan -n 5000 -2 n zzzseqs_glam.out zzzseqs.all.fna", intern=T )
## Then search for 2nd strongest motif by masking:
##  system( "./progs/glam2mask -o zzzseqs_masked1.fna zzzseqs_glam.out zzzseqs.fna" )
##  system( "./progs/glam2 -2 -a 6 -b 25 -q 1 -o zzzseqs_glam2.out n zzzseqs_masked1.fna" )
##  qqq2 <- system( "./progs/glam2scan -n 5000 -2 n zzzseqs_glam2.out zzzseqs.all.fna", intern=T )
glam.one.cluster <- function( k, seq.type="upstream", min.seqs=cluster.rows.allowed[ 1 ],
                               max.seqs=cluster.rows.allowed[ 2 ], verbose=T, unlink=T, ... ) {
  if ( is.numeric( k ) ) rows <- get.rows( k )
  else rows <- k
  seqs <- get.sequences( rows, seq.type=seq.type, ... ) ##mask=mask, blast=blast, 
  if ( is.null( seqs ) || length( seqs ) < min.seqs ) return( list( k=k ) ) 
  if ( length( seqs ) < min.seqs || length( seqs ) > max.seqs ) return( list( k=k ) ) 

  ##if ( is.null( seqs ) || length( seqs ) < min.meme.seqs ) return( NULL )
  ##if ( uniq ) seqs <- seqs[ ! duplicated( seqs ) ] ## unique() doesnt preserve names
  ##if ( length( seqs ) < min.meme.seqs || length( seqs ) > max.meme.seqs ) return( NULL )
  all.seqs <- genome.info$all.upstream.seqs[[ seq.type ]]
  ##if ( verbose ) cat( length( seqs ), "SEQUENCES.\n" )

  tmp <- strsplit( meme.cmd, " " )[[ 1 ]]
  w.min <- tmp[ which( tmp == "-minw" ) + 1 ]
  w.max <- tmp[ which( tmp == "-maxw" ) + 1 ]

  fna.file <- my.tempfile( "glamseqs.fna." )
  fna1.file <- my.tempfile( "glamseqs1.fna." )
  all.fna.file <- my.tempfile( "glamseqs.all.fna." )
  out.file <- my.tempfile( "glam.out." )
  
  glam.out <- list()
  cat( paste( ">", names( seqs ), "\n", seqs, sep="" ), file=fna.file, sep="\n", append=F )
  file.copy( fna.file, fna1.file, overwrite=T )
  cat( paste( ">", names( all.seqs ), "\n", all.seqs, sep="" ), file=all.fna.file, sep="\n", append=F )
  for ( i in 1:n.motifs[[ seq.type ]][ iter ] ) {
    cmd <- paste( sprintf( "%s/glam2 -2 -n 5000 -a", progs.dir ), w.min, "-b", w.max, "-q 1 -O", out.file, "n",
                 fna1.file )
    if ( verbose ) cat( i, cmd, "\n" )
    system( cmd )
    glam.out[[ i ]] <- list()
    glam.out[[ i ]]$txt <- readLines( paste( out.file, "/glam2.txt", sep="" ) )
    glam.out[[ i ]]$meme <- readLines( paste( out.file, "/glam2.meme", sep="" ) )
    cmd <- paste( sprintf( "%s/glam2scan -n 5000 -2 n ", progs.dir ), out.file, "/glam2.txt ", all.fna.file, sep="" )
    if ( verbose ) cat( i, cmd, "\n" )
    glam.out[[ i ]]$scan <- system( cmd, intern=T )
    if ( i < n.motifs[[ seq.type ]][ iter ] ) {
      cmd <- paste( sprintf( "%s/glam2mask -o ", progs.dir ), fna1.file, " ", out.file, "/glam2.txt ",
                   fna.file, sep="" ) ##zzzseqs1.fna zzzseqs_glam.out/glam2.txt zzzseqs.fna"
      if ( verbose ) cat( i, cmd, "\n" )
      system( cmd )
      file.copy( fna1.file, fna.file, overwrite=T )
    }
  }
  if ( unlink ) system( paste( 'rm -rf', fna.file, fna1.file, all.fna.file, out.file ) )

  ## Parse out the scores, and start/stop/strand
  for ( i in 1:length( glam.out ) ) {
    tmp <- glam.out[[ i ]]$scan
    tmp <- tmp[ tmp != "" ]
    tmp <- tmp[ ! 1:length( tmp ) %in% grep( "^\\s+", tmp, perl=T ) ]
    tmp <- tmp[ -(1:3) ]
    tmp <- do.call( rbind, strsplit( tmp, "\\s+" ) )
    scores <- as.numeric( tmp[ ,6 ] ); names( scores ) <- tmp[ ,1 ]
    ##glam.out[[ i ]]$scores <- scores
    glam.out[[ i ]]$posns <- data.frame( gene=tmp[ ,1 ], start=as.integer( tmp[ ,2 ] ), end=as.integer( tmp[ ,4 ] ),
                                        strand=tmp[ ,5 ], score=scores )
  }

  ## Parse out the pssms
  for ( i in 1:length( glam.out ) ) {
    tmp <- glam.out[[ i ]]$txt
    score <- as.numeric( strsplit( grep( "^Score:\\s+", tmp, perl=T, val=1 )[ 1 ], "\\s+", perl=T )[[ 1 ]][ 2 ] )
    start <- grep( "a  c  g  t Del Ins Score", tmp, fixed=T )[ 1 ] + 1
    end <- grep( "^Score:\\s+", tmp, perl=T )[ 2 ] - 2
    tmp <- tmp[ start:end ]
    tmp <- gsub( "^\\s+", "", tmp[ seq( 1, length( tmp ), by=2 ) ], perl=T )
    tmp <- do.call( rbind, strsplit( tmp, "\\s+" ) )[ ,1:4 ]
    tmp <- t( apply( tmp, 1, as.numeric ) )
    end <- grep( "^Score:\\s+", tmp, perl=T )[ 2 ] - 2
    attr( tmp, "score" ) <- score
    colnames( tmp ) <- c( "A", "C", "G", "T" ); glam.out[[ i ]]$pssm <- tmp
  }

  for ( i in 1:length( glam.out ) ) {
    attr( glam.out[[ i ]], "meme.out" ) <- glam.out[[ i ]]$meme
    glam.out[[ i ]]$txt <- glam.out[[ i ]]$scan <- glam.out[[ i ]]$meme <- NULL
  }
  invisible( glam.out )
}

## Write out a single bicluster's (usually 2) motifs to MEME format
cluster.meme.motif.lines <- function( k, seq.type=names( meme.scores )[ 1 ], logodds=F ) { 
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )
  memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
  if ( is.null( memeOut ) ) return( lines )
  for ( i in 1:length( memeOut ) ) {
    pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
    mat.type <- "letter-probability matrix"
    if ( logodds ) mat.type <- "log-odds matrix"
    lines <- c( lines, "",
               sprintf( "MOTIF bic_%03d_%02d", k, i ),
               sprintf( "BL   MOTIF bic_%03d_%02d width=0 seqs=0", k, i ),
               sprintf( "%s: alength= 4 w= %d nsites= %d E= %.3e", mat.type, nrow( pssm ),
                       memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
    if ( ! logodds ) lines <- c( lines, apply( pssm, 1, function( i )
                                              sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    else lines <- c( lines, apply( round( log( pssm + 0.01 ) ), 1, function( i )
                                  sprintf( "%6d %6d %6d %6d", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
  }
  lines
}

all.motifs.to.mast.file <- function( ks=1:k.clust, seq.type=names(mot.weights)[1],
                                    e.value.cutoff=100, resid.cutoff=0.8 ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  cluster.motif.lines <- function( k ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    lines <- character()
    memeOut <- meme.scores[[ seq.type ]][[ k ]]
    if ( is.null( memeOut ) || memeOut == "" ) return( lines )
    memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
    if ( is.null( memeOut ) ) return( lines )
    if ( clusterStack[[ k ]]$resid > resid.cutoff ) return( lines )
    ##max.motifs <- max( max.motifs, length( memeOut ) )
    for ( i in 1:length( memeOut ) ) {
      if ( memeOut[[ i ]]$e.value > e.value.cutoff ) next
      pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
      lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                         memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
      lines <- c( lines, apply( pssm, 1, function( i )
                               sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    }
    lines
  }

  lines.t <- c( lines, do.call( c, lapply( ks, cluster.motif.lines ) ) )
  tfile <- my.tempfile( "tomtom_t_", )
  cat( lines.t, file=tfile, sep="\n" ) ## Write out all motifs
  tfile
}  

## All vs. all comparison of motifs in every cluster using MEME package's tomtom program
## Note that tomtom can only handle up to 5000 motifs at a time!!!
## But this can be changed in the tomtom.c code -- I have done it on pinnacle to go up to 15000 motifs
motif.similarities.tomtom <- function( query=1:k.clust, target=1:k.clust, query.mot=NA, target.mot=NA,
                                      seq.type=names(mot.weights)[1],
                                      e.value.cutoff=100, resid.cutoff=0.8, dist.meth="ed", q.thresh=0.5,
                                      min.overlap=4, q.pseudo=0, t.pseudo=0, min.gene.overlap=NA,
                                      desymmetrize=T, unlink=T, files.only=F, verbose=T, ... ) {
  meme.let <- c( "A", "C", "G", "T" )
  lines <- c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "",
             "Background letter frequencies (from dataset with add-one prior applied):" )
  lines <- c( lines, paste( names( unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ),
                       sprintf( "%.3f", unlist( genome.info$bg.list[[ seq.type ]][ meme.let ] ) ), collapse=" " ) )

  query <- query[ ! is.na( query ) ]
  target <- target[ ! is.na( target ) ]
  if ( is.na( query.mot ) ) query.mot <- rep( NA, length( query ) )
  if ( is.na( target.mot ) ) target.mot <- rep( NA, length( target ) )
  ##max.motifs <- 0
  
  cluster.motif.lines <- function( k, mot ) { ## Write out a single bicluster's (usually 2) motifs to MEME format
    lines <- character()
    memeOut <- meme.scores[[ seq.type ]][[ k ]]
    if ( is.null( memeOut ) || memeOut == "" ) return( lines )
    memeOut <- meme.scores[[ seq.type ]][[ k ]]$meme.out
    if ( is.null( memeOut ) ) return( lines )
    if ( clusterStack[[ k ]]$resid > resid.cutoff ) return( lines )
    ##max.motifs <- max( max.motifs, length( memeOut ) )
    for ( i in 1:length( memeOut ) ) {
      if ( ! is.na( mot ) && i != mot ) next
      if ( memeOut[[ i ]]$e.value > e.value.cutoff ) next
      pssm <- memeOut[[ i ]]$pssm ##; colnames( pssm ) <- meme.let ##col.let; pssm <- pssm[ ,meme.let ]
      lines <- c( lines, "",
                 sprintf( "MOTIF bic_%03d_%02d_%.3f_%.3e", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "BL   MOTIF bic_%03d_%02d_%.3f_%.3e width=0 seqs=0", k, i, clusterStack[[ k ]]$resid, memeOut[[ i ]]$e.value ),
                 sprintf( "letter-probability matrix: alength= 4 w= %d nsites= %d E= %.3e", nrow( pssm ),
                         memeOut[[ i ]]$sites, memeOut[[ i ]]$e.value ) )
      lines <- c( lines, apply( pssm, 1, function( i )
                               sprintf( "%5.3f %5.3f %5.3f %5.3f", i[ 1 ], i[ 2 ], i[ 3 ], i[ 4 ] ) ) )
    }
    lines
  }

  cmd <- "%s/tomtom -verbosity 1 -q-thresh %.3f -dist %s -min-overlap %d -text -query-pseudo %.3f -target-pseudo %.3f -target %s"
  
  ##lines.q <- lines.t <- lines
  if ( is.na( min.gene.overlap ) ) { ## If no cluster-specific targets, just create target file once for all clusters
    lines.t <- c( lines, do.call( c, lapply( target, function( k )
                                            cluster.motif.lines( k, target.mot[ which( target == k ) ] ) ) ) )
    if ( verbose ) cat( "TARGET MOTIFS:", length( grep( "MOTIF", lines.t ) ) / 2, range( target ), "\n" )
    tfile <- my.tempfile( "tomtom_t_", ); cat( lines.t, file=tfile, sep="\n" ) ## Write out all target motifs
    cmd <- sprintf( cmd, progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile )
  }

  if ( ! is.na( min.gene.overlap ) ) c.rows <- lapply( target, get.rows )
  mc <- get.parallel( length( query ) ) ## Parallelize running query motifs against all target motifs
  if ( is.na( mc$par ) ) mc$par <- 1 ## check on this!
  if ( mc$par > length( query ) ) mc$par <- length( query )
  if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) ) {
    mc$apply <- lapply; mc$par <- 1 }
  tout <- do.call( rbind, ##mc$apply( query, function( k ) {
                  mc$apply( 1:mc$par, function( par ) {
    ks <- query[ seq( par, length( query ), by=mc$par ) ] ##query[ query %in% seq( par, max( query ), by=mc$par ) ]
    ##ks <- k
    lines.q <- c( lines, do.call( c, lapply( ks, function( k )
                                            cluster.motif.lines( k, query.mot[ which( query == k ) ] ) ) ) )
    if ( is.na( query.mot ) ) n.query.motifs <- length( grep( "MOTIF", lines.q ) ) / 2
    else n.query.motifs <- length( grep( "MOTIF", lines.q ) )
    if ( verbose ) cat( "QUERY MOTIFS:", n.query.motifs, range( ks ), "\n" )
    if ( n.query.motifs <= 0 ) return( NULL )
    if ( ! is.na( min.gene.overlap ) ) {
      for ( k in target ) {
        rows <- get.rows( k )
        t.ok <- sapply( c.rows, function( r ) sum( r %in% rows ) >= min.gene.overlap )
        t.ok[ k ] <- FALSE
      }
      lines.t <- c( lines, do.call( c, lapply( target[ t.ok ], function( k )
                                              cluster.motif.lines( k, target.mot[ which( target == k ) ] ) ) ) )
      if ( verbose ) cat( "TARGET MOTIFS:", length( grep( "MOTIF", lines.t ) ) / 2, range( target ), "\n" )
      tfile <- my.tempfile( "tomtom_t_", ); cat( lines.t, file=tfile, sep="\n" ) ## Write out all target motifs
      cmd <- sprintf( cmd, progs.dir, q.thresh, dist.meth, min.overlap, q.pseudo, t.pseudo, tfile )
    }
    ##qfile <- my.tempfile( "tomtom_q_" )
    qfile <- paste( gsub( "tomtom_t_", "tomtom_q_", tfile ), "_", min( ks, na.rm=T ), "_",
                   max( ks, na.rm=T ), sep="" )
    cat( lines.q, file=qfile, sep="\n" )
    cmd <- paste( cmd, "-query", qfile )
    if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) ) {
      if ( is.character( files.only ) ) cat( paste( cmd, " > ", qfile, ".out\n", sep="" ), file=files.only, append=T )
      return( paste( cmd, " > ", qfile, ".out", sep="" ) )
    }
    if ( verbose ) print( cmd )
    tout <- system( cmd, intern=T )
    if ( unlink ) unlink( qfile )
    tout <- do.call( rbind, strsplit( tout, "\t" ) )
    if ( ! is.null( tout ) ) {
      colnames( tout ) <- tout[ 1, ,drop=F ]; tout <- tout[ -1, ,drop=F ]
      tout <- as.data.frame( tout[ tout[ ,1 ] != tout[ ,2 ], ,drop=F ] ) ## < ##de-symmetrize the output )upper-tri)
    }
    print(dim(tout))
    tout
  } ) ) ##, mc.preschedule=F ) )
  if ( files.only == TRUE || ( ! is.logical( files.only ) && ! is.na( files.only ) && ! is.null( files.only ) ) )
    return( tout )

  cat( "GOT", nrow( tout ), "motif alignments.\n" )
  ## Summarize into data frame
  q.id <- t( sapply( strsplit( as.character( tout[ ,1 ] ), "_" ), function( i ) as.integer( i[ 2:3 ] ) ) )
  q.res.ev <- t( sapply( strsplit( as.character( tout[ ,1 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) ) )
  t.id <- t( sapply( strsplit( as.character( tout[ ,2 ] ), "_" ), function( i ) as.integer( i[ 2:3 ] ) ) )
  t.res.ev <- t( sapply( strsplit( as.character( tout[ ,2 ] ), "_" ), function( i ) as.numeric( i[ 4:5 ] ) ) )
  tout2 <- data.frame( biclust1=as.integer( q.id[ ,1 ] ), motif1=as.integer( q.id[ ,2 ] ),
                      resid1=as.numeric( q.res.ev[ ,1 ] ), e.value1=as.numeric( q.res.ev[ ,2 ] ),
                      biclust2=as.integer( t.id[ ,1 ] ), motif2=as.integer( t.id[ ,2 ] ),
                      resid2=as.numeric( t.res.ev[ ,1 ] ), e.value2=as.numeric( t.res.ev[ ,2 ] ),
                      offset=as.integer( as.character( tout$`Optimal offset` ) ),
                      p.value=as.numeric( as.character( tout$`p-value` ) ),
                      q.value=as.numeric( as.character( tout$`q-value` ) ),
                      overlap=as.integer( as.character( tout$`Overlap` ) ),
                      consensus1=as.factor( as.character( tout$`Query consensus` ) ),
                      consensus2=as.factor( ifelse( tout$`Orientation` == "-",
                        rev.comp( as.character( tout$`Target consensus` ) ),
                        as.character( tout$`Target consensus` ) ) ),
                      orientation=tout$`Orientation` ) ##,
  if ( exists( cmd ) ) attr( tout2, "tomtom.cmd" ) <- cmd
  tout2 <- tout2[ order( tout2$q.value, tout2$p.value ), ]
  rm( tout )
  
  if ( desymmetrize ) tout2 <- desymmetrize.tomtom.results( tout2 )
  ##else tout <- tout2
  tout2
}

## lower tri != upper tri because of asymmetries in background - so lets use the min of a vs b and b vs a:
## do it matrix-y by min-ing the upper- and lower- tri of a p-value matrix
desymmetrize.tomtom.results <- function( tt.out ) {
  mot.names <- unique( c( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                         paste( tt.out$biclust2, tt.out$motif2, sep="_" ) ) )
  mot.names.2 <- cbind( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                       paste( tt.out$biclust2, tt.out$motif2, sep="_" ) )
  mot.names.2a <- cbind( mot.names.2[ ,2 ], mot.names.2[ ,1 ] )
  tmp <- matrix( NA, nrow=length( mot.names ), ncol=length( mot.names ) )
  rownames( tmp ) <- colnames( tmp ) <- mot.names
  lt <- lower.tri( tmp )
  tmp[ mot.names.2 ] <- tt.out$p.value
  tmp2 <- cbind( tmp[ lt ], t( tmp )[ lt ] )
  tmp2[ is.na( tmp2 ) ] <- Inf
  tmp[ lt ] <- apply( tmp2, 1, min, na.rm=T )
  tmp[ is.infinite( tmp ) ] <- NA
  rm( tmp2 )
  tmp.good <- ! is.na( tmp[ lt ] )

  tmp2 <- tmp * NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$biclust1
  out <- data.frame( biclust1=tmp2[ lt ][ tmp.good ] )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$motif1
  out <- cbind( out, data.frame( motif1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$resid1
  out <- cbind( out, data.frame( resid1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$e.value1
  out <- cbind( out, data.frame( e.value1=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$biclust2
  out <- cbind( out, data.frame( biclust2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$motif2
  out <- cbind( out, data.frame( motif2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$resid2
  out <- cbind( out, data.frame( resid2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$e.value2
  out <- cbind( out, data.frame( e.value2=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$offset
  out <- cbind( out, data.frame( offset=tmp2[ lt ][ tmp.good ] ) )
  out <- cbind( out, data.frame( p.value=tmp[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$q.value
  out <- cbind( out, data.frame( q.value=tmp2[ lt ][ tmp.good ] ) )
  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- tt.out$overlap
  out <- cbind( out, data.frame( overlap=tmp2[ lt ][ tmp.good ] ) )
  ##if ( is.factor( tt.out$orientation ) ) {
  ##  tmp2[,] <- NA; tmp2[ mot.names.2 ] <- tmp2[ mot.names.2a ] <- as.numeric( tt.out$orientation )
  ##  out <- cbind( out, data.frame( orientation=tmp2[ lt ][ tmp.good ] ) )
  ##}
  
  rm( tmp2 )
  if ( "consensus1" %in% colnames( tt.out ) || "consensus2" %in% colnames( tt.out ) ||
      is.character( tt.out$orientation ) ) {
    tmp.c <- matrix( "", nrow=length( mot.names ), ncol=length( mot.names ) )
    rownames( tmp.c ) <- colnames( tmp.c ) <- mot.names
    if ( "consensus1" %in% colnames( tt.out ) ) {
      tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$consensus1 )
      out <- cbind( out, data.frame( consensus1=tmp.c[ lt ][ tmp.good ] ) )
    }
    if ( "consensus2" %in% colnames( tt.out ) ) {
      tmp.c[,] <- ""; tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$consensus2 )
      out <- cbind( out, data.frame( consensus2=tmp.c[ lt ][ tmp.good ] ) )
    }
    if ( is.factor( tt.out$orientation ) ) tt.out$orientation <- as.character( tt.out$orientation )
    if ( is.character( tt.out$orientation ) ) {
      tmp.c[,] <- ""; tmp.c[ mot.names.2 ] <- tmp.c[ mot.names.2a ] <- as.character( tt.out$orientation )
      out <- cbind( out, data.frame( orientation=tmp.c[ lt ][ tmp.good ] ) )
    }
    rm( tmp.c )
  }
  
  out <- subset( out, ! is.na( p.value ) )
  out <- out[ order( out$p.value, out$q.value ), ]
  out
}

## agglomerative clustering with complete linkage (thanks Chris!) of tomtom similarities to get
##   groups of distinct motifs
## Idea - use # genes to weight each motif's contribution to the output "combined" pssms
cluster.tomtom.results <- function( tt.out, seq.type=names(mot.weights)[1],
                                   e.value.cutoff=Inf, p.value.cutoff=0.05, resid.cutoff=Inf, n.cutoff=3,
                                   make.pssms=T, min.size=3, k.cut=0.5, return.aligned.pssms=F,
                                   n.gene.weight=F, ... ) { ##, plot=F ) {
  if ( is.na( e.value.cutoff ) ) e.value.cutoff <- Inf
  if ( is.na( resid.cutoff ) ) resid.cutoff <- Inf
  if ( is.na( p.value.cutoff ) ) p.value.cutoff <- Inf
  ##if ( ! is.na( e.cutoff ) ) tt.out <- subset( tt.out, e.value1 <= e.cutoff & e.value2 <= e.cutoff )
  ##if ( ! is.na( resid.cutoff ) ) tt.out <- subset( tt.out, resid1 <= resid.cutoff & resid2 <= resid.cutoff )
  ##if ( ! is.na( p.cutoff ) ) tt.out <- subset( tt.out, p.value <= p.cutoff )
  if ( ! is.infinite( e.value.cutoff ) || ! is.infinite( resid.cutoff ) || ! is.infinite( p.value.cutoff ) )
    tt.out <- subset( tt.out, e.value1 <= e.value.cutoff & e.value2 <= e.value.cutoff &
                     resid1 <= resid.cutoff & resid2 <= resid.cutoff & p.value <= p.value.cutoff )
  mot.names <- c( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                 paste( tt.out$biclust2, abs( tt.out$motif2 ), sep="_" ) )
  mot.tab <- sort( table( mot.names ), decreasing=T )
  mot.names <- names( mot.tab )
  mot.names.2 <- cbind( paste( tt.out$biclust1, tt.out$motif1, sep="_" ),
                       paste( tt.out$biclust2, abs( tt.out$motif2 ), sep="_" ) )
  tmp <- matrix( 1, nrow=length( mot.names ), ncol=length( mot.names ) )
  rownames( tmp ) <- colnames( tmp ) <- mot.names
  tmp[ mot.names.2 ] <- tt.out$p.value
  tmp[ cbind( mot.names.2[ ,2 ], mot.names.2[ ,1 ] ) ] <- tt.out$p.value
  tmp[ tmp > p.value.cutoff ] <- 1
  hc <- hclust( as.dist( tmp ), method="complete" )
  if ( k.cut < 1 ) cs <- cutree( hc, h=k.cut )
  else cs <- cutree( hc, k=k.cut )
  rm( tmp ); gc()
  ##out.tt <- list()
  ##for ( i in as.integer( names( sort( table( cs ), decreasing=T ) ) ) ) {
  ##require( multicore )
  cat( "HERE:", range( as.integer( names( sort( table( cs ), decreasing=T ) ) ) ), "\n" )
  meme.let <- c( "A", "C", "G", "T" )
  mc <- get.parallel( length( table( cs ) ) ) ## Parallelize running query motifs against all target motifs
  out.tt <- mc$apply( as.integer( names( sort( table( cs ), decreasing=T ) ) ), function( i ) {
    ##out.tt <- list()
    wh <- t( apply( do.call( rbind, strsplit( names( which( cs == i ) ), "_" ) ), 1, as.integer ) )
    if ( nrow( wh ) < n.cutoff ) return( NULL ) ##next
    cat( i, nrow( wh ), "\n" )
    wh.tmp <- apply( wh, 1, paste, collapse="_" )
    b1 <- paste( tt.out$biclust1, abs( tt.out$motif1 ), sep="_" )
    b2 <- paste( tt.out$biclust2, abs( tt.out$motif2 ), sep="_" )
    tt <- subset( tt.out, b1 %in% wh.tmp & b2 %in% wh.tmp )
    if ( nrow( tt ) <= 0 ) return( NULL ) ##next
    if ( "q.value" %in% names( tt ) ) tt <- tt[ order( tt$p.value, tt$q.value ), ]
    else tt <- tt[ order( tt$p.value ), ]
    b1 <- paste( tt$biclust1, abs( tt$motif1 ), sep="_" )
    b2 <- paste( tt$biclust2, abs( tt$motif2 ), sep="_" )
    ##if ( "e.value1" %in% names( tt ) ) tt <- tt[ order( tt$e.value1, tt$e.value2 ), ]
    ##else tt <- tt[ order( tt$p.value ), ]
    if ( make.pssms ) {
      tto <- tt
      mot.names <- unique( c( b1, b2 ) )
      ##rm( b1, b2, wh.tmp )
      ##if ( length( mot.names ) <= 3 ) next
      tmp1 <- as.integer( strsplit( mot.names[ 1 ], "_" )[[ 1 ]] )
      tmp.mot <- meme.scores[[ seq.type ]][[ tmp1[ 1 ] ]]$meme.out[[ tmp1[ 2 ] ]]
      pssm <- orig.pssm <- tmp.mot$pssm
      colnames( pssm ) <- meme.let
      if ( n.gene.weight ) pssm <- pssm * tmp.mot$sites
      if ( return.aligned.pssms ) aligned.pssms <- list()

      if ( length( mot.names ) >= min.size ) {
        ## Get first motif, then take all alignments to it, and construct combined pssm      
        if ( "width1" %in% names( tto ) ) max.width <- max( c( tto$width1, tto$width2 ) )
        else max.width <- max( nchar( c( as.character( tto$consensus1 ), as.character( tto$consensus2 ) ) ) )
        for ( jj in 1:max.width ) pssm <- rbind( rep( 0, 4 ), pssm, rep( 0, 4 ) )
        if ( return.aligned.pssms ) aligned.pssms[[ mot.names[ 1 ] ]] <- pssm
        first.ind <- which( apply( pssm, 1, function( j ) any( j != 0 ) ) )[ 1 ]
        ##last.ind <- first.ind + nrow( meme.scores[[ seq.type ]][[ tmp1[ 1 ] ]]$meme.out[[ tmp1[ 2 ] ]]$pssm ) - 1
        orig.width <- nrow( orig.pssm ) ##meme.scores[[ seq.type ]][[ tmp1[ 1 ] ]]$meme.out[[ tmp1[ 2 ] ]]$pssm )
        ##colnames( pssm ) <- meme.let
        tttt <- subset( tto, ( biclust1 == tmp1[ 1 ] & motif1 == tmp1[ 2 ] ) |
                       ( biclust2 == tmp1[ 1 ] & abs( motif2 ) == tmp1[ 2 ] ) )
        rm( tto )
        
        for ( m in mot.names[ 2:length( mot.names ) ] ) {
          ##cat(i,m,length(mot.names),"\n")
          tmp <- as.integer( strsplit( m, "_" )[[ 1 ]] )
          ##if ( all( tmp == tmp1 ) ) next
          ttt <- unique( subset( tttt, ( biclust1 == tmp[ 1 ] & abs( motif1 ) == tmp[ 2 ] ) | 
                        ( biclust2 == tmp[ 1 ] & abs( motif2 ) == tmp[ 2 ] ) ) )
          if ( nrow( ttt ) <= 0 ) next
          else if ( nrow( ttt ) > 1 ) ttt <- ttt[ 1, ]
          tmp.mot <- meme.scores[[ seq.type ]][[ tmp[ 1 ] ]]$meme.out[[ tmp[ 2 ] ]]
          pssm2 <- tmp.mot$pssm
          if ( n.gene.weight ) pssm2 <- pssm2 * tmp.mot$sites
          ##if(m=="128145_1")stop()
          if ( ( "orientation" %in% colnames( ttt ) && ttt$orientation == "-" ) || ttt$motif2 < 0 ) {
            if ( ttt$biclust2 == tmp1[ 1 ] && abs( ttt$motif2 ) == tmp1[ 2 ] )
              offset <- first.ind + orig.width - ttt$offset - nrow( pssm2 )
            else offset <- first.ind - ttt$offset ##- nrow( pssm2 ) ##+ orig.width 
            pssm2 <- pssm2[ ,4:1 ][ nrow( pssm2 ):1, ]
          } else {
            if ( ttt$biclust1 == tmp1[ 1 ] && abs( ttt$motif1 ) == tmp1[ 2 ] ) offset <- first.ind - ttt$offset
            else offset <- first.ind + ttt$offset
          }
          pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] <- pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] + pssm2
          if ( return.aligned.pssms ) {
            tmp.pssm <- pssm * 0
            tmp.pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] <-
              tmp.pssm[ offset:( offset + nrow( pssm2 ) - 1 ), ] + pssm2
            aligned.pssms[[ m ]] <- tmp.pssm
          }
        }
      }
      first.ind2 <- which( apply( pssm, 1, function( j ) any( j != 0 ) ) )[ 1 ]
      pssm <- pssm[ -( 1:( first.ind2 - 1 ) ), ]
      last.ind2 <- which( apply( pssm, 1, function( j ) all( j == 0 ) ) )[ 1 ]
      if ( ! is.na( last.ind2 ) ) pssm <- pssm[ -( ( last.ind2 - 1 ):nrow( pssm ) ), ]
      pssm <- ( pssm + 1e-9 ) / ( max( pssm, na.rm=T ) + 1e-9 )
      ##if ( plot ) viewPssm( pssm, main.title=paste( i, nrow( wh ), length( mot.names ) ) )
      attr( tt, "combined.pssm" ) <- pssm
      if ( return.aligned.pssms ) {
        for ( m in names( aligned.pssms ) ) {
          pssm <- aligned.pssms[[ m ]]
          pssm <- pssm[ -( 1:( first.ind2 - 1 ) ), ]
          pssm <- pssm[ -( ( last.ind2 - 1 ):nrow( pssm ) ), ]
          aligned.pssms[[ m ]] <- pssm
        }
        attr( tt, "aligned.pssms" ) <- aligned.pssms
      }
      attr( tt, "wh" ) <- wh
      attr( tt, "mot.names" ) <- mot.names
    }
    tt
  } )
  ## if ( plot ) {
  ##   for ( i in 1:length( out.tt ) ) {
  ##     if ( is.null( out.tt[[ i ]] ) ) next
  ##     wh <- attr( out.tt[[ i ]], "wh" )
  ##     mot.names <- attr( out.tt[[ i ]], "mot.names" )
  ##     viewPssm( pssm, main.title=paste( i, nrow( wh ), length( mot.names ) ) )
  ##   }
  ## }
  attr( out.tt, "hc.out" ) <- hc
  out.tt
}

motif.similarities.custom <- function( query=1:(k.clust-1), target=query, seq.type=names(mot.weights)[1],
                                      e.value.cutoff=100, resid.cutoff=0.99, dist.meth="ed", min.overlap=4,
                                      p.values=T, p.cutoff=0.01, p.correct=T, out.p.only=F, verbose=F,
                                      consensus=F, filter=T, add.to=NULL ) {
  mc <- get.parallel( k.clust ) ## Parallelize running query motifs against all target motifs
  mc$apply <- lapply

  ##out <- list()
  query <- sort( query ); target <- sort( target )
  ##for ( i in query ) { ##1:( k.clust - 1 ) ) {
  outA <- do.call( c, mc$apply( query, function( i ) {
    out <- list()
    if ( verbose ) cat( i, "\n" )
    if ( length( get.rows( i ) ) <= 1 ) return( out )
    if ( clusterStack[[ i ]]$resid > resid.cutoff ) return( out ) ##next
    memeOut <- meme.scores[[ seq.type ]][[ i ]]$meme.out
    if ( is.null( memeOut ) ) return( out ) ##next
    for ( ii in 1:length( memeOut ) ) {
      if ( is.na( memeOut[[ ii ]]$e.value ) || memeOut[[ ii ]]$e.value > e.value.cutoff ) next
      pssm1 <- memeOut[[ ii ]]$pssm
      if ( is.null( pssm1 ) ) next
      for ( j in target[ target >= i ] ) { ##(i+1):k.clust ) {
        ##if ( verbose && ii == 1 ) cat( i, j, "\n" )
        if ( length( get.rows( j ) ) <= 1 ) next
        if ( clusterStack[[ j ]]$resid > resid.cutoff ) next
        memeOut2 <- meme.scores[[ seq.type ]][[ j ]]$meme.out
        if ( is.null( memeOut2 ) ) next
        for ( jj in 1:length( memeOut2 ) ) {
          if ( is.na( memeOut2[[ jj ]]$e.value ) || memeOut2[[ jj ]]$e.value > e.value.cutoff ) next
          pssm2 <- memeOut2[[ jj ]]$pssm
          if ( is.null( pssm2 ) ) next
          if ( ! is.null( add.to ) && nrow( subset( add.to, biclust1 == i & motif1 == ii & biclust2 == j &
                                                   motif2 == jj ) ) >= 1 ) next
          ##if ( verbose == 2 ) cat( i, ii, j, jj, "\n" )
          ## Last column is # of tests (increases for each shift, and 2x again for reverse test)
          tmp <- compare.pssms( pssm1, pssm2, rev=F, weight=F, score=dist.meth, min.ov=min.overlap )
          out[[ paste( i, ii, j, jj ) ]] <- cbind( i, ii, clusterStack[[ i ]]$resid, memeOut[[ ii ]]$e.value,
                                                  memeOut[[ ii ]]$width,
                                                  ##if ( consensus ) pssm.to.string( pssm1 ) else "",
                                                  j, jj, clusterStack[[ j ]]$resid,
                                                  memeOut[[ jj ]]$e.value, memeOut[[ jj ]]$width,
                                                  ##if ( consensus ) pssm.to.string( pssm2 ) else "",
                                                  tmp, nrow( tmp ) ) 
          tmp <- compare.pssms( pssm1, pssm2, rev=T, weight=F, score=dist.meth, min.ov=min.overlap )
          out[[ paste( i, ii, j, -jj ) ]] <- cbind( i, ii, clusterStack[[ i ]]$resid, memeOut[[ ii ]]$e.value,
                                                   memeOut[[ ii ]]$width,
                                                   ##if ( consensus ) pssm.to.string( pssm1 ) else "",
                                                   j, -jj, clusterStack[[ j ]]$resid,
                                                   memeOut[[ jj ]]$e.value, memeOut[[ ii ]]$width,
                                                   ##if ( consensus ) pssm.to.string( pssm2 ) else "",
                                                   tmp, nrow( tmp ) ) ##abs( tmp[ ,1 ] ) + 1 ) ##( abs( tmp[ ,1 ] ) * 2 + 1 ) * 2 )
        }
      }
    }
    out
  } ) ) ##, mc.preschedule=F ) )

  outA <- do.call( rbind, outA ); rownames( outA ) <- NULL
  if ( is.null( outA ) ) return( outA )
  ##outA <- as.data.frame( outA )
  colnames( outA ) <- c( "biclust1", "motif1", "resid1", "e.value1", "width1", ##"consensus1",
                       "biclust2", "motif2", "resid2", "e.value2", "width2", ##"consensus2",
                       "offset", "overlap", dist.meth, "n.tests" )
  ##for ( i in c( 1, 2, 5, 6, 7, 10, 11, 12, 14 ) ) outA[[ i ]] <- as.factor( outA[[ i ]] ) ## Make it smaller (RAM)
  quantiles <- cdfs <- NULL
  out2 <- NULL
  ## Compute for each overlap value, the score corresponding to p.value=0.01
  ovs <- table( outA[ ,"overlap" ] )
  if ( ! is.na( p.cutoff ) ) {
    ##all.out <- rbind( do.call( rbind, fwd ), do.call( rbind, rev ) )
    quantiles <- do.call( rbind, mc$apply( as.integer( names( ovs ) ), function( ov )
                                          c( ov, quantile( outA[ outA[ ,"overlap" ] %in%
                             ( ov + ( if ( ovs[ as.character( ov ) ] < 1000 ) c( -1, 0, 1 ) else 0 ) ), dist.meth ],
                                                          1-p.cutoff ) ) ) )
  }

  ## Compute for each overlap value and score, the corresponding p-value using the CDFs of the scores distribution
  ##    for the given overlap
  if ( p.values ) { 
    ##cdfs <- list()
    ##for ( ov in as.integer( names( ovs ) ) ) cdfs[[ ov ]] <- ecdf( outA[ outA[ ,6 ] %in%
    cdfz <- mc$apply( as.integer( names( ovs ) ), function( ov ) ecdf( outA[ outA[ ,"overlap" ] %in%
                          ( ov + ( if ( ovs[ as.character( ov ) ] < 1000 ) c( -1, 0, 1 ) else 0 ) ), dist.meth ] ) )
    cdfs <- list(); for ( i in 1:length( ovs ) ) cdfs[[ as.integer( names( ovs ) )[ i ] ]] <- cdfz[[ i ]]; rm( cdfz )
    qq <- numeric(); qq[ quantiles[ ,1 ] ] <- quantiles[ ,2 ]
    print( dim( outA ) )
    if ( filter ) { out2 <- outA[ outA[ ,dist.meth ] > qq[ outA[ ,"overlap" ] ], ]; rm( qq ) }
    else out2 <- outA
    if ( out.p.only ) rm( outA )##; gc() }
    ##out2 <- out2[ order( out2[ ,7 ], decreasing=T ), ]
    print( dim( out2 ) )
    ##out2 <- cbind( out2, p.value=1 - apply( out2, 1, function( j ) cdfs[[ j[ 6 ] ]]( j[ 7 ] ) ) )
    pvs <- rep( NA, nrow( out2 ) )
    tmp <- sapply( as.integer( names( ovs ) ), function( i ) {
      pvs[ which( out2[ ,"overlap" ] == i ) ] <- cdfs[[ i ]]( out2[ which( out2[ ,"overlap" ] == i ), dist.meth ] )
      pvs } )
    if ( out.p.only ) rm( cdfs )##; gc() }
    print( dim( out2 ) )
    pvs <- 1 - apply( tmp, 1, sum, na.rm=T ); rm( tmp )##; gc()
    pvs <- pvs * out2[ ,"n.tests" ] / 2 / 2 ## Correct for multiple tests and a HACK to make it close to tomtom's pvs... does this factor of 4 make any sense?
    out2 <- cbind( out2, p.value=pvs )
    print( dim( out2 ) )
    if ( filter && ! is.na( p.cutoff ) ) { out2 <- out2[ pvs <= p.cutoff, ]; pvs <- pvs[ pvs <= p.cutoff ] }
    out2 <- out2[ order( pvs ), ]
    print( dim( out2 ) )
    rownames( out2 ) <- NULL
    colnames( out2 )[ 1:14 ] <- c( "biclust1", "motif1", "resid1", "e.value1", "width1", ##"consensus1",
                                  "biclust2", "motif2", "resid2", "e.value2", "width2", ##"consensus2",
                                  "offset", "overlap", dist.meth, "n.tests" )
    ##out2 <- out2[ ! duplicated( abs( out2[ ,c( 1:2, 5:6 ) ] ) ), ]
  }
  if ( consensus ) {
    out2 <- as.data.frame( out2 )
    for ( ind in c( "biclust1", "motif1", "width1", "biclust2", "motif2", "width2", "offset", "overlap", "n.tests" ) )
      out2[[ ind ]] <- as.integer( out2[[ ind ]] )
    out2 <- cbind( out2, as.data.frame( t( apply( out2, 1, function( i ) c( pssm.to.string( meme.scores[[ seq.type ]][[ i[ "biclust1" ] ]]$meme.out[[ i[ "motif1" ] ]]$pssm ), pssm.to.string( meme.scores[[ seq.type ]][[ i[ "biclust2" ] ]]$meme.out[[ i[ "motif2" ] ]]$pssm ) ) ) ) ) )
      colnames( out2 )[ c( -1, 0 ) + ncol( out2 ) ] <- c( "consensus1", "consensus2" )
  }
  ##}
  ##if ( ! consensus && exists( "out" ) ) out <- out[ ,! colnames( out ) %in% c( "consensus1", "consensus2" ) ]
  ##list( fwd=fwd, rev=rev, quantiles=quantiles, cdfs=cdfs )
    print( dim( out2 ) )
  if ( ! out.p.only ) return( list( out=as.data.frame( outA ), out.p=as.data.frame( out2 ),
                                   quantiles=quantiles, cdfs=cdfs ) )
  else return( as.data.frame( out2 ) )
}

compare.pssms <- function( pssm1, pssm2, rev.comp=F, weight=F, min.ov=6, score="cor" ) {
  getEntropy = function( pssm ) {
    pssm[ pssm == 0 ] <- 0.00001
    entropy <- apply( pssm, 1, function( i ) -sum( i * log2( i ) ) )
    return( entropy )
  }

  if ( rev.comp ) pssm2 <- pssm2[ nrow( pssm2 ):1, 4:1 ]
  if ( score == "cor" || score == "pearson" ) cors <- cor( t( pssm1 ), t( pssm2 ), method="pearson" )
  else if ( score == "spearman" ) cors <- cor( t( pssm1 ), t( pssm2 ), method="spearman" )
  else if ( score == "ed" ) cors <-
    ##-apply( pssm2, 1, function( i ) apply( pssm1, 1, function( j ) sqrt( sum( (j-i)^2 ) ) ) )
    -as.matrix( dist( rbind( pssm1, pssm2 ), "euclidean" ) )[ 1:nrow( pssm1 ), nrow( pssm1 ) + 1:nrow( pssm2 ) ]
    ## Note package "proxy" has a 2-matrix 'dist' method but it is about 5x slower than using stats::dist!
    ##-as.matrix( proxy::dist( pssm1, pssm2, "Euclidean" ) )
  ## cors is pairwise distance between each row in each pssm
  if ( weight ) { ## Weight column match scores by (minimum of) information content of columns
    score1 <- ( 2 - getEntropy( pssm1 ) ) ## / 2
    score2 <- ( 2 - getEntropy( pssm2 ) ) ## / 2
    scores <- outer( score1, score2, fun="min" ) / 2
    cors <- cors * scores
  }

  rev <- FALSE
  if ( ncol( cors ) > nrow( cors ) ) { rev <- TRUE; cors <- t( cors ) }
  min.ind <- ncol( cors ) ##min( dim( cors ) )
  max.ind <- nrow( cors ) ##max( dim( cors ) )
  vec <- 1:min.ind
  covs <- do.call( rbind, lapply( (-max.ind+min.ov):(max.ind-min.ov), function( i ) {
    vec2 <- vec + i
    vec2a <- vec2 > 0 & vec2 <= max.ind
    vec2 <- vec2[ vec2a ]
    if ( length( vec2 ) < min.ov ) return( NULL )
    c( i, length( vec2 ), mean( cors[ cbind( vec2, vec[ vec2a ] ) ] ) )
  } ) )
  ##covs <- covs[ covs[ ,2 ] >= min.ov, ,drop=F ]
  if ( ! rev ) covs[ ,1 ] <- -covs[ ,1 ] ## Make sign of shift same as tomtom's
  ## Output shift (name) is amt. that pssm1 is shifted to get this match score w/ pssm2
  ##names( covs ) <- as.character( (-max.ind/2+1):(max.ind/2-1) )
  return( covs ) ##[ which.max( covs ) ] )
}


if ( FALSE ) { ## Test 'motif.similarities.custom' against tomtom in speed and p-values returned
  load("~/Sites/cMonkey/output/cmonkey_4.4.3_hpy_819x57_10_Feb_26_11:47:29/cm_session.RData")
  sys.source("~/scratch/biclust/cmonkey-motif-other.R",envir=e)
  sys.source("~/scratch/biclust/cmonkey-motif.R",envir=e)
  system.time(tmp1<-e$motif.similarities.tomtom(seq.type='upstream',min.ov=6,q.thresh=1,desym=F)) ## 15.1 sec
  system.time(tmp2<-e$motif.similarities.custom(seq.type='upstream',min.ov=6,verbose=T,p.cutoff=0.99)) ## 9.7 sec
  system.time(e$motif.similarities.tomtom(seq.type='upstream',min.ov=6,q.thresh=0.1,desym=T)) ## 14.6 sec
  system.time(e$motif.similarities.custom(seq.type='upstream',min.ov=6,verbose=T,p.cutoff=0.1)) ## 5.0 sec
  pv=NULL
  out2=tmp2$out.p
  ##out2$p.value=out2$p.value*out2$n.tests/2/2
  ##out2$p.value=1-(1-out2$p.value)^out2$n.tests
  for( i in 1:nrow(out2)){
    pv1=out2$p.value[i]
    pv2=subset(tmp1,biclust1==out2[i,1]&motif1==out2[i,2]&biclust2==out2[i,3]&motif2==abs(out2[i,4]))
    pv=rbind(pv,c(pv1,pv2$p.value))
    ##if(pv1>0.01&&pv2$p.value<=0.01)print(pv2)
  }
  par(mfrow=c(2,2))
  plot(log10(pv),ylim=c(-10,1),pch=20,cex=0.5,col=ifelse(pv[,1]>0.05&pv[,2]<=0.05,"red",
                                                ifelse(pv[,2]>0.05&pv[,1]<=0.05,"green",
                                                       ifelse(pv[,1]<=0.05&pv[,2]<=0.05,"blue","black") ) ) )
  plot(rank(pv[,1]),rank(pv[,2]),pch=20,cex=0.5,col=ifelse(pv[,1]>0.05&pv[,2]<=0.05,"red",
                                                ifelse(pv[,2]>0.05&pv[,1]<=0.05,"green",
                                                       ifelse(pv[,1]<=0.05&pv[,2]<=0.05,"blue","black") ) ) )
  hist(pv[1,]/(pv[,2]+1e-9),breaks=50)

  par(mfrow=c(5,5))
  tmp=e$cluster.tomtom.results(tmp2,resid.cutoff=0.5,plot=T)  
}
