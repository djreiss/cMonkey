###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

#ifndef PACKAGE

## Experimental and unfinished functions (mostly for data loading and post-proc analysis of clusters)
try( source( "KEGG_scripts/load.kegg.raw.new.R" ) )

cluster.GO.annotations <- function( env, k, which.pfc=c("P","F","C")[1] ) {
  if ( is.null( genome.info$goa ) ) {
    cat( "Loading GO ontology database ... this may take a while (it will be cached for subsequent calls).\n" )
    ## This func is in KEGG_scripts/load.kegg.raw.new.R:
    goa <- load.GO.annotations( org.id=genome.info$org.id$V1[ 1 ] )
    goa <- goa[ ! duplicated( goa[ ,-(5:6) ] ), ]
    env$genome.info$goa <- goa
    ##assign( "genome.info$goa", goa )
  }
  ##goa <- get( "genome.info$goa" )
  goa <- env$genome.info$goa
  if ( is.numeric( k[ 1 ] ) ) k <- get.rows( k )
  syns <- unique( unlist( get.synonyms( k ) ) )
  unique( goa[ goa$names %in% syns & goa$GO.PFC %in% which.pfc, c( "names", "GO.id", "GO.name", "GO.PFC" ) ] )
}

cluster.GO.pvalues <- function( env, ks=1:k.clust, ... ) {
  out <- data.frame()
  go.total <- NA
  for ( k in ks ) {
    if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
    else { rows <- ks; k <- '' }
    tab <- cluster.GO.annotations( env, rows, ... )
    if ( is.na( go.total ) ) go.total <- length( unique( genome.info$goa$names ) )
    syns <- unique( unlist( get.synonyms( rows ) ) )
    syns <- syns[ syns %in% genome.info$goa$names ]; names( syns ) <- NULL
    tt <- table( tab$GO.id )
    tab <- subset( tab, GO.id %in% names( tt[ tt > 1 ] ) )
    for ( i in names( tt[ tt > 1 ] ) ) { 
      ttt <- subset( tab, GO.id == i )
      gns <- paste( ttt$names, collapse=" ", sep=" " )
      cat.total <- nrow( subset( genome.info$goa, GO.id == i ) )
      out <- rbind( out, data.frame( k=k, genes=gns, n.in=length( syns ), n.hit=nrow( ttt ), n.tot=go.total,
                                    n.hit.tot=cat.total, p.value=phyper( nrow( ttt ), cat.total,
                                                           go.total - cat.total, length( syns ), lower=F ),
                                    GO.name=ttt$GO.name[ 1 ], GO.id=i ) )
    }
    if ( all( is.character( ks ) ) ) break
  }
  if ( nrow( out ) > 1 ) return( out[ order( out$p.value ), ] )
  out
}

cluster.KEGG.pvalues <- function( ks=1:k.clust ) {
  cluster.KEGG.annotations <- function( k ) {
    if ( ! exists( "genome.info$koa" ) ) {
      cat( "Loading KEGG annotations database ... this may take a while (it will be cached for subsequent calls).\n" )
      koa <- load.KEGG.annotations( org.id=organism )
      assign( "genome.info$koa", koa )
    }
    koa <- get( "genome.info$koa" )
    if ( is.numeric( k[ 1 ] ) ) k <- get.rows( k )
    syns <- unlist( get.synonyms( k ) )
    unique( koa[ koa$gene %in% syns, ] )
  }
  
  out <- data.frame()
  ko.total <- NA
  for ( k in ks ) {
    if ( is.numeric( k[ 1 ] ) ) rows <- get.rows( k )
    else { rows <- ks; k <- "" }
    tab <- cluster.KEGG.annotations( rows )
    ##print(tab)
    if ( is.na( ko.total ) ) ko.total <- length( unique( `genome.info$koa`$gene ) )
    syns <- unique( c( rows, unlist( get.synonyms( rows ) ) ) )
    syns <- syns[ syns %in% `genome.info$koa`$gene ]; names( syns ) <- NULL
    tt <- table( as.character( tab$path ) )
    tab <- subset( tab, path %in% names( tt[ tt > 1 ] ) )
    for ( i in names( tt[ tt > 1 ] ) ) { 
      ttt <- subset( tab, path == i )
      gns <- paste( ttt$gene, collapse=" ", sep=" " )
      cat.total <- nrow( subset( `genome.info$koa`, path == i ) )
      ##print(subset( `genome.info$koa`, path == i ) )
      ##cat(i,gns,cat.total,"\n")
      out <- rbind( out, data.frame( k=k, genes=gns, n.in=length( syns ), n.hit=nrow( ttt ), n.tot=ko.total,
                                    n.hit.tot=cat.total, p.value=phyper( nrow( ttt ), cat.total,
                                                           ko.total - cat.total, length( syns ), lower=F ),
                                    KEGG.path=ttt$path.name[ 1 ], KEGG.path.id=i ) )
    }
    if ( all( is.character( ks ) ) ) break
  }
  if ( nrow( out ) > 1 ) return( out[ order( out$p.value ), ] )
  out
}

bicluster.pvalues <- function( ks=1:k.clust, scores=r.scores ) {
  mc <- get.parallel( length( ks ) )
  do.call( c, mc$apply( ks, function( k ) {
    rows <- get.rows( k )
    sc.in <- scores[ rows, k ]
    sc.out <- scores[ ! rownames( scores ) %in% rows, k ]
    t.test( sc.out, sc.in, alt='g' )$p.value
  } ) )
}

cluster.score.pvalues <- function( ks=1:k.clust, score=c("resid","network")[1], rats=ratios[[ 1 ]],
                                  nets=names( net.weights ), n.samp=1000, out.resids=F, verbose=F ) {  
  residual.submatrix <- function( rats ) {
    d.rows <- rowMeans( rats, na.rm=T )
    d.cols <- colMeans( rats, na.rm=T )
    d.all <- mean( d.rows, na.rm=T )

    rij <- rats + d.all
    rij[,] <- rij[,] - matrix( d.cols, nrow=nrow( rij ), ncol=ncol( rij ), byrow=T )
    rij[,] <- rij[,] - matrix( d.rows, nrow=nrow( rij ), ncol=ncol( rij ), byrow=F )
    ##rij[,] <- rij[,] - outer( d.rows, d.cols, '+' ) ## more elegant but slower!
    mean( abs( rij ), na.rm=TRUE )
  }

  subnetwork.score <- function( nets, rows ) {
    n.scores <- sapply( nets, function( n ) {
      net <- networks[[ n ]]
      cons <- net[ as.character( net$protein1 ) %in% rows, c( "protein2", "combined_score" ), drop=F ]
      if ( is.null( cons ) || nrow( cons ) <= 0 ) return( 0 )
      cons <- cons[ as.character( cons$protein2 ) %in% rows, , drop=F ]
      if ( is.null( cons ) || nrow( cons ) <= 0 ) return( 0 )
      sum( cons$combined_score )
    } )
    -weighted.mean( n.scores, net.weights[ nets ], na.rm=T )
  }

  tmp <- t( sapply( clusterStack[ ks ], function( k ) c( k$nrows, k$ncols ) ) )
  out <- rep( NA, max( ks ) )

  mc <- get.parallel( n.samp )
  
  for ( k in ks ) {
    if ( ! is.na( out[ k ] ) ) next
    rows <- get.rows( k ); nrows <- length( rows )
    cols <- get.cols( k ); ncols <- length( cols )
    if ( score == "resid" ) {
      resids.sample <- unlist( mc$apply( 1:n.samp, function( i )
                                        residual.submatrix( rats[ sample( 1:nrow( rats ), nrows ),
                                                                 sample( 1:ncol( rats ), ncols ) ] ) ) )
                                        ##mean( e$r.scores[ sample( 1:nrow( rats ), nrows ), k ], na.rm=T ) ) )
      dupes <- which( is.na( out ) & tmp[ ,1 ] == nrows & tmp[ ,2 ] %in% ( (-3):3 + ncols ) )
    } else if ( score == "network" ) {
      resids.sample <- unlist( mc$apply( 1:n.samp, function( i )
                                        subnetwork.score( nets, sample( rownames( rats ), nrows ) ) ) )
      dupes <- which( is.na( out ) & tmp[ ,1 ] == nrows )
    }
    for ( kk in dupes[ dupes %in% ks ] ) {
      rows <- get.rows( kk ); cols <- get.cols( kk )
      if ( score == "resid" ) resid <- residual.submatrix( rats[ rows, cols ] )
                                       ## mean( e$r.scores[ rows, kk ], na.rm=T )
      else if ( score == "network" ) resid <- subnetwork.score( nets, rows )
      out[ kk ] <- mean( resids.sample <= resid )
      if ( verbose ) cat( kk, nrows, ncols, resid, out[ kk ], "\n" )
    }
  }
  ##if ( out.resids ) { attr( resids.sample, "resid" ) <- resid; return( resids.sample ) }
  out[ ks ]
}

## My own code in case I had a function such as this:
## history <- list(); cm.func.each.iter <- function() history[[ as.character( iter ) ]] <<- row.membership
process.history <- function( e, history ) {
  out1 <- out2 <- NULL
  row.memb <- NULL
  for ( i in 1:length( history ) ) {
    old.row.memb <- row.memb 
    row.memb <- t( apply( history[[ i ]], 1, function( i ) 1:e$k.clust %in% i ) )
    if ( is.null( old.row.memb ) ) next
    changed <- sum( row.memb != old.row.memb, na.rm=T )
    ##cat( i, changed, "\n" )
    out1 <- rbind( out1, c( i, changed ) )
  }
  par( mfrow=c( 2, 1 ) ); plot( out1 )
  row.memb <- NULL
  for ( i in seq( 1, length( history ), by=10 ) ) {
    old.row.memb <- row.memb 
    row.memb <- t( apply( history[[ i ]], 1, function( i ) 1:e$k.clust %in% i ) )
    if ( is.null( old.row.memb ) ) next
    changed <- sum( row.memb != old.row.memb, na.rm=T )
    ##cat( i, changed, "\n" )
    out2 <- rbind( out2, c( i, changed ) )
  }
  plot( out2 )
  for(i in 1:length(history))print(c(i,e$get.rows(1,history[[1]])%in%e$get.rows(1,history[[i]])))
}

## Get orthologs from e.g. http://rsat.ulb.ac.be/rsat/data/genomes/Halobacterium_sp/blast_hits/ ???
## Using 2 environments - based on cmonkey.init() for 2 different (or same?) organisms
## e.g. q_Halobacterium_sp_db_Pyrococcus_furiosus_ranks.tab.gz
## another option is to get the fasta files via e.g.
##  http://rsat.ulb.ac.be/rsat/data/genomes/Halobacterium_sp/genome/Halobacterium_sp_aa.fasta
## and run blast to get best reciprocal hits (better b/c some org's aren't in the RSAT orthologs tabs)
get.best.blast.reciprocal.hits <- function( env1, env2, do.blast=T, e.cutoff=1e-5, fraction.cutoff=0.5,
                                           filter=T ) {
  rsat.url <- "http://rsat.ulb.ac.be/rsat/data/genomes/"
  org1 <- env1$genome.info$species
  org2 <- env2$genome.info$species

  if ( ! do.blast ) { ## Get the table from rsat?
    fname <- paste( "data/", org1, "/q_", org1, "_db_", org2, "_ranks.tab.gz", sep="" )
    url <- paste( rsat.url, "/Halobacterium_sp/blast_hits/q_", org1, "_db_", org2, "_ranks.tab.gz", sep="" )
    cMonkey:::dlf( fname, url )
    tab <- read.delim( gzfile( fname ) )
    tab <- subset( tab, q_rank == 1 & s_rank == 1 )
    tab <- subset( tab, ! duplicated( query ) & ! duplicated( subject ) )
    return( tab )
  }

  fname1 <- paste( "data/", org1, "/", org1, "_aa.fasta", sep="" )
  cMonkey:::dlf( fname1, paste( rsat.url, org1, "/genome/", org1, "_aa.fasta", sep="" ) )
  fname2 <- paste( "data/", org2, "/", org2, "_aa.fasta", sep="" )
  cMonkey:::dlf( fname2, paste( rsat.url, org2, "/genome/", org2, "_aa.fasta", sep="" ) )
  
##   cm.attach( "env1" )
##   seqs1 <- get.sequences( "all", seq.opt="refseq" )
##   cm.detach( "env1" )
##   cm.attach( "env2" )
##   seqs2 <- get.sequences( "all", seq.opt="refseq" )
##   cm.detach( "env2" )

  ## See recommendation from Moreno-Hagelsieb and Latimer "Choosing BLAST options for better detection of
  ##    orthologs as reciprocal best hits" Bioinformatics 24:319-324, 2008 for best BLASTP options. 
  ##bl2seq.cmd <- "./progs/bl2seq"
  formatdb.cmd <- "./progs/formatdb -p T -o F -i %s -l %s"
  
  tmp.log.file <- tempfile( "formatdb.log" )
  cmd <- sprintf( formatdb.cmd, fname1, tmp.log.file )
  print( cmd ); system( cmd )

  blast.cmd <- "./progs/blastall -d %s -i %s -p blastp -F \"m S\" -m 8 -b 200000" ##-s T 
  cmd <- sprintf( blast.cmd, fname1, fname2 )
  print( cmd ); blast.out <- system( cmd, intern=TRUE, ignore=TRUE )

  if ( substr( blast.out[ 1 ], 1, 12 ) == "# BLASTN 2.2" ) blast.out <- blast.out[ -(1:3) ]
  out <- t( sapply( strsplit( blast.out, "\t" ), cbind ) )
  ##out <- do.call( rbind, strsplit( blast,out, "\t" ) )
  ##colnames( out ) <- c( "Query id", "Subject id", "% identity", "alignment length", "mismatches", "gap openings",
  ##                     "q. start", "q. end", "s. start", "s. end", "e-value", "bit score" )
  out <- data.frame( `Query id`=out[ ,1 ], `Subject id`=out[ ,2 ], `% identity`=as.numeric( out[ ,3 ] ),
                    `alignment length`=as.integer( out[ ,4 ] ), mismatches=as.integer( out[ ,5 ] ),
                    `gap openings`=as.integer( out[ ,6 ] ), `q. start`=as.integer( out[ ,7 ] ),
                    `q. end`=as.integer( out[ ,8 ] ), `s. start`=as.integer( out[ ,9 ] ),
                    `s. end`=as.integer( out[ ,10 ] ), `e-value`=as.numeric( out[ ,11 ] ),
                    `bit score`=as.numeric( out[ ,12 ] ) )

  read.fasta <- function( fname ) {
    lines <- readLines( gzfile( fname ) )
    lines <- lines[ lines != "" ]
    starts <- grep( "^>", lines, perl=T )
    stops <- c( starts[ 2:length( starts ) ], length( lines ) + 1 )
    seqs <- sapply( 1:length( starts ), function( i ) paste( lines[ ( starts[ i ] + 1 ):( stops[ i ] - 1 ) ],
                                                            collapse="", sep="" ) )
    names( seqs ) <- gsub( "^>", "", lines[ starts ], perl=T )
    seqs
  }

  if ( ! filter ) return( out )
  ## Filter - e-value cutoff
  out <- subset( out, e.value <= e.cutoff )
  
  ## get rid of alignments that are < 50% of either of the 2 sequences
  lens1 <- nchar( read.fasta( fname1 ) )
  lens2 <- nchar( read.fasta( fname2 ) )
  out <- cbind( out, fraction.x=out$alignment.length / lens1[ as.character( out$Subject.id ) ],
               fraction.y=out$alignment.length / lens2[ as.character( out$Query.id ) ] )
  out <- subset( out, fraction.x >= fraction.cutoff & fraction.y >= fraction.cutoff )
  
  ## Note - Query.id comes from env2 species; Subject.id comes from env1 species
  names1 <- unique( subset( env1$genome.info$feature.names, names %in% attr( env1$ratios, "rnames" ) |
                           names %in% env1$genome.info$transl.table )[ , c( "id", "names" ) ] )
  names2 <- unique( subset( env2$genome.info$feature.names, names %in% attr( env2$ratios, "rnames" ) |
                           names %in% env2$genome.info$transl.table )[ , c( "id", "names" ) ] )
  out <- merge( out, names2, by.x='Query.id', by.y='id', all.x=T, all.y=F )
  out <- merge( out, names1, by.x='Subject.id', by.y='id', all.x=T, all.y=F )
  
  ## Sort it so first one (which is kept) is the best
  out <- out[ order( -out$bit.score, out$e.value ), ]
  ## Look for reciprocal best hits
  ##out <- out[ ! duplicated( out[ ,1:2 ] ), ] ## First get rid of dupes (keeping best hit b/c it is sorted)
  out <- subset( out, ! duplicated( Subject.id ) & ! duplicated( Query.id ) )
  out <- subset( out, ! is.na( names.y ) & ! is.na( names.x ) )
  out
}

## ## Map orfs of organism in env1 to genome sequence in env2
## get.orf.mappings.to.genome <- function( env1, env2 ) {
##   rsat.url <- "http://rsat.ulb.ac.be/rsat/data/genomes/"
##   org1 <- env1$genome.info$species
##   org2 <- env2$genome.info$species

## ##   fname1 <- paste( "data/", org1, "/", org1, "_aa.fasta", sep="" )
## ##   dlf( fname1, paste( rsat.url, org1, "/genome/", org1, "_aa.fasta", sep="" ) )
## ##   fname2 <- paste( "data/", org2, "/", org2, "_aa.fasta", sep="" )
## ##   dlf( fname2, paste( rsat.url, org2, "/genome/", org2, "_aa.fasta", sep="" ) )

##   cm.attach( "env1" )
##   orf.seqs1 <- get.sequences( rownames( env1$ratios ), seq.opt="gene" )
##   cm.detach( "env1" )
##   orf.seqs1.fname <- tempfile( "orf1.fst." )
##   cat( paste( ">", names( orf.seqs1 ), "\n", orf.seqs1, sep="" ), file=orf.seqs1.fname, sep="\n" )

##   chr.seqs <- env2$genome.info$genome.seqs
##   chr.fname <- tempfile( "chr.fst." )
##   cat( paste( ">", names( chr.seqs ), "\n", chr.seqs, sep="" ), file=chr.fname, sep="\n" )

##   formatdb.cmd <- "./progs/formatdb -p T -o F -i %s -l %s"
  
##   tmp.log.file <- tempfile( "formatdb.log" )
##   cmd <- sprintf( formatdb.cmd, chr.fname, tmp.log.file )
##   print( cmd ); system( cmd )

##   blast.cmd <- "./progs/blastall -d %s -i %s -p blastn -n T -W 56 -m 8 -b 200000" ##-s T 
##   cmd <- sprintf( blast.cmd, chr.fname, orf.seqs1.fname )
##   print( cmd ); blast.out <- system( cmd, intern=TRUE, ignore=TRUE )
## }

## TODO: get orthologs from OrthoMCL database (groups or reciprocal best hits?)
## Using 2 environments - based on cmonkey.init() for 2 different (or same?) organisms
## THIS FUNCTION IS NOT FINISHED
## get.orthologs <- function( env1, env2, e.cutoff=1e-5, prefix1="VNG", prefix2="SSO" ) {
##   ##if ( ! file.exists( "data/orthomcl" ) ) dir.create( "data/orthomcl/" )
##   dlf( "data/orthomcl/seqs_id_map_orthomcl-2.txt.gz",
##       "http://www.orthomcl.org/common/downloads/2/seqs_id_map_orthomcl-2.txt.gz",
##       "Downloading OrthoMCL database files to local cache..." )
##   dlf( "data/orthomcl/groups_orthomcl-2.txt.gz",
##       "http://www.orthomcl.org/common/downloads/2/groups_orthomcl-2.txt.gz",
##       "Downloading OrthoMCL database files to local cache..." )
##   dlf( "data/orthomcl/recip_best_hits_orthomcl-2.txt.gz",
##       "http://www.orthomcl.org/common/downloads/2/recip_best_hits_orthomcl-2.txt.gz",
##       "Downloading OrthoMCL database files to local cache..." )

##   seq.ids <- read.delim( gzfile( "data/orthomcl/seqs_id_map_orthomcl-2.txt.gz" ), head=F ) ## This file is BIG
##   seq.ids <- cbind( seq.ids, refseq=sapply( strsplit( as.character( seq.ids$V2 ), "|", fixed=T ), "[", 2 ) )
##   seq.ids1 <- subset( seq.ids, refseq %in% env1$genome.info$feature.names$names )
##   seq.ids2 <- subset( seq.ids, refseq %in% env2$genome.info$feature.names$names )
##   ##seq.ids1 <- seq.ids[ grep( paste( "^", env1$cmonkey.params$organism, "\\|", sep="" ), seq.ids$V2 ), ]
##   ##seq.ids2 <- seq.ids[ grep( paste( "^", env2$cmonkey.params$organism, "\\|", sep="" ), seq.ids$V2 ), ]
##   rm( seq.ids ); gc()

##   ##seq.ids1 <- cbind( seq.ids1, refseq=sapply( strsplit( as.character( seq.ids1$V2 ), "|", fixed=T ), "[", 2 ) )
##   ##seq.ids2 <- cbind( seq.ids2, refseq=sapply( strsplit( as.character( seq.ids2$V2 ), "|", fixed=T ), "[", 2 ) )

## ##   if ( ! is.null( genome.info$transl.table ) ) 
## ##     names <- rbind( names, unique( subset( genome.info$feature.names, names %in% genome.info$transl.table )[ ,c( "id", "names" ) ] ) )
## ##   merged <- merge( names, seq.ids, by.x="id", by.y="refseq" )[ ,c( "V1", "names" ) ]
## ##   if ( ! is.null( genome.info$transl.table ) ) {
## ##     tmp <- merged$names[ ! merged$names %in% attr( ratios, "rnames" ) & merged$names %in% genome.info$transl.table ]
## ##     tmp <- genome.info$transl.table[ genome.info$transl.table %in% tmp ];
## ##     tmp2 <- names( tmp ); names( tmp2 ) <- tmp; tmp <- tmp2; rm( tmp2 )
## ##     merged$names[ ! merged$names %in% attr( ratios, "rnames" ) & merged$names %in% genome.info$transl.table ] <-
## ##       tmp[ merged$names[ ! merged$names %in% attr( ratios, "rnames" ) & merged$names %in% genome.info$transl.table ] ]
## ##   }
  
##   best.hits <- read.delim( gzfile( "data/orthomcl/recip_best_hits_orthomcl-2.txt.gz" ), head=F ) ## This file is VERY BIG
##   best.hits <- subset( best.hits, V1 %in% c( as.character( seq.ids1$V1 ), as.character( seq.ids2$V1 ) ) &
##                       V2 %in% c( as.character( seq.ids1$V1 ), as.character( seq.ids2$V1 ) ) ); gc()
##   if ( ! is.na( e.cutoff ) ) best.hits <- subset( best.hits, V4 <= e.cutoff & V4 <= e.cutoff )
##   colnames( best.hits ) <- paste( "X", 1:ncol( best.hits ), sep="" )
##   best.hits <- merge( best.hits, seq.ids1, by.x="X1", by.y="V1", all.x=T )
##   best.hits <- merge( best.hits, seq.ids1, by.x="X2", by.y="V1", all.x=T )
##   best.hits <- merge( best.hits, seq.ids2, by.x="X1", by.y="V1", all.x=T )
##   best.hits <- merge( best.hits, seq.ids2, by.x="X2", by.y="V1", all.x=T )
##   ## Note this table contains org1->org1 and org2->org2 homologs too -- the next 2 lines gets rid of those.
##   best.hits <- best.hits[ ,c( "X4", "X5", "V2.y", "V2.x.1" ) ]
##   best.hits <- subset( best.hits, ! is.na( V2.y ) & ! is.na( V2.x.1 ) )
## ##   best.hits <- merge( best.hits, merged, by.x="V2", by.y="V1", all=T )

##   best.hits <- cbind( best.hits[ ,c( "X4", "X5" ) ],
##                      org1=sapply( strsplit( as.character( best.hits$V2.y ), "|", fixed=T ), "[", 1 ),
##                      ref1=sapply( strsplit( as.character( best.hits$V2.y ), "|", fixed=T ), "[", 2 ),
##                      org2=sapply( strsplit( as.character( best.hits$V2.x.1 ), "|", fixed=T ), "[", 1 ),
##                      ref2=sapply( strsplit( as.character( best.hits$V2.x.1 ), "|", fixed=T ), "[", 2 ) )

## ##   n1 <- unique( subset( env1$genome.info$feature.names, names %in% env1$networks$string$protein1 |
## ##                names %in% env1$networks$string$protein2 )[ ,1:2 ] )
## ##   n2 <- unique( subset( env2$genome.info$feature.names, names %in% env2$networks$string$protein1 |
## ##                names %in% env2$networks$string$protein2 )[ ,1:2 ] )
## ##   best.hits <- merge( best.hits, n1, by.x="ref1", by.y="id", all.x=T )
## ##   best.hits <- merge( best.hits, n1, by.x="ref2", by.y="id", all.x=T )
## ##   best.hits <- merge( best.hits, n2, by.x="ref1", by.y="id", all.x=T )
## ##   best.hits <- merge( best.hits, n2, by.x="ref2", by.y="id", all.x=T )

##   envs <- list( env1, env2 )
##   env.orgs <- sapply( envs, function( i ) i$cmonkey.params$organism )
##   string <- list( unique( c( as.character( env1$networks$string$protein1 ),
##                             as.character( env1$networks$string$protein2 ) ) ),
##                  unique( c( as.character( env2$networks$string$protein1 ),
##                            as.character( env2$networks$string$protein2 ) ) ) )
##   syns.tab <- NULL
##   for ( i in 1:nrow( best.hits ) ) {
##     e1 <- which( env.orgs == best.hits$org1[ i ] )
##     syns1 <- get.synonyms( as.character( best.hits$ref1[ i ] ), ft=envs[[ e1 ]]$genome.info$feature.names )
##     syns1 <- syns1[ syns1 %in% rownames( envs[[ e1 ]]$ratios ) | syns1 %in% string[[ e1 ]] ]
##     if ( length( syns1 ) > 1 ) syns1 <- syns1[ which.max( nchar( syns1 ) ) ]
##     else if ( length( syns1 ) <= 0 ) {
##       syns1 <- get.synonyms( as.character( best.hits$ref1[ i ] ), ft=envs[[ e1 ]]$genome.info$feature.names )
##       if ( ! is.na( prefix1 ) ) syns1 <- grep( paste( "^", prefix1, sep="" ), syns1, perl=T, val=T )
##       else syns1 <- NA
##     }
##     if ( length( syns1 ) <= 0 ) syns1 <- NA
    
##     e2 <- which( env.orgs == best.hits$org2[ i ] )
##     syns2 <- get.synonyms( as.character( best.hits$ref2[ i ] ), ft=envs[[ e2 ]]$genome.info$feature.names )
##     syns2 <- syns2[ syns2 %in% rownames( envs[[ e2 ]]$ratios ) | syns2 %in% string[[ e2 ]] ]
##     if ( length( syns2 ) > 1 ) syns2 <- syns2[ which.max( nchar( syns2 ) ) ]
##     else if ( length( syns2 ) <= 0 ) {
##       syns2 <- get.synonyms( as.character( best.hits$ref2[ i ] ), ft=envs[[ e2 ]]$genome.info$feature.names )
##       if ( ! is.na( prefix2 ) ) syns2 <- grep( paste( "^", prefix2, sep="" ), syns2, perl=T, val=T )
##       else syns2 <- NA
##     }
##     if ( length( syns2 ) <= 0 ) syns2 <- NA

##     syns.tab <- rbind( syns.tab, data.frame( syn1=syns1, syn2=syns2 ) )
##   }
##   best.hits <- cbind( best.hits, syns.tab )
## }

## Get the total fraction of iterations that each gene is in the same cluster as every other gene
## Must have had "save.history=TRUE"; will use this history.
## Can use this as a 1-distance matrix for hierarchical clustering?
cluster.cloud <- function( history ) {
  r.mat <- matrix( 0, nrow=attr( ratios, "nrow" ), ncol=attr( ratios, "nrow" ) )
  rownames( r.mat ) <- colnames( r.mat ) <- attr( ratios, "rnames" )
  c.mat <- matrix( 0, nrow=attr( ratios, "ncol" ), ncol=attr( ratios, "ncol" ) )
  rownames( c.mat ) <- colnames( c.mat ) <- attr( ratios, "cnames" )
  ##for ( i in 1:length( history ) ) {
  require( multicore )
  out <- mclapply( 0:3, function( ii ) { ##length( history ), function( i ) {
    for ( i in ( 1:length( history ) )[ 1:length( history ) %% 4 == ii ] ) {
      if ( i %% 100 == ii ) cat( ii, i, Sys.getpid(), "\n" )
      rm <- history[[ i ]]$row.membership
      cm <- history[[ i ]]$col.membership
      for ( k in 1:k.clust ) {
        rows <- unique( which( rm == k, arr=T )[ ,1 ] )
        r.mat[ rows, rows ] <- r.mat[ rows, rows ] + 1
        cols <- unique( which( cm == k, arr=T )[ ,1 ] )
        c.mat[ cols, cols ] <- c.mat[ cols, cols ] + 1
      }
    }
    list( r=r.mat, c=c.mat )
  } )
  for ( i in 1:length( out ) ) {
    r.mat <- r.mat + out[[ i ]]$r
    c.mat <- c.mat + out[[ i ]]$c
  }
  diag( r.mat ) <- diag( c.mat ) <- NA
  list( r=r.mat, c=c.mat )
}

## if ( FALSE ) {
##   qqq <- cluster.cloud( history )
##   dists <- max( qqq$r, na.rm=T ) - qqq$r
##   hc <- hclust( as.dist( dists ) )
##   memb <- cutree( hc, k=50 )
## }

## if ( FALSE ) {
## for(k in 1:k.clust){
##   qqq<-adjust.clust(k,motif=T,expand=length(get.rows(k))>50)
##   row.membership<-qqq$rm
##   if ( ! is.null( qqq$ms ) ) meme.scores[[k]]<-qqq$ms
##   plotClust(k)
## }
## }

## test <- function() {
##   rm <- t( apply( rr.scores[,], 1, order, decreasing=T )[ 1:n.clust.per.row, ,drop=F ] ) ##[ iter ]
##   rm <- t( apply( rm, 1, sort ) ); if ( n.clust.per.row == 1 ) rm <- t( rm ) ##[ iter ]

##   for ( k in 1:k.clust ) {
##     r <- get.rows( k )
##     r.out <- rownames( rr.scores )[ ! rownames( rr.scores ) %in% r ]
##     ll <- cluster.loglik( k )
##     for ( g in rownames( rr.scores ) ) {
##       if ( g %in% r ) {
##         l.in <- rr.scores[ r[ r != g ], k ]
##         l.out <- rr.scores[ c( r.out, g ), k ]
##       } else {
##         l.in <- rr.scores[ c( r, g ), k ]
##         l.out <- rr.scores[ r.out[ r.out != g ], k ]
##       }
##       ll.new <- sum( log( c( l.in, 1 - l.out ) + 1e-99 ), na.rm=T )
##       if ( ll.new < ll ) cat( g, g %in% r, ll.new - ll, "\n" )
##     }
##   }
## }

#endif
