###################################################################################
## cMonkey - version 4, Copyright (C) David J Reiss, Institute for Systems Biology
##                      dreiss@systemsbiology.org
## This software is provided AS IS with no warranty expressed or implied. Neither 
## the authors of this software nor the Institute for Systems Biology shall be held
## liable for anything that happens as a result of using this software
###################################################################################

#' Initialize a cMonkey environment for biclustering
#' NOTE: to run cMonkey as simple (slow) kmeans on ratios data, set:
#'   n.clust.per.row <- 1; n.clust.per.col <- k.clust
#'   no.genome.info <- TRUE; post.adjust <- FALSE; net.weights <- mot.weights <- numeric()
#'  
#' @param env  The input ratios matrix or the file name for a tab delinated file containging one
#' @param ...  Used to set many cMonkey parameters.  See "set.param" lines for details
#' @export
#' @usage e <- cmonkey.init( NULL, organism="sce" )
cmonkey.init <- function( env=NULL, ... ) {
  if ( ! exists( "cmonkey.params" ) ) {
    cmonkey.params <- new.env( hash=T ) ##, parent=globalenv() )
    ##if ( ! is.null( env ) && is.environment( env ) ) cmonkey.params <- new.env( hash=T, parent=env )
  }
    tmp.e <- environment( cMonkey:::cmonkey ) ## Packaged - get the env. that the "cmonkey" function is stored in
  
  if ( ! is.null( env ) && ( is.list( env ) || is.environment( env ) ) ) { ## if env is an input data list or env (e.g. "halo") copy its items here.
    ## Allow override of variables in "env" by ones that came in by commandline.
    for ( i in names( env ) ) if ( ! i %in% names( list( ... ) ) ) assign( i, env[[ i ]] )
    if ( is.list( env ) ) env <- NULL
  }

  for ( i in ls( tmp.e ) ) {
    f2 <- NULL
    if ( ( ! is.null( env ) && exists( i, envir=env, inherit=F ) ) ) { ## User-override of function (pre-declared), in env
      f <- try( get( i, envir=env ) ) ## User-defined override function
      f2 <- try( get( i, envir=tmp.e ) ) ## Original (overridden cmonkey package) function
    } else if ( exists( i, envir=.GlobalEnv, inherit=F ) ) { ## User-override of function (pre-declared), in global env
      f <- try( get( i, envir=.GlobalEnv ) ) ## User-defined override function
      f2 <- try( get( i, envir=tmp.e ) ) ## Original (overridden cmonkey package) function
    } else if ( exists( i ) ) { ## User-override of function (pre-declared), somewhere (e.g. in attached data)
      f <- try( get( i ) ) ## User-defined override function
      f2 <- try( get( i, envir=tmp.e ) ) ## Original (overridden cmonkey package) function
    } else {
      f <- try( get( i, envir=tmp.e ) ) ## Copy original func from cmonkey package env.
    }
    if ( class( f ) == "function" ) 
      environment( f ) <- sys.frames()[[ length( sys.frames() ) ]] ## Make each func's env be this func's frame
    assign( i, f ) ## Copy over all cmonkey functions & objs into this frame - eventually get copied into env
    if ( ! is.null( f2 ) && class( f2 ) == "function" && object.size( f2 ) != object.size( f ) ) {
      environment( f2 ) <- sys.frames()[[ length( sys.frames() ) ]] ## Make each func's env be this func's frame
      assign( paste( "super", i, sep="." ), f2 ) ## Copy over overridden cmonkey function w/ "super." at end of name
      if ( ! is.null( env ) ) {
        assign( paste( "super", i, sep="." ), f2, envir=env )
        environment( env[[ paste( "super", i, sep="." ) ]] ) <- env
      }
    }
  }
  rm( f, f2, tmp.e )

  if ( ! is.null( env ) ) for ( i in ls( env ) ) assign( i, get( i, env ) ) ## Copy over params/data pre-set in env

  args <- mget( names( formals() ), env=as.environment( -1 ) )
  for ( i in names( args ) ) if ( ! i %in% c( "...", "env" ) ) set.param( i, args[[ i ]] )
  for ( i in names( list( ... ) ) ) if ( i != "env" ) set.param( i, list( ... )[[ i ]] )
  rm( args )

  if ( sink.number() > 0 ) for ( i in 1:sink.number() ) try( sink(), silent=T ) ## Sometimes R gets into a funky sunky state
  set.param( "save.logfile", FALSE ) ## Divert output to a logfile as well as the terminal
  if ( save.logfile != FALSE ) sink( save.logfile, split=T,
         append=(exists("dont.init")&&dont.init)||(exists("is.inited")&&!is.inited) ) 

  if ( ! exists( "organism" ) ) {
    ##cat( "\33[31mWARNING: No organism was set; using \"hpy\".\33[0m\n" )
    message( "WARNING: No organism was set; using \"hpy\"." )
    organism <- "hpy"
    Sys.sleep( 3 )
  }
  set.param( "organism", organism )  

  ## If ratios doesnt exist but row.weights does, create a list w/ the objects named by names(row.weights)
  if ( ! exists( "ratios" ) && exists( "row.weights" ) ) {
    try( { ratios <- lapply( names( row.weights ), get ) ## default is length-1 vector with name "ratios" so this works.
    names( ratios ) <- names( row.weights ) } )
  }

  ## Ratios can be a single matrix/df/filename or a list of those; if a single matrix/df, make it a list w/that
  if ( ( exists( "ratios" ) && ! is.null( ratios ) ) ) {
    if ( is.matrix( ratios ) || is.data.frame( ratios ) ) ratios <- list( ratios=load.ratios( ratios ) )
    else if ( is.character( ratios ) ) ratios <- lapply( ratios, load.ratios )
    else ratios <- lapply( ratios, function( r ) as.matrix( load.ratios( r ) ) )
    ratios <- ratios[ sapply( ratios, function( r ) all( dim( r ) > 0 ) ) ]
    attr( ratios, "rnames" ) <- sort( unique( unlist( lapply( ratios, rownames ) ) ) )
    attr( ratios, "cnames" ) <- sort( unique( unlist( lapply( ratios, colnames ) ) ) )
    attr( ratios, "nrow" ) <- length( attr( ratios, "rnames" ) )
    attr( ratios, "ncol" ) <- length( attr( ratios, "cnames" ) ) ## Summary attributes (for multiple ratios matrices)
    if ( is.null( names( ratios ) ) ) {
      names( ratios ) <- paste( "ratios", 1:length( ratios ), sep='.' )
      if ( exists( 'row.weights' ) && ! all( names( ratios ) %in% names( row.weights ) ) )
        for ( i in names( ratios ) ) row.weights[ i ] <- row.weights[ 1 ]
    }
    for ( n in names( ratios ) ) {
      if ( ncol( ratios[[ n ]] ) > 1 ) {
        attr( ratios[[ n ]], "maxRowVar" ) <- mean( apply( ratios[[ n ]][,], 1, var, use="pair" ), na.rm=T ) ##* 1.2
        attr( ratios[[ n ]], "all.colVars" ) <- apply( ratios[[ n ]][,], 2, var, use="pair", na.rm=T )
      }
    }
    rm( n )
  }
  if ( exists( "ratios" ) && is.null( names( ratios ) ) )
    names( ratios ) <- paste( "ratios", 1:length( ratios ), sep='.' )
  if ( ! is.null( env ) && exists( "ratios" ) ) assign( "ratios", ratios, envir=env )

  ## if ( exists( "is.eukaryotic" ) && is.eukaryotic ) {
  ##   set.param( "operon.shift", FALSE )
  ##   set.param( "remove.low.complexity.subseqs", TRUE )
  ## }

  ## Default param settings (if already set, they will not be changed):
  set.param( "cog.org", "?" ) ## If "?" or NA, will try "Org" (for organism=="org")
  set.param( "rsat.species", "?" ) ## If "?" or NA, will try getting it from data/KEGG/KEGG_all_species.tab
  set.param( "n.iter", 2000 ) ##1000
  set.param( "n.clust.per.row", 2 ) ##n.clust.per.row <- 2
  if ( exists( "ratios" ) && ! is.null( ratios ) ) {
    set.param( "k.clust", round( attr( ratios, "nrow" ) * n.clust.per.row / 20 ) ) ## dflt avg clust size of 20
  } else {
    set.param( "k.clust", 100 )
  }
  set.param( "n.clust.per.col", if ( exists( "ratios" ) && attr( ratios, "ncol" ) >= 60 ) round( k.clust / 2 ) else round( k.clust * 2 / 3 ) ) ## dflt to 1/2 of conds in each clust (avg)
  set.param( "row.iters", seq( 1, n.iter, by=2 ) )
  set.param( "col.iters", seq( 1, n.iter, by=5 ) )
  ##set.param( "meme.iters", seq( 399, n.iter, by=100 ) ) ## Which iters to re-run meme?
  ##set.param( "meme.iters", c( seq( 600, 1200, by=100 ), seq( 1250, 1500, by=50 ), seq( 1525, 1800, by=25 ),
  ##                           seq( 1810, max( n.iter, 1820 ), by=10 ) ) )
  set.param( "meme.iters", seq( 100, n.iter, by=100 ) )
  ##set.param( "mot.iters", seq( 100, n.iter, by=10 ) ) ## Which iters to use results of most recent meme run in scores
  ##set.param( "mot.iters", seq( 601, max( n.iter, 605 ), by=3 ) ) ## Which iters to use results of most recent meme run in scores
  set.param( "mot.iters", seq( 100, n.iter, by=10 ) ) ## Which iters to use results of most recent meme run in scores
  set.param( "net.iters", seq( 1, n.iter, by=7 ) ) ## Which iters to re-calculate network scores?
  set.param( "row.scaling", 1 ) ##6 )  ## Seems to work best for Mpn, works good for Halo (6 is with quantile.normalize turned on)
  set.param( "row.weights", c( ratios=1 ) ) ## Optionally load multiple ratios files and set relative weights
  set.param( "mot.scaling", seq( 0, 1, length=n.iter*3/4 ) ) ##* 0.5
  set.param( "mot.weights", c( `upstream meme`=1 ) ) ##, `upstream weeder`=0.5, `upstream spacer`=1, `upstream memepal`=1 ) ) ## Sequence and algorithm for for motif search: Optionally use different automatically computed sequences (e.g. `downstream meme`=1) or an input file (e.g. `fstfile=zzz.fst meme`=1) (csvfile too!) and motif algos (e.g. weeder, spacer, prism, meme, memepal)
  set.param( "net.scaling", seq( 0, 0.5, length=n.iter*3/4 ) ) ##0.1 0.25
  ## Net weights and grouping weights - names must correspond to full file paths (sifs) that are to be read in.
  set.param( "net.weights", c( string=0.5, operons=0.5 ) ) ## prolinks=0.5 Relative scaling(s) of each network
  ## Can use pre-set nets: "operons"; "prolinks.(GN/GC/PP/RS)"; "predictome.(chromo/comp/fusion/phylogenetic)";
  ##           "string.(combined/neighborhood/fusion/cooccurence/coexpression/experimental/database/textmining)"
  set.param( "grouping.weights", numeric() ) ## Vector of weights of 3-column tsvs to read (p1, score, group)
  set.param( "plot.iters", seq( 2, n.iter, by=25 ) ) ## set to NA or 0 for no plotting
  set.param( "post.adjust", TRUE ) ## Post-processing of clusters to make sure that each gene that belongs in a given cluster but is not in it (because another cluster out-competed for it) is added
  set.param( "parallel.cores", TRUE ) ## 3 ## NA -> no parallelization; TRUE - detect # of cores and use all
  set.param( "parallel.cores.motif", TRUE ) ## can be different for motifing (uses less RAM so can use more cores?)
  
  ## Probably no need to change these:
  set.param( "max.changes", c( rows=0.5, cols=5 ) ) ##c( rows=0.2, cols=5 ) ) ## 0.2 => 1 in 5 chance of each row getting updated per iter
  set.param( "cluster.rows.allowed", c( 3, 70 ) ) ##200 ) ) ## Min/max number of rows to allow in a bicluster
  set.param( "merge.cutoffs", c( n=0.3, cor=0.975 ) ) ## n=0.3 => merge 1 pair of clusters every 3 iters; if n>1 then merge that number of pairs of clusters every iter; cor is correlation cutoff
  ## Note: to use seeded clusters use the "list=" row seeding method and set "fuzzy.index" to close to 0... seems to work (note setting it to 0 is probably a bad idea).
  set.param( "fuzzy.index", 0.75 * exp( -( 1:n.iter ) / (n.iter/4) ) ) ## (n.iter/6) ## hack to add stochasticity
  ##set.param( "fuzzy.index", 0.7 * exp( -( 1:n.iter ) / (n.iter/3) ) + 0.05 ) ## (n.iter/6) ## hack to add stochasticity
  set.param( "translation.tab", NULL ) ## custom 2-column translation table to be used for additional synonyms
  set.param( "seed.method", c( rows="kmeans", cols="best" ) ) ## "net=string:5" "rnd" "kmeans" "trimkmeans=TRIM" "rnd" "list=FILENAME" "rnd=NG" "cor=NG" "net=netname:NG" "netcor=netname:NG" "custom" -- NG is # of genes per seeded cluster; "best" or "rnd" is option for cols
  set.param( "maintain.seed", NULL ) ## List of lists of vectors of rows to maintain for each k: force seeded rows or cols in each cluster to STAY there! e.g. maintain.seed=list(rows=list(`3`=c(gene1,gene2,gene3))) ; This should be used in conjunection with seed.method["rows"]=="custom" or "list=..."
  ##set.param( "string.links.url", "http://string82.embl.de/newstring_download/protein.links.v8.2.txt.gz" ) ## Need to update this when they update their version number

  ##set.param( "n.motifs", c( rep( 1, n.iter/2 ), rep( 2, n.iter/4 ), 3 ) ) ##rep( 2, n.iter/3 ) ) ) ##, rep( 2, n.iter/4 + 50 ) ) )
  set.param( "n.motifs", c( rep( 1, n.iter/3 ), rep( 2, n.iter/3 ) ) ) ##rep( 2, n.iter/3 ) ) ) ##, rep( 2, n.iter/4 + 50 ) ) )
  ##set.param( "motif.width.range", c( 6, 24 ) ) ## Can be an iter-based param
  if ( file.exists( "./progs" ) && file.exists( "./progs/meme" ) ) {
    set.param( "progs.dir", "./progs/" )
  } else if ( "package:cMonkey" %in% search() && file.exists( sprintf( "%s/progs/",
                                                                      system.file( package="cMonkey" ) ) ) ) {
    set.param( "progs.dir", sprintf( "%s/progs/", system.file( package="cMonkey" ) ) )
  } else if ( any( mot.scaling > 0 ) && ( ! exists( "no.genome.info" ) || ! no.genome.info ) ) {
    message( "WARNING: You do not have meme/mast/dust installed in the correct location.\nTrying to install it now.\n" )
    install.binaries()
    set.param( "progs.dir", sprintf( "%s/progs/", system.file( package="cMonkey" ) ) )
    if ( "package:cMonkey" %in% search() && ! file.exists( sprintf( "%s/progs/", system.file( package="cMonkey" ) ) ) )
      message( "WARNING: Could not install meme. Please see the website for installation instructions." )
  }
  ##set.param( "meme.cmd", paste( progs.dir, "meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $compute -pal=non", sep="/" ) ) ## -nomatrim -prior addone -spfuzz 1 -pal" ) ## Can use -pal=non (or exlude the option), -pal=pal (or just -pal), or -pal=both
  ## Parameters for meme.cmd:
  ##   "-bfile $bgFname" -- if this is omitted, then no background file is submitted; the default meme bg. is used (i.e. derived from the input sequences)
  ##   "-psp $pspFname" -- if this is omitted, no position-specific prior is used (default)
  ##   "-cons $none" -- if this is changed to "-cons $compute" then the consensus from previous meme run on this cluster is used as seed for this meme run (if the previous motif had a good E-value)
  ##   "-pal=non" -- if this is changed to "-pal=pal" then force palindrome search; "-pal=both": try both pal and non-pal and use the result with the lowest E-value. THIS IS DEFUNCT... now use memepal and memeboth motifing option instead
  set.param( "meme.cmd", paste( progs.dir, "meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 24 -mod zoops -nostatus -text -cons $compute", sep="/" ) ) ##-allw -pal=non -cons $none
  set.param( "mast.cmd", sprintf( "%s/mast $memeOutFname -d $fname -bfile $bgFname -nostatus -stdout -text -brief -ev 99999 -mev 99999 -mt 0.99 -seqp -remcorr", progs.dir ) )
  set.param( "dust.cmd", sprintf( "%s/dust $fname", progs.dir ) )
  ##set.param( "meme.addl.args", "-time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw %2$d -maxw %3$d -mod zoops" ) ## -nomatrim -prior addone -spfuzz 1" )
  ##set.param( "mast.addl.args", "" ) ##-ev 99999 -mev 99999 -mt 0.99 -seqp -remcorr" ) ##"-ev 10 -mev 10 -mt 0.1
  ##set.param( "meme.consensus", "compute" ) ## "" to not use the best e-val consens from prev meme run on clust
  ##set.param( "meme.consensus.e.val.limit", 0.1 ) ## Use consensus of previous meme run only if E-value was < this
  ##set.param( "motif.palindrome.option", "non" ) ## "non", "pal", or "both" (try both and use one w/ best E-value)
  ##set.param( "uniquify.seqs", TRUE )
  ##set.param( "remove.low.complexity.subseqs", TRUE ) ##FALSE ) ##TRUE
  set.param( "operon.shift", TRUE )
  ##set.param( "meme.seqs.allowed", cluster.rows.allowed ) ## Min/max number of seqs to allow to be fed to meme
  set.param( "bg.order", 3 ) ##bg.order <- NA ##0 ##3   ## NA -> no global background; use the input sequences
  set.param( "recalc.bg", TRUE ) ## if recalc.bg==TRUE, recalc bg for each MEME run using only seqs for genes that are NOT in cluster. Ideally, would be TRUE always, but could be slow for big genomes.
  set.param( "motif.upstream.search", c( -20, 150 ) ) ##-50, 250 ) ##-50, 200 )
  set.param( "motif.upstream.scan", c( -30, 250 ) ) ##-30, 150 ) ##-50, 200 )
##  if ( any( mot.scaling > 0 ) && ( ! file.exists( meme.cmd ) || ! file.exists( mast.cmd ) ) )
##    stop( paste( "Motif finding is requested but", meme.cmd, "and/or", mast.cmd, "is not installed!" ) )

  set.param( "rsat.urls", c( "http://rsat.ccb.sickkids.ca/", "http://rsat.ulb.ac.be/rsat/", 
                            "http://embnet.ccg.unam.mx/rsa-tools" ) ) ## Right now only first one is used
  set.param( "stats.iters", c( 1, seq( 5, n.iter, by=5 ) ) )
  set.param( "cm.script.each.iter", "cm.script.each.iter.R" ) ## a vector of R script file names to source every iter!

  set.param( "date.run", format( Sys.time(), "%y %b %d %H:%M:%S" ) ) ##date()
  set.param( "cmonkey.version", cm.version ) ## New versions: 4m -> 4.1.3 (4l -> 4.1.2; 4a -> 4.0.1, etc)
  set.param( "session.info", unlist( list( R.version, Sys.info(), Sys.getenv(), sessionInfo() ) ), quiet=T ) ## overkill?
  set.param( "time.started", date() )
  if ( exists( "ratios" ) && ! is.null( ratios ) ) {
    set.param( "cmonkey.filename", paste( "cmonkey", cmonkey.version, organism,
                                         paste( sapply( ratios, dim ), collapse="x" ),
                                         gsub( " ", "_", date.run ), sep="_" ) )
  } else {
    set.param( "cmonkey.filename", paste( "cmonkey", cmonkey.version, organism, "0x0",
                                         gsub( " ", "_", date.run ), sep="_" ) )
  }
  ##set.param( "rnd.seed", as.integer( Sys.time() ) )
  ## Guaranteed unique up to 1e-6 seconds! (True?)
  op <- options( digits.secs=10 )
  set.param( "rnd.seed", as.integer( substr( gsub( '[-:. ]', "", as.character( Sys.time() ) ), 12, 20 ) ) )
  options( op ); rm( op )
  set.seed( rnd.seed )
  set.param( "big.memory", FALSE ) ##50 * 2^20 ) ## Matrices that are > 50 MB are stored as file-backed big.memory
  set.param( "big.memory.verbose", FALSE )

  if ( organism == "hsa" ) rsat.urls[ 1 ] <- rsat.urls[ 2 ] ## Special case for hsa - not hosted on mirrors.
  
  if ( ! exists( "rsat.species" ) || rsat.species == "?" || is.na( rsat.species ) ) {
    err <- dlf( "data/KEGG/KEGG_taxonomy.txt", 'http://baliga.systemsbiology.net/cmonkey/taxonomy.txt' ) ##"ftp://ftp.genome.jp/pub/kegg/genes/taxonomy" )
    if ( class( err ) != "try-error" ) {
      tab <- read.delim( "data/KEGG/KEGG_taxonomy.txt", sep='\t', comment='#', head=F, as.is=T )
      rsat.spec <- as.character( subset( tab, V2 == organism, select="V4", drop=T ) )[ 1 ]; rm( tab )
      if ( any( strsplit( rsat.spec, "" )[[ 1 ]] == "(" ) ) rsat.spec <- gsub( '\\s\\(.*\\)', "", rsat.spec )
    } else {
      rsat.spec <- '?'
    }
    rsat.spec <- gsub( " ", "_", rsat.spec, fixed=T )
    kegg.spec <- rsat.spec
    
    ## Check to see if species is in RSAT dir... if not, guess or ask.
    if ( ! file.exists( "data/RSAT_genomes_listing.txt" ) ) {
      require( RCurl )
      tmp <- strsplit( getURL( paste( rsat.urls[ 1 ], "/data/genomes/", sep="" ) ), "\n" )[[ 1 ]]
      writeLines( tmp, con="data/RSAT_genomes_listing.txt" ) ## Cache the listing -- it takes some time to dld
    } ##!else { ##if ( ! require( RCurl ) ) {
      ##stop( "Please install 'RCurl' package.\n" )
    ##}

    vals <- character()
    if ( file.exists( "data/RSAT_genomes_listing.txt" ) ) {
      tmp <- readLines( "data/RSAT_genomes_listing.txt" )
      vals <- grep( rsat.spec, tmp, fixed=T, val=T )
    }

    if ( ##! file.exists( "data/RSAT_genomes_listing.txt" ) ||
        length( vals ) <= 0 ) {
      message( "Could not find correct organism for RSAT... will try to guess..." )

      ## err <- dlf( "data/KEGG/KEGG_all_species.tab", "ftp://ftp.genome.jp/pub/kegg/genes/etc/all_species.tab" )
      ## if ( class( err ) != "try-error" ) {
      ##   tab <- read.delim( "data/KEGG/KEGG_all_species.tab", sep='\t', comment='#', head=F, as.is=T )
      ##   vals <- grep( rsat.spec, tab$V9, val=T )
      ##   if ( length( vals ) >= 1 ) rsat.spec <- vals
      ## } else {
      
      max.dist <- 0.5; vals <- rep( "", 3 )
      while( length( vals ) > 1 ) {
        vals <- agrep( rsat.spec, tmp, ignore=T, max.dist=max.dist, val=T ) ## tailored for Halo...
        max.dist <- max.dist - 0.01
        if ( length( vals ) <= 0 ) {
          max.dist <- max.dist + 0.02
          vals <- agrep( rsat.spec, tmp, ignore=T, max.dist=max.dist, val=T )
          break
        }
      }
      if ( length( vals ) > 1 ) {
        rsat.spec <- sapply( strsplit( vals, "[<>/]" ), "[", 8 )
        message( "Found ", length( rsat.spec ), " matches..." )
        rsat.spec <- rsat.spec[ menu( rsat.spec, graphics=F, title="Please choose one." ) ]
      }
      if ( length( vals ) == 1 ) {
        rsat.spec <- strsplit( vals, "[<>/]" )[[ 1 ]][ 8 ]
        message( "Found one match: ", rsat.spec, " ..." )
        message( "If this is not correct, you're not quite out of luck -- set the 'rsat.species' parameter manually." )
      }
      ##}
    }
    set.param( "rsat.species", rsat.spec, override=T )
    ##dlf( paste( "data/STRING/species.", string.version, ".txt", sep="" ),
    ##    paste( "http://string.embl.de/newstring_download/species.", string.version, ".txt", sep="" ) )
    ## url <- string.links.url
    ## fname <- strsplit( url, "/" )[[ 1 ]]; fname <- sprintf( "data/STRING/%s", fname[ length( fname ) ] )
    ## dlf( gsub( ".gz", "", gsub( "protein.links", "species", fname ) ),
    ##     gsub( ".gz", "", gsub( "protein.links", "species", url ) ) )
    rm( tmp, rsat.spec, err, vals ) ##, url, fname )
  } else {
    set.param( "rsat.species", rsat.species )
  }

  if ( ! exists( "taxon.id" ) || taxon.id == "?" || is.na( taxon.id ) || length( taxon.id ) <= 0 ) {
    fname <- dlf( "data/GO/proteome2taxid", "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid" )
    tab <- read.delim( gzfile( "data/GO/proteome2taxid" ), head=F )
    taxon.id <- subset( tab, V1 == gsub( "_", " ", rsat.species ) )$V2
    if ( length( taxon.id ) <= 0 ) taxon.id <- subset( tab, grepl( gsub( "_", " ", rsat.species ), V1 ) )$V2[ 1 ]
    set.param( "taxon.id", taxon.id, override=T )
    rm( tab, fname )
  }

  if ( ! exists( "cog.org" ) || cog.org == '?' || is.na( cog.org ) ) {
    tmp <- strsplit( organism, "" )[[ 1 ]]; tmp[ 1 ] <- toupper( tmp[ 1 ] )
    cog.o <- paste( tmp, sep='', collapse='' )
    if ( cog.o == '' ) cog.o <- '?'
    set.param( 'cog.org', cog.o, override=T )
    rm( cog.o, tmp )
  } else {
    set.param( "cog.org", cog.org )
  }

  message( "Organism is ", organism, " ", cog.org, " ", rsat.species, " ", taxon.id )

  genome.loc <- paste( rsat.urls[ 1 ], "/data/genomes/", rsat.species, "/genome/", sep="" )
  fname <- paste( "data/", rsat.species, "/organism.tab", sep="" )
  err <- dlf( fname, paste( genome.loc, "/organism.tab", sep="" ) )
  org.tab <- readLines( fname )
  org.tab <- strsplit( org.tab[ length( org.tab ) ], "\t" )[[ 1 ]]
  is.eukaryotic <- any( grepl( "Eukaryota", org.tab ) )
  cat( "Is eukaryote:", is.eukaryotic, "\n" )
  rm( err, org.tab, genome.loc, fname )
  if ( is.eukaryotic ) {
    message( "Organism is a eukaryote; presuming there are no operons." )
    set.param( "is.eukaryotic", TRUE, override=T )
    set.param( "operon.shift", FALSE, override=T )
    set.param( "discard.genome", TRUE, override=T ) ## big genomes and no operon shifting, so we can discard it.
    set.param( "recalc.bg", FALSE, override=T ) ## Takes too much RAM to do this once per motif run
    ##set.param( "remove.low.complexity.subseqs", TRUE, override=T )
    if ( "operons" %in% names( net.weights ) ) {
      net.weights <- net.weights[ names( net.weights ) != "operons" ]
      set.param( "net.weights", net.weights, override=T )
    }
  }
  
  if ( get.parallel( 100, verbose=T )$mc ) on.exit( try( kill( children(), SIGKILL ) ), add=T ) ## Make sure no zombies are lying around
  on.exit( { if ( sink.number() > 0 ) for ( i in 1:sink.number() ) try( sink(), silent=T ) } ) ## Sometimes R gets into a funky sunky state

  ## Make sure the total "weights" sum to one - the "scaling"s set the total scaling
  if ( sum( net.weights, na.rm=T ) > 0 ) net.weights <- net.weights / sum( net.weights, na.rm=T )
  if ( sum( row.weights, na.rm=T ) > 0 ) row.weights <- row.weights / sum( row.weights, na.rm=T )
  if ( sum( mot.weights, na.rm=T ) > 0 ) mot.weights <- mot.weights / sum( mot.weights, na.rm=T )

  ## Update all motif-related parameters so they're seq.type specific.
  ## If it's a >1 length vector, make it a list of vectors; else make it a vector
  for ( i in c( "n.motifs", "meme.cmd", "mast.cmd", ##"meme.consensus", ##"meme.addl.args", "mast.addl.args", 
               ##"meme.consensus.e.val.limit", "motif.palindrome.option", "motif.width.range", 
               ##"uniquify.seqs",
               ##"remove.low.complexity.subseqs","meme.seqs.allowed",
               "meme.iters", 
               "operon.shift", "bg.order", 
               "motif.upstream.search", "motif.upstream.scan" ) ) {
    v <- get( i )
    if ( all( names( mot.weights ) %in% names( v ) ) ) next
    if ( is.vector( v ) && length( v ) > 1 ) v <- list( `1`=v ) ##`upstream meme`=v )
    ##!else names( v )[ 1 ] <- "upstream meme"
    names( v ) <- names( mot.weights )[ 1 ] 
    for ( n in names( mot.weights )[ ! names( mot.weights ) %in% names( v ) ] ) {
      if ( is.list( v ) ) v[[ n ]] <- v[[ 1 ]]
      else if ( is.vector( v ) ) v[ n ] <- v[ 1 ]
      names( v )[ length( v ) ] <- names( mot.weights )[ length( v ) ]
    }
    assign( i, v )
  }
  rm( v )

  if ( ! is.null( env ) ) for ( i in ls() ) {
    if ( i %in% c( "i", "env" ) ) next
    tmp <- get( i )
    if ( class( tmp ) == "function" ) environment( tmp ) <- env
    assign( i, tmp, envir=env )
  }
  
  ## INIT
  if ( ! is.na( rsat.species ) && ( ! exists( "genome.info" ) || genome.info$species != rsat.species ) ) {
    cat( "Initializing genome info for organism", organism, "\n" )
    
    set.param( "no.genome.info", FALSE )
    genome.info <- get.genome.info()
    if ( ! is.null( env ) ) assign( "genome.info", genome.info, envir=env )
    
    if ( is.na( taxon.id ) || length( taxon.id ) <= 0 ) {
      taxon.id <- genome.info$taxon.id
      set.param( "taxon.id", taxon.id, override=T )
      message( "Organism is ", organism, " ", cog.org, " ", rsat.species, " ", taxon.id )
    }

    ## Get common prefix from feature.names and use those genes (assume >40% of ORF names have this suffix)
    if ( exists( 'ratios' ) && ! is.null( ratios ) ) tmp <- toupper( attr( ratios, "rnames" ) )
    else if ( exists( 'genome.info' ) && ! is.null( genome.info$feature.names ) ) {
      tmp <- toupper( subset( genome.info$feature.names, type == "primary", select="names", drop=T ) )
      if ( exists( 'ratios' ) && ! is.null( ratios ) )
        tmp <- tmp[ toupper( tmp ) %in% toupper( attr( ratios, "rnames" ) ) ]
    }
    qqq <- sapply( 1:4, function( nch ) max( table( substr( tmp, 1, nch ) ) ) / length( tmp ) ); nch <- 0
    if ( any( qqq > 0.9 ) ) { nch <- which( qqq > 0.9 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >60% of names
    else if ( any( qqq > 0.6 ) ) { nch <- which( qqq > 0.6 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >60% of names
    else if ( any( qqq > 0.4 ) ) { nch <- which( qqq > 0.4 ); nch <- nch[ length( nch ) ] } ## Longest prefix in >40% of names
    ##nch <- 0; while( max( table( substr( tmp, 1, nch + 1 ) ) ) / length( tmp ) > 0.4 ) nch <- nch + 1
    prefix <- NA
    if ( nch > 0 ) {
      prefix <- names( which.max( table( substr( tmp, 1, nch ) ) ) )
      message( "Assuming gene/probe names have common prefix '", prefix, "'." )
      genome.info$gene.prefix <- prefix
    } else {
      message( "Could not find a common gene/probe identifier prefix. This only matters if there's no expression matrix." )
      prefix <- genome.info$gene.prefix <- NA
    }

    if ( TRUE ) {
      tmp2 <- tmp
      if ( length( unique( nchar( tmp2 ) ) ) > 1 ) {
        nc <- max( nchar( tmp2 ) )
        for ( i in 1:length( tmp2 ) ) {
          tmp2[ i ] <- paste( tmp2[ i ], paste( rep( ' ', nc - nchar( tmp2[ i ] ) ), sep='', collapse='' ), sep='' )
        }
      }
      tmp2 <- do.call( rbind, strsplit( tmp2, '' ) )
      regex <- apply( tmp2, 2, function( i ) sort( unique( i ) ) )
      for ( i in 1:length( regex ) ) {
        ii <- as.integer( regex[[ i ]] )
        if ( ! any( is.na( ii ) ) ) {
          if ( length( ii ) == length( ii[ 1 ]:ii[ length( ii ) ] ) && all( sort( ii ) == ii[1]:ii[length(i)] ) )
            regex[[ i ]] <- paste( '[', paste( ii[ 1 ], ii[ length( ii ) ], sep='-' ), ']', sep='' )
        }
        if ( length( regex[[ i ]][ regex[[ i ]] != ' ' ] ) > 1 ) regex[[ i ]] <- c( '[', regex[[ i ]], ']' )
        if ( any( regex[[ i ]] == '' | regex[[ i ]] == ' ' | is.na( regex[[ i ]] ) ) )
          regex[[ i ]] <- c( regex[[ i ]][ regex[[ i ]] != ' ' ], '?' )
      }
      regex <- paste( unlist( lapply( regex, paste, sep='', collapse='' ) ), sep='', collapse='' )
      regex <- gsub( '.', '\\.', regex, fixed=T )
      message( "Assuming gene/probe names have regex '", regex, "'." )
      genome.info$gene.regex <- regex
    }    

    genome.info$all.gene.names <- unique( as.character( subset( genome.info$feature.names,
                                                        grepl( paste( "^", genome.info$gene.regex, sep="" ), names,
                                                                     ignore=T, perl=T ), select="names", drop=T ) ) )
    if ( length( genome.info$all.gene.names ) ) { ## regex is still in testing phase!
      genome.info$all.gene.names <- unique( as.character( subset( genome.info$feature.names,
                                                        grepl( paste( "^", genome.info$gene.prefix, sep="" ), names,
                                                                     ignore=T, perl=T ), select="names", drop=T ) ) )
    }

    if ( ! is.null( env ) ) assign( "genome.info", genome.info, envir=env )

    ## Get operon predictions; do this before getting sequences so we can op-shift if desired
    genome.info$operons <- NULL
    if ( ( operon.shift || "operons" %in% names( net.weights ) ) && ! no.genome.info ) {
      tmp.operons <- try( get.operon.predictions( "microbes.online" ) )
      if ( class( tmp.operons ) == "try-error" ) {
        message( "Could not fetch operons file. Assuming it doesn't exist (eukaryote?)" )
        set.param( "is.eukaryotic", TRUE, override=T )
        set.param( "operon.shift", FALSE, override=T ); operon.shift[ 1:length( operon.shift ) ] <- FALSE
        ##set.param( "remove.low.complexity.subseqs", TRUE, override=T )
        if ( "operons" %in% names( net.weights ) ) {
          net.weights <- net.weights[ names( net.weights ) != "operons" ]
          set.param( "net.weights", net.weights, override=T )
        }
      } else {
        genome.info$operons <- tmp.operons
      }
      rm( tmp.operons )
      if ( ! is.null( env ) ) assign( "genome.info", genome.info, envir=env )
    }
    
    if ( ! exists( 'ratios' ) || is.null( ratios ) ) {
      message( "WARNING: No ratios matrix -- will generate an 'empty' one with all annotated ORFs for 'probes'." )
      if ( ! is.null( genome.info$gene.regex ) ) rows <- unique( as.character( subset( genome.info$feature.names,
                                                                    grepl( paste( "^", genome.info$gene.regex, sep="" ), names,
                                                                    ignore=T, perl=T ), select="names", drop=T ) ) )
      else rows <- unique( as.character( subset( genome.info$feature.names, type=="primary", select="names", drop=T ) ) )

      ratios <- list( ratios=t( t( rep( NA, length( rows ) ) ) ) ); rownames( ratios$ratios ) <- rows
      attr( ratios, "rnames" ) <- sort( unique( rows ) ); rm( rows )
      attr( ratios, "nrow" ) <- length( attr( ratios, "rnames" ) )
      attr( ratios, "ncol" ) <- 1
      cat( "Ratios: ", attr( ratios, "nrow" ), "x", 1, "\n" )
    }
    rm( nch, prefix, regex, tmp, qqq )

    ## Get upstream/downstream seqs, compute bg model, optionally discard genome seqs (for memory)
    if ( ! no.genome.info && length( mot.weights ) > 0 ) { ##is.null( genome.info$genome.seqs ) ) {
      genome.info$all.upstream.seqs <- genome.info$bg.list <- list()
      genome.info$bg.fname <- character()

      ##if ( ! all( is.na( bg.order ) ) ) {
      for ( i in names( mot.weights ) ) {
        cat( "Pre-computing all '", i, "' seqs (", paste( motif.upstream.scan[[ i ]], collapse=", " ), ")...\n", sep="" )
        ## Note we don't filter all seqs (used for background) - even removing ATGs; is this okay?
        genome.info$all.upstream.seqs[[ i ]] <- get.sequences( ##attr( ratios, "rnames" ),
                                                              genome.info$all.gene.names, seq.type=i,
                                                              distance=motif.upstream.scan[[ i ]], filter=F )
        if ( ! is.null( env ) ) assign( "genome.info", genome.info, envir=env )
        message( sum( ! attr( ratios, "rnames" ) %in% names( genome.info$all.upstream.seqs[[ i ]] ) ),
                " probes have no '", i, "' sequence." )
        if ( ! is.na( bg.order[ i ] ) ) {
          cat( "Pre-computing '", i, "' residue bg distrib (order=", bg.order[ i ], ")...\n", sep="" )
          tmp.seqs <- if ( ! is.null( genome.info$all.upstream.seqs[[ i ]] ) ) genome.info$all.upstream.seqs[[ i ]]
          else get.sequences( ##attr( ratios, "rnames" ),
                             genome.info$all.gene.names, distance=motif.upstream.search[[ i ]], seq.type=i, filter=F )
          genome.info$bg.fname[ i ] <- my.tempfile( "meme.tmp", suf=".bg" ) 
          capture.output(
                         genome.info$bg.list[[ i ]] <- mkBgFile( tmp.seqs, order=bg.order[ i ], verbose=T,
                                                                 bgfname=genome.info$bg.fname[ i ],
                                                                 use.rev.comp=grepl( "-revcomp", meme.cmd[ i ] ) ) )
          rm( tmp.seqs )
        } else {
          message( "NOT USING a global sequence background distribution!" )
        }
        if ( ! is.null( env ) ) assign( "genome.info", genome.info, envir=env )
      }
      ##}
##!ifndef 
##      if ( big.memory ) genome.info$all.upstream.seqs <-
##        list.reference( genome.info$all.upstream.seqs, sprintf( "%s/all.genome.seqs", cmonkey.filename ) )
##!endif
    }
    
    networks <- list()
    if ( ! is.na( net.iters ) && any( net.iters %in% 1:n.iter ) ) {
      ## Note vague info on STRING scores at http://string.embl.de/newstring_cgi/show_info_page.pl ..
      ## low confidence - 20% (or better), medium confidence - 50%, high confidence - 75%, highest confidence - 95%
      ##if ( file.exists( 'data/STRING/string.csv' ) ) { ## HACK WARNING (Chris P.)
      ##  networks[[ "string" ]] <- read.csv( 'data/STRING/string.csv', row.names=1, header=TRUE )
      ##} else
      if ( length( grep( "string", names( net.weights ) ) ) > 0 ) {
        if ( "string" %in% names( net.weights ) ) { ##|| "string.combined" %in% names( net.weights ) ) {
          ##if ( "string.combined" %in% names( net.weights ) )
          ##  names( net.weights )[ names( net.weights ) == "string.combined" ] <- "string"
          if ( exists( "string.links" ) ) { ## Let a user pre-load it as a 3-column data frame
            string <- string.links
          } else {
            ##if ( exists( "get.STRING.links.NEW" ) ) string <- get.STRING.links.NEW( genome.info$org.id$V1[ 1 ] )
            ##!else
            cat( "Loading STRING network.\n" )
            string <- get.STRING.links( genome.info$org.id$V1[ 1 ] ) ##, detailed=F )
            ##string <- subset( string, combined_score >= 500 )
            ##cat( "Read in", nrow( string ), "STRING edges that pass cutoff (500); weight =", net.weights[ "string" ], "\n" )
          }
          if ( ! is.null( string ) && nrow( string ) > 0 ) {
            cat( "Read in", nrow( string ), "STRING edges; weight =", net.weights[ "string" ], "\n" )
            string$combined_score <- string$combined_score / max( string$combined_score, na.rm=T ) * 1000
            string$combined_score <- 1000 * exp( string$combined_score / 1000 ) / exp( 1 )
            networks[[ "string" ]] <- string
          } else {
            warning( "Could not load STRING network. Either", organism, "is not there or your gene names are not standard." )
          }
          rm( string )
        }
        ## if ( length( grep( "string.", names( net.weights ) ) ) > 0 ) {
        ##   string <- get.STRING.links( genome.info$org.id$V1[ 1 ], detailed=T )
        ##   for ( n in grep( "string.", names( net.weights ), val=T ) ) {
        ##     tp <- strsplit( n, ".", fixed=T )[[ 1 ]][ 2 ]
        ##     if ( tp %in% colnames( string ) ) {
        ##       str <- string[ ,c( "protein1", "protein2", tp ) ]
        ##       colnames( str )[ 3 ] <- "combined_score"
        ##       str <- subset( str, combined_score >= 500 )
        ##       cat( "Read in", nrow( str ), n, "edges that pass cutoff (500); weight =", net.weights[ n ], "\n" )
        ##       str$combined_score <- 1000 * exp( string$combined_score / 1000 ) / exp( 1 )
        ##       networks[[ n ]] <- str
        ##     }
        ##   }
        ##   rm( string, tp, str )
        ## }
      }
      if ( ! is.null( env ) ) assign( "networks", networks, envir=env )

      if ( "operons" %in% names( net.weights ) && ! is.null( genome.info$operons ) ) {
        cat( "Converting operon predictions into a network...\n" )
        tmp <- tapply( genome.info$operons$gene, genome.info$operons$head ); names( tmp ) <- genome.info$operons$gene
        mc <- get.parallel( length( unique( tmp ) ) )
        out.sif <- do.call( rbind, mc$apply( unique( tmp ), function( j ) {
          whch <- which( tmp == j )
          gs <- names( whch )
          if ( length( gs ) <= 1 || length( gs ) > attr( ratios, "nrow" ) / 20 ) return( NULL )
          tmp.sif <- t( combn( gs, 2 ) )
          tmp.sif <- tmp.sif[ tmp.sif[ ,1 ] != tmp.sif[ ,2 ], ,drop=F ]
          data.frame( protein1=tmp.sif[ ,1 ], protein2=tmp.sif[ ,2 ], combined_score=rep( 1000, nrow( tmp.sif ) ) )
        } ) )
        ##out.sif <- do.call( rbind, out.sif )
        if ( ! is.null( out.sif ) && nrow( out.sif ) > 0 ) {
          out.sif$combined_score <- rep( 1000, nrow( out.sif ) )
          colnames( out.sif ) <- c( "protein1", "protein2", "combined_score" )
          networks[[ "operons" ]] <- out.sif
        }
        rm( tmp, mc, out.sif )
      }
      if ( ! is.null( env ) ) assign( "networks", networks, envir=env )
      
      ## Read in prolinks interactions from the web (including operon prediction edges)
      
      ## Read in additional networks from sifs (3-column file: p1, s, p2). If s is a string
      ##   (e.g. "pp") then edge weight is assumed to be 1; otherwise s can be numeric and
      ##   an edge weight
      ## First check to see if network is already existing in memory as a 2 or 3-column matrix
      if ( exists( "net.weights" ) && length( net.weights ) > 0 && ! is.null( names( net.weights ) ) ) {
        for ( i in names( net.weights ) ) {
          if ( i %in% names( networks ) ) next ## already done (e.g. operons)
          if ( file.exists( i ) ) { 
            cat( "Loading sif interactions from file:", i, "; weight =", net.weights[ i ], "\n" ) 
            sif <- load.sif.interactions( i ) 
          } else if ( exists( i ) && ! is.null( ncol( get( i ) ) ) && ncol( get( i ) ) >= 2 ) { ## Exists as name of network variable in memory
            cat( "Using network '", i, "' that exists in memory already; weight = ", net.weights[ i ], "\n", sep="" )
            sif <- get( i ) ## Must have column names 'protein1', 'protein2', 'combined_score'
            if ( ncol( sif ) == 2 ) sif <- cbind( sif, rep( 1, nrow( sif ) ) ) ## Add a weights column
            colnames( sif ) <- c( "protein1", "protein2", "combined_score" )
          } else {
            next ## "operons", "string", "prolinks", "predictome" are dealt with elsewhere
          }
          networks[[ basename( i ) ]] <- sif; rm( sif )
        }
      }
      if ( ! is.null( env ) ) assign( "networks", networks, envir=env )

      ## Read in additional gene "groupings" from sif (3-column file: g, s, p). If s is a string
      ##   (e.g. "pd") then grouping weight is assumed to be 1; otherwise s can be numeric and
      ##   a grouping weight. "g" denotes a grouping, e.g. what TF binds to its upstream region,
      ##   or a COG group or GO function, etc.
      ## Groups are converted into a set of (weighted) "interactions" between pairs of genes in the same group
      if ( exists( "grouping.weights" ) && length( grouping.weights ) > 0 ) {
        if ( exists( "net.weights" ) ) net.weights <- c( net.weights, grouping.weights )
        else net.weights <- grouping.weights
        for ( i in names( grouping.weights ) ) {
          if ( i %in% names( networks ) ) next ## already done (e.g. operons)
          if ( file.exists( i ) ) {
            cat( "Loading groupings from file:", i, "; weight =", grouping.weights[ i ], "\n" ) 
            sif <- load.sif.interactions( i ) 
          } else if ( exists( i ) && ! is.null( ncol( get( i ) ) ) && ncol( get( i ) ) >= 2 ) { ## Exists as name of network variable in memory
          ##} else {
            cat( "Using groupings from '", i, "' that exists in memory already; weight = ", grouping.weights[ i ], "\n", sep="" )
            sif <- get( i )
            if ( ncol( sif ) == 2 ) sif <- cbind( sif, combined_score=rep( 1, nrow( sif ) ) ) ## Add a weights column
          }
          ## We will assume that column with fewer unique names is the "groups" column.
          colnames( sif ) <- c( "group", "protein", "combined_score" )
          if ( sum( unique( as.character( sif$protein ) ) %in% attr( ratios, "rnames" ) ) <
              sum( unique( as.character( sif$group ) ) %in% attr( ratios, "rnames" ) ) ) {
            sif <- sif[ ,c( 2, 1, 3 ) ]
            colnames( sif ) <- c( "group", "protein", "combined_score" )
          }
          sif <- sif[ order( sif$group ), ]
          tmp <- tapply( sif$protein, sif$group ); names( tmp ) <- as.character( sif$protein )
          cat( "Converting", length( unique( tmp ) ), "groupings to a network (this may take a while for big grouping files)..." )
          mc <- get.parallel( length( unique( tmp ) ) )
          out.sif <- mc$apply( unique( tmp ), function( j ) { 
            whch <- which( tmp == j )
            gs <- names( whch )
            if ( length( gs ) <= 1 || length( gs ) > attr( ratios, "nrow" ) / 20 ) return( NULL )
            ws <- sif$combined_score[ whch ]; names( ws ) <- gs
            tmp.sif <- t( combn( gs, 2 ) ) 
            tmp.sif <- tmp.sif[ tmp.sif[ ,1 ] != tmp.sif[ ,2 ], ,drop=F ] 
            tmp.sif <- data.frame( protein1=tmp.sif[ ,1 ], protein2=tmp.sif[ ,2 ], combined_score=( ws[ tmp.sif[ ,1 ] ] + ws[ tmp.sif[ ,2 ] ] ) / 2 )
            rownames( tmp.sif ) <- NULL
            if ( j %% 100 == 0 ) cat( j ); cat( "." )
            tmp.sif
          } ) 
          cat( length( unique( tmp ) ), "... " )
          out.sif <- do.call( rbind, out.sif )
          colnames( out.sif ) <- c( "protein1", "protein2", "combined_score" )
          networks[[ basename( i ) ]] <- out.sif ##basename( grouping.files[ i ] ) ]] <- out.sif
          cat( "DONE\n" )
        }
        rm( sif, tmp, out.sif, i, mc )
      }
      if ( ! is.null( env ) ) assign( "networks", networks, envir=env )

      ## Post-processing step for all networks to make them kosher
      ##if ( exists( "networks" ) ) {
      ## Need to make sure names of genes in the networks are mapped correctly to probe names
      ## NOTE: dont need to use genome.info$transl.table for Halo for the STRING entries (for some reason)...
      ## Convert protein names to integers (corresponding to row number in ratios) for speed... NOTE that
      ##   I also throw out edges for genes that are not in the ratios
      for ( n in names( networks ) ) {
        nn <- networks[[ n ]]
        if ( nrow( nn ) <= 0 ) {
          message( "WARNING: no edges in network", n, "... skipping." );
          if ( length( grep( n, seed.method[ 1 ] ) ) > 0 ) {
            message( "ALSO, we have to change the row seeding method from", seed.method, "to 'kmeans'." )
            seed.method[ "rows" ] <- "kmeans"
            set.param( "seed.method", seed.method, override=T )
          }
          next
        }
        nodes <- unique( c( as.character( nn$protein1 ), as.character( nn$protein2 ) ) )
        cat( nrow( nn ), "edges,", length( nodes ), "nodes in network", n, "\n" )
        
        ## Remove self-interactions because we don't care about them here
        nn <- subset( nn, as.character( protein1 ) != as.character( protein2 ) ) 
        
        ## Merge duplicate edges (summing up their weights)  (will die miserably if too many nodes so...)
        ## Actually, we probably dont need to do this for the clustering (see get.network.scores() --
        ##    but it will make the plotting look bad, so it might be worth it (even though it's slow)
        dupes <- duplicated( nn[ ,c( "protein1", "protein2" ) ] )
        if ( sum( dupes ) > 0 ) { ## This creates an nxn table so for really big networks this may be a problem...
          cat( "Merging", sum( dupes ), "duplicate edges in network", n,
              "; this could take a while for networks with lots of nodes...\n" )
          tmp.nn <- subset( nn, dupes )
          dupe.nodes <- unique( c( as.character( tmp.nn$protein1 ), as.character( tmp.nn$protein2 ) ) )
          if ( length( dupe.nodes ) < 6000 ) { ## arbitrary size
            tmp <- tapply( tmp.nn$combined_score, tmp.nn[ ,c( "protein1", "protein2" ) ], sum, na.rm=T )
            tmp2 <- which( ! is.na( tmp ), arr=T )
            nn.new <- data.frame( protein1=rownames( tmp )[ tmp2[ ,1 ] ], protein2=colnames( tmp )[ tmp2[ ,2 ] ],
                                 combined_score=tmp[ tmp2 ] )
            rm( tmp, tmp2 )
            nn <- rbind( nn.new, nn )
            rm( nn.new )
            nn <- nn[ ! duplicated( nn[ ,c( "protein1", "protein2" ) ] ), ]
          }
          rm( tmp.nn, dupe.nodes )
        }
        
        if ( exists( "ratios" ) && ! is.null( ratios ) && ! any( nodes %in% attr( ratios, "rnames" ) ) ) {
          ## ath in STRING has genes like AT1G09180.1 whereas I only know about AT1G09180 so lets double check this
          if ( median( nchar( nodes ) ) > median( nchar( attr( ratios, "rnames" ) ) ) &&
              any( substr( nodes, 1, median( nchar( attr( ratios, "rnames" ) ) ) ) %in% attr( ratios, "rnames" ) ) ) {
            nn$protein1 <- substr( as.character( nn$protein1 ), 1, median( nchar( attr( ratios, "rnames" ) ) ) )
            nn$protein2 <- substr( as.character( nn$protein2 ), 1, median( nchar( attr( ratios, "rnames" ) ) ) )
            nodes <- unique( c( as.character( nn$protein1 ), as.character( nn$protein2 ) ) )
          }

          if ( ! is.null( genome.info$synonyms ) ) { ## Reconcile node names against synonyms for those not in ratios
            rr <- attr( ratios, "rnames" )[ ! attr( ratios, "rnames" ) %in% nodes ]
            if ( length( rr ) > 0 ) {
              cat( "Reconciling network", n, length( rr ), "node names with probe names...\n" )
              syns <- get.synonyms( rr )
              mc <- get.parallel( length( syns ) )
              is.there <- unlist( mc$apply( syns, function( i ) any( i %in% nodes ) ) )
              
              syns <- syns[ is.there ]
              nnc1 <- as.character( nn$protein1 ); nnc2 <- as.character( nn$protein2 )
              nnc1.t <- ! nnc1 %in% attr( ratios, "rnames" )
              nnc2.t <- ! nnc2 %in% attr( ratios, "rnames" )
              ##                 for ( i in names( syns ) ) {
              ##                   nnc1[ nnc1.t & nnc1 %in% syns[[ i ]] ] <- i
              ##                   nnc2[ nnc2.t & nnc2 %in% syns[[ i ]] ] <- i
              ##                 }
              mc <- get.parallel( 2 ) ## Can we parallelize this better?
              tmp <- mc$apply( 1:2, function( ii ) {
                for ( i in names( syns ) ) {
                  if ( ii == 1 ) nnc1[ nnc1.t & nnc1 %in% syns[[ i ]] ] <- i
                  else nnc2[ nnc2.t & nnc2 %in% syns[[ i ]] ] <- i
                }
                if ( ii == 1 ) return( nnc1 ) else return( nnc2 )
              } )
              nnc1 <- tmp[[ 1 ]]; nnc2 <- tmp[[ 2 ]]; rm( tmp, nnc1.t, nnc2.t )
              cat( sum( ! is.there ), "probes have no nodes in", n, "network (but",
                  sum( attr( ratios, "rnames" ) %in% nodes, na.rm=T ) + sum( is.there ), "do)\n" )
              nn$protein1 <- nnc1; nn$protein2 <- nnc2
              tmp <- nnc1 %in% attr( ratios, "rnames" ) & nnc2 %in% attr( ratios, "rnames" )
              nn <- subset( nn, tmp == TRUE )
              rm( tmp, syns, is.there, nnc1, nnc2, nnc1.t, nnc2.t, tmp, rr, i )
            }
          }
        } else {
          cat( sum( ! attr( ratios, "rnames" ) %in% nodes ), "probes have no nodes in", n, "network (but",
              sum( attr( ratios, "rnames" ) %in% nodes, na.rm=T ), "do)\n" )
        }

        ## symmetrize the network!
        ttmp <- nn[ ,c( 2, 1, 3 ) ]; colnames( ttmp ) <- colnames( nn )
        nn <- rbind( nn, ttmp ); rm( ttmp )
        
        nn <- nn[ ! duplicated( nn[ ,c( "protein1", "protein2" ) ] ), ]
        cat( n, "network filtered, symmetrized and uniquified:", nrow( nn ), "edges.\n" )
        networks[[ n ]] <- nn
        if ( ! is.null( env ) ) assign( "networks", networks, envir=env )
      }
      rm( n, nn, nodes, dupes )
      
      ## Need to make sure sum of edge weights over all networks is the same so that bigger networks
      ##    dont get undue influence just because they're big.
      if ( length( networks ) > 1 ) {
        sums <- sapply( networks, function( n ) sum( n$combined_score, na.rm=T ) )
        ms <- max( sums[ sums > 0 ], na.rm=T )
        if ( length( sums ) > 0 && ! is.na( ms ) ) for ( n in names( networks ) )
          if ( sums[ n ] > 0 ) networks[[ n ]]$combined_score <- networks[[ n ]]$combined_score / sums[ n ] * ms
        rm( n, sums, ms )
      }
      if ( ! is.null( env ) ) assign( "networks", networks, envir=env )
      ##}
      
      if ( ! is.null( names( net.weights ) ) ) names( net.weights ) <- basename( names( net.weights ) )
      ##set.param( "net.weights", net.weights, override=T )
      ##names( networks ) <- basename( names( networks ) )
      ##if ( exists( "net.weights", envir=.GlobalEnv ) ) {
      ##  assign( "net.weights.old", get( "net.weights", envir=.GlobalEnv ), envir=.GlobalEnv ) ## If it was set here
      ##  rm( net.weights, envir=.GlobalEnv ) ##         originally, it will override this updated version, so rename it.
      ##}
    }    

    
    if ( ! no.genome.info && cog.org != '' && cog.org != '?' && ! is.null( cog.org ) && all(is.na(plot.iters) | plot.iters==0) ) {
      ## COG code from NCBI whog file
      cat( "Loading COG functional codes (for plotting), org. code", cog.org, ": trying NCBI whog file...\n" )
      genome.info$cog.code <- get.COG.code( cog.org ) ##, genome.info$feature.names ) ##genome.info$transl.table,
      ##     if ( is.null( genome.info$cog.code ) || all( is.na( genome.info$cog.code ) ) ) {
      ##       if ( file.exists( "cmonkey-experimental.R" ) ) {
      ##         ## COG code from KEGG KO (usually fewer genes assigned)
      ##         ## Experimental:
      ##         source( "cmonkey-experimental.R", local=T )
      ##         cat( "Loading COG functional codes (for plotting): using KEGG orthology files...\n" )
      ## WHERE DID THIS FUNCTION GO???
      ##         genome.info$cog.code <- get.COG.code.from.KEGG( organism ) ##, genome.info$feature.names ) ##genome.info$transl.table,
      ##       }
      ##     }
      cat( sum( ! is.na( genome.info$cog.code ) ), "genes have a COG code (",
          if ( is.null( genome.info$cog.code ) ) attr( ratios, "nrow" ) else sum( is.na( genome.info$cog.code ) ),
          "do not)\n" )
    }

  }

  iter <- 0
  meme.scores <- clusterStack <- list()
  for ( i in names( mot.weights ) ) {
    meme.scores[[ i ]] <- list()
    meme.scores[[ i ]][[ k.clust + 1 ]] <- ""
  }
  stats <- row.scores <- col.scores <- mot.scores <- net.scores <- NULL

  if ( ! exists( "favorite.cluster" ) )
    favorite.cluster <- function() min( which( sapply( 1:k.clust, function( k ) length( get.rows( k ) ) ) >
                                              cluster.rows.allowed[ 1 ] * 2 ) )

  ##if ( n.iter == 3000 ) n.iter <- 2000 ## 2000 is just as good; but this way maintain the schedules
  row.scaling <- extend.vec( row.scaling )
  mot.scaling <- extend.vec( mot.scaling )
  net.scaling <- extend.vec( net.scaling ) ##, envir=cmonkey.env )
  n.motifs <- lapply( n.motifs, extend.vec ) ##, envir=cmonkey.env )
  ##min.motif.width <- lapply( motif.width.range, function( i ) extend.vec( i[ 1 ] ) ) ##, envir=cmonkey.env )
  ##max.motif.width <- lapply( motif.width.range, function( i ) extend.vec( i[ 2 ] ) ) ##, envir=cmonkey.env )
  ##n.clust.per.row <- extend.vec( n.clust.per.row ) ##, envir=cmonkey.env )
  ##n.clust.per.col <- extend.vec( n.clust.per.col ) ##, envir=cmonkey.env )
  fuzzy.index <- extend.vec( fuzzy.index ) ##, envir=cmonkey.env )

  is.inited <- TRUE
  if ( is.null( env ) ) env <- new.env( hash=T, parent=globalenv() )
  else parent.env( env ) <- globalenv()
  parent.env( cmonkey.params ) <- env
  attr( env, "class" ) <- c( "environment", "cmonkey" ) ## For future use??
  for ( i in ls() ) {
    if ( i %in% c( "i", "env" ) ) ##, "DATE", "VERSION", "cm.version", "time.started", "cmonkey.init", 
##                   "cmonkey.re.init", "cog.org", "date.run", "dlf", "extend.vec", "get.COG.code",
##                   "get.STRING.links", "session.info", "load.ratios", "load.sif.interactions", 
##                   "get.condition.groups", "get.genome.info", "get.operon.predictions", "get.predictome.links",
##                   "get.prolinks.links", "getMastPValuesAndEValues", "getMemeMotifInfo", "getMemeMotifPssm",
##                   "pssm.to.string", "rnd.seed", "rev.comp", "viewPssm", "mkBgFile", "get.mast.pvals",
##                   "remove.low.complexity", "residual.submatrix", "runMast", "runMeme", "system.time.limit", 
##                   "rsat.species", "rsat.urls", "save.logfile", "seed.clusters", "seed.method",
##                   "set.param", "string.version", "taxon.id", "update.cmonkey.env" ) )
      next
    tmp <- get( i )
    if ( class( tmp ) == "function" ) environment( tmp ) <- env
    assign( i, tmp, envir=env )
  }

  if ( exists( "favorite.cluster" ) ) ##{
    env$favorite.cluster <- favorite.cluster; environment( env$favorite.cluster ) <- env
    ##print( env$favorite.cluster ) }
  if ( exists( "cm.func.each.iter" ) ) {
    env$cm.func.each.iter <- cm.func.each.iter; environment( env$cm.func.each.iter ) <- env
    ##print( env$cm.func.each.iter )
    try( env$cm.func.each.iter(), silent=T ) ## User-defined func. to run each iteration (here, at end of initialization)
  }

  
  cat( "INITIALIZATION IS COMPLETE.\n" )
  env$iter <- env$iter + 1
  invisible( env )
}
