#ifndef PACKAGE

source( "cmonkey.R" )

cm.package <- function( data=F, bigdata=F, install=T, update.web=F, check=F, version=cm.version, beta=F ) {
  ## Can get halo ratios and translation.tab and other relevant data via data(halo);attach(halo)
  source.files <- c( list.files( patt="^cmonkey.*\\.R$" ), 'cmonkey.R', 'cmonkey-optim1.R' )
  source.files <- source.files[ ! source.files %in% c( "cmonkey.code.update.R", "cmonkey-experimental.R",
                                                      "cmonkey-motif-other.R", "cmonkey-utils.R",
                                                      "cmonkey-ensemble.R", "cmonkey-ensemble-funcs.R",
                                                      "cmonkey-ensemble-funcs2.R", "cmonkey-optim2.R",
                                                      "cmonkey-data-load-mo.R" ) ]
  print( source.files )
  
  ##cm.name <- "cMonkey"
  if ( beta ) {
    ##cm.name <- paste( cm.name, "beta", sep="." ) ## Include experimental functions
    ##version <- paste( version, "99", sep="." ) ## Include experimental functions
    update.web <- FALSE ## DO NOT PUT ON THE WEB! bigdata <- 
  }

  onLoad <- function( libname, pkgname ) { ##.onAttach
    packageStartupMessage( "Loading ", pkgname, " version ", VERSION, " (", DATE, ")" )
    packageStartupMessage( "Copyright (C) David J Reiss, Institute for Systems Biology; dreiss@systemsbiology.org." )
    packageStartupMessage( "https://github.com/dreiss-isb/cmonkey" )
    if ( grepl( "beta", VERSION ) ) return()
    vers <- try( readLines( "http://dreiss-isb.github.io/cMonkey/VERSION" ), silent=T )
    if ( class( vers ) != "try-error" ) {
      vers <- gsub( " ", "", vers )
      if ( vers != VERSION ) packageStartupMessage( "\nYou are not using the most current version of cMonkey.\nPlease consider upgrading to v", vers, " via:\n\n> install.packages('devtools', dep=T)\n> require(devtools)\n> install_github('cmonkey', 'dreiss-isb', subdir='cMonkey')" )
      else packageStartupMessage( "Congratulations! You are using the latest version of cMonkey.\n" )
    } else {
      packageStartupMessage( "Could not check to see if you are using the latest version of cMonkey." )
    }
  }
  
  source( "../generic_scripts/construct.package.R" )

  ## Note: eventually, want to include "ref" package, turn off the "gsubs" option, make "matrix.reference"
  ##    return ref(m), etc.. These are currently all turned off for the official package.
  pkgname <- 'cMonkey'; if ( beta ) pkgname <- 'cMonkey.beta'
  construct.package( pkgname, version=version, source.files=source.files, nocpp=beta,
                    functions.visible=c( "cmonkey", "cmonkey.init" ), ##, "plotClust", "cluster.summary", "plotStats",
                      ##"get.rows", "get.cols", "clusters.w.func", "clusters.w.conds", "clusters.w.genes", 
                      ##"write.project", "update.cmonkey.env", "save.cmonkey.env" ),
                    functions.excluded=if ( beta ) "NONE" else
                    c( "cm.package", "construct.package", "pssm.to.consensus",
                      ### FUNCTIONS NOT YET NEEDED IN OFFICIAL PACKAGE:
                      "get.prolinks.links", "get.predictome.links", "load.ratios.GEO", "load.ratios.MicrobesOnline",
                      ### OLD FUNCTIONS:
                      "get.STRING.links.OLD", "pareto.adjust.weights.OLD",
                      ### EXPERIMENTAL FUNCTIONS:
                      "get.condition.groups", "get.col.weights", "get.row.weights", 
                      "filter.updated.memberships", "pareto.adjust.weights", 
                      "consolidate.duplicate.clusters", "re.seed.empty.clusters",
                      ##"row.col.membership.from.clusterStack",
                      ##"matrix.reference", 
                      "list.reference", "ffify.env", "un.ffify.env",
                      ## "get.STRING.links.NEW",
                      ### EXPERIMENTAL POST PROCESSING FUNCTIONS:
                      "cluster.GO.annotations", "cluster.GO.pvalues", "cluster.KEGG.pvalues",
                      "motif.similarities.tomtom", "desymmetrize.tomtom.results", "cluster.tomtom.results",
                      "motif.similarities.custom", "cluster.custom.results", "compare.pssms",
                      "process.history", "compare.clusters",  "adjust.clust.2", ##"update.cmonkey.env",
                      "residual.pvalue", "pssm.motif.lines", "cluster.score.pvalues", "bicluster.pvalues",
                      ### NEW (EXPERIMENTAL) MOTIF FUNCTIONS:
                      "get.sequence.psps", "get.row.scores.NEW",
                      "weeder.one.cluster", "spacer.one.cluster",
                      "cosmo.one.cluster", "glam.one.cluster", "gadem.one.cluster", 
                      "blast.align", "parse.blast.out", "blast.match.seqs", "all.dna.seqs",
                      ### NEW OUTPUT FUNCTIONS:
                      "write.bicluster.network"
                      ),
                    ##data=list( halo="halo", hpy="hpy", mpn="mpn" ), ##ecoli="ecoli", yeast="yeast", ath="ath",
                    ##objects.included=c( cm.init="cm.init", cm.main="cm.main" ),
                    required=if ( beta ) c( "ref", "bigmemory", "filehash" ) else "",
                    suggested=c( "RCurl", "doMC", "foreach", "igraph0", "RSVGTipsDevice", "hwriter", "parallel" ), 
                      ##"ref", "bigmemory", "filehash" ), ##"fUtilities", , "Cairo", "trimcluster"
                    short.desc="cMonkey integrated biclustering algorithm",
                    long.desc="cMonkey integrated biclustering algorithm for learning co-regulated gene modules",
                    url="https://github.com/dreiss-isb/cmonkey",
                    reference='"Integrated biclustering of heterogeneous genome-wide datasets for the inference of global regulatory networks",\nby David J Reiss, Nitin S Baliga and Richard Bonneau: \\url{http://www.biomedcentral.com/1471-2105/7/280}',
                    onLoad=onLoad, gsubs=if ( beta ) list( c( "[,]", "" ), c( "[]", "" ) ) else NULL )

  pkg.name <- sprintf( "lib/%s_%s.tar.gz", pkgname, version )
  ##if ( beta ) {
  ##system( sprintf( "mv %s lib/cMonkey_%s_beta.tar.gz", pkg.name, version ) )
  ##  pkg.name <- sprintf( "lib/cMonkey_%s_beta.tar.gz", version )
  ##}
  if ( install ) system( sprintf( "R CMD INSTALL %s", pkg.name ) )

  if ( check ) {
    cwd <- setwd( "lib" )
    system( sprintf( "R CMD CHECK %s", pkg.name ) )
    setwd( cwd )
  }

  if ( beta ) return()
  
  if ( data ) {
    ## Halo data goes in data package
    cat( "Packaging Halo data...\n" )
    halo <- list( organism="hal", cog.org="Hbs", rsat.species="Halobacterium_sp", taxon.id=64091, k.clust=250 )
    load( "data/halo/halo.rats.RData" ); halo$ratios <- ratios
    halo$translation.tab <- read.delim( "data/halo/vng7000_to_vng5000.tsv" )
    colnames( halo$translation.tab ) <- c( "V1", "V2" )
    halo$bat.clust.genes <- c("VNG0654C","VNG0656H","VNG0828H","VNG0829G","VNG0830G","VNG0831G","VNG0832C","VNG1039H",
                              "VNG1370G","VNG1405C","VNG1458G","VNG1459H","VNG1461H","VNG1462G","VNG1463G","VNG1464G",
                              "VNG1465G","VNG1467G","VNG1468H","VNG1628G","VNG1630H","VNG1655H","VNG1656H","VNG1657H",
                              "VNG1755G","VNG1882G","VNG1883G","VNG1884G","VNG2137G","VNG2535H")
    halo$cm.func.each.iter <- function(){if(!iter%in%plot.iters)return();fc<-clusters.w.genes(bat.clust.genes);
                                         cat("HERE:",which.max(fc),max(fc),"\n")}
    halo$favorite.cluster <- function() which.max( clusters.w.genes( bat.clust.genes ) )
    halo$net.weights <- c( string=1 ) ## prolinks=1 ##0.5, operons=0.5 )
    halo <<- halo

    ## Hpy data goes in data package
    cat( "Packaging Hpy data...\n" )
    hpy <- list( organism="hpy", k.clust=75 )
    load( "data/hpy/hpy.rats.RData" ); hpy$ratios <- ratios
    source( "cmonkey-data-load.R" ); hpy$pp.ints <- load.sif.interactions( "data/hpy/HPY-PPINTS.sif.gz" )
    hpy$cm.func.each.iter=function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func("flagell");##flag.clusters();
                                     cat("FLAG CLUSTER:",which.max(fc),max(fc),"\n")}
    hpy$favorite.cluster <- function() which.max( clusters.w.func( "flagell" ) )
    hpy$net.weights <- c( string=1, pp.ints=1 ) ## operons=0.5, ##, `data/hpy/HPY-PPINTS.sif.gz`=1 )
    hpy <<- hpy

    ## Mycoplasma data goes in data package
    cat( "Packaging Myc data...\n" )
    mpn <- list( organism="mpn", k.clust=75 )
    ratios <- read.delim( gzfile( "data/Mycoplasma_pneumoniae/GSE14015_series_matrix.txt.gz" ), comment='!' )
    cnames <- grep( "^!Sample_title", readLines( gzfile( "data/Mycoplasma_pneumoniae/GSE14015_series_matrix.txt.gz" ),
                                                n=100 ), val=T )
    cnames <- gsub( "\302\272", "deg", gsub( "\"", "", strsplit( cnames, "\t" )[[ 1 ]] ) )
    rownames( ratios ) <- sapply( strsplit( as.character( ratios[[ 1 ]] ), "|", fixed=T ), "[", 1 )
    colnames( ratios ) <- cnames; rm( cnames )
    ratios <- as.matrix( ratios[ ,-1 ] )
    mpn$ratios <- ratios
    mpn <<- mpn
    
    source( "../generic_scripts/construct.package.R" )
    construct.package( "cMonkey.data", version=version, source.files=NULL, ##source.files,
                      data=list( halo="halo", hpy="hpy", mpn="mpn" ), ##required="cMonkey",
                      short.desc="Additional example data sets for use with the cMonkey integrated biclustering algorithm",
                      long.desc="Additional example data sets for use with the cMonkey integrated biclustering algorithm" )

    if ( install ) system( sprintf( "R CMD INSTALL lib/cMonkey.data_%s.tar.gz", version ) )
  }
  
  if ( bigdata ) { 
    cat( "Packaging Eco data...\n" ) ## Ecoli data package goes in "bigdata" package
    ecoli <- list( organism="eco", k.clust=300, rsat.species="Escherichia_coli_K12" )
    ecoli$ratios.old <- read.delim( file=gzfile( "data/ecoli/ratios.table.txt-cmpaper.gz" ), sep="\t", as.is=T, header=T )
    ecoli$ratios <- as.matrix( read.delim( file=gzfile( "data/ecoli/ratios.table.txt.gz" ), sep="\t", as.is=T, header=T ) )
    ecoli$cm.func.each.iter=function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func("flagell");##flag.clusters();
                                    cat("FLAG CLUSTER:",which.max(fc),max(fc),"\n")}
    ecoli$favorite.cluster <- function() which.max( clusters.w.func( "flagell" ) )
    ecoli <<- ecoli

    cat( "Packaging Sce data...\n" ) ## Yeast data package goes in "bigdata" package
    yeast <- list( organism="sce", is.eukaryotic=T, k.clust=450, motif.upstream.search=c( -30, 250 ),
                  motif.upstream.scan=c( -30, 500 ) )
    yeast$meme.cmd <- "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 12 -mod zoops -nostatus -text -cons $none"
    yeast$ratios <- as.matrix( read.delim( file=gzfile( "data/yeast/ratios.table.txt-GOOD.gz" ), sep="\t", as.is=T, header=T ) )
    yeast$preprocess.ratios <- function( ratios ) { ## remove rows with <10% of conditions changing
      cat( "SCE: Filtering out small change (<10% with >2-fold) rows from ratios matrix...\n" )
      nochange <- apply( 2^ratios, 1, function( i ) mean( ! is.na( i ) & ( i < 0.5 | i > 2 ), na.rm=T ) ) < 0.1
      ratios <- ratios[ ! nochange, ]
      super.preprocess.ratios( ratios ) ## Default preprocessing and mean/variance scaling
    }
    source( "cmonkey-data-load.R" )
    yeast$pp.ints <-load.sif.interactions( "data/yeast/yeast-ppInt.sif.gz" )
    yeast$bind.ints <-load.sif.interactions( "data/yeast/BIND_4932.sif.gz" )
    yeast$chip.chip <-load.sif.interactions( "data/yeast/yeast-chipChip.sif.gz" )
    yeast$net.weights <- c( string=0.25, pp.ints=0.5, bind.ints=0.1 ) ##string=0.25
    yeast$grouping.weights <- c( chip.chip=1 )
    yeast$cm.func.each.iter <- function(){if(!iter%in%plot.iters)return();fc<-clusters.w.func('proteasom');
                                       cat("HERE:",which.max(fc),max(fc,na.rm=T),"\n")}
    yeast$favorite.cluster <- function() which.max( clusters.w.func( "proteasom" ) )
    yeast <<- yeast

    cat( "Packaging Ath data...\n" ) ## arabidobsis data "ath" package goes in "bigdata" package
    ath <- list( organism="ath", rsat.species="Arabidopsis_thaliana", is.eukaryotic=T, k.clust=600,
                motif.upstream.search=c( -20, 500 ), motif.upstream.scan=c( -30, 1000 ) )
    ath$meme.cmd <- "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -maxw 12 -mod zoops -nostatus -text -cons $none"
    ratios <- "data/ath/AtGE_Abiotic_SR_2DRnew.txt.gz" ##AtGE_Abiotic_SR_2DR.txt.gz"
    ratios <- read.delim( file=gzfile( ratios ), sep="\t", as.is=T, header=T )
    ##attr( ratios, "rnames" ) <- toupper( ratios[ ,1 ] ); ath$ratios <- as.matrix( ratios[ ,-1 ] ); rm( ratios )
    rownames( ratios ) <- toupper( ratios[ ,1 ] ); ath$ratios <- as.matrix( ratios[ ,-1 ] ); rm( ratios )
    ath$preprocess.ratios <- function( ratios ) { ## remove rows with <50% of conditions changing
      cat( "ATH: Filtering out small change (<10% with >2-fold) rows from ratios matrix...\n" )
      nochange <- apply( 2^ratios, 1, function( i ) mean( ! is.na( i ) & ( i < 0.5 | i > 2 ), na.rm=T ) ) < 0.1
      ratios <- ratios[ ! nochange, ]
      super.preprocess.ratios( ratios ) ## Default preprocessing and mean/variance scaling
    }
    ath <<- ath
    
    source( "../generic_scripts/construct.package.R" )
    construct.package( "cMonkey.bigdata", version=version, source.files=NULL, ##source.files,
                      data=list( ecoli="ecoli", yeast="yeast", ath="ath" ), ##required="cMonkey",
                      short.desc="Additional example (big)data for use with the cMonkey integrated biclustering algorithm",
                      long.desc="Additional example (big)data for use with the cMonkey integrated biclustering algorithm" )

    if ( install ) system( sprintf( "R CMD INSTALL lib/cMonkey.bigdata_%s.tar.gz", version ) )
  }
  
  if ( update.web ) {
    system( sprintf( "cp -fv lib/index.html lib/cMonkey*_%s.tar.gz ~/Sites/cMonkey/", version ) )
    system( sprintf( "rpl VERSION \"%s\" ~/Sites/cMonkey/index.html", version ) )
    system( "cp -fv ~/Sites/cMonkey/index.html ~/Sites/cMonkey/cmonkey.html" )
    if ( install ) {
      cwd <- setwd( "~/Library/R/packages" ) ## This will change - works for pinnacle!
      system( sprintf( "zip -r cMonkey_%s.zip cMonkey", version ) )
      system( sprintf( "zip -r cMonkey.data_%s.zip cMonkey.data", version ) )
      if ( bigdata ) system( sprintf( "zip -r cMonkey.bigdata_%s.zip cMonkey.bigdata", version ) )
      system( sprintf( "mv -v cMonkey*.zip %s/lib/", cwd ) )
      setwd( cwd )
    }
    md5sums <- system( sprintf( "md5sum lib/cMonkey*_%s*", version, cwd ), intern=T )
    cat( sprintf( "VERSION %s", version ), md5sums, "\n", sep="\n", file="lib/md5sums.txt" )
    cat( version, "\n", file="lib/VERSION" )
    system( sprintf( "scp lib/VERSION lib/md5sums.txt ~/Sites/cMonkey/cmonkey.html lib/cMonkey*_%s.tar.gz lib/cMonkey*_%s.zip bragi:/local/apache2/htdocs/cmonkey/", version, version ) )
    system( sprintf( "scp lib/cMonkey_%s.tar.gz bragi:/local/apache2/htdocs/cmonkey/cMonkey_latest.tar.gz", version ) )
  }
}

cmonkey.unit.test <- function( plot=T ) {
  rm( list=ls()[ ! ls() %in% c( 'plot' ) ] )
  rm( list=ls()[ ! ls() %in% c( 'plot' ) ], envir=.GlobalEnv )
  require( cMonkey )
  setwd( tempdir() )
  print( getwd() )
  data( hpy, package='cMonkey.data' )
  e <- cmonkey.init( hpy )
  if ( plot ) cmonkey( e, dont.init=T )
  else cmonkey( e, dont.init=T, plot.iters=0 )
  
  ## require( RUnit ) ## dont know how to use this yet
}

cmonkey.init.ec2 <- function( inst, ... ) {
  source( "~/scratch/halo/generic_scripts/ec2-utils.R" )
  ec2.setenv()
  ##inst <- ec2.spawn.spot.instances( ... )
  ##ec2.setup( inst, ec2.tools=F, full=F )
  f <- tempfile()
  cran.repos <- 'cran.cnr.Berkeley.edu'
  cat( "sudo killall -9 console-kit-daemon\n",
      "sudo mv /usr/sbin/console-kit-daemon /usr/sbin/console-kit-daemon.bkup\n",
      "sudo cp /bin/true /usr/sbin/console-kit-daemon\n",
      sprintf( "tar -xvzf %s\n", basename( options( "ec2.tools.tarball" )$ec2.tools.tarball ) ),
      "mkdir biclust; cd biclust; tar -xvzf ~/z_for_ec2.tgz; cd ~/\n",
      ##"sudo dd if=/dev/zero of=/mnt/2Gb.swap bs=1M count=2048; sudo mkswap /mnt/2Gb.swap\n",
      ##"sudo swapon /mnt/2Gb.swap; sudo sysctl vm.swappiness=10\n",
      "sudo apt-get -f -y --force-yes update\n",
      "sudo apt-get -f -y --force-yes install r-base-core r-base-dev tcsh r-recommended rpl s3cmd libcurl4-openssl-dev\n", ##r-cran-mass r-cran-vr r-cran-matrix ruby subversion libopenssl-ruby littler python-rpy pdftk openjdk-6-jdk openjdk-6-jre emacs23-nox pdftk r-cran-doSNOW heirloom-mailx sendmail
      "echo \"*       soft    nofile  1024\" | sudo tee -a /etc/security/limits.conf\n",
      "echo \"*       hard    nofile  65535\" | sudo tee -a /etc/security/limits.conf\n",
      "echo \"ulimit -n 65535\" >> ~/.bashrc\n", ## set the limit at each login
      "R CMD INSTALL lib/cMonkey_4.8.4.tar.gz\n",
      "mkdir ~/R; echo \".libPaths( c( \\\"~/R/\\\", .libPaths() ) )\" >~/.Rprofile\n",
      "echo \"require( utils ); require( graphics ); require( stats )\" >>~/.Rprofile\n",
      "echo \"my.utils <- attach( NULL, name=\\\"my.utils\\\" )\n",
      "sys.source( \\\"~/my.util.R\\\", env=my.utils )\" >>~/.Rprofile\n",
      sprintf( "Rscript -e \"Sys.setenv(MAKE=\\\"make -j 8\\\");install.packages( c(%s), dep=F, lib=\\\"~/R\\\", repos=\\\"http://%s/\\\" )\"\n", "\\\"cMonkey\\\",\\\"parallel\\\",\\\"igraph0\\\",\\\"RSVGTipsDevice\\\",\\\"glmnet\\\",\\\"lars\\\"", cran.repos ),
      sprintf( "Rscript -e \"Sys.setenv(MAKE=\\\"make -j 8\\\");install.packages( c(%s), dep=T, lib=\\\"~/R\\\", repos=\\\"http://%s/\\\" )\"\n", "\\\"doMC\\\",\\\"RCurl\\\"", cran.repos ),
      "mkdir ~/progs; cd ~/progs; wget \"http://meme.nbcr.net/downloads/old_versions/meme_4.3.0.tar.gz\"\n",
      "tar -xvzf meme_4.3.0.tar.gz; cd ~/progs/meme_4.3.0; mkdir local\n",
      "./configure --prefix=`pwd`/local/ --enable-dependency-tracking --enable-opt --disable-shared --disable-fast-install --enable-serial\n",
      "make -j 8; make install; cd ~/progs/; ln -s meme_4.3.0/local/bin/meme\n",
      "ln -s meme_4.3.0/local/bin/mast; ln -s meme_4.3.0/local/bin/dust; ln -s meme_4.3.0/local/bin/tomtom\n",
      "cd ~/progs; tar -xvzf ~/weeder1.4.2.tar.gz; cd Weeder1.4.2; ./compileall\n",
      "cd ~/progs; ln -s Weeder1.4.2/weederlauncher.out; ln -s Weeder1.4.2/weederTFBS.out; ln -s Weeder1.4.2/adviser.out\n",
      sep="", file=f )
  system( "tar cvf - `find . -name '*.R'` lib/cMonkey_4.8.4.tar.gz |gzip -c >z_for_ec2.tgz" )
  ec2.upload( inst, sprintf( "z_for_ec2.tgz '%s' '%s' '%s' '%s'", options( "ec2.tools.tarball" )$ec2.tools.tarball, f,
                            '/Users/dreiss/scratch/halo/generic_scripts/my.util.R', './progs/weeder1.4.2.tar.gz' ),
             wait=T, via.s3=F )
  ec2.exec( inst, sprintf( "source %s", basename( f ) ) )
  ec2.get.from.s3( inst, file="cmonkey_data.tgz", bucket="s3://dreiss-data/cmonkey/" )
  ec2.exec( inst, "cd ~/biclust; tar -xvzf ~/cmonkey_data.tgz; ln -s ~/progs" )
  inst
}
#endif
