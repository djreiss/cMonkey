
load.ratios.MicrobesOnline <- function( id=taxon.id ) { ## mostly b/c ecoli taxon.id=83333, but need to use id=511145
  f <- sprintf( "data/%s/microbesOnlineExprData_%d.tsv", rsat.species, id )
  if ( ! file.exists( sprintf( "%s.gz", f ) ) ) {
    timeout <- options( timeout=600 ) ## they can be very slow...
    try( dlf( f, sprintf( "http://www.microbesonline.org/cgi-bin/microarray/getData.cgi?taxId=%d", id ) ) )
    system( sprintf( "gzip -fv %s", f ) )
    options( timeout=timeout$timeout )
  }
  cat( "Loading:", f, "\n" )
  rats <- read.delim( gzfile( sprintf( "%s.gz", f ) ), row.names=1, sep="\t", as.is=T, head=T )
  as.matrix( rats )
}

load.genome.info.MicrobesOnline <- function( id=taxon.id ) {
  f <- sprintf( "data/%s/microbesOnlineGenomeInfo_%d.tsv", rsat.species, id )
  if ( ! file.exists( sprintf( "%s.gz", f ) ) ) {
    try( dlf( f, sprintf( "http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab", id ) ) )
    system( sprintf( "gzip -fv %s", f ) )
  }
  cat( "Loading:", f, "\n" )
  out <- read.delim( gzfile( sprintf( "%s.gz", f ) ), row.names=1, sep="\t", as.is=T, head=T )
  out
} 

## org.code can be string, e.g. 'Escherichi coli' or number e.g. 511145
load.data <- function( org.code=NULL ) { 
  try( dlf( 'data/NCBI/taxdump.tar.gz', sprintf( "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" ) ) )
  system( 'cd data/NCBI; tar -xvzf taxdump.tar.gz names.dmp; gzip -v names.dmp' )
  if ( ! is.null( org.code ) ) {
    lines <- system( sprintf( 'gunzip -c data/NCBI/names.dmp.gz | grep "%s"', as.character( org.code ) ), intern=T )
  }
  lines <- do.call( rbind, lapply( strsplit( lines, '\t' ), function( i ) i[ i != '' & i != '|' ] ) )[ ,1:2 ]
  while ( length( unique( lines[ ,1 ] ) ) > 1 ) {
    ## Add a menu to let the user select the right code
  }
  taxon.id <- as.integer( unique( lines[ ,1 ] ) )
  org.name <- unique( lines[ ,2 ] )

  ## Genome sequences
  try( dlf( sprintf( '%s/genome.txt', data.dir ),
           sprintf( 'http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%s;export=genome', taxon.id ) ) )
  system( sprintf( 'gzip -v %s/genome.txt', data.dir ) )
  genome <- read.fasta( gzfile( sprintf( '%s/genome.txt.gz', data.dir ) ) )
  
  ## Gene annotations
  data.dir <- sprintf( 'data/%s/', gsub( '[.,; ]', '_', org.name[ 1 ] ) )
  try( dlf( sprintf( '%s/genome_info.tsv', data.dir ),
           sprintf( 'http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=%d;export=tab', taxon.id ) ) )
  system( sprintf( 'gzip -v %s/genome_info.tsv', data.dir ) )
  genome.info <- read.delim( gzfile( sprintf( '%s/genome_info.tsv.gz', data.dir ) ), head=T )

  invisible( list( taxon.id=taxon.id, org.name=org.name, genome=genome, genome.info=genome.info ) )
}

## Use go annotations from microbes online annotation table
cluster.GO.annotations <- function( genes, org.code=taxon.id, which.pfc=c( 'BP', 'MF', 'CC' ) ) {
  anno <- load.genome.info.MicrobesOnline( id=org.code )
  anno <- subset( anno, sysName != '' & GO != '' )
  onto <- which.pfc[ 1 ]

  require( topGO )
  gos <- strsplit( as.character( anno$GO ), ',' )
  names( gos ) <- as.character( anno$sysName )
  allGenes <- as.factor( 1 - ( names( gos ) %in% genes ) ); names( allGenes ) <- names( gos )
  
  GOdata <- new( "topGOdata", description='test', ontology=onto, allGenes=1-allGenes,
                annot=annFUN.gene2GO, gene2GO=gos, nodeSize=2 )
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups( GOdata, test.stat )
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata, test.stat)    
  test.stat <- new("classicScore", testStatistic = GOtTest, name = "T tests")
  resultTT <- getSigGroups(GOdata, test.stat)
  allRes1 <- GenTable( GOdata, KS=resultKS, T=resultTT, F=resultFisher, orderBy='KS', topNodes=20 )
  if ( plot ) showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')
  allRes1
}
