e$write.project <- function( ks=sapply( as.list( clusterStack ), "[[", "k" ), ##save.session=T, ##pdfs=T, ##dev="SVG", 
                          out.dir=NULL, gaggle=T, seq.type="upstream meme", gzip=T,
                          output=c("svg","pdf","png","html","main","rdata"), ... ) { ##network=F, 
  if ( is.null( out.dir ) ) out.dir <- cmonkey.filename
  cat( "Outputing to", out.dir, "\n" )
  ##if ( ! is.null( dev ) && dev == "SVG" )
  ##require( RSVGTipsDevice )
  if ( ! file.exists( out.dir ) ) ##&& ! is.null( dev ) && dev == "SVG" )
    dir.create( out.dir, recursive=T, showWarnings=F )
  ##require( fUtilities )

  ##clusterStack <- get.clusterStack( ks=1:k.clust )
  ##ks <- sapply( clusterStack, "[[", "k" )
  clusterStack <- clusterStack[ ks ]
  mc <- get.parallel( length( ks ) )
  has.pdftk <- length( system( "which pdftk", intern=T ) ) > 0 ## Use for compressing pdfs

  if ( ! file.exists( paste( out.dir, "/svgs", sep="" ) ) )
    dir.create( paste( out.dir, "/svgs", sep="" ), showWarnings=F )
  if ( ! file.exists( paste( out.dir, "/pdfs", sep="" ) ) ) ##pdfs && 
    dir.create( paste( out.dir, "/pdfs", sep="" ), showWarnings=F )
  if ( ! file.exists( paste( out.dir, "/htmls", sep="" ) ) )
    dir.create( paste( out.dir, "/htmls", sep="" ), showWarnings=F )
  
  ##if ( ! is.null( dev ) && dev == "SVG" ) {
  if ( "svg" %in% output ) {
    require( RSVGTipsDevice )
    if ( ! file.exists( sprintf( "%s/svgs/stats.svg", out.dir ) ) &&
        ! file.exists( sprintf( "%s/svgs/stats.svgz", out.dir ) ) ) {
      ##cat( out.dir, "/svgs/stats.svg", "\n", sep="" )
      cat( "STATS...\n" )
      devSVGTips( sprintf( "%s/svgs/stats.svg", out.dir ), toolTipMode=2, title="Biclustering statistics",
                 xmlHeader=T )
      par( family="Arial" )
      plotStats( new.dev=F )
      dev.off()
    }

    require( igraph )
    cat( "SVGS: " )
    for ( qqq in 1:3 ) {
    mc$apply( ks, function( k ) { ## will this work in parallel? seems to. No, it doesn't.
      ##k <- ks[ i ]
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/svgs/cluster%04d.svg", out.dir, k ) ) ||
          file.exists( sprintf( "%s/svgs/cluster%04d.svgz", out.dir, k ) ) ) return( NULL )
      devSVGTips( sprintf( "%s/svgs/cluster%04d.svg", out.dir, k ), toolTipMode=2,
                 title=sprintf( "Bicluster %04d", k ), xmlHeader=T )
      ##try(
      plotClust( k, w.motifs=T, seq.type=seq.type, ... ) ##)
      dev.off()
    } )
    }
    cat( "\n" )
  }
  
  ##if ( pdfs ) {
  if ( "pdf" %in% output ) {
    require( igraph ) ## load it here
    cat( "PDFS: " )
    mc$apply( ks, function( k ) { ## will this work in parallel? seems to.
      ##k <- ks[ i ]
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/pdfs/cluster%04d.pdf", out.dir, k ) ) ) return( NULL )
      if ( has.pdftk ) pdf( sprintf( "%s/pdfs/cluster%04d.pdf", out.dir, k ) ) ## Will be compressed later
      else cairo_pdf( sprintf( "%s/pdfs/cluster%04d.pdf", out.dir, k ) ) ## Writes out smaller PDFs
      try( plotClust( k, w.motifs=T, seq.type=seq.type, ... ), silent=T )
      dev.off()
    } )
    cat( "\n" )
  }
  
  if ( gaggle && "html" %in% output ) { ## Embed the SVG in some gaggleized html for use w/ firegoose --
    require( hwriter ) ## see http://www.ebi.ac.uk/~gpau/hwriter/
    ## see http://gaggle.systemsbiology.net/docs/geese/firegoose/microformat/
    ## Note that "embed" might not be the best way to embed the SVG --
    ##    see http://www.w3schools.com/svg/svg_inhtml.asp for other options.
    ##'%^%' <- function( a, b ) paste( a, b, sep="\n" )
    cat( "HTMLS: " )
    mc$apply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d.html", out.dir, k ) ) ) return()
      rows <- sort( get.rows( k ) )
      if ( length( rows ) <= 0 ) return()
      short.names <- get.long.names( rows, short=T )
      short.names <- cbind( rows, short.names )
      rownames( short.names ) <- colnames( short.names ) <- NULL
      long.names <- get.long.names( rows, short=F )
      long.names <- cbind( rows, long.names )
      rownames( long.names ) <- colnames( long.names ) <- NULL
      refseq.names <- unique( unlist( get.synonyms( rows ) ) )
      refseq.names <- grep( "^NP_", refseq.names, val=T )
      upstream.seqs <- try( get.sequences( k, filter=F, uniq=F ), silent=T ) ##uniquify=F, silent=T )
      if ( class( upstream.seqs ) == "try-error" || is.null( upstream.seqs ) ) {
        upstream.seqs <- rep( "", length( rows ) ); names( upstream.seqs ) <- rows }
      upstream.seqs <- cbind( names( upstream.seqs ), upstream.seqs )
      rownames( upstream.seqs ) <- colnames( upstream.seqs ) <- NULL

      htmltext <-
        paste( c( "<html><head><title>Bicluster %K (%FILE)</title>", ##, k, cmonkey.filename ),
                 "<style type=\"text/css\">",
                 "  .hidden {", "     display: none;", "   }",
                 "  .gaggle-data {", "     color: green;", "     font-size: xx-small;", "   }",
                 "   p {", "     color: red;", "     font-size: x-small;", "   }",
                 "</style>",
                 "<script type=\"text/javascript\">",
                 "   function toggleVisible(id){",
                 "      if (document.getElementById){",
                 "         obj = document.getElementById(id);",
                 "         if (obj) {",
                 "            if (obj.style.display == 'none'){",
                 "               obj.style.display = 'block';",
                 "            } else {",
                 "               obj.style.display = 'none';",
                 "            }",
                 "         }",
                 "      }",
                 "   }",
                 "</script>",
                 "</head>",
                 "<table><tr><td>",
                 "<iframe src=\"../svgs/cluster%K03d%K.svg\" width=\"600\" height=\"520\" frameborder=\"0\"></iframe>",
                 "</td><td>",
                 
                 "<p><a href=\"#bicluster%K03d%K\" onclick=\"toggleVisible('bicluster%K03d%K'); return false;\">[+]</a>",
                 "Show/hide bicluster #%K rows and columns.</p>", 
                 "<div id=\"bicluster%K03d%K\" style=\"display:none;\" class=\"gaggle-data bicluster\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>", length( rows ), length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-cluster\">",
                 "      <ol class=\"gaggle-rowNames\">",
                 paste( "<li>", sort( rows ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   <ol class=\"gaggle-columnNames\">",
                 paste( "<li>", sort( get.cols( k ) ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_genes\" onclick=\"toggleVisible('bicluster%K03d%K_genes'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (genes).</p>", 
                 "<div id=\"bicluster%K03d%K_genes\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K genes</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( rows ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( rows ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_short_names\" onclick=\"toggleVisible('bicluster%K03d%K_short_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (short gene names).</p>", 
                 "<div id=\"bicluster%K03d%K_short_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K short names</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( short.names ) ),
                 "   <span class=\"gaggle-namelist-tag hidden\">short_name</span>",
                 hwrite( short.names, table.class="toc", col.class=list( NA, "short_name" ), border=1,
                        table.style="font-family: monospace; font-size: xx-small; color: green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_long_names\" onclick=\"toggleVisible('bicluster%K03d%K_long_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (long gene names).</p>", 
                 "<div id=\"bicluster%K03d%K_long_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K long names</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( long.names ) ),
                 "   <span class=\"gaggle-namelist-tag hidden\">long_name</span>",
                 hwrite( long.names, table.class="toc", col.class=list( NA, "long_name" ), border=1,
                        table.style="font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_refseq_names\" onclick=\"toggleVisible('bicluster%K03d%K_refseq_names'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K rows (NCBI RefSeq gene IDs).</p>", 
                 "<div id=\"bicluster%K03d%K_refseq_names\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K NCBI RefSeq IDs</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( refseq.names ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( refseq.names ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_upstream_seqs\" onclick=\"toggleVisible('bicluster%K03d%K_upstream_seqs'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K gene upstream sequences.</p>", 
                 "<div id=\"bicluster%K03d%K_upstream_seqs\" style=\"display:none;\" class=\"gaggle-data genes\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K upstream sequences</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", nrow( upstream.seqs ) ), 
                 "   <span class=\"gaggle-namelist-tag hidden\">upstream</span>",
                 hwrite( upstream.seqs, table.class="toc", col.class=list( NA, "upstream" ), border=1, 
                        table.style="font-family: monospace; font-size: xx-small; color:green; border-collapse: collapse" ),
                 "   </div>",
                 
                 "<p><a href=\"#bicluster%K03d%K_arrays\" onclick=\"toggleVisible('bicluster%K03d%K_arrays'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K columns (arrays; conditions).</p>", 
                 "<div id=\"bicluster%K03d%K_arrays\" style=\"display:none;\" class=\"gaggle-data arrays\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K arrays</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%d</span>", length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-namelist\">",
                 "      <ol>",
                 paste( "<li>", sort( get.cols( k ) ), "</li>", sep="", collapse="" ),
                 "      </ol>",
                 "   </div>",
                 "</div>",

                 "<p><a href=\"#bicluster%K03d%K_ratios\" onclick=\"toggleVisible('bicluster%K03d%K_ratios'); return false;\">[+]</a>", 
                 "Show/hide bicluster #%K ratios.</p>", 
                 "<div id=\"bicluster%K03d%K_ratios\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                 "   <span class=\"gaggle-name hidden\">bicluster %K ratios</span>", 
                 "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                 sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>", length( rows ), length( get.cols( k ) ) ),
                 "   <div class=\"gaggle-matrix-tsv\">",
                 "        RATIOS",
                 "   </div>",
                 "</div>",

                 if ( ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out ) &&
                     ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]] ) )
                 paste( "<p><a href=\"#bicluster%K03d%K_pssm1\" onclick=\"toggleVisible('bicluster%K03d%K_pssm1'); return false;\">[+]</a>", 
                       "Show/hide bicluster #%K motif PSSM #1.</p>", 
                       "<div id=\"bicluster%K03d%K_pssm1\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                       "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #1</span>", 
                       "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                       sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>",
                               nrow( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm ),
                               ncol( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm ) ),
                       "   <div class=\"gaggle-matrix-tsv\">",
                       "           MOTIF1",
                       "   </div>",
                       "</div>" ) else "",

                 if ( length( meme.scores[[ seq.type ]][[ k ]]$meme.out ) >= 2 &&
                     ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]] ) )
                 paste( "<p><a href=\"#bicluster%K03d%K_pssm2\" onclick=\"toggleVisible('bicluster%K03d%K_pssm2'); return false;\">[+]</a>", 
                       "Show/hide bicluster #%K motif PSSM #2.</p>", 
                       "<div id=\"bicluster%K03d%K_pssm2\" style=\"display:none;\" class=\"gaggle-data ratios\">", 
                       "   <span class=\"gaggle-name hidden\">bicluster %K motif PSSM #2</span>", 
                       "   <span class=\"gaggle-species hidden\">%SPECIES</span>", 
                       sprintf( "   <span class=\"gaggle-size hidden\">%dx%d</span>",
                               nrow( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm ),
                               ncol( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm ) ), 
                       "   <div class=\"gaggle-matrix-tsv\">",
                       "           MOTIF2",
                       "   </div>",
                       "</div>" ) else "",
                 
                 "</td></table>",
                 if ( "pdf" %in% output ) sprintf( "<a href=\"../pdfs/cluster%04d.pdf\">View PDF version</a>", k ) else "",
                 "</html>" ), collapse="\n" )

      rm( short.names, long.names, refseq.names, upstream.seqs )

      htmltext <- gsub( "%K03d%K", sprintf( "%04d", k ), htmltext )
      htmltext <- gsub( "%K", k, htmltext )
      htmltext <- gsub( "%FILE", cmonkey.filename, htmltext )
      htmltext <- gsub( "%SPECIES", gsub( "_", " ", rsat.species ), htmltext )

      tmp <- as.data.frame( get.cluster.matrix( rows, get.cols( k ) ) )
      tmp <- cbind( GENES=rownames( tmp ), tmp )
      tf <- tempfile()
      write.table( tmp, file=tf, sep="\t", quote=F, row.names=F ); rm( tmp )
      htmltext <- sub( "RATIOS", paste( readLines( tf ), collapse="\n" ), htmltext )
      unlink( tf )

      if ( ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out ) ) {
        if ( ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]] ) ) {
          tmp <- as.data.frame( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 1 ]]$pssm )
          if ( ! is.null( tmp ) && nrow( tmp ) > 0 ) {
            tmp <- cbind( 1:nrow( tmp ), tmp )
            colnames( tmp ) <- c( "POSITION", "A", "C", "G", "T" )
            write.table( tmp, file=tf, sep="\t", quote=F, row.names=F )
            htmltext <- sub( "MOTIF1", paste( readLines( tf ), collapse="\n" ), htmltext )
            unlink( tf )
          }
          rm( tmp )
        }
        
        if ( length( meme.scores[[ seq.type ]][[ k ]]$meme.out ) >= 2 &&
            ! is.null( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]] ) ) {
          tmp <- as.data.frame( meme.scores[[ seq.type ]][[ k ]]$meme.out[[ 2 ]]$pssm )
          if ( ! is.null( tmp ) && nrow( tmp ) > 0 ) {
            tmp <- cbind( 1:nrow( tmp ), tmp )
            colnames( tmp ) <- c( "POSITION", "A", "C", "G", "T" )
            write.table( tmp, file=tf, sep="\t", quote=F, row.names=F )
            htmltext <- sub( "MOTIF2", paste( readLines( tf ), collapse="\n" ), htmltext )
            unlink( tf )
          }
          rm( tmp )
        }
      }
      rm( tf )

      cat( htmltext, file=sprintf( "%s/htmls/cluster%04d.html", out.dir, k ), sep="\n" )
      rm( htmltext )
    }, mc.preschedule=F )
    cat( "\n" )
  }

  ##require( Cairo )
  if ( "png" %in% output ) {
    ##parallel.cores <- 1
    mc <- get.parallel( length( ks ), para=1 )
    cat( "PROFILES: " )
    ##for ( k in ks ) {
    mc$apply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d_profile.png", out.dir, k ) ) ) return() ##next
      try( {
        c <- get.clust( k )
        png( sprintf( "%s/htmls/cluster%04d_profile.png", out.dir, k ), width=128, height=64, antialias="subpixel" )
        par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
        plotCluster( c, main="", no.par=T, ... )
        dev.off() }, silent=T )
    } )
    cat( "\n" )
    
    cat( "NETWORKS: " )
    require( igraph )
    mc$apply( ks, function( k, ... ) {
    ##for ( k in ks ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d_network.png", out.dir, k ) ) ) return() ##next
      try( {
        png( sprintf( "%s/htmls/cluster%04d_network.png", out.dir, k ), width=64, height=64, antialias="subpixel" )
        par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
        c <- get.clust( k )
        plotCluster.network( c, cex=0.3, no.legend=T, ... )
        dev.off() }, silent=T )
    } )
    cat( "\n" )
    
    cat( "MOTIFS: " )
    ##for ( k in ks ) {
    mc$apply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      e.vals <- lapply( meme.scores[[ seq.type ]][[ k ]]$meme.out, "[[", "e.value" )
      pssms <- lapply( meme.scores[[ seq.type ]][[ k ]]$meme.out, "[[", "pssm" )
      if ( length( pssms ) < 2 ) {
        for ( i in ( length( pssms ) + 1 ):2 ) {
          pssms[[ i ]] <- matrix( 0.25, nrow=6, ncol=4 )
          e.vals[[ i ]] <- Inf
        }
      }
      for ( pp in 1:length( pssms ) ) {
        if ( file.exists( sprintf( "%s/htmls/cluster%04d_pssm%d.png", out.dir, k, pp ) ) ) return() ##next
        try( { 
          png( sprintf( "%s/htmls/cluster%04d_pssm%d.png", out.dir, k, pp ), width=128, height=64,
              antialias="subpixel" )
          if ( is.matrix( pssms[[ pp ]] ) )
            try( viewPssm( pssms[[ pp ]], e.val=NA, mot.ind=pp, main.title=sprintf( "e=%.3g", e.vals[[ pp ]] ),
                          cex.main=0.7 ), silent=T ) ##e.val=NA, mot.ind=pp
          dev.off() }, silent=T )
      }
    } )
    cat( "\n" )
    
    cat( "MOTIF POSITIONS: " )
    ##for ( k in ks ) {
    mc$apply( ks, function( k, ... ) {
      if ( k %% 25 == 0 ) cat( k ) else cat( "." )
      if ( file.exists( sprintf( "%s/htmls/cluster%04d_mot_posns.png", out.dir, k ) ) ) return() ##next
      try( {
        png( sprintf( "%s/htmls/cluster%04d_mot_posns.png", out.dir, k ), width=128, height=12+6*length(get.rows(k)),
            antialias="subpixel" )
        par( mar=rep(0.5,4)+0.1, mgp=c(3,1,0)*0.75 )
        c <- plotClust( k, dont.plot=T, ... ) ##seq.type=seq.type, 
        plotClusterMotifPositions( c, cex=0.4, no.key=T, ... )
        dev.off() }, silent=F )
    } )
    cat( "\n" )
  }

  if ( "main" %in% output ) {
    mc <- get.parallel( length( ks ) )
    cat( "WRITING MAIN HTML TABLE..." )
    require( hwriter ) ## see http://www.ebi.ac.uk/~gpau/hwriter/
    dlf( paste( out.dir, "hwriter.css", sep="/" ), "http://www.ebi.ac.uk/~gpau/hwriter/hwriter.css" )
    dlf( paste( out.dir, "sorttable.js", sep="/" ), "http://www.kryogenix.org/code/browser/sorttable/sorttable.js" )
    cat( "..." )
    ##cluster.summ <- cluster.summary( ... )
    ##if ( nrow( cluster.summ ) <= 0 )
    cluster.summ <- cluster.summary( e.cutoff=NA, nrow.cutoff=NA )
    write.table( cluster.summ, file=paste( out.dir, "/cluster.summary.tsv", sep="" ), quote=F, sep="\t" )
    cat( "..." )
    html <- openPage( paste( out.dir, "/index.html", sep="" ), link.javascript="sorttable.js",
                     title=paste( "cMonkey bicluster summary for run", cmonkey.filename ), link.css='hwriter.css' )
    hwrite( paste( "<h2>cMonkey bicluster summary for run", cmonkey.filename, "</h2>" ), html )
    hwrite( '<ul><li>Download a tab-delimited version of this table', html, link="cluster.summary.tsv",
           style='font-size:75%' )
    hwrite( "<li>Download a list of each bicluster's gene members", html, link="cluster.members.genes.txt",
           style='font-size:75%' )
    hwrite( "<li>Download a list of each bicluster's array/condition members", html,
           link="cluster.members.arrays.txt", style='font-size:75%' )
    cat( "..." )
    ##hwrite( "<li>Network summary of biclusters", html, link="svgs/main.svg", style='font-size:75%' )
    ##if ( "svg" %in% output )
    hwrite( "<li>Plots of summary statistics of biclustering run", html,
                                    link="svgs/stats.svg", style='font-size:75%' )
    ##if ( "rdata" %in% output )
    hwrite( "<li>Saved cMonkey R session file", html, link="cm_session.RData",
                                      style='font-size:75%' )
    hwrite( "<li>Summary of cMonkey input parameters</ul>", html, link="cm.params.txt", style='font-size:75%' )
    hwrite( "<br><center><b>Bicluster summary</b></center><br>", html )
    hwrite( "<br><center><b>Sort the table by a given column by clicking on the column's header.<br>Click on bicluster link in first column for more info.</b></center><br>", html,
           style='font-size:60%' )
    cat( "..." )
    himg0 <- hwriteImage( sprintf( "htmls/cluster%04d_profile.png", as.integer( rownames( cluster.summ ) ) ),
                         table=F )
    himg0 <- hwrite( paste( himg0, sprintf( "Residual = %.3f", cluster.summ$resid ), sep="<br>" ), center=TRUE,
                    table=F )
    himg0a <- hwriteImage( sprintf( "htmls/cluster%04d_network.png", as.integer( rownames( cluster.summ ) ) ),
                          table=F )
    e.val.1 <- lapply( meme.scores[[ seq.type ]][ as.integer( rownames( cluster.summ ) ) ],
                      function( i ) i$meme.out[[ 1 ]]$e.value )
    for ( i in 1:length( e.val.1 ) ) if ( is.null( e.val.1[[ i ]] ) ) e.val.1[[ i ]] <- NA
    himg1 <- hwriteImage( sprintf( "htmls/cluster%04d_pssm1.png", as.integer( rownames( cluster.summ ) ) ),
                         table=F, title=sprintf( "E-val = %.3g", unlist( e.val.1 ) ) )
    ##himg1a <- hwrite( as.character( cluster.summ$consensus1 ), table=T )
    himg1 <- hwrite( paste( himg1, as.character( cluster.summ$consensus1 ), sep="<br>" ), center=TRUE, table=F )
    e.val.2 <- lapply( meme.scores[[ seq.type ]][ as.integer( rownames( cluster.summ ) ) ],
                      function( i ) i$meme.out[[ 2 ]]$e.value )
    for ( i in 1:length( e.val.2 ) ) if ( is.null( e.val.2[[ i ]] ) ) e.val.2[[ i ]] <- NA
    himg2 <- hwriteImage( sprintf( "htmls/cluster%04d_pssm2.png", as.integer( rownames( cluster.summ ) ) ),
                         table=F, title=sprintf( "E-val = %.3g", unlist( e.val.2 ) ) )
    himg2 <- hwrite( paste( himg2, as.character( cluster.summ$consensus2 ), sep="<br>" ), center=TRUE, table=F )
    himg2a <- hwriteImage( sprintf( "htmls/cluster%04d_mot_posns.png", as.integer( rownames( cluster.summ ) ) ),
                          table=F )
    cluster.summ$score <- sprintf( "%.3f", cluster.summ$score )
    ##cluster.summ$resid <- sprintf( "%.3f", cluster.summ$resid ) 
    rn <- rownames( cluster.summ )
    cat( "..." )
    cluster.summ.orig <- cluster.summ
    cluster.summ <- cbind( bicluster=cluster.summ$k, n.genes=cluster.summ$nrow,
                          n.arrays=sapply( as.integer( rownames( cluster.summ ) ),
                            function( i ) length( get.cols( i ) ) ), score=cluster.summ$score,
                          residual=sprintf( "%.3f", cluster.summ$resid ) )
    if ( "score.norm" %in% colnames( cluster.summ.orig ) )
      cluster.summ <- cbind( cluster.summ, score.norm=sprintf( "%.3f", cluster.summ.orig$score.norm ) ) ##
    rownames( cluster.summ ) <- rn
    rows <- list(); for ( k in as.integer( rn ) ) rows[[ k ]] <- sort( get.rows( k ) )
    himg3 <- hwrite( sapply( as.integer( rn ), function( k ) paste( rows[[ k ]], collapse=" " ) ), table=F )
    cat( "...\n" )
    himg4 <- hwrite( unlist( mc$apply( as.integer( rn ), function( k ) { if ( k %% 25 == 0 ) cat( k ) else cat( "." );
      tmp <- get.long.names( rows[[ k ]], short=T ); tmp <- unique( tmp[ ! tmp %in% rows[[ k ]] & tmp != "" ] )
      paste( tmp, collapse=" " ) } ) ), table=F ); cat( "\n" )
    himg5 <- hwrite( unlist( mc$apply( as.integer( rn ), function( k ) { if ( k %% 25 == 0 ) cat( k ) else cat( "." );
      tmp <- get.long.names( rows[[ k ]], short=F ); tmp <- unique( tmp[ ! tmp %in% rows[[ k ]] & tmp != "" ] )
      paste( tmp, collapse=" | " ) } ) ), table=F ); cat( "\n" )
    ## Let's make the table sortable using code from http://www.kryogenix.org/code/browser/sorttable/
    e.val.1[ is.na( e.val.1 ) ] <- 9e9; e.val.2[ is.na( e.val.2 ) ] <- 9e9
    nas <- rep( NA, nrow( cluster.summ ) )
    hwrite( cbind( cluster.summ[ ,1:min(6,ncol(cluster.summ)) ], profile=himg0, network=himg0a, motif1=himg1, motif2=himg2,
                  motif.posns=himg2a, probe.names=himg3, short.names=himg4, 
                  long.names=himg5 ), html, row.names=F, table.style='text-align:center;font-size:70%;font-family:Arial', table.class='sortable',
           row.style=list( 'font-weight:bold;text-align:center;font-size:70' ), 
           col.style=list( probe.names='font-size:70%', orf.names='font-size:50%', short.names='font-size:50%', long.names='font-size:50%', motif1='font-size:50%', motif2='font-size:50%' ),
           col.sorttable_customkey=list( residual=sprintf( "%.3f", cluster.summ.orig$residual ),
             score.norm=if ( "score.norm" %in% colnames( cluster.summ.orig ) )
             sprintf( "%.3f", cluster.summ.orig$score.norm ) else NULL,
             profile=sprintf( "%.3f", cluster.summ.orig$resid ),
             motif1=sprintf( "%.30f", unlist( e.val.1 ) ), e.val1=sprintf( "%.30f", unlist( e.val.1 ) ),
             motif2=sprintf( "%.30f", unlist( e.val.2 ) ), e.val2=sprintf( "%.30f", unlist( e.val.2 ) ) ),
           col.class=list( network=c( "sorttable_nosort", nas ), ##motif1=c( "sorttable_nosort", nas ), motif2=c( "sorttable_nosort", nas ),
             motif.posns=c( "sorttable_nosort", nas ) ),
           col.link=list( sprintf( "htmls/cluster%04d.html", as.integer( rownames( cluster.summ ) ) ) ) )
    closePage( html, splash=F )
    ##}
    
    for ( i in sapply( 1:k.clust, function( k ) c( k, sort( get.rows( k ) ) ) ) )
      cat( i, "\n", file=paste( out.dir, "/cluster.members.genes.txt", sep="" ), append=T )
    for ( i in sapply( 1:k.clust, function( k ) c( k, sort( get.cols( k ) ) ) ) )
      cat( i, "\n", file=paste( out.dir, "/cluster.members.arrays.txt", sep="" ), append=T )

    tmp <- capture.output( for ( name in ls( cmonkey.params ) ) {
      cat( name, "= " ); str( get( name, envir=cmonkey.params ), no.list=T ) } ) 
    cat( tmp, file=paste( out.dir, "/cm.params.txt", sep="" ), sep="\n", collapse="\n" )

    ## This is to make it display from the ISB central web server (but set xmlHeader=T works instead)
    ##if ( file.exists( "/sw/bin/rpl" ) || file.exists( "/usr/bin/rpl" ) )
    ##  system( paste( "rpl '<svg version' '<?xml version=\"1.0\"?><svg version' ", out.dir, "/svgs/*.svg" ) )
  }
    
  ## TODO: compress the svg's to svgz's (using gzip) and replace '.svg' with '.svgz' in all htmls.
  ## Optional because some web servers (ahem, Microsoft?) can't handle svgz's
  ## TODO: use pdftk to compress the pdfs, too
  if ( gzip ) { ##&& length( system( "which rpl", intern=T ) ) > 0 ) {
    rpl <- function( find, replace, file, ... ) {
      f <- readLines( file ); f <- gsub( find, replace, f, ... ); writeLines( f, con=file ) }
    ##if ( "svg" %in% output ) {
    system( sprintf( "gzip -v %s/svgs/*.svg", out.dir ) )
    for ( f in list.files( paste( out.dir, "/svgs", sep="" ), full=T ) ) if ( grepl( ".svg.gz", f, fixed=T ) )
      system( sprintf( "mv -v %s %s", f, sub( ".svg.gz", ".svgz", f, fixed=T ) ) )
    ##system( sprintf( "rpl .svg .svgz %s/*.html %s/htmls/*.html", out.dir, out.dir ) )
    ##for ( f in
    mc$apply( c( list.files( sprintf( "%s/htmls", out.dir ), pattern=glob2rx( "*.html" ), full=T ),
                list.files( out.dir, pattern=glob2rx( "*.html" ), full=T ) ), function( f ) {
                  cat( f, "\n" ); rpl( '.svg"', '.svgz"', f, fixed=T ) } )
    ##for ( f in list.files( out.dir, pattern=glob2rx( "*.html" ), full=T ) ) {
    ##  cat( f, "\n" ); rpl( '.svg"', '.svgz"', f, fixed=T ) }
    ##}
    if ( has.pdftk )
      mc$apply( list.files( paste( out.dir, "/pdfs", sep="" ), full=T ), function( f )
               if ( grepl( ".pdf", f, fixed=T ) ) {
                 ##system( sprintf( "ls -al %s", f ) )
                 system( sprintf( "pdftk %s output %s.tmp compress", f, f ) )
                 system( sprintf( "/bin/mv -fv %s.tmp %s", f, f ) )
                 ##system( sprintf( "ls -al %s", f ) )
               } )
  }
  if ( "rdata" %in% output ) save.cmonkey.env( file=paste( out.dir, "/cm_session.RData", sep="" ) )
  out.dir
}

environment( e$write.project ) <- e