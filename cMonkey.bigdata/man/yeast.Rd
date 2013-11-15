\name{yeast}
\alias{yeast}
\docType{data}
\title{
##   ## data name/kind ... ##
}
\description{
##  ## A concise (1-5 lines) description of the dataset. ##
}
\usage{data(yeast)}
\format{
  The format is:
List of 15
 $ organism             : chr "sce"
 $ is.eukaryotic        : logi TRUE
 $ k.clust              : num 450
 $ motif.upstream.search: num [1:2] -30 250
 $ motif.upstream.scan  : num [1:2] -30 500
 $ meme.cmd             : chr "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -max"| __truncated__
 $ ratios               : num [1:6160, 1:667] -0.859 -0.314 0.155 0.609 0.544 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:6160] "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
  .. ..$ : chr [1:667] "X29.alpha.factor.release.sample016" "X35.alpha.factor.release.sample015" "X37.alpha.factor.release.sample018" "X39.alpha.factor.release.sample027" ...
 $ preprocess.ratios    :function (ratios)  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 164 32 169 5 32 5 164 169
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ pp.ints              :'data.frame':	5614 obs. of  3 variables:
  ..$ protein1      : Factor w/ 2342 levels "Q0045","Q0250",..: 1 2 3 3 4 5 5 5 6 6 ...
  ..$ protein2      : Factor w/ 2731 levels "Q0045","Q0085",..: 941 785 941 785 584 2675 141 1558 2142 2426 ...
  ..$ combined_score: num [1:5614] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 ...
 $ bind.ints            :'data.frame':	68193 obs. of  3 variables:
  ..$ protein1      : Factor w/ 5497 levels "(Ace)aapa(Mcm)",..: 3050 3093 3092 3093 3049 3049 3092 3052 3047 2982 ...
  ..$ protein2      : Factor w/ 21421 levels "15S_rRNA2","18S_rRNA",..: 8329 8329 8331 8331 8329 8331 8324 8324 8316 8334 ...
  ..$ combined_score: num [1:68193] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 ...
 $ chip.chip            :'data.frame':	4358 obs. of  3 variables:
  ..$ protein1      : Factor w/ 106 levels "YBL008W","YBL021C",..: 33 33 28 28 31 98 80 37 35 12 ...
  ..$ protein2      : Factor w/ 2363 levels "YAL003W","YAL005C",..: 1 2 3 4 5 6 6 7 7 7 ...
  ..$ combined_score: num [1:4358] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 ...
 $ net.weights          : Named num [1:3] 0.25 0.5 0.1
  ..- attr(*, "names")= chr [1:3] "string" "pp.ints" "bind.ints"
 $ grouping.weights     : Named num 1
  ..- attr(*, "names")= chr "chip.chip"
 $ cm.func.each.iter    :function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 176 32 177 87 32 87 176 177
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ favorite.cluster     :function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 178 31 178 84 31 84 178 178
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
}
\details{
##  ## If necessary, more details than the __description__ above ##
}
\source{
##  ## reference to a publication or URL from which the data were obtained ##
}
\references{
##  ## possibly secondary sources and usages ##
}
\examples{
data(yeast)
## maybe str(yeast) ; plot(yeast) ...
}
\keyword{datasets}
