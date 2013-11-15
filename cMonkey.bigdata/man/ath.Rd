\name{ath}
\alias{ath}
\docType{data}
\title{
##   ## data name/kind ... ##
}
\description{
##  ## A concise (1-5 lines) description of the dataset. ##
}
\usage{data(ath)}
\format{
  The format is:
List of 9
 $ organism             : chr "ath"
 $ rsat.species         : chr "Arabidopsis_thaliana"
 $ is.eukaryotic        : logi TRUE
 $ k.clust              : num 600
 $ motif.upstream.search: num [1:2] -20 500
 $ motif.upstream.scan  : num [1:2] -30 1000
 $ meme.cmd             : chr "./progs/meme $fname -bfile $bgFname -psp $pspFname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %1$d -evt 1e9 -minw 6 -max"| __truncated__
 $ ratios               : num [1:22591, 1:108] -0.1586 0.3295 -0.0196 0.1403 0.3878 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22591] "AT1G01010" "AT1G01030" "AT1G01040" "AT1G01050" ...
  .. ..$ : chr [1:108] "Cold_05_S" "Cold_1_S" "Cold_3_S" "Cold_6_S" ...
 $ preprocess.ratios    :function (ratios)  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 189 30 194 5 30 5 189 194
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
data(ath)
## maybe str(ath) ; plot(ath) ...
}
\keyword{datasets}
