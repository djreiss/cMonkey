\name{hpy}
\alias{hpy}
\docType{data}
\title{
##   ## data name/kind ... ##
}
\description{
##  ## A concise (1-5 lines) description of the dataset. ##
}
\usage{data(hpy)}
\format{
  The format is:
List of 7
 $ organism         : chr "hpy"
 $ k.clust          : num 75
 $ ratios           : num [1:819, 1:57] -1.083 -0.838 2.228 -0.822 -1.112 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:819] "HP0002" "HP0003" "HP0004" "HP0005" ...
  .. ..$ : chr [1:57] "TC1-6HR" "TC1-12HR" "TC1-18HR" "TC1-24HR" ...
 $ pp.ints          :'data.frame':	1465 obs. of  3 variables:
  ..$ protein1      : Factor w/ 244 levels "HP0001","HP0002",..: 1 1 1 2 2 2 2 3 3 3 ...
  ..$ protein2      : Factor w/ 615 levels "HP0001","HP0002",..: 24 145 320 2 470 549 565 41 54 71 ...
  ..$ combined_score: num [1:1465] 1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 ...
 $ cm.func.each.iter:function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 121 27 122 85 27 85 121 122
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ favorite.cluster :function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 123 29 123 80 29 80 123 123
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ net.weights      : Named num [1:2] 1 1
  ..- attr(*, "names")= chr [1:2] "string" "pp.ints"
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
data(hpy)
## maybe str(hpy) ; plot(hpy) ...
}
\keyword{datasets}
