\name{halo}
\alias{halo}
\docType{data}
\title{
##   ## data name/kind ... ##
}
\description{
##  ## A concise (1-5 lines) description of the dataset. ##
}
\usage{data(halo)}
\format{
  The format is:
List of 11
 $ organism         : chr "hal"
 $ cog.org          : chr "Hbs"
 $ rsat.species     : chr "Halobacterium_sp"
 $ taxon.id         : num 64091
 $ k.clust          : num 250
 $ ratios           : num [1:2072, 1:268] 0.271 0.564 0.37 -1.046 -0.486 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2072] "VNG0001H" "VNG0002G" "VNG0003C" "VNG0005H" ...
  .. ..$ : chr [1:268] "NRC-1__HO_D_vs_NRC-1" "NRC-1__HO_L_vs_NRC-1" "NRC-1__LO_D_vs_NRC-1" "NRC-1__LO_L_vs_NRC-1" ...
 $ translation.tab  :'data.frame':	147 obs. of  2 variables:
  ..$ V1: Factor w/ 147 levels "VNG5003H","VNG5005H",..: 1 2 3 4 5 6 7 8 9 10 ...
  ..$ V2: Factor w/ 147 levels "VNG7002","VNG7003",..: 1 2 3 4 5 6 7 8 9 10 ...
 $ bat.clust.genes  : chr [1:30] "VNG0654C" "VNG0656H" "VNG0828H" "VNG0829G" ...
 $ cm.func.each.iter:function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 110 31 111 81 31 81 110 111
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ favorite.cluster :function ()  
  ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 112 30 112 88 30 88 112 112
  .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x103a376e8> 
 $ net.weights      : Named num 1
  ..- attr(*, "names")= chr "string"
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
data(halo)
## maybe str(halo) ; plot(halo) ...
}
\keyword{datasets}
