If you came here looking for the cMonkey R package, you're at the right place.

<!--https://github.com/dreiss-isb/cmonkey/blob/master/cMonkey_4.9.10.tar.gz?raw=true-->

To install (in R):

```
install.packages('devtools', dep=T)
require(devtools)
install_github('cmonkey', 'dreiss-isb', subdir='cMonkey')
```

To install the data (or bigdata) packages:

```
install_github('cmonkey', 'dreiss-isb', subdir='cMonkey.data')
install_github('cmonkey', 'dreiss-isb', subdir='cMonkey.bigdata')
```

... and then go here 

<!--http://baliga.systemsbiology.net/cmonkey -->
http://dreiss-isb.github.io/cMonkey

for some additional documentation.

---

If you want to run cMonkey from raw R scripts, you've come to the right place. All of the R files that are directly part of the cMonkey package are in the Rcode directory. (Note this code has been preprocessed to remove all of my testing and experimental functions.) The *original* files that contain my un-pre-processed, experimental code is in the ORIGINAL directory. 

So, you can run cMonkey locally (without the package), after doing the following (in R):

```python
source( "ORIGINAL/cmonkey.R", chdir=T )    ## This should be all you need

## And This is only necessary if you don't have the 'progs/' dir in your current dir (and will only work if
##     you have previously installed the cMonkey package):
progs.dir <- sprintf( "%s/progs/", system.file( package="cMonkey" ) )

cmonkey( ... ) ## run as you would if you had loaded the cMonkey package

## or:

e <- cMonkey.init(...); cmonkey( e, ..., dont.init=T )

## will work as well.
```

Visit http://dreiss-isb.github.io/cMonkey for more instructions on using the package for your organism.


