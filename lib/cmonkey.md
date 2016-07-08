# The cMonkey _R_ package
<table style="text-align:center;font-size:70%" border="1">
<p>You have reached the **updated** [_cMonkey_](http://baliga.systemsbiology.net/cmonkey) site with a new [_R_](http://www.r-project.org) package, installation and running instructions, and additional information. 
<p>If you are looking for the original _cMonkey_ site, you can [find it here](http://baliga.systemsbiology.net/drupal/content/cmonkey). The original _cMonkey_ algorithm was published with the manuscript
["Integrated biclustering of heterogeneous genome-wide datasets for the inference of global regulatory networks"](http://www.biomedcentral.com/1471-2105/7/280), by David J Reiss, Nitin S Baliga and Richard Bonneau.
<!--<p><font color="#ff0000">*NEW*</font>&nbsp The _cMonkey_ [R](http://www.r-project.org) package is [now on CRAN](http://cran.r-project.org/web/packages/cMonkey/index.html). <!--[CRAN](http://cran.r-project.org/web/packages/).-->
<!--<p>_cMonkey_ packages are [now hosted on R-forge](https://r-forge.r-project.org/projects/cmonkey/).-->
<p><font color="#ff0000">*NEW*</font>&nbsp _cMonkey_ source code is [now on Github](https://github.com/djreiss/cmonkey).
<tr>
<td>organism</td><td>n.genes</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>n.arrays</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>profile</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>network</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>motif1</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>motif2</td><td style="font-weight:bold;text-align:center;font-size:70%"><a>motif.posns</td></tr>
<td><a href="output/cmonkey_4.1.5_sce/htmls/cluster136.html">_S. cerevisiae_</a></td><td>12</td><td>445</td></td><td>![output/cmonkey_4.1.5_sce/htmls/cluster136_profile.png](output/cmonkey_4.1.5_sce/htmls/cluster136_profile.png)</img></td><td>![output/cmonkey_4.1.5_sce/htmls/cluster136_network.png](output/cmonkey_4.1.5_sce/htmls/cluster136_network.png)</img></td><td>![output/cmonkey_4.1.5_sce/htmls/cluster136_pssm1.png](output/cmonkey_4.1.5_sce/htmls/cluster136_pssm1.png)</img></td><td>![output/cmonkey_4.1.5_sce/htmls/cluster136_pssm2.png](output/cmonkey_4.1.5_sce/htmls/cluster136_pssm2.png)</img></td><td>![output/cmonkey_4.1.5_sce/htmls/cluster136_mot_posns.png](output/cmonkey_4.1.5_sce/htmls/cluster136_mot_posns.png)</img></td></tr>
<td>[_E. coli_](output/cmonkey_4.1.5_eco/htmls/cluster104.html)</td><td>24</td><td>329</td><td>![output/cmonkey_4.1.5_eco/htmls/cluster104_profile.png](output/cmonkey_4.1.5_eco/htmls/cluster104_profile.png)</img></td><td>![output/cmonkey_4.1.5_eco/htmls/cluster104_network.png](output/cmonkey_4.1.5_eco/htmls/cluster104_network.png)</img></td><td>![output/cmonkey_4.1.5_eco/htmls/cluster104_pssm1.png](output/cmonkey_4.1.5_eco/htmls/cluster104_pssm1.png)</img></td><td>![output/cmonkey_4.1.5_eco/htmls/cluster104_pssm2.png](output/cmonkey_4.1.5_eco/htmls/cluster104_pssm2.png)</img></td><td>![output/cmonkey_4.1.5_eco/htmls/cluster104_mot_posns.png](output/cmonkey_4.1.5_eco/htmls/cluster104_mot_posns.png)</img></td></tr>
<td>[_H. pylori_](output/cmonkey_4.1.5_hpy/htmls/cluster054.html)</td><td><a>12</td><td><a>36</td><td>![output/cmonkey_4.1.5_hpy/htmls/cluster054_profile.png](output/cmonkey_4.1.5_hpy/htmls/cluster054_profile.png)</img></td><td>![output/cmonkey_4.1.5_hpy/htmls/cluster054_network.png](output/cmonkey_4.1.5_hpy/htmls/cluster054_network.png)</img></td><td>![output/cmonkey_4.1.5_hpy/htmls/cluster054_pssm1.png](output/cmonkey_4.1.5_hpy/htmls/cluster054_pssm1.png)</img></td><td>![output/cmonkey_4.1.5_hpy/htmls/cluster054_pssm2.png](output/cmonkey_4.1.5_hpy/htmls/cluster054_pssm2.png)</img></td><td>![output/cmonkey_4.1.5_hpy/htmls/cluster054_mot_posns.png](output/cmonkey_4.1.5_hpy/htmls/cluster054_mot_posns.png)</img></td></tr>
<td><a href="output/cmonkey_4.1.5_hal/htmls/cluster037.html">_H. salinarum_</a></td><td>21</td><td>175</td><td>![output/cmonkey_4.1.5_hal/htmls/cluster037_profile.png](output/cmonkey_4.1.5_hal/htmls/cluster037_profile.png)</img></td><td>![output/cmonkey_4.1.5_hal/htmls/cluster037_network.png](output/cmonkey_4.1.5_hal/htmls/cluster037_network.png)</img></td><td>![output/cmonkey_4.1.5_hal/htmls/cluster037_pssm1.png](output/cmonkey_4.1.5_hal/htmls/cluster037_pssm1.png)</img></td><td>![output/cmonkey_4.1.5_hal/htmls/cluster037_pssm2.png](output/cmonkey_4.1.5_hal/htmls/cluster037_pssm2.png)</img></td><td>![output/cmonkey_4.1.5_hal/htmls/cluster037_mot_posns.png](output/cmonkey_4.1.5_hal/htmls/cluster037_mot_posns.png)</img></td></tr>
<td>[_A. thaliana](output/cmonkey_4.1.5_ath/htmls/cluster494.html)</a></td><td>19</td><td>66</td><td>![output/cmonkey_4.1.5_ath/htmls/cluster494_profile.png](output/cmonkey_4.1.5_ath/htmls/cluster494_profile.png)</img></td><td>![output/cmonkey_4.1.5_ath/htmls/cluster494_network.png](output/cmonkey_4.1.5_ath/htmls/cluster494_network.png)</img></td><td>![output/cmonkey_4.1.5_ath/htmls/cluster494_pssm1.png](output/cmonkey_4.1.5_ath/htmls/cluster494_pssm1.png)</img></td><td>![output/cmonkey_4.1.5_ath/htmls/cluster494_pssm2.png](output/cmonkey_4.1.5_ath/htmls/cluster494_pssm2.png)</img></td><td>![output/cmonkey_4.1.5_ath/htmls/cluster494_mot_posns.png](output/cmonkey_4.1.5_ath/htmls/cluster494_mot_posns.png)</img></td></tr>
</table>

* * *

# SUMMARY

<p>The latest version of <i>cMonkey_ is still under active development and has now been applied successfully to many systems, including plants (_Arabidopsis thaliana_) and mammals (_Homo sapiens_).
<p>The _cMonkey_ package will enable a user to run the integraged biclustering algorithm on their own microarray data, for their own organism of interest. During initialization, it will **automatically** download and integrate additional data for that organism, including:

<font color="#ff0000">IMPORTANT NOTE:</font> _cMonkey_ currently does not completely support _H. sapiens_. We are working on a version for that specific task. In the meantime, [see below](#hsa) for instructions on using _cMonkey_ to run _only_ on your Hsa expression data (i.e. without motifs or networks).

# INSTALLATION

<p>To run [_cMonkey_](http://baliga.systemsbiology.net/cmonkey) on your own expression data, <!--we recommend using--> you will need to use a UNIX-y operating system (e.g., Mac OS-X or Linux). _cMonkey_ <!--will run-->is not currently supported on Windows<!--, but parallelization will be more difficult to set up (see [below](#parallel) and [below](#windows))-->. For Windows users, [Cygwin](http://www.cygwin.com), [VirtualBox](https://www.virtualbox.org/), or [Amazon AWS](http://aws.amazon.com/) are inexpensive options to obtain access to a UNIX system.
<p>You will only need to do the following steps once:

<ol>
<li> Install the latest version of [_R_](http://www.r-project.org) (I am currently using version 2.11.0). It should work for versions 2.9.x and up.
<li> Install the following _R_ packages and their dependencies (all are helpful, but none are absolutely required): <tt>RCurl, doMC, igraph0, RSVGTipsDevice, and hwriter</tt> by typing in _R_: <pre>install.packages( c( 'RCurl', 'doMC', 'igraph0', 'RSVGTipsDevice', 'hwriter' ) )</pre> 

<!--<li> The latest version of _cMonkey_ may now be directly installed via the usual method, by typing in _R_: <pre>install.packages( "cMonkey" )</pre>
<ul>
<li> Alternatively, d--><li>Download and install [the latest (currently version 4.9.7) _cMonkey_ _R_ source package](cMonkey_4.9.7.tar.gz) on your system ([MD5](http://baliga.systemsbiology.net/cmonkey/md5sums.txt)). In _R_, type: 
<pre>download.file( 'http://baliga.systemsbiology.net/cmonkey/cMonkey_4.9.7.tar.gz', 'cMonkey_4.9.7.tar.gz' )
install.packages( 'cMonkey_4.9.7.tar.gz', repos=NULL, type='source' )
</pre>

<!--<li> If you plan to do motif detection (and of course you do!), download and compile/install the following:

<!--... now you're all set!

-->

[Email me](mailto:dreiss.isb@gmail.com.org) if you have any problems with any of these instructions.

* * *

# EXAMPLES: RUNNING _cMonkey_

<!--Now for the examples (including loading the _cMonkey_ package):-->

<a name="results"></a>

# EXAMPLES: EXPLORATION OF RESULTS

Once a _cMonkey_ run is complete, you may use the following to explore the results:

<!--... and more to come!-->
More functions and parameters will be documented on an as-asked-about basis.

* * *
<a name="params"></a>

# _cMonkey_ PARAMETERS

_cMonkey_ has a multitude of **internal parameters** which affect various aspects of its performance and data integration. 

For example, additional 

 may be included by a simple tweaking of the parameters. Most of these are currently undocumented but please contact me if you are interested in such possibilities.

Input parameters and data may be pre-set in one of several different ways, including:

... or any combination thereof.

**Some parameters which may be of general interest:**

More functions and parameters will be documented on an as-asked-about basis.
