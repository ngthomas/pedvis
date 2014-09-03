# pedvis

A tiny package to take a data frame that lists a pedigree in three columns 
(kid, pa, ma) and make a .dot file for a marriage node diagram of the thing.
It also will write it out as a factor graph.  And it will run dot on it if
you have dot in your system path.


This is sort of a hairball at the moment, and I would do it all differently
if I had the time to do it over, but oh well.  it works for what I want to do
at the moment.

The documentation isn't complete, but I mostly want to get this thing up on
GitHub so Clemento, my homie at the lab, can draw some pedigrees.

## Preliminaries

1. You should have `dot` from the `Graphviz` package installed.  If not, get it 
from http://www.graphviz.org/ for your system and install it and make sure that 
`dot` is in your PATH.

2. You should have `epstodpdf` on your system.  i.e. you ought to have ImageMagick
or something of the sort.

3. The program will try to open a PDF file with a system call to the "open" command,
which might be a Mac thing, and not portable.


## Getting the package and some examples.

The easiest way to get and install this package is with the devtools package from 
within R:
```r
devtools::install_github("pedvis", "eriqande")
```

Then, if that was successful, here you can see what the input is all about and do some examples
```r
library(pedvis)  # load the package

# here is a simple input pedigree
simple_test_ped()
```
Note that the names of the individuals have to be 
*strings* NOT factors or integers.  So, use `stringsAsFactors = FALSE`
if reading in a data frame with `read.table`.

Here we draw that simple pedigree as a marriage node diagram and we put labels on everyone and
we make some of the p-node factors invisible (that should sound cryptic...).
And we make a few of the individuals shaded by telling the function we have observed data on them.
```r
stp <- simple_test_ped()
all_names <- unique(unlist(stp))
out_list <- ped2dot(stp, 
                    ShowLabelNodes = all_names, 
                    pfactorNodeStyle = "invis", 
                    pfactorEdgeStyle = "invis", 
                    ObsNodes = c(letters[1:8], 1:5, 10:12), 
                    outf = "testy"
                    )
```

Note that running this function will create a file called testy.dot in the current working directory.  It will 
then run dot on it and make testy.ps, and then convert that to a pdf with epstopdf.  Then it will open that
PDF file with the default PDF viewer.  It also makes an svg version so you can edit it with Inkscape.