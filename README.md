# pedvis

A tiny package to take a data frame that lists a pedigree in three columns 
(kid, pa, ma) and make a .dot file for a marriage node diagram of the thing.
It also will write it out as a factor graph.  And it will run dot on it if
you have dot in your system path.


This is sort of a hairball at the moment, and I would do it all differently
if I had the time to do it over, but oh well.  it works for what I want to do.  

The documentation isn't complete, but I mostly want to get this thing up on
GitHub so Clemento, my homie at the lab, can draw some pedigrees.

## Preliminaries

You should have `dot` from the `Graphviz` package installed.  If not, get it 
from http://www.graphviz.org/ for your system and install it and make sure that 
`dot` is in your PATH.

## Getting the package and some examples.