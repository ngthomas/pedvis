
# this is here for testing while developing
x <- data.frame(
  Kid = c(letters[1:8], "A", "B", "C", "U", "X"),
  Pa  = c("A", "A", "A", "B", "B", "C", "D", "E", "1", "2", "3", "4", "5"),
  Ma  = c("Z", "Z", "Y", "Z", "X", "U", "U", "V", "10", "10", "11", "12", "12"),
  
  stringsAsFactors = F
)




#' write a dot file from a pedigree specified in a data frame
#' 
#' @param x the data frame holding the pedigree.  There must be at 
#' least three columns (kid, pa, and ma) which should be character
#' vectors, NOT factors and NOT integers.  Coerce them to characters
#' if they are actually integers.
#' @param pa  string giving the name of the column holding the fathers
#' in each row of the pedigree.  (Or can be thought of as "parent 1" in the 
#' case in which sex is not treated.) Default is "Pa"
#' @param ma the name of the column holding the mother in each row (or "parent 2"). Default is "Ma"
#' @param kid the name of the column the kid
#' @param ObsNodes a character vector of the identifiers of the kids, mas, and pas that have observed
#' genotype data associated with them.
#' @param ShowLabelNodes a character vector of the identifiers of the kids, mas, and pas that should have their
#' labels printed inside their nodes.
#' @param ObsStyle a named list of the extra style values that you want the observed nodes to have
#' @param outf Name of the output file (so files will be outf.dot, outf.ps, outf.pdf, etc.)
#' @export
ped2dot <- function(x, pa = "Pa", ma = "Ma", kid = "Kid", 
                    ObsNodes = character(0),
                    ShowLabelNodes = character(0),
                    ObsStyle = list(style="filled", fillcolor="gray"),
                    outf = "pedvis-ped") {
  stopifnot(is.data.frame(x)) 
  stopifnot(all(c(pa, ma, kid) %in% names(x)))
  stopifnot(is.character(x[[kid]]))
  stopifnot(is.character(x[[pa]]))
  stopifnot(is.character(x[[ma]]))
  
  # add another column which are the names of the marriage nodes
  x$mn <- paste(x[[pa]], x[[ma]], sep="x")
  
  
  # now make a big list of output strings and things
  ret <- list()
  # this is for the dot preamble.  Could modify eventually to 
  # allow function input to modify this.
  ret$pream <- c("digraph xxx {",
             "label =\"  \"",
             "ranksep=0.4;",
             "nodesep=0.4")
                          
  
  # here we set the properties of the marriage nodes
  mn_props <- lapply(x$mn, function(z) list(shape = "circle",
                                                style = "filled",
                                                label = "\"\"",
                                                height = .06,
                                                fillcolor = "black"))
  names(mn_props) <- x$mn
  
  # here we set the properties of the nodes that represent individuals in the pedigree
  unique_labels <- unique(c(x[[kid]], x[[pa]], x[[ma]]))
  node_props <- vector("list", length(unique_labels))
  names(node_props) <- unique_labels
  # set some default node values
  node_props <- lapply(node_props, function(z) list(shape = "circle", 
                                                            regular = 1,
                                                            label="\"\"",
                                                            height = .35) )
  
  # now, make any one listed as a dad a square (Note, I am not checking sex
  # consistency at this point.)
  node_props[x[[pa]]] <- lapply(node_props[x[[pa]]], function(z) {z$shape = "box"; z$height = .295; z})

  
  # now, anyone listed as observed gets the ObsStyle applied:
  node_props[ObsNodes] <- lapply(node_props[ObsNodes], function(x) c(x, ObsStyle))
  
  
  # and anyone listed as ShowLabel gets their label printed
  tmp  <- lapply(ShowLabelNodes, function(x) {y <- node_props[[x]]; y$label = x; y})
  node_props[ShowLabelNodes] <- tmp
  
  
  # now, make ObsFactor nodes.  These are little squares that sit at the same rank as 
  # any observed node.
  ofactor_nodes <- paste("of", ObsNodes, sep="_")
  
  # here we set the properties of the ofactor nodes
  ofactor_props <- lapply(ofactor_nodes, function(z) list(shape = "square",
                                            style = "filled",
                                            label = "\"\"",
                                            height = .10,
                                            fillcolor = "black"))
  names(ofactor_props) <- ofactor_nodes
  
  # and here we set the text to make the ranks of the ofactors equal to their observed nodes
  ofactor_ranks <- paste("{ rank = same; ", "\"", ofactor_nodes, "\"; ", "\"", ObsNodes, "\"; }",  sep="")
  
  
  # and finally, lets make p-factor nodes.  i.e. those nodes representing the allele frequencies
  # in the population (the prior on genotypes of founders...)
  # first, the founders are any nodes that have no parents listed.  
  # THIS IS WHERE I AM WORKING ON THINGS AT THE MOMENT.
  #Founders <- 
  
  # now store all the lines for the styles of nodes
  ret$node_props <- prop2string(node_props)
  ret$mn_props <- prop2string(mn_props)
  ret$ofactor_props <- prop2string(ofactor_props)
  
  # store the ranks
  #ret$ofactor_ranks <- ofactor_ranks
  
  # and here we store all the text specifying edges from ofactors to observed nodes
  # note that we pump their weight way up!
  ret$of2obs_nodes <- paste("\"", ObsNodes, "\"", " -> ", "\"", ofactor_nodes, "\"", "  [dir=none, weight=10000];", sep="")
  
  # here we store all the lines for specifying edges from and to the marriage nodes
  ret$mn2kid_arrows <- paste("\"", x$mn, "\"", " -> ", "\"", x[[kid]], "\"", "  [dir=none];", sep="")
  ret$par2mn_arrows <- c(unique(paste("\"", x[[pa]], "\"", " -> ", "\"", x$mn, "\"", "  [dir=none];", sep="")),
                         unique(paste("\"", x[[ma]], "\"", " -> ", "\"", x$mn, "\"", "  [dir=none];", sep="")) )
  
  
  
  ret$closer <- "}"
  
  cat(unlist(ret), sep="\n", file = paste(outf, ".dot", sep=""))
  
  # I have to modify this to use outf
  system("dot -Tps pedvis-ped.dot -o pedvis-ped.ps; epstopdf pedvis-ped.ps; open pedvis-ped.pdf")
  
  # down here return more than just ret...(later...)
  ret
  
}



#' given a named list of properties, this turns each one into a string.
#' 
#' It returns a vector of such strings.
#' 
#' @export
prop2string <- function(y) {
  sapply(names(y), function(q) {
    z <- y[[q]]
    paste("\"", q, "\"", " [",
      paste(paste(names(z), z, sep="="), collapse = ", "),
      "];", collapse="", sep="")
  })
}