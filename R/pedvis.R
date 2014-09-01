
#' A simple pedigree to use for testing and demonstration
#' 
#' This returns a data frame suitable for passing to ped2dot
#' @param No parameters for now.
#' @export
simple_test_ped <- function() {
  data.frame(
    Kid = c(letters[1:8], "A", "B", "C", "U", "X"),
    Pa  = c("A", "A", "A", "B", "B", "C", "D", "E", "1", "2", "3", "4", "5"),
    Ma  = c("Z", "Z", "Y", "Z", "X", "U", "U", "V", "10", "10", "11", "12", "12"),
    stringsAsFactors = F
  )
}




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
#' @param ProngNodes a character vector of the identifiers of the kids, mas, and pas that should be considered
#' prong nodes, and therefore be printed a little differently.  Additionally, the edges incident to a prong node
#' will be drawn with ProngEdgeStyle.
#' @param ObsStyle a named list of the extra style values that you want the observed nodes to have
#' @param ProngStyle a named list of the extra style values that you want the prong nodes to have
#' @param ProngEdgeStyle a named list of the extra style values that you want the edges connected to prong nodes to have
#' @param outf Name of the output file without any suffixes (so files will be outf.dot, outf.ps, outf.pdf, etc.)
#' @examples
#' # simple pedigree as marriage node diagram
#' junk1 <- ped2dot(simple_test_ped(), outf = "ped2dot_ex1", pfactorNodeStyle = "invis", pfactorEdgeStyle = "invis")
#' 
#' # color in some of the individuals as having observed data
#' junk2 <- ped2dot(simple_test_ped(), ObsNodes = c("a", "b", "f", "C", "D", "Z", "U", "10"), outf = "ped2dot_ex2")
#' 
#' # draw the same thing, but make the allele freq factors invisible
#' junk3 <- ped2dot(simple_test_ped(), ObsNodes = c("a", "b", "f", "C", "D", "Z", "U", "10"), 
#'         pfactorNodeStyle = "invis", pfactorEdgeStyle = "invis",
#'         outf = "ped2dot_ex3")
#'  
#' # like example 2 but with labels on the variable nodes (big circles and squares)
#' junk4 <- ped2dot(simple_test_ped(), ObsNodes = c("a", "b", "f", "C", "D", "Z", "U", "10"), 
#'          outf = "ped2dot_ex4", ShowLabelNodes = unique(unlist(simple_test_ped())))
#'          
#' # put some prong nodes on there
#' junk5 <- ped2dot(simple_test_ped(), ObsNodes = c("a", "b", "f", "C", "D", "Z", "U", "10"), 
#'              outf = "ped2dot_ex5", ShowLabelNodes = unique(unlist(simple_test_ped())), ProngNodes = c("c", "d", "e", "E", "V", "3", "11", "12", "4"))
#' @export
ped2dot <- function(x, pa = "Pa", ma = "Ma", kid = "Kid", 
                    ObsNodes = character(0),
                    ShowLabelNodes = character(0),
                    ProngNodes = character(0),
                    ObsStyle = list(style="filled", fillcolor="gray"),
                    ProngStyle = list(style="dashed"),
                    ProngEdgeStyle = list(style="dashed"),
                    outf = "pedvis-ped",
                    pfactorNodeStyle = "filled",
                    pfactorEdgeStyle = "solid",
                    Draw_O_factors = FALSE,
                    RankSep = 1.0,
                    NodeSep = 1.0
                    ) {
  stopifnot(is.data.frame(x)) 
  stopifnot(all(c(pa, ma, kid) %in% names(x)))
  stopifnot(is.character(x[[kid]]))
  stopifnot(is.character(x[[pa]]))
  stopifnot(is.character(x[[ma]]))
  
  po_intersect <- intersect(ObsStyle, ProngStyle)
  if(length(po_intersect) > 0) {
    po_inter_strs <- paste(po_intersect, collapse = ", ")
    stop(paste("Nodes cannot be ObsNodes and ProngNodes at the same time.  Problems with:", po_inter_strs))
  }
  
  # add another column which are the names of the marriage nodes
  x$mn <- paste(x[[pa]], x[[ma]], sep="x")
  
  
  # now make a big list of output strings and things
  ret <- list()
  # this is for the dot preamble.  Could modify eventually to 
  # allow function input to modify this.
  ret$pream <- c("digraph xxx {",
             "label =\"  \"",
             paste("ranksep=", RankSep, ";", sep=""),
             paste("nodesep=", NodeSep, sep=""),
             "compress=false")
                          
  
  # here we set the properties of the marriage nodes
  mn_props <- lapply(x$mn, function(z) list(shape = "circle",
                                                style = "filled",
                                                label = "\"\"",
                                                height = 0.225,
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
                                                            height = 0.86,
                                                            fixedsize = "true") )
  
  # now, make any one listed as a dad a square (Note, I am not checking sex
  # consistency at this point.)
  node_props[x[[pa]]] <- lapply(node_props[x[[pa]]], function(z) {z$shape = "box"; z$height = .73; z})

  
  # now, anyone listed as observed gets the ObsStyle applied:
  node_props[ObsNodes] <- lapply(node_props[ObsNodes], function(x) c(x, ObsStyle))
  
  # anyone listed as a prong gets ProngStyle applied
  node_props[ProngNodes] <- lapply(node_props[ProngNodes], function(x) c(x, ProngStyle))
  
  # and anyone listed as ShowLabel gets their label printed
  tmp  <- lapply(ShowLabelNodes, function(x) {y <- node_props[[x]]; y$label = x; y})
  node_props[ShowLabelNodes] <- tmp
  
  
  # now, make ObsFactor nodes.  These are little squares that sit at one rank below  
  # any observed node. 
  ofactor_nodes <- paste("of", ObsNodes, sep="_")
  
  # here we set the properties of the ofactor nodes
  ofactor_props <- lapply(ofactor_nodes, function(z) list(shape = "square",
                                            style = "filled",
                                            label = "\"\"",
                                            height = 0.225,
                                            fillcolor = "black"))
  names(ofactor_props) <- ofactor_nodes
   
  
  # and finally, lets make p-factor nodes.  i.e. those nodes representing the allele frequencies
  # in the population (the prior on genotypes of founders...)
  # first, the founders are any nodes that have no parents listed.  
  mapa <- unique(c(x[[ma]], x[[pa]]))
  Founders <- mapa[!(mapa %in% x[[kid]])]
  pfactor_nodes <- paste("pf", Founders, sep="_")
  pfactor_props <- lapply(pfactor_nodes, function(z) list(shape = "diamond",
                                                          style = pfactorNodeStyle,
                                                          regular = 1,
                                                          label = "\"\"",
                                                          height = 0.25,
                                                          fillcolor = "black"))
  names(pfactor_props) <- pfactor_nodes
  
  # now store all the commands in the dot file for the styles of nodes
  ret$node_props <- prop2string(node_props)
  ret$mn_props <- prop2string(mn_props)
  if(length(ObsNodes) > 0 && Draw_O_factors == TRUE) {ret$ofactor_props <- prop2string(ofactor_props)}
  ret$pfactor_props <- prop2string(pfactor_props)
  
  
  # and here we store all the text specifying edges from ofactors to observed nodes
  # note that we pump their weight way up!
  if(length(ObsNodes) > 0 && Draw_O_factors == TRUE) {
    obsf_edgestyle <- rep("", length(ObsNodes))
    obsf_edgestyle[ObsNodes %in% ProngNodes] <- "style = dashed,";  # this lets us prong out the obs-nodes if need be
    ret$of2obs_nodes <- paste("\"", ObsNodes, "\"", " -> ", "\"", ofactor_nodes, "\"", "  [", obsf_edgestyle, " dir=none, weight=10000];", sep="")}
  
  # and here is text specifying the edges connecting the pfactor nodes to the founders
  # default edge style
  pf2founders_style <- rep(paste("  [dir=none, weight=10000, style=", pfactorEdgeStyle,"]"), length=length(Founders))
  # now change those for founders that are prong nodes
  pf2founders_style[ Founders %in% ProngNodes] <- paste(" [dir = none, weight=10000,", list2keyval(ProngEdgeStyle), "]")
  ret$pf2founders <- paste("\"", pfactor_nodes, "\"", " -> ", "\"", Founders, "\"", pf2founders_style,";", sep="")
  
  
  # here we store all the lines for specifying edges from and to the marriage nodes
  x$mn2kid_style = " [dir = none];"  # set the default style for this edge
  x$par2mn_style = " [dir = none];"  # set the default style for this kind of edge
  
  # now, any edge connected to a prong should get additional styling
  x$mn2kid_style[x[[kid]] %in% ProngNodes] <- paste(" [dir = none,", list2keyval(ProngEdgeStyle), "]")
  x$par2mn_style[ (x[[ma]] %in% ProngNodes) | (x[[pa]] %in% ProngNodes)] <- paste(" [dir = none,", list2keyval(ProngEdgeStyle), "]")
  
  
  
  
  ret$mn2kid_arrows <- paste("\"", x$mn, "\"", " -> ", "\"", x[[kid]], "\"", x$mn2kid_style, sep="")
  ret$par2mn_arrows <- c(unique(paste("\"", x[[pa]], "\"", " -> ", "\"", x$mn, "\"", x$par2mn_style, sep="")),
                         unique(paste("\"", x[[ma]], "\"", " -> ", "\"", x$mn, "\"", x$par2mn_style, sep="")) )
  
  
  
  ret$closer <- "}"
  
  cat(unlist(ret), sep="\n", file = paste(outf, ".dot", sep=""))
  
  # I have to modify this to use outf
  DOTFile <- paste(outf, "dot", sep = ".")
  PSFile <- paste(outf, "ps", sep = ".")
  PDFFile <- paste(outf, "pdf", sep = ".")
  SVGFile <- paste(outf, "svg", sep = ".")
  SVG_CALL <- paste("dot -Tsvg", DOTFile, "-o",  SVGFile)
  CALL <- paste("dot -Tps", DOTFile, "-o",  PSFile, ";",  "epstopdf",  PSFile, ";",  "open",  PDFFile)
  
  system(SVG_CALL)
  system(CALL)
  
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


#' given a named list, this give a string with comma-separated keyval pairs
#'
#' @param z the named list
#' @export
#' @examples
#' list2keyval(list(color="red", width="fat", number=2))
list2keyval <- function(z) {
  paste(sapply(names(z), function(y) paste(y,z[y], sep="=")), collapse = ", ")
}