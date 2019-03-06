devtools::install_github("eriqande/pedvis")
library(pedvis)  # load the package

# here is a simple input pedigree

stp <- simple_test_ped()
all_names <- unique(unlist(stp))
out_list <- ped2dot(stp, 
                    ShowLabelNodes = all_names, 
                    pfactorNodeStyle = "invis", 
                    pfactorEdgeStyle = "invis", 
                    ObsNodes = c(letters[1:8], 1:5, 10:12), 
                    outf = "test/testy.svg"
)


library(tidyverse)
library(xml2)

#' Add incoming and outcoming message arrows along all edges in the given SVG pedigree file 
#' 
#' The product 
#' @param infile.svg Name of the input svg file
#' @param outfile.svg Name of the output svg file
#' @export
#' 
addMsgArrow <- function(infile.svg, outfile.svg){
x <- read_xml(infile.svg)

# a way to grab class information: whether it is an edge or node
xml.class.ls <- xml_children(x)[1][[1]]  %>% xml_children(.) %>% sapply(., function(i) xml_attr(i,"class"))

edge.ls <- which(xml.class.ls=="edge")
# grabbing any edge information
path_params <- xml_children(x)[1][[1]]  %>% xml_children(.) %>% .[edge.ls] %>%
  sapply(., function(i)  xml_children(i) %>% .[2] %>% xml_attr(.,"d"))


bezier.pt.ls <- str_split(path_params, pattern="[MC, ]+") %>% sapply(., function(i) {
  len <- length(i)
  split.pt <- (len-5)*0.5 #15 / 9 for reg
  
  first.pt <- ceiling(split.pt/3)-1 
  second.pt <- split.pt - 2*first.pt
  
  first.pivot.pt <- 4 + first.pt*2
  second.pivot.pt <- ifelse(first.pivot.pt == 4, 6, 4 + second.pt*2)
  
  
  x4<- i[len-1]
  y4<- i[len]
  return(c(i[2],i[3],
           i[first.pivot.pt],i[first.pivot.pt+1],
           i[second.pivot.pt],i[second.pivot.pt+1],
           i[len-1],i[len]))
}) %>% t %>% as_tibble()

colnames(bezier.pt.ls) <- c("p1x", "p1y",
                            "p2x", "p2y",
                            "p3x", "p3y",
                            "p4x", "p4y")


time.list <- rescale.curve(bezier.pt.ls, n.breaks = 80)
new.bezier.pt <- bezier_segment_pts(bezier.pt.ls, time.list[,1],time.list[,2])
up.offset.bezier.pt <- offset.bezier(new.bezier.pt, direction = 1)  %>% 
  apply(., 1, function(i) paste0( "M",i[1],",",i[2],
                                  "C",i[3],",",i[4],
                                  " ",i[5],",",i[6],
                                  " ",i[7],",",i[8]))

down.offset.bezier.pt <- offset.bezier(new.bezier.pt, direction = -1) %>% 
  apply(., 1, function(i) paste0( "M",i[1],",",i[2],
                                  "C",i[3],",",i[4],
                                  " ",i[5],",",i[6],
                                  " ",i[7],",",i[8]))


parentNodes <- xml_children(x)[1][[1]]  


lapply(1:length(edge.ls), function(i) {
  
  edge.indx <- edge.ls[i]
  
  
  orig.id <- parentNodes %>% xml_children(.) %>% .[edge.indx] %>% xml_attr(.,"id")
  orig.title <- parentNodes %>% xml_children(.) %>% .[edge.indx] %>% xml_children(.) %>% .[1] %>% xml_contents()
  
  add.In <- read_xml(paste0('<g id="',orig.id,'a" class="edge"><title>',orig.title,'_up</title>
                   <path marker-end="url(#head)" stroke-width="1.5" fill="none" stroke="black" d="',up.offset.bezier.pt[i],'"/>
                     </g>'))
  
  xml_add_child(xml_children(x)[1][[1]]  ,
                .value=add.In,
                .where = length(xml_children(xml_children(x)[1][[1]])),
                .copy = TRUE)
  
  add.In <- read_xml(paste0('<g id="',orig.id,'b" class="edge"><title>',orig.title,'_down</title>
                   <path marker-end="url(#head)" fill="none" stroke-width="1.6" stroke="black" d="',down.offset.bezier.pt[i],'"/>
                            </g>'))
  
  xml_add_child(xml_children(x)[1][[1]] ,
                .value=add.In,
                .where = length(xml_children(xml_children(x)[1][[1]])),
                .copy = TRUE)
})

xml_add_sibling(x %>% xml_child(),
                .value=read_xml('<defs> <marker id="head" orient="auto"
                                markerWidth="2" markerHeight="4"
                                refX="1" refY="1.5">
                                <!-- triangle pointing right (+x) -->
                                <path d="M0,0 V2 L2,2 Z" fill="black"/>
                                </marker>
                                </defs>'),
                .where = "before")

write_xml(x, file=outfile.svg)
}






# find the norm slope of the point and transform the points
# dB/dt | t=0 := 3*(P1-P0)
# let's assume we will create a parallel curve separate with a distance of 1. 


# find t1 and t2 to restrict a fixed arrow length 

# rescale the curve based on t1 and t2



