# contains all the necessary fns to perform bezier curve manipulation

#' return x and y point given all of the bezier curve points and time points
#' @param bezier.pt tibble data with all four bezier x y points
#' @param t time from 0 to 1
#' 
# bezier.pt.ls <- dplyr::tibble(p1x = 514.7528, p1y = -261.7987,
# p2x = 561.1639, p2y = -235.9111,
# p3x = 653.0751, p3y = -184.6441,
# p4x = 679.3386, p4y = -169.9946)
# get.curve.pt(bezier.pt.ls, 1)

get.curve.pt <- function (bezier.pt, t){
  
  p.x <- bezier.pt %>%
    dplyr::select(p1x, p2x, p3x, p4x) %>%
    as.matrix() %>%
    as.numeric() %>%
    matrix(., ncol=4, byrow=F)
  
  p.y <-
    bezier.pt %>%
    dplyr::select(p1y,p2y, p3y, p4y) %>%
    as.matrix() %>%
    as.numeric() %>%
    matrix(., ncol=4, byrow=F)
  
  t.set <- c(1, t**2, t**3, t**4)
  prefix.set <- matrix(c(1,0,0,0,
                         -3,3,0,0,
                         3,-6,3,0,
                         -1,3,-3,1), ncol=4, nrow=4, byrow = T)
  
  rbind(t.set %*% prefix.set %*% t(p.x),
        t.set %*% prefix.set %*% t(p.y)) %>%
    t() %>%
    as.tibble() %>%
    dplyr::rename(x="V1", y="V2")
}



#'  returns a new beginning and ending time point (t_begin and t_end) to trim out a curve path to
#' @param bezier.pt a n by 8 data table contains original bezier points
#' @param length.arrow length of the unit arrow
#' @param n.breaks number of break points to approximate arc length

rescale.curve <- function(bezier.pt, length.arrow = 10, n.breaks = 20) {
  
  n.curve <- dim(bezier.pt)[1]
  t.breakpt <- seq(0,1, length.out =n.breaks)
  
  xy.timept <- lapply(t.breakpt, function(i) get.curve.pt(bezier.pt, i)) %>%
    bind_rows()
  
  n.entries <- dim(xy.timept)[1]
  
  # approximate length of the curve
  seg.curve <- (xy.timept[(n.curve+1):n.entries,]-xy.timept[1:(n.entries-n.curve),])**2 %>%
    rowSums() %>%
    sqrt(.) %>%
    matrix(.,nrow=n.curve)
  
  len.curve <- rowSums(seg.curve)
  
  # get t_begin and t_end breakpt
  ideal.intv <- cbind(len.curve/2 - length.arrow/2, len.curve/2 + length.arrow/2)
  # get cumulative length
  cum.len <- apply(seg.curve, 1, cumsum) %>% t
  
  left.indx <- apply(cum.len >= ideal.intv[,1],1, function(x) min(which(x))) -1
  right.indx <- apply(cum.len <= ideal.intv[,2],1, function(x) max(which(x))) + 1
  
  cbind((t.breakpt[-1])[left.indx],
        (t.breakpt[-1])[right.indx])
}


#' return bezier control pts (for either x or y) for internal segment or original curve
#' @param t0 how far along from 0 to 1 to start the segment
#' @param t1 where to end the segment
#' @param bezier.pt a n by 8 data table contains original bezier points
bezier_segment_pts <- function(bezier.pt, t0, t1) {
  #if(t0 > t1) stop("t0 must be less than t1")
  u0 <- 1 - t0
  u1 <- 1 - t1
  
  p1 <- bezier.pt %>% dplyr::select(p1x, p1y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p2 <- bezier.pt %>% dplyr::select(p2x, p2y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p3 <- bezier.pt %>% dplyr::select(p3x, p3y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p4 <- bezier.pt %>% dplyr::select(p4x, p4y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  
  # compute the new ones:
  P1 <- (u0*u0*u0*p1 + (t0*u0*u0 + u0*t0*u0 + u0*u0*t0)*p2 + (t0*t0*u0 + u0*t0*t0 + t0*u0*t0)*p3 + t0*t0*t0*p4)
  P2 <- (u0*u0*u1*p1 + (t0*u0*u1 + u0*t0*u1 + u0*u0*t1)*p2 + (t0*t0*u1 + u0*t0*t1 + t0*u0*t1)*p3 + t0*t0*t1*p4)
  P3 <- (u0*u1*u1*p1 + (t0*u1*u1 + u0*t1*u1 + u0*u1*t1)*p2 + (t0*t1*u1 + u0*t1*t1 + t0*u1*t1)*p3 + t0*t1*t1*p4)
  P4 <- (u1*u1*u1*p1 + (t1*u1*u1 + u1*t1*u1 + u1*u1*t1)*p2 + (t1*t1*u1 + u1*t1*t1 + t1*u1*t1)*p3 + t1*t1*t1*p4)
  
  cbind(P1, P2, P3, P4)
}

# offset bezier curve using Tiller and Hanson method
## https://math.stackexchange.com/questions/465782/control-points-of-offset-bezier-curve/467038#467038


#' returns a dataframe of offset curve (for both direction with id)
#' @param bezier.pt a n by 8 data table contains original bezier points
#' @param offset.distance unit distance offsetting limits
#' @param direction on top or below of the appointed edge: assigned downstream pointing arrow on top
#' 
offset.bezier <- function(bezier.pt, offset.distance = 10, direction = 1) {
  
  # compare p0 to p1
  dy <- bezier.pt[,4] - bezier.pt[,2]
  dx <- bezier.pt[,3] - bezier.pt[,1]
  
  m <- dx/dy
  
    
  new.x1 <- offset.distance*direction/(sqrt(1+m^2))
  new.y1 <- abs(m)*new.x1 
  
  #if(direction == -1 & m <= 0 ) new.x1 <- new.x1 * -1
  intervene <- which(m != 0 & m > 0 ) 
  new.y1[intervene] <- new.y1[intervene] * -1
  
  # compare p3 to p2
  dy <- bezier.pt[,8] - bezier.pt[,6]
  dx <- bezier.pt[,7] - bezier.pt[,5]
  
  m <- dx/dy 
  
  new.x2 <- offset.distance*direction/(sqrt(1+m^2))
  new.y2 <- abs(m)*new.x2 
  
  intervene <- which(m != 0 & m > 0 ) 
  new.y2[intervene] <- new.y2[intervene] * -1
  
  #if(direction == -1 &  m <= 0 ) new.x2 <- new.x2 * -1
  #if(direction == -1 & m != 0 & m > 0 ) new.y2 <- new.y2 * -1
  
  recoord <- cbind(bezier.pt[,1] + new.x1, bezier.pt[,2] + new.y1,
    bezier.pt[,3] + new.x1, bezier.pt[,4] + new.y1,
    bezier.pt[,5] + new.x2, bezier.pt[,6] + new.y2,
    bezier.pt[,7] + new.x2, bezier.pt[,8] + new.y2)
  
  if(direction == -1) return(recoord[,c(7,8,5,6,3,4,1,2),drop=FALSE])
  
  return(recoord)
}

#' Add incoming and outcoming message arrows along all edges in the given SVG pedigree file 
#' 
#' The output svg figure contains all added message arrows
#' @param infile.svg String. Directory path for the input svg file. Required
#' @param outfile.svg String. Directory path for the output svg file. Required
#' @param edge.sel  String or list. Selected which edge to labeled. if left emptied (as default), all edges are marked with bidirectional arrows. The edge labeling uses individual id and its relationship to its connected m-node (m-node's id uses "Pa's id"x"mom's id"; e.g. 8x9). The edge label between m-node and its associated parent has the parent's id comes first following with m-node label. e.g "8->8x9". The edge label between m-node and its associated offspring has the m-node id following the offspring id. e.g. "8x9->10".
#' @param dir.sel String or list. Specifies the direction of the arrows correspond to edge.sel. 1=pointing toward first mentioned node of the edge label. 2=point toward the second node of the edge label. 3=draws both arrows, where each points to each othe way
#' @param color.sel String. Specifies the color (limited to 'black', 'blue' or 'red' - case sensitive) for the arrow. If unspecified, black will be used. Else, the color list should be the same length as edge.sel
#' @param len.arrow length of the arrow
#' @param stroke.width width of the arrow. 1.6 as default.
#' @export
#' 
addMsgArrow <- function(infile.svg, outfile.svg,
                        edge.sel = "",
                        dir.sel = "",
                        color.sel = "",
                        len.arrow = 10,
                        stroke.width = 1.6) {
  
  if(length(edge.sel) != length(dir.sel)) {
    stop("The edge list (edge.sel) and direction list (dir.sel) need to be match in length.")
  }
  
  if(color.sel[1] != "" ) {
    
    if (sum(!color.sel %in% c("black","red","blue")) >0) stop ("The color list (color.sel) has to be either black, red or blue.")
    if(length(color.sel) != length(dir.sel)) stop("The edge list (edge.sel) and color list (color.sel) need to be match in length.")
  }
  
  x <- xml2::read_xml(infile.svg)
  
  # a way to grab class information: whether it is an edge or node
  xml.class.ls <- xml2::xml_children(x)[1][[1]] %>% xml2::xml_children(.) %>% sapply(., function(i) xml2::xml_attr(i,"class"))
  
  # grabbing any edge information
  edge.ls <- which(xml.class.ls=="edge")
  path_params <- xml2::xml_children(x)[1][[1]] %>% xml2::xml_children(.) %>% .[edge.ls] %>%
    sapply(., function(i)  xml2::xml_children(i) %>% .[2] %>% xml2::xml_attr(.,"d"))
  
  path.title <- xml2::xml_children(x)[1][[1]] %>% xml2::xml_children(.) %>% .[edge.ls] %>% sapply(., function(i) xml2::xml_children(i) %>% xml2::xml_text(.,"title")) %>% .[1,]
  
  # clean out any edge explicitly defined
  if(edge.sel[1] != "") {
    edge.pickup <- sapply(path.title, function(i) i %in% edge.sel)
    path_params <- path_params[edge.pickup]
    edge.ls <- edge.ls[edge.pickup]
  }
  
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
  }) %>% t %>% dplyr::as_tibble()
  
  colnames(bezier.pt.ls) <- c("p1x", "p1y",
                              "p2x", "p2y",
                              "p3x", "p3y",
                              "p4x", "p4y")
  
  time.list <- rescale.curve(bezier.pt.ls, length.arrow=len.arrow,  n.breaks = 80)
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
  
  parentNodes <- xml2::xml_children(x)[1][[1]]  
  
  lapply(1:length(edge.ls), function(i) {
    
    edge.indx <- edge.ls[i]
    
    orig.id <- parentNodes %>% xml2::xml_children(.) %>% .[edge.indx] %>% xml2::xml_attr(.,"id")
    orig.title <- parentNodes %>% xml2::xml_children(.) %>% .[edge.indx] %>% xml2::xml_children(.) %>% .[1] %>% xml2::xml_text(.,"title")
    
    dir.spec <- 3
    color.def <- "black"
    
    
    if(dir.sel[1] != "") {
      dir.spec <- dir.sel[which(orig.title == edge.sel)[1]]
    }
    
    if(color.sel[1] != ""){
      color.def <- color.sel[which(orig.title == edge.sel)[1]]
    }
    
    if (dir.spec == 2 || dir.spec == 3 ) {
      
      add.In <- xml2::read_xml(paste0('<g id="',orig.id,'a" class="edge"><title>',orig.title,'_up</title>
                                      <path marker-end="url(#head',color.def,')" stroke-width="',stroke.width,'" fill="none" stroke="',color.def,'" d="',up.offset.bezier.pt[i],'"/>
                                      </g>'))
      
      xml2::xml_add_child(xml2::xml_children(x)[1][[1]]  ,
                          .value=add.In,
                          .where = length(xml2::xml_children(xml2::xml_children(x)[1][[1]])),
                          .copy = TRUE)
    }
    
    if (dir.spec == 1 || dir.spec == 3 ) {
      add.In <- xml2::read_xml(paste0('<g id="',orig.id,'b" class="edge"><title>',orig.title,'_down</title>
                                      <path marker-end="url(#head',color.def,')" fill="none" stroke-width="',stroke.width,'" stroke="',color.def,'" d="',down.offset.bezier.pt[i],'"/>
                                      </g>'))
      
      xml2::xml_add_child(xml2::xml_children(x)[1][[1]] ,
                          .value= add.In,
                          .where = length(xml2::xml_children(xml2::xml_children(x)[1][[1]])),
                          .copy = TRUE)
    }
  })
  
  xml2::xml_add_sibling(x %>% xml2::xml_child(),
                        .value=xml2::read_xml('<defs> <marker id="headblue" orient="auto"
                                              markerWidth="2" markerHeight="4"
                                              refX="1" refY="1.5">
                                              <!-- triangle pointing right (+x) -->
                                              <path d="M0,0 V2 L2,2 Z" fill="blue"/>
                                              </marker>
                                              <marker id="headred" orient="auto"
                                              markerWidth="2" markerHeight="4"
                                              refX="1" refY="1.5">
                                              <!-- triangle pointing right (+x) -->
                                              <path d="M0,0 V2 L2,2 Z" fill="red"/>
                                              </marker>
                                              <marker id="headblack" orient="auto"
                                              markerWidth="2" markerHeight="4"
                                              refX="1" refY="1.5">
                                              <!-- triangle pointing right (+x) -->
                                              <path d="M0,0 V2 L2,2 Z" fill="black"/>
                                              </marker>
                                              </defs>'),
                        .where = "before")
  
  xml2::write_xml(x, file=outfile.svg)
}

#' Assigning the time step for which edge messages are allowed to propagate 
#' 
#' If the pedigree has a cycle, the process stops once it couldn't find any possible llegal updates
#' 
#' @param marriage.tbl data frame. 3 column marriage dataframe contains kid's id, parent 0's id, parent 1's id
#' @param p.node.on boolean. account for p-node edge propagation. Default: false
#' @param g.node.on boolean. account for g-node edge propagation. Default: false
#' @export

CountFGSteps <- function(marriage.tbl, p.node.on=FALSE, g.node.on = FALSE) {
  
  marriage.tbl <- dplyr::tbl_df(marriage.tbl)
  
  node.tbl <- dplyr::bind_rows(
    marriage.tbl %>% dplyr::group_by(Pa, Ma) %>% dplyr::summarise(indiv = Pa[1], m_node = paste0(Pa[1], "x", Ma[1]), is_parent = 1),
    marriage.tbl %>% dplyr::group_by(Pa, Ma) %>% dplyr::summarise(indiv = Ma[1], m_node = paste0(Pa[1], "x", Ma[1]), is_parent = 1),
    marriage.tbl %>% dplyr::group_by(Pa, Ma) %>% dplyr::mutate(indiv = Kid, m_node = paste0(Pa, "x", Ma), is_parent = 0) %>%
      dplyr::select(-Kid)
  ) %>% dplyr::ungroup() %>% dplyr::select(-Pa, -Ma) %>%
    dplyr::mutate() %>%
    dplyr::group_by(indiv) %>%
    dplyr::mutate(n.mNode = n(),
           msg.in.indiv = 0,
           msg.out.indiv = 0,
           step.in = 0,
           step.out = 0) %>%
    dplyr::ungroup() %>% dplyr::group_by(m_node) %>% dplyr::mutate(n.indiv = n())
  
  n.indiv <- length(unique(node.tbl$indiv))
  
  ## need to check for cyclic structure
  
  last.updates <- -1
  while(sum(c(node.tbl$step.in,node.tbl$step.out) ==0) > 0 && sum(c(node.tbl$step.in,node.tbl$step.out) ==0) != last.updates) {
    
    last.updates <- sum(c(node.tbl$step.in,node.tbl$step.out) ==0)
    #passing message from var to m-node
    step.ct <- max(c(node.tbl$step.in,node.tbl$step.out))+1
    
    node.tbl <- node.tbl %>%
      dplyr::group_by(indiv) %>%
      dplyr::mutate(is.ready = (n.mNode-1 == sum(msg.in.indiv)-msg.in.indiv)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(step.out = ifelse(step.out == 0 & is.ready, step.ct, step.out),
             msg.out.indiv = ifelse(is.ready, 1, 0))
    
    
    step.ct <- max(c(node.tbl$step.in,node.tbl$step.out))+1
    
    node.tbl <- node.tbl %>%
      dplyr::group_by(m_node) %>%
      dplyr::mutate(is.ready = (n.indiv-1 == sum(msg.out.indiv)-msg.out.indiv)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(step.in = ifelse(step.in == 0 & is.ready, step.ct, step.in),
             msg.in.indiv = ifelse(is.ready, 1, 0))
  }
  
  ## works on p-node propagation
  is.complete <- sum(c(node.tbl$step.in,node.tbl$step.out) ==0) == 0 
  n.steps <- max(c(node.tbl$step.in,node.tbl$step.out))+1
  if(p.node.on) {
    founder.ls <- setdiff(c(node.tbl$indiv),c(marriage.tbl$Kid))
    
    add.p.node <- dplyr::tibble(step = rep(1,length(founder.ls)),
                                edge.lab = paste0("pf_",founder.ls, "->", founder.ls),
                                dir.sel =2,
                                color.sel = "black")
    
    if(is.complete) {
      add.p.node <- dplyr::bind_rows(add.p.node,
                                     dplyr::tibble(step = rep(n.steps+1,length(founder.ls)),
                                  edge.lab = paste0("pf_",founder.ls, "->", founder.ls),
                                  dir.sel =1,
                                  color.sel = "black"))
    }
    
  }
  
  if (g.node.on) {
    indiv.ls <- unique(node.tbl$indiv)
    add.g.node <- dplyr::tibble(step = rep(1,n.indiv),
                                edge.lab = paste0(indiv.ls, "->of_", indiv.ls),
                                dir.sel =1,
                                color.sel = "black")
    
    if(is.complete) {
      add.g.node <- dplyr::bind_rows(add.g.node, dplyr::tibble(step = rep(n.steps+1, n.indiv),
                                  edge.lab = paste0(indiv.ls, "->of_", indiv.ls),
                                  dir.sel =2,
                                  color.sel = "black"))
    }
    
  }
  
  
  timed.tbl <- tidyr::gather(node.tbl, "msg.dir", "step", 7:8) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(edge.lab = ifelse(is_parent ==1,
                             paste0(indiv,"->",m_node),
                             paste0(m_node,"->",indiv)),
           dir.sel = ifelse((is_parent ==1 & msg.dir == "step.out") |
                              (is_parent ==0 & msg.dir == "step.in"),2,1),
           color.sel = ifelse(msg.dir == "step.out", "blue", "red")) %>%
    dplyr::select(step, edge.lab, dir.sel, color.sel)
  
  
  if(p.node.on || g.node.on) {
    timed.tbl <- timed.tbl %>%
      dplyr::mutate(step = ifelse(step >0, step + 1, step ))
    
    if (p.node.on) timed.tbl <- dplyr::bind_rows(timed.tbl, add.p.node)
    if (g.node.on) timed.tbl <- dplyr::bind_rows(timed.tbl, add.g.node)
  }
  
  return(timed.tbl)
}

