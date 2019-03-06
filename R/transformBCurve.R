# contains all the necessary fns to perform bezier curve manipulation

#' return x and y point given all of the bezier curve points and time points
#' @param bezier.pt tibble data with all four bezier x y points
#' @param t time from 0 to 1
#' @example get.curve.pt(bezier.pt.ls, 1)

get.curve.pt <- function (bezier.pt, t){
  
  p.x <- bezier.pt %>%
    select(p1x, p2x, p3x, p4x) %>%
    as.matrix() %>%
    as.numeric() %>%
    matrix(., ncol=4, byrow=F)
  
  p.y <-
    bezier.pt %>%
    select(p1y,p2y, p3y, p4y) %>%
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
    rename(x="V1", y="V2")
}



#  returns a new beginning and ending time point (t_begin and t_end) to trim out a curve path to
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
  
  p1 <- bezier.pt %>% select(p1x, p1y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p2 <- bezier.pt %>% select(p2x, p2y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p3 <- bezier.pt %>% select(p3x, p3y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  p4 <- bezier.pt %>% select(p4x, p4y) %>% as.matrix() %>% as.numeric() %>% matrix(., ncol=2, byrow=F)
  
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
  
  if(direction == -1) return(recoord[,c(7,8,5,6,3,4,1,2)])
  
  return(recoord)
}

