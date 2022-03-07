#' Generate numbers from a Halton Sequence with a random start
#'
#' For efficiency, this function can generate points along a random start Halton Sequence for
#' predefined Halton 
#'
#' @param n Number of points required
#' @param seeds Random starting point in each dimension
#' @param bases Co-prime base for the Halton Sequence
#' @param boxes Halton boxes that points are required to be generated in
#' @param J Defines the Halton frame, and relates to the number of boxes used.
#'
#' @return Matrix with the columns, order of point, x in [0,1) and y in [0,1)
#'
#' @examples
#' \dontrun{
#' # First 10 points in the Halton Sequence for base 2,3
#' pts <- RSHalton(n = 10)
#' # First 10 points in the Halton Sequence for base 2,3 with starting point at the 15th and 22nd index.
#' pts <- RSHalton(n = 10, seeds = c(14, 21))
#' }
#' @export
RSHalton <- function(n = 10, seeds = c(0,0),bases = c(2,3), boxes = 0, J = c(0,0)) {
  ##
  ## Generate n points from a random start d dimensional Halton sequence.
  ##
  ## Inputs:
  ##
  ## n 			sample size
  ## bases    	coprime bases e.g. c(2,3) (Halton Sequence)
  ## seeds  	random seeds  e.g. c(0,0) (Halton Sequence)
  ## boxes 		Index of the Halton Sequence that the Box falls in
  ## B  		Number of boxes to divide Halton Sequence into


  ########### Initialize #########################################
  d <- length(bases);	pts <- mat.or.vec(n, d)
  if (length(seeds) != d){
    seeds <- rep(seeds[1],d)
  }

  boxes <- where2Start(boxes = boxes, J = J[1:2], seeds = seeds[1:2], bases = bases[1:2])
  B <- prod(bases[1:2]^J)

  ########### Main Loop #########################################
  for (i in 1:d) {
    b <- bases[i];   	u <- seeds[i];
    k <- (rep(u + boxes, ceiling(n/length(boxes))) + rep(seq(0, ceiling(n/length(boxes)) - 1)*B, each = length(boxes)))[1:n]
    xk <- (k %% b)/b;
    for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
      xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
    }
    pts[,i] <- cbind(xk)
  }
  pts <- cbind(k+1-u, pts)
  return(pts)
}

#' Create a Halton Frame raster, defined in sf based on the bounding box master sample.
#'
#' Make a Halton Frame based on B = 2^J[1]*3^J[2] grid cells. If rotation is required, will return rotated.
#' This function is an internal function simply to select sub BAS points without having to do spatial clipping at the point level.
#'
#' @param base Co-prime base for BAS, do not change from 2,3.
#' @param J Definition for the number of grid cells of Halton frame.
#' @param bb Bounding box shapefile with centroid, random seed, rotation.
#' @param rotate Boolean if you want to rotate shape before exporting.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' bb <- getBB()
#' haltonFrame <- makeFrame(J = c(8,4), bb = bb)
#' }
#' @export
makeFrame <- function(base = c(2,3), J = c(2,2), bb, rotate = FALSE)
{
  b.bounds <- st_bbox(bb)
  B <- prod(base^J)
  halt.grid <- raster(extent(matrix( b.bounds, 2, 2 )), nrow=base[2]^J[2], ncol=base[1]^J[1])
  halt.grid <- rasterToPolygons(halt.grid)
  projection(halt.grid) <- st_crs(bb)$proj4string
  if(rotate) return(rotate.shp(st_as_sf(halt.grid), bb))
  return(st_as_sf(halt.grid))
}

#' Clip a Halton Frame based on the bounding box to the current shape.
#'
#' Take a shapefile as a sf object and clips boxes from a Halton frame around it. Size of those boxes is chosen
#' by choosing J, the number of base 2,3 powers to subdivide. Intended for internal use but can be useful in other
#' context.
#'
#' @param shp shape as spatial features object to wrap into Halton frame.
#' @param bb Master Sample bounding box.
#' @param base Co-prime base for BAS, do not change from 2,3.
#' @param J Definition for the number of grid cells of Halton frame.
#' @param projstring Projection that the master sample is in, can be passed as part of the bounding box.
#' @param rotate Boolean if you want to rotate shape before exporting.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' bb <- getBB()
#' data(NS_bioregion)
#' haltonBoxes <- shape2Frame(shp = NS_biogregion, J = c(6,4), bb = bb)
#' }
#' @export
shape2Frame <- function(shp, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL, rotate = FALSE)
{
  if( !is.null( bb))
  {
	bb.bounds <- st_bbox(bb)
	scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
	shift.bas <- bb.bounds[1:2]
	theta <- attr(bb, "rotation")
	cntrd <- attr(bb, "centroid")
	projstring <- st_crs(bb)$proj4string
  }else{ return("Define Bounding Box Please.")}

  if( is.null( projstring)) {
	projstring <- getProj()
	cat("Assuming Projection\n")
	}
  if(st_crs(shp) != st_crs(projstring)) shp <- st_transform(shp, projstring)

  #Stretch bounding box to Halton Frame Size:
  shp2 <- rotate.shp(shp, bb, back = FALSE)
  bb2 <- st_bbox(shp)
  xy <- (bb2 - shift.bas[c(1,2,1,2)])/scale.bas[c(1,2,1,2)]
  lx <- floor(xy[1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- floor(xy[2] / (1/base[2]^J[2]))/(base[2]^J[2])
  ux <- ceiling(xy[3] /(1/base[1]^J[1]))/(base[1]^J[1])
  uy <- ceiling(xy[4] /(1/base[2]^J[2]))/(base[2]^J[2])
  nx <- (ux-lx)*base[1]^J[1]
  ny <- (uy-ly)*base[2]^J[2]

  bb.new <- c(lx,ly, ux, uy)*scale.bas[c(1,2,1,2)] + shift.bas[c(1,2,1,2)]
  halt.frame <- raster(extent(matrix( bb.new , 2, 2)), nrow=ny, ncol=nx)
  projection(halt.frame) <- projstring
  halt.poly <- rasterToPolygons(halt.frame)
  if(rotate) return(rotate.shp(st_as_sf(halt.poly), bb))
  return(st_as_sfc(halt.poly))
}


#' Internal function to find the ordering of the first box according to the random seed.
#'
#' This is a function to find which Halton Box the initial BAS point from the Master Sample falls into and
#' thus use it to order the remaining boxes based on the initial. It also helps us tracks
#' the master sample index as we skip boxes that have no resource.
#'
#' @param J Definition for the number of grid cells of Halton frame.
#' @param seeds Master Sample random seed.
#' @param bases Co-prime bases should really always be 2,3
#' @param boxes ordering of boxes that have been clipped to be reordered according to the master sample seed.
#'
#' @return vector of reordered Halton indices.
#'
#' @export
where2Start <- function(J = c(1,1), seeds = c(0,0), bases = c(2,3), boxes = NULL)
{
  B <- prod(bases^J)
  L <- seeds %% bases^J
  boxInit <- SolveCongruence(matrix(L, ncol = 2, nrow = 1), bases, J)

  if(is.null(boxes)) return()
  boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
  return(sort(boxes))
}