########################
#Master Sample Code
#Paul van Dam-Bates
########################
#Take in a polygon and generate points from the master sample
#Step 1: Determine Island that the polygon falls in
#Step 2: Use random seed from that island to start Halton Sequence
#Step 3: Ouptut number of points required clipped for that region.
#--------------------------------------------------------------------

#' @import raster
#' @import sp
#' @import sf
#' @import rgeos
#' @import Rcpp
#' @import rgdal
#' @import raster
NULL

## usethis namespace: start
#' @useDynLib BASMasterSample, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL


#' @export
#Halton Sequence:
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

  boxes <- where2Start(boxes = boxes, J = J, seeds = seeds, bases = bases)
  B <- prod(bases^J)

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


#' @export
rot <- function(a) 
{
		matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
}

#' @export
rotate.bb <- function(shp, theta)
{	
	if (class(shp)[1] != "sf") 
	{
		shp <-  st_as_sf(shp)
	}

	bb <- st_as_sfc(st_bbox(shp))
	
	cntrd <- st_centroid(bb)
	bb.rot <- (bb - cntrd) * rot(theta) + cntrd
	bb.new <- st_as_sfc(st_bbox(bb.rot))

	attr(bb.new, "rotation") = theta
	attr(bb.new, "centroid") = st_coordinates(cntrd)
	return(bb.new)
}

# Give random points on the unit square and map to rotated bounding box.
#' @export
rotate.scale.coords <- function(coords, bb, back = TRUE)
{	
	coords <- pts[,2:3]
	theta <- ifelse(back, -1, 1) * attr(bb, "rotation")	# Rotate backwards
	cntrd <- attr(bb, "centroid")

	bb.bounds <- st_bbox(bb)
	bb.scale <- diag(2) * (bb.bounds[3:4] - bb.bounds[1:2])
	
	coords <- st_multipoint(coords, dim = "XY") %>% st_geometry()
	coords.scale <- coords*bb.scale + bb.bounds[1:2]
	coords.rot <- (coords.scale - cntrd) * rot(theta) + cntrd

	st_crs(coords.rot) <- st_crs(bb)
	return(coords.rot)
}

# Give random points on the unit square and map to rotated bounding box.
#' @export
rotate.poly <- function(shp, bb, back = TRUE)
{	
	theta <- ifelse(back, -1, 1) * attr(bb, "rotation")	# Rotate backwards
	cntrd <- attr(bb, "centroid")

	shp <- st_transform(shp, st_crs(bb))
	shp <- st_geometry(shp)
	shp.rot <- (shp - cntrd) * rot(theta) + cntrd
	st_crs(shp.rot) <- st_crs(bb)
	return(shp.rot)
}


#' @export
#Solve Congruence to get order of boxes:
systemCong <- function(L = c(1/4, 1/3), J = c(2,2), base = c(2,3))
{
  x <- 0:(base[1]^J[1]-1)/1/base[1]^J[1]	#Fix rounding errors!
  y <- 0:(base[2]^J[2]-1)/1/base[2]^J[2]

  L <- c(x[which.min(abs(x - L[1]))], y [which.min(abs(y - L[2]))])

  a1 <- sum( ( floor((L[1] + .Machine$double.eps)*base[1]^(1:J[1])) %% base[1]) * base[1]^( 1:J[1]-1))
  a2 <- sum( ( floor((L[2] + .Machine$double.eps) *base[2]^(1:J[2])) %% base[2]) * base[2]^( 1:J[2]-1))
  mod <- base^J
  B <- prod(mod)
  possible <- 0:(B-1)
  sol1 <- possible %% mod[1] == a1
  sol2 <- possible %% mod[2] == a2
  return(possible[which(sol1 + sol2 == 2)])
}

#' @name makeFrame
#' @title Make a Halton grid over the bounding box
#' @export
#Create a Halton Grid over the Bounding Box
makeFrame <- function(base = c(2,3), J = c(2,2), bb, rotate = NULL)
{
  b.bounds <- st_bbox(bb)
  B <- prod(base^J)
  halt.grid <- raster(extent(matrix( b.bounds, 2, 2 )), nrow=base[2]^J[2], ncol=base[1]^J[1])
  halt.grid <- rasterToPolygons(halt.grid)
  projection(halt.grid) <- st_crs(bb)$proj4string
  if(!is.null(rotate)) return(rotate.poly(st_as_sf(halt.grid), bb))
  return(st_as_sf(halt.grid))
}

# Wrap a Halton Frame over the sample Shape.
#' @name shape2Frame
#' @title Wrap a Halton Frame over the sample shape
#' @export
shape2Frame <- function(shp, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL, rotate = NULL)
{
  if( !is.null( bb))
  {
	bb.bounds <- st_bbox(bb)
	scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
	shift.bas <- bb.bounds[1:2]
	theta <- attr(bb, "rotation")
	cntrd <- attr(bb, "centroid")
  }else{ return("Define Bounding Box Please.")}

  if( is.null( projstring)) {
	projstring <- getProj()
	cat("Assuming Projection\n")
	}
  if(st_crs(shp) != st_crs(projstring)) shp <- st_transform(shp, projstring)

  #Stretch bounding box to Halton Frame Size:
  shp2 <- rotate.poly(shp, bb, back = FALSE)
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
  if(!is.null(rotate)) return(rotate.poly(st_as_sf(halt.poly), bb))
  return(st_as_sfc(halt.poly))
}


#' @export
#Where to start the Halton Sequence
where2Start <- function(J = c(1,1), seeds = c(0,0), bases = c(2,3), boxes = NULL)
{
  B <- prod(bases^J)
  L <- seeds %% bases^J
  boxInit <- SolveCongruence(matrix(L, ncol = 2, nrow = 1), bases, J)

  if(is.null(boxes)) return()
  boxes <- ifelse(boxes < boxInit, B + (boxes - boxInit), boxes - boxInit)
  return(sort(boxes))
}

buildMS <- function(shp, d = 2, showOutput = TRUE)
{
  
  seed <- floor(runif(d, 0, 10000))
  # We always use base 2,3
  base <- c(2,3,5)[1:d]

  # Just work with sf:
  if (class(shp)[1] != "sf") 
  {
    shp <-  st_as_sf(shp)
  }
  
  # Create a Random Rotation:
  theta <- runif(1, -pi, pi)
  build.bb <- rotate.bb(shp, theta = theta)
  st_crs(build.bb) <- st_crs(shp)
  attr(build.bb, "seed") <- seed
  
  if(showOutput){
	cat("Seed:", seed, "\n")
	cat("Rotation:", theta, "Radians\n")
  }
  return(build.bb)
}

#' @name getProj
#' @title Define spatial objects in projection of the master sample.
#' @description Default projection of the master sample. Needed for consistency for the entire bounding box.
#' @export
getProj <- function()
{	#BC Albers
    #http://spatialreference.org/ref/epsg/3005/
    msproj <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
	return(msproj)
}


#' @name getBB
#' @title Get the bounding box for other functions.
#' @description Bounding box under NAD83/BC Albers.
#' @export
getBB <- function()
{
  bb.df <- c("xmin" = 85148, "ymin" = 33745, "xmax" = 1280999, "ymax" = 1351981)
  bb <- st_as_sfc(st_bbox(bb.df))

  attr(bb, "rotation") = 0
  attr(bb, "centroid") = st_centroid(bb)  
  
  st_crs(bb) <- st_crs(getProj())
  return(bb)
}

#' @name getSeed
#' @title Random seed definition for Western Canada Marine Master Sample.
#' @description Defines the random seed specific to the Western Canada Marine Master Sample. 
#' @export
getSeed <- function()
{
  seed <- c(37916, 85846)
  return(seed)
}

getRotation <- function()
{
  return(0)
}


#' @name masterSample
#' @title Generate sample points using BAS master sample
#' @description Generates BAS sample points in a specified sample frame based on New Zealand terrestrial mastersample. Users need to specify 'island' the sample frame is in and how many points are required.
#' @export
masterSample <- function(shp, N = 100, bb = NULL, nExtra = 5000){
  
  # We always use base 2,3
  base <- c(2,3)

  # Updating to work for sf only. Start here...
  if (class(shp)[1] != "sf") 
  {
    shp <-  st_as_sf(shp)
  }

  # Set up Western Canada Marine Master Sample as default, general for any.
  # bb now includes its rotation as well.
  if(is.null(bb))
  {  
	bb <- getBB()
	msproj <- getProj()
	seed <- getSeed()
  }else{
	msproj <- st_crs(bb)$proj4string
	seed <- attr(bb, "seed")[1:2]	# 3rd dimension is not yet supported...
  }
	
  orig.crs <- NULL
  if(st_crs(shp) != st_crs(msproj)) 
  {
	orig.crs <- st_crs(shp)$proj4string
	shp <- st_transform(shp, msproj)
  }  
  
  #Scale and shift Halton to fit into bounding box
  bb.bounds <- st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  cntrd <- attr(bb, "centroid")
  theta <- attr(bb, "rotation")

  #We can use Halton Boxes to speed up code when the polygons are small and all over the place.
  #Kind of like magic!
  draw <- N + nExtra
  
  J <- c(0, 0)
  hal.frame <- shape2Frame(shp, J = J, bb = bb, projstring = msproj)
  area.shp <- as.numeric(sum(st_area(shp)))
  while(area.shp < 0.25*as.numeric(st_area(hal.frame))[1])	# Subset again:
  {
	if(base[2]^J[2] > base[1]^J[1]){ 
		J[1] <- J[1] + 1
	}else{
		J[2] <- J[2] + 1
	}
	hal.frame <- shape2Frame(shp, J = J, bb = bb, projstring = msproj)
  }
	
	hal.fr2 <- rotate.poly(hal.frame, bb)
	hal.indx <- which(rowSums(st_intersects(hal.fr2, shp, sparse = FALSE)) > 0)
	hal.pts <- st_centroid(hal.frame)[hal.indx,] %>% st_coordinates
	
	# Find the corner Halton Pts
	box.lower <- t(apply(hal.pts, 1, FUN = function(x){(x - shift.bas)/scale.bas}))
	A <- GetBoxIndices(box.lower, base, J)
	halt.rep <- SolveCongruence(A, base, J)	
	B <- prod(c(2,3)^J)

	# I like to know how many divisions we had to make...
	print(J)

  getSample <- function(k = 0, endPoint = 0){
    if(k == 0){ seedshift <- seed
    }else seedshift <- endPoint + seed
    pts <- RSHalton(n = draw, seeds = seedshift, bases = c(2,3), boxes = halt.rep, J = J)
	xy <- cbind(pts[,2]*scale.bas[1] + shift.bas[1], pts[,3]*scale.bas[2] + shift.bas[2])
	if(theta != 0) xy <- sweep ( sweep(xy, 2,  cntrd, FUN = "-") %*% rot(-theta), 2,  cntrd, FUN = "+")
	pts.coord <- st_as_sf(data.frame(SiteID = pts[,1] + endPoint, xy), coords = c(2,3))
	st_crs(pts.coord) <- st_crs(bb)
    pts.coord <- pts.coord[shp,]
    return(pts.coord)
  }
  
  pts.sample <- getSample()
  while(nrow(pts.sample) == 0) {
    draw <- draw * 2
    pts.sample <- getSample()
  }

  di <- 1
  while(nrow(pts.sample) < N){
    last.pt <- pts.sample$SiteID[nrow(pts.sample)]
    new.pts <- getSample(k = di, endPoint = last.pt)
    if(nrow(new.pts) > 0) pts.sample <- rbind(pts.sample, new.pts)
    di <- di + 1
  }

  smp <- pts.sample[1:N,]
  if(!is.null(orig.crs)) 
  {
	smp <- st_transform(smp, orig.crs)
  }  
  return(smp)
}

#############
# Take BAS Point and Make Halton Frame around it.
#############
#' @name point2Frame
#' @title Make a halton frame around a BAS point
#' @export
point2Frame <- function(pt, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL)
{
  if(!is.null(bb))
  {
    scale.bas <- bb[,2] - bb[,1]
    shift.bas <- bb[,1]
  }else{ return("Define Bounding Box Please.")}

  if(is.null(projstring)) projstring <- proj4string(pt)

  xy <- coordinates(pt)
  xy <- (xy - shift.bas)/scale.bas
  lx <- floor(xy[1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- floor(xy[2] / (1/base[2]^J[2]))/(base[2]^J[2])
  frame.order <- systemCong(c(lx, ly), base = base, J = J)
  lx <- lx*scale.bas[1] + shift.bas[1]
  ly <- ly*scale.bas[2] + shift.bas[2]
  framei <- Polygons(list(Polygon(cbind(c(lx,lx,lx + scale.bas[1]/base[1]^J[1], lx + scale.bas[1]/base[1]^J[1],lx), c(ly,ly + scale.bas[2]/base[2]^J[2],ly + scale.bas[2]/base[2]^J[2], ly, ly)))),
                     ID = frame.order)
  framei <- SpatialPolygonsDataFrame(SpatialPolygons(list(framei), proj4string = CRS(projstring)), data = data.frame(Order = frame.order, row.names = frame.order))
  return(framei)
}

#' This is an internal Master Sample function to assign
#' indices to points based on the discrete Halton Box overlay.
#'
#' @param input An sp or sf spatial points. Accepts either. Or even a data frame with X, Y names.
#' @param J Integer for number of Halton Boxes to make. If not set it defaults to 100m roughly.
#' @param bb Master Sample bounding box.
#' @param base Generally 2,3. If you change it read the literature.
#' @param seed Master Sample Seed for two bases.
#' @param s1 Permutation of x orderings (0, 1)
#' @param s2 Permutation of y orderings (0, 1, 2)
#'
#' @export
getIndividualBoxIndices <- function(pts, J = NULL, bb)
{
	# We always use base 2,3
	base <- c(2,3)

	# Updating to work for sf only. Start here...
	if (class(pts)[1] != "sf") 
	{
		pts <-  st_as_sf(pts)
	}

	# Set up Western Canada Marine Master Sample as default, general for any.
	# bb now includes its rotation as well.
	if(is.null(bb))
	{  
		bb <- getBB()
		msproj <- getProj()
		seed <- getSeed()
	}else{
		msproj <- st_crs(bb)$proj4string
		seed <- attr(bb, "seed")[1:2]	# 3rd dimension is not yet supported...
	}

	orig.crs <- NULL
	if(st_crs(pts) != st_crs(msproj)) 
	{
		orig.crs <- st_crs(pts)$proj4string
		pts <- st_transform(pts, msproj)
	}  

	#Scale and shift Halton to fit into bounding box
	bb.bounds <- st_bbox(bb)
	scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
	shift.bas <- bb.bounds[1:2]

	cntrd <- attr(bb, "centroid")
	theta <- attr(bb, "rotation")	
	
	if(is.null(J)) J <- round(log(scale.bas/100)/log(base), 0)

	B <- prod(base^J)
	
	Bx <- base[1]^J[1]
	By <- base[2]^J[2]

	xy <- st_coordinates(pts)
	# Rotate pts to the bounding box:
	if(theta != 0) xy <- sweep ( sweep(xy, 2,  cntrd, FUN = "-") %*% rot(theta), 2,  cntrd, FUN = "+")
	# Scale to 0-1
	xy <- cbind((xy[,1] - shift.bas[1])/scale.bas[1], (xy[,2] - shift.bas[2])/scale.bas[2])
	Axy <- cbind(floor((xy[,1] + 2*.Machine$double.eps)*Bx), floor((xy[,2] + 2*.Machine$double.eps)*By))

	haltonIndex <- SolveCongruence(Axy, base = base, J = J)

	# Adjust everything for the Master Sample random seed.
	a1 <- seed[1:2] %% base^J
	boxInit <- SolveCongruence(matrix(a1, ncol = 2), base = base, J = J)
	
	# Adjusted index:
	haltonIndex <- ifelse(haltonIndex < boxInit, B + (haltonIndex - boxInit), haltonIndex - boxInit)
	# Return the Halton Index for all "pts" in dat that are passed to this function.
	
	pts$HaltonIndex <- haltonIndex
	return(pts)
}