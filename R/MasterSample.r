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
makeFrame <- function(base = c(2,3), J = c(2,2), bb)
{
  B <- prod(base^J)
  halt.grid <- raster(extent(as.matrix( bb )), nrow=base[2]^J[2], ncol=base[1]^J[1])
  halt.grid <- rasterToPolygons(halt.grid)
  return(halt.grid)
}

# Wrap a Halton Frame over the sample Shape.
#' @name shape2Frame
#' @title Wrap a Halton Frame over the sample shape
#' @export
shape2Frame <- function(shp, bb = NULL, base = c(2,3), J = c(2,2), projstring = NULL)
{
  if( !is.null( bb))
  {
    scale.bas <- bb[,2] - bb[,1]
    shift.bas <- bb[,1]
  }else{ return("Define Bounding Box Please.")}

  if( is.null( projstring)) {
	projstring <- getProj()
	cat("Assuming Projection\n")
	}
  if(crs(shp)@projargs != projstring) shp <- spTransform(shp, projstring)

  #Stretch bounding box to Halton Frame Size:
  bb2 <- bbox(shp)
  xy <- (bb2 - shift.bas)/scale.bas
  lx <- floor(xy[1,1] / (1/base[1]^J[1]))/(base[1]^J[1])
  ly <- floor(xy[2,1] / (1/base[2]^J[2]))/(base[2]^J[2])
  ux <- ceiling(xy[1,2] /(1/base[1]^J[1]))/(base[1]^J[1])
  uy <- ceiling(xy[2,2] /(1/base[2]^J[2]))/(base[2]^J[2])
  nx <- (ux-lx)*base[1]^J[1]
  ny <- (uy-ly)*base[2]^J[2]

  bb.new <- data.frame(min = c(lx, ly), max = c(ux, uy), row.names = c("x","y"))
  halt.frame <- raster(extent(as.matrix( bb.new*scale.bas + shift.bas )), nrow=ny, ncol=nx)
  projection(halt.frame) <- projstring
  halt.poly <- rasterToPolygons(halt.frame)
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
  bb <- data.frame(min = c(85148,33745), max = c(1280999,1351981), row.names = c("x","y")) 
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

#' @name masterSample
#' @title Generate sample points using BAS master sample
#' @description Generates BAS sample points in a specified sample frame based on New Zealand terrestrial mastersample. Users need to specify 'island' the sample frame is in and how many points are required.
#' @export
masterSample <- function(shp, N = 100, msproj = NULL, bb = NULL, seed = NULL){
  
  # We always use base 2,3
  base <- c(2,3)

  # Set up Western Canada Marine Master Sample as default, general for any.
  if(is.null(msproj)) msproj <- getProj()
  if(is.null(bb)) bb <- getBB()
  if(is.null(seed)) seed <- getSeed()
  
  if(crs(shp)@projargs != msproj) 
  {
	orig.crs <- crs(shp)@projargs
	shp <- spTransform(shp, msproj)
  }
  
  #Scale and shift Halton to fit into bounding box
  scale.bas <- bb[,2] - bb[,1]
  shift.bas <- bb[,1]

  #We can use Halton Boxes to speed up code when the polygons are small and all over the place.
  #Kind of like magic!
  draw <- N + 5000
  
  J <- c(0, 0)
  hal.frame <- shape2Frame(shp, J = J, bb = bb, projstring = msproj)
  area.shp <- sum(raster::area(shp))
  while(area.shp < 0.25*raster::area(hal.frame)[1])	# Subset again:
  {
	if(base[2]^J[2] > base[1]^J[1]){ 
		J[1] <- J[1] + 1
	}else{
		J[2] <- J[2] + 1
	}
	hal.frame <- shape2Frame(shp, J = J, bb = bb, projstring = msproj)	
  }
	
	boxes <- which(rowSums(gIntersects(shp, hal.frame, byid = TRUE)) > 0)
	hal.polys <- hal.frame[boxes,]@polygons
	
	# Find the corner Halton Pts
	box.lower <- do.call("rbind", lapply(hal.polys, FUN = function(x){data.frame(t(x@labpt))}))
	box.lower <- t(apply(box.lower, 1, FUN = function(x){(x - shift.bas)/scale.bas}))
	A <- GetBoxIndices(box.lower, base, J)
	halt.rep <- SolveCongruence(A, base, J)	
	B <- prod(c(2,3)^J)

	# I like to know how many divisions we had to make...
	print(J)

  getSample <- function(k = 0, endPoint = 0){
    if(k == 0){ seedshift <- seed
    }else seedshift <- endPoint + seed
    pts <- RSHalton(n = draw, seeds = seedshift, bases = c(2,3), boxes = halt.rep, J = J)
    pts[,2] <- pts[,2]*scale.bas[1] + shift.bas[1]
    pts[,3] <- pts[,3]*scale.bas[2] + shift.bas[2]
    pts.coord <- SpatialPointsDataFrame(cbind(pts[,2],pts[,3]),proj4string=CRS(msproj), data.frame(SiteID = paste0(pts[,1] + endPoint)))
    indx <- gIntersects(shp, pts.coord, byid = TRUE)
    pts.coord <- pts.coord[rowSums(indx) > 0,]
    return(pts.coord)
  }

  pts.sample <- getSample()
  while(nrow(pts.sample) == 0) {
    draw <- draw * 2
    pts.sample <- getSample()
  }

  di <- 1
  while(nrow(pts.sample) < N){
    last.pt <- pts.sample@data$Count[nrow(pts.sample)]
    new.pts <- getSample(k = di, endPoint = last.pt)
    if(nrow(new.pts) > 0) pts.sample <- rbind(pts.sample, new.pts)
    di <- di + 1
  }
	
  smp <- pts.sample[1:N,]
  if(orig.crs != msproj) 
  {
	smp <- spTransform(smp, orig.crs)
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

#' @name lineSamp
#' @title Generate samples in linear features based on BAS mastersample
#' @export
# Sampling a linear feature:
lineSamp <- function (n = 10, x, seed = 0, halt = TRUE) 
{	
	cc <- do.call("c", coordinates(x))
	# Trick here is to figure out where the "breaks" in the lines occur
	# By index cumulative sum we can identify them
	brks <- sapply(cc, nrow)	
	brk <- cumsum(brks[-length(brks)])	
    cc.df <- do.call("rbind", cc)
    cc.mat <- as.matrix(cc.df)
    lengths = LineLength(cc.mat, longlat = FALSE, sum = FALSE)
	lengths[brk] <- 0	# Remove the length for discontinuities.
    csl = c(0, cumsum(lengths))
    maxl = csl[length(csl)]
    if (halt == TRUE) {
        pts = lineHalton(n, u = seed) * maxl
    }
    else {
        pts = runif(n) * maxl
    }
    int = findInterval(pts, csl, all.inside = TRUE)
    where = (pts - csl[int])/diff(csl)[int]
    xy = cc.mat[int, , drop = FALSE] + where * (cc.mat[int + 
        1, , drop = FALSE] - cc.mat[int, , drop = FALSE])
    samp <- SpatialPoints(xy, proj4string = CRS(proj4string(x)))
	return(samp)
}

#' @export
lineHalton <- function(n = 10, u = 0, b = 5) {
  k <- u:(u+n-1);    xk <- (k %% b)/b;
  for (j in 1:(ceiling(logb(u+n,b)) + 2)) {
    xk <- xk + (floor(k/(b^j)) %% b)/(b^(j+1));
  }
  return(xk)
}
