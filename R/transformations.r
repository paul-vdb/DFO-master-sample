#' Rotation Matrix
#'
#' Generate a rotation matrix for rotating objects later 
#'
#' @param a radians of rotation.
#'
#' @return Matrix
#'
#' @examples
#' \dontrun{
#' rot(pi)
#' }
#' @export
rot <- function(a) 
{
		matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)
}

#' Rotate Bounding box by theta radians
#'
#' Given some shp defined as the boundary of interest, rotate it around the centroid
#' and return the rotation and the centroid as attributes. This is used for
#' defining a Master Sample bounding box that has random rotation while ensuring that
#' the new rotated bounding box fits the shp.
#'
#' @param shp A spatial file with the spatial boundary of the sample.
#' @param theta Radians of rotation. Positive to the right of pi/2, negative to the left.
#'
#' @return Matrix
#'
#' @examples
#' \dontrun{
#' data(NS_bioregion)
#' bb.new <- rotate.bb(NS_bioregion, -pi/3)
#' }
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

#' Scale and rotate points from the unit square to a defined projection.
#'
#' Given some coordinates on [0,1)x[0,1), shift and scale them to the bounding box, and then rotate
#' them given the bounding box rotation defined by the Master Sample.
#'
#' @param coords Output from RSHalton() to be converted to the spatial surface of interest.
#' @param bb Special shape file defining the bounding box with attributes for centroid and rotation.
#' @param back Boolean for whether or not the rotation is back to the original rotated bounding box.
#'
#' @return sf spatial points with projection defined in bb.
#'
#' @examples
#' \dontrun{
#' pts <- RSHalton(n = 10)
#' bb <- getBB()
#' pts.shp <- rotate.scale.coords(coords = pts, bb)
#' }
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

#' Rotate a polygon around the centroid of a Master Sample bounding box.
#'
#' Given some polygon within the bounding box of a Master Sample rotate it by theta defined
#' by that bounding box either backwards or forwards.
#'
#' @param shp Any polygon within the sample frame defined as a spatial features object.
#' @param bb Special shape file defining the bounding box with attributes for centroid and rotation.
#' @param back Boolean for whether or not the rotation is back to the original rotated bounding box.
#'
#' @return rotated sf spatial object.
#'
#' @examples
#' \dontrun{
#' data(NS_bioregion)
#' bb <- getBB()
#' pts.shp <- rotate.shp(shp = NS_bioregion, bb = bb)
#' }
#' @export
rotate.shp <- function(shp, bb, back = TRUE)
{	
	theta <- ifelse(back, -1, 1) * attr(bb, "rotation")	# Rotate backwards
	cntrd <- attr(bb, "centroid")

	shp <- st_transform(shp, st_crs(bb))
	shp <- st_geometry(shp)
	shp.rot <- (shp - cntrd) * rot(theta) + cntrd
	st_crs(shp.rot) <- st_crs(bb)
	return(shp.rot)
}
