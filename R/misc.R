
remove_small <- function(intr, d) {

	intr <- terra::disagg(intr)
	intr$area <- terra::expanse(intr, unit="km")
	x <- intr[intr$area> d,]
	y <- intr[intr$area<= d,]
	y <- terra::disagg(terra::aggregate(y, dissolve=TRUE))
	terra::values(x) <- NULL
	terra::values(y) <- NULL
	if (nrow(y) > 0 & nrow(x) > 0) {
		intr <- terra::combineGeoms(x, y, boundary=TRUE, distance=TRUE)
		intr$area <- terra::expanse(intr, unit="km")
	}
	intr
}


adjust_range <- function(x, sp, land, CAmin=50000, CAmax=250000) { 
	ca_add <- terra::buffer(sp, CAmin, quadsegs=12) #|> terra::aggregate()
	ca_remove <- terra::buffer(sp, CAmax) # |> terra::aggregate()
	x <- terra::mask(x, ca_remove, updatevalue=NA)
	x <- terra::rasterize(ca_add, x, update=TRUE)
	terra::mask(x, land)
}

