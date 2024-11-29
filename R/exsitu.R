n_zones <- function(x, min_area, m=3) {
	round(pmax(1, pmin(x, m*sqrt(x))))
}


make_zones <- function(x, range, n, spread=TRUE) {

	if (spread) {
		rr <- terra::rasterize(range, rast(x), 1:nrow(range))
		p <- terra::as.polygons(rr)
		p$area <- expanse(p, "km") / 1000 
		avga <- sum(p$area) / n
		p$n <- round(p$area / avga)
		totn <- sum(p$n)
		while (totn > n) {
			p$error[i] <- p$area[i] - p$n[i] * avga
			i <- which.max(p$error)
			p$n[i] <- p$n[i] - 1
			totn <- sum(p$n)
		}
		while (totn < n) {
			p$error[i] <- p$area[i] - p$n[i] * avga
			i <- which.min(p$error)
			p$n[i] <- p$n[i] + 1
			totn <- sum(p$n)
		}
		p <- p[p$n > 0, ]
		seeds <- lapply(1:nrow(p), \(i) spatSample(p[i,], p$n[i])) |> vect() |> crds()
		km <- terra::k_means(mask(x, rr), seeds, iter.max = 25)
	} else {

		xm <- terra::mask(x, range)
		km <- terra::k_means(xm, n, iter.max = 25)
	}

	terra::as.polygons(km)
}


get_samplesize <- function(range, fun=n_zones, ...) {
	if (inherits(range, "SpatRaster")) {
		range <- terra::as.polygons(range)
		range <- range[range[,1,drop=TRUE]==1]
	} else if (!inherits(range, "SpatVector")) {
		stop("range should be a SpatVector")
	}
	a <- terra::expanse(range, unit="km")
	n <- fun(a, ...)
	#n <- max(nmin, n)
	#z <- max(1, min(n, round(a/min_area)))
	return(list(range=range, area=a, n=n))
}


branch_length <- function(tree, samp, adjust, adjfun=log10 ) {

## TODO take into account that you may need multiple obs per leaf

	if (adjust) {
		tab <- table(samp)
		tab <- adjfun(tab) + 1
		# or? 
		# tab <- pmax(1, adjfun(tab))

	}
	sample <- as.character(unique(samp))
    if (is.null(tree$edge.length)) {stop("Tree has no branch lengths, cannot compute pd") }
    
	absent <- tree$tip.label[!(tree$tip.label %in% sample)]
    if (length(sample) == 0) {
        GD <- 0
    } else if (length(sample) == 1) {
		# also adjust for sample size
        GD <- max(tree$edge.length)
    } else if (length(absent) == 0) {
		if (adjust) {
			i <- match(names(tab), tree$tip.label)
			tree$edge.length[i] <- tree$edge.length[i] * tab
		}
        GD <- sum(tree$edge.length)
    } else {
        tree <- ape::drop.tip(tree, absent)
		if (adjust) {
			i <- match(names(tab), tree$tip.label)
			tree$edge.length[i] <- tree$edge.length[i] * tab
        }
		GD <- sum(tree$edge.length)
    }
    return(GD)
}


small_ssize_penalty <- function(ssize, score, minssize=10) {
# linear, make curvilinear
	ifelse(ssize > minssize, score, score * ssize / minssize)
}


ex_cvs_old <- function(regions, sample, env=NULL, adjust=FALSE, minssize=10) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	if (minssize <= 0) {
		return(0)
	}
		
	if (nrow(regions) == 1) {
		# cannot make a tree 
		s <- terra::extract(regions, sample)
		if (length(na.omit(terra::extract(regions, s)[,2])) > 1) { 
			return(small_ssize_penalty(length(sample), 1, minssize))
		} else {
			return(0)		
		}
	}
		
	xy <- terra::centroids(regions)
	d <- terra::distance(xy, unit="km")
	if (!is.null(env)) {
		e <- terra::extract(env, regions, fun=mean, na.rm=TRUE)
		ed1 <- dist(e[,2])
		ed2 <- dist(e[,3])		
	} 
	
	a <- cbind(g=as.vector(d), tmp=as.vector(ed1), prc=as.vector(ed2))
	
	
	x <- stats::hclust(d)
	x <- ape::as.phylo(x)
	terra::values(regions) <- NULL
	
	# if nsamples for a region > threshold (10) one neighbor that is empty can get an observation
	##
	
	s <- terra::extract(regions, sample)

	actual_pd <- branch_length(x, s[,2], adjust=adjust)	
	potential_pd <- branch_length(x, 1:nrow(regions), adjust=adjust)
	score <- min(1, actual_pd/potential_pd)
	
	small_ssize_penalty(length(sample), score, minssize)
	
}


EnvDist <- function(env, regions, envfun) {
	e <- terra::extract(env, regions, fun=mean, na.rm=TRUE, ID=FALSE)
	ed1 <- dist(e[,1])
	ed2 <- dist(e[,2])	
	ed <- data.frame(ed1, ed2)
	names(ed) <- names(env)
	# express environmental distance expressed as geographic distance
	envd <- envfun(ed)
	structure(envd, class = 'dist', Size=attr(ed1, "Size"))
}



XS_dendro <- function(regions, sample, env=NULL, envfun=NULL, clustmethod="complete", adjust=FALSE, minssize=10) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	stopifnot(minssize > 0)

	if (nrow(sample) <= 0) {
		list(score=0)
	}
		
	if (nrow(regions) == 1) {
		# cannot make a tree 
		return(small_ssize_penalty(nrow(sample), 1, minssize))
	}
		
	xy <- terra::centroids(regions)
	d <- terra::distance(xy, unit="km")
	if (!is.null(env)) {
		envd <- EnvDist(env, regions, envfun)
		d <- (d + envd)
	} 

	x <- stats::hclust(d, method=clustmethod)
	x <- ape::as.phylo(x)
	terra::values(regions) <- NULL
	
	# if nsamples for a region > threshold (10) one neighbor that is empty can get an observation
	##
	
	s <- terra::extract(regions, sample)

	actual_pd <- branch_length(x, s[,2], adjust=adjust)	
	potential_pd <- branch_length(x, 1:nrow(regions), adjust=adjust)
	score <- min(1, actual_pd/potential_pd)
	
	score <- small_ssize_penalty(length(sample), score, minssize)
	list(score=score)
}


get_network1 <- function(regions, sample, land=NULL, maxlink=750, weed=TRUE) {

	terra::values(regions) <- data.frame(id=1:nrow(regions))
	patches <- terra::disagg(terra::aggregate(regions))
	patches$pid <- 1:nrow(patches) # pid = patch id 

	x <- terra::centroids(regions, inside=TRUE)

	x$pid <- terra::extract(patches, x)$pid
	patches <- patches[unique(x$pid), ]

	x <- terra::round(x, 5)  # to help merge

## make graph of delauny edges 
	d <- terra::disagg(terra::as.lines(terra::delaunay(x)), TRUE)


## should not cross the sea
#	if (!is.null(land)) {
#		d$id <- 1:nrow(d)
#		eland <- terra::as.polygons(terra::ext(land))
#		eland <- terra::erase(eland, land)
#		rem <- terra::intersect(d, eland)$id
#		d <- d[!(d$id %in% rem), ]
#	}
	

	g <- terra::geom(d)[, c(1, 3:4)]
	p <- terra::geom(x)[, c(1, 3:4)]
	colnames(p)[1] <- "att"
	m <- merge(g, p)
	m <- m[order(m$geom), ]

	rr <- data.frame(
		from = m$att[seq(1, nrow(m), 2)],
		to = m$att[seq(2, nrow(m), 2)]
	) 
	rr[,1:2] <- apply(rr[,1:2], 1, sort) |> t()
	rr <- rr[order(rr[,1], rr[,2]), ] |> unique()

	rr <- cbind(rr, w=1, dst=terra::distance(x[rr[,1], ], x[rr[,2], ], pairwise=TRUE, unit="m")/1000)

	if (!weed) return(r)




	#vx <- values(x)

## only the shortest between not connected ranges, unless distance is very short
	rr$fpid <- x$pid[match(rr$from, x$id)]
	rr$tpid <- x$pid[match(rr$to, x$id)]

	# in different patch, and far apart
#	i <- (rr$fpid != rr$tpid) & (rr$dst > maxlink)  
#	prr <- rr[i,]

	# same patch
#	rr <- rr[!i, ]
#	rr <- rr[rr$dst < maxlink, , drop=FALSE]
	

	i <- rr$dst > maxlink
	prr <- rr[i,]
	rr <- rr[!i, ]

	# keep shortest distance between two patches
	prr <- suppressWarnings(prr[ prr$dst == ave(prr$dst, prr$fpid, prr$tpid, FUN=min), ])
	tst <- rr[, c("fpid", "tpid")]
	tst$test <- 1
	m <- merge(prr, tst, by=c("fpid", "tpid"), all.x=TRUE)
	m <- m[is.na(m$test), ]	
	m$test <- NULL
	if (nrow(m) > 0) {
		nx <- nrow(x)
		m <- m[order(m$dst), ]
		# if you cannot get there otherwise
		g <- igraph::graph_from_data_frame(rr, directed = FALSE) |> igraph::decompose()
		gg <- lapply(g, \(x)igraph::V(x)$name)
		for (i in 1:nrow(m)) {
			if ((length(gg)==1) & (length(gg[[1]]) == nx)) break
			a <- sapply(gg, \(x) m$from[i] %in% x) 
			b <- sapply(gg, \(x) m$to[i] %in% x) 
			if (any((a+b)==1)) {
				rr <- rbind(rr, m[i, ])
				g <- igraph::graph_from_data_frame(rr, directed = FALSE) |> igraph::decompose()
				gg <- lapply(g, \(x)igraph::V(x)$name)
			}
		}
	}
	
	xy <- terra::crds(x)
	colnames(xy) <- c("xf", "yf")
    rr <- cbind(rr, xy[rr$from, ])
	colnames(xy) <- c("xt", "yt")
    rr <- cbind(rr, xy[rr$to, ])
	
#	v = exsitu:::make_spatvect(rr)
#	v$DID <- 1:nrow(d)
#	v$flen <- terra::perim(v) / 1000
#	eland <- terra::as.polygons(terra::ext(regions))
#	eland <- rbind(terra::erase(eland, regions), regions)
#	rem <- terra::intersect(v, eland)
#	rem$clen <- terra::perim(rem) / 1000
#	rem$frac <- rem$clen / rem$flen
#	rem <- rem[rem$id == rem$from | rem$id == rem$to, ] 
#	vrem <- values(rem)
#	a <- aggregate(vrem["frac"], vrem[, c("from", "to")], sum)
#	a <- a[a$frac < 0.66, ]
	

#	plot(regions, border="gray")
#	lines(v, lwd=.5, col="blue")

	rr

}


make_spatvect <- function(rr) {
	x <- as.matrix(rr[, c("xf", "yf", "xt", "yt")])
	a <- lapply(1:nrow(x), \(i)  cbind(i, matrix(x[i,], nrow=2, byrow=2)))
	b <- do.call(rbind, a)
	v <- terra::vect(b, "lines")
	values(v) <- rr
	crs(v) <- "+proj=longlat"
	v
}



get_network2 <- function(regions, sample, maxdist=1500) {

	terra::values(regions) <- data.frame(id=1:nrow(regions))
	patches <- terra::disagg(terra::aggregate(regions))
	patches$pid <- 1:nrow(patches) # pid = patch id 

	x <- terra::centroids(regions, inside=TRUE)

	x$pid <- terra::extract(patches, x)$pid
	patches <- patches[unique(x$pid), ]
	np <- nrow(patches)

	x <- terra::round(x, 5)  # to help merge
	xy <- terra::crds(x)

	adj <- adjacent(regions)
	adj <- data.frame(unique(t(apply(adj, 1, sort))))
	colnames(adj) <- c("from", "to")
	adj$adj <- 1
	if (np > 1) {
		up <- sort(unique(x$pid))
		dx <- as.matrix(distance(x))
		diag(dx) <- NA
		i <- match(colnames(dx), x$id)
		pid <- x$pid[i]
		for (p in up) {
			s <- dx[pid != p, pid == p, drop=FALSE]
			rid <- pid[pid != p]
			upp <- sort(unique(rid))
			for (pp in upp) {
				ss <- s[rid == pp, ,drop=FALSE]
				j <- which.min(apply(ss, 1, min))
				k <- which.min(ss[j, ])
				adj <- rbind(adj, c(sort(as.integer(c(rownames(ss)[j], colnames(ss)[k]))), 0))
			}
		}
		adj <- unique(adj)
	}

	colnames(xy) <- c("xf", "yf")
	adj <- cbind(adj, xy[adj$from, ])
	colnames(xy) <- c("xt", "yt")
	adj <- cbind(adj, xy[adj$to, ])

	adj <- cbind(adj, w=1, dst=terra::distance(x[adj$from, ], x[adj$to, ], pairwise=TRUE, unit="m")/1000)

	if (np > 2) {
		padj <- adj[adj$adj == 0, ]
		adj <- adj[adj$adj != 0, ]
		nx <- nrow(x)
		padj <- padj[order(padj$dst), ]
		g <- igraph::graph_from_data_frame(adj, directed = FALSE) |> igraph::components()
		for (i in 1:nrow(padj)) {
			adj2 <- rbind(adj, padj[i,])
			gg <- igraph::graph_from_data_frame(adj2, directed = FALSE) |> igraph::components()
			if ((gg$no < g$no) | (length(gg$membership) > length(g$membership))) {
				g <- gg
				adj <- adj2
			}
		}
	}

	adj[adj$dst < maxdist, ]

}



XS_net <- function(regions, sample, env=NULL, envfun=NULL, minssize=10, maxlink=1500, return_network=FALSE) {

## TODO  RH
# fix the adjust effect such that when you have many observations in one zones
# they can only contribute to their neighbors. Do not increase branch length to avoid that
# one region does not compensate for another

	stopifnot(minssize > 0)

	if (nrow(sample) <= 0) {
		list(score=0)
	}
#
# weed=TRUE
# land=NULL
#	rr <- get_network1(regions, sample, land=land, maxlink=maxlink, weed=weed)
	rr <- get_network2(regions, sample, maxlink)
	if (return_network) {
		return(make_spatvect(rr))
	}


	if (!is.null(env)) {
		envd <- as.matrix(EnvDist(env, regions, envfun))
		rr$envdst <- envd[as.matrix(rr[,1:2])]
		rr$geodst <- rr$dst
		# sum geo and env dist
		rr$dst <- rr$dst + rr$envdst
	} 
	
	gg <- igraph::graph_from_data_frame(rr, directed = FALSE)
	igraph::E(gg)$weight <- rr$dst * rr$w
	y <- unique(terra::extract(regions, sample)[,2])
	rr$w2 <- rowSums(!matrix(as.matrix(rr[,1:2]) %in% y, ncol=2)) / 2
	igraph::E(gg)$weight2 <- rr$dst * rr$w2

	n <- igraph::count_components(gg)
	score <- nodes <- rep(NA, n)
	dg <- igraph::decompose(gg)
	for (k in 1:n) {
		g <- dg[[k]]
		dst <- igraph::distances(g)
		d1 <- sum(dst) 

		igraph::E(g)$weight <- igraph::E(g)$weight2

		if (length(y) > 1) {
			b <- combn(as.character(y), 2)
			nms <- igraph::V(g)$name
			haveb <- apply(matrix(b %in% nms, nrow=2), 2, all)
			b <- b[,haveb,drop=FALSE]
			if (ncol(b) > 0) {
				for (i in 1:ncol(b)) {
					if (!igraph::are_adjacent(g, b[1,i], b[2,i])) {
						g <- igraph::add_edges(g, b[,i], weight=0)
					}
				}
			}
		}
		d2 <- igraph::distances(g)
		score[k] <- 1 - (sum(d2) / d1)
		nodes[k] <- length(g)
	}
	score <- weighted.mean(score, nodes)
	
	score <- small_ssize_penalty(length(sample), score, minssize)
	list(score=score, dist=rr)
}

