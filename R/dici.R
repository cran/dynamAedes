#' Estimate of of mosquito dispersal
#'
#' Estimates of dispersal (in km^2) for the simulated mosquito population when \code{scale = "lc"}.
#' @param input_sim matrix. \code{dynamAedes} compressed output matrix (\code{compressed=TRUE}).
#' @param coords matrix. A matrix reporting the spatial coordinates of the temperature observations.
#' @param eval_date numeric. Define the day of evaluation; it refers to the column number of the input temperature matrix. 
#' @param breaks numeric vector. Quantile breaks, default interquartile range: \code{c(0.25,0.5,0.75)}.
#' @param space See below for more details.
#' @return if space=FALSE then it returns a dataframe with quantiles of the distribution of dispersal distances; if space=TRUE (experimental) then it returns the invaded cells on the last day of model simulations for each of the iterations.
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

dici <- function(input_sim=NULL, coords=NULL, eval_date=NULL, breaks=c(0.25,0.5,0.75), space=FALSE) {
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > than number of simulated days...")
	}
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-dispersal model, set scale='lc' in dynamAedes()")
	} else {
		# Invaded cells
		inv_cells <- lapply(input_sim, function(x) {
			z1 <- lapply(x, function(d) {
				z <- data.frame(d[,which(colSums(d)>0)]); 
				colnames(z) <- which(colSums(d)>0)
			})
			z2 <- rbind(sapply(z1,unique))
		})
		# Distance from intro to each invaded cell per day
		if(!space){
			cdist <- lapply(1:length(inv_cells), function(x) {
				.meuc(c1=inv_cells[[x]], c2=inv_cells[[x]], coords=coords)
			})
			return(cbind.data.frame(.returndis(distl=cdist,days=eval_date, breaks=breaks),day=eval_date))
		}else{
			maxdate <- max(eval_date)
			outs <- lapply(inv_cells, function(x) {
				coords1 <- as.data.frame(coords)
				coords1$inv <- 0
				coords1$inv[x[[maxdate]]] <- 1
				names(coords1)=c("X", "Y", paste0("day",maxdate,"_inv_cells_Iteration"))
				coords1 <- rasterFromXYZ(coords1)
				return(coords1)
			})
			return(stack(outs))
		}
	}
}

# Ad-hoc functions
.returndis <- function(distl=NA, days=NA, breaks=breaks) {
	out <- do.call(rbind.data.frame, 
		lapply(days, function(x) {
			round(quantile(unlist(sapply(1:length(distl),
				function(y) {mean(if(x<=length(distl[[y]])) as.integer(distl[[y]][[x]]) else NA, na.rm=T)})), probs=breaks, na.rm=T),1)
		}))
	names(out) <- c("25%","50%","75%")
	return(out)
}

.euc <- function(xs, ys) { sqrt((xs[1]-xs[2])^2 + (ys[1]-ys[2])^2) }
.meuc <- function(c1, c2, coords) {
	c3 <- coords[unlist(c1[,1])[1],]
	outd <- list(NA)
	for ( rw in 1:length(c2) ) {
		outd[[rw]] <- sapply(unlist(c2[rw]), function(x) {
			.euc(c(as.numeric(c3[1]), as.numeric(coords[x,1])), c(as.numeric(c3[2]), as.numeric(coords[x,2])))
		})
	}
	return(outd)
}
