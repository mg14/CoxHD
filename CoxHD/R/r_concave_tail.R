r.TailProbs <- function(eta, B, r) {
# TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
# We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
# s is -1/r
	MAXa <- 100000
	MINa <- 0.0001
	s <- -1/r
	etaB <- eta * B
	k_start <- (ceiling(2 * etaB) + 1)
	if(k_start > B) stop("eta is too large")
	
	a_vec <- rep(MAXa,B)
	
	Find.a <- function(prev_a) uniroot(Calc.a, lower = MINa, upper = prev_a, tol = .Machine$double.eps^0.75)$root
	Calc.a <- function(a) {
		denom <- sum((a + 0:k)^(-s))
		num <- sum((0:k) * (a + 0:k)^(-s))
		num / denom - etaB
	}
	
	for(k in k_start:B) a_vec[k] <- Find.a(a_vec[k-1])
	
	OptimInt <- function(a) {
		num <- (k + 1 - etaB) * sum((a + 0:(t-1))^(-s))
		denom <- sum((k + 1 - (0:k)) * (a + 0:k)^(-s))
		1 - num / denom
	}
	
	output <- rep(1, B)
	
	prev_k <- k_start
	for(t in k_start:B) {
		cur_optim <- rep(0, B)
		for (k in prev_k:(B-1)) cur_optim[k] <- optimize(f=OptimInt, lower = a_vec[k+1], upper = a_vec[k], maximum  = TRUE)$objective
		output[t] <- max(cur_optim)
		prev_k <- which.max(cur_optim)
	}
	return(output)
}

minD <- function(theta, B, r = c(-1/2, -1/4)) {
  pmin(c(rep(1, B), r.TailProbs(theta^2, B, r[1])), r.TailProbs(theta, 2*B, r[2]))
}
