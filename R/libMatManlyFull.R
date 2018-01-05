options(show.error.messages = FALSE)
eucl.dist <- function(X, Y){

	sqrt(sum(X - Y)^2)

}

Manly.bic <- function(logl, n, M){
	return(-2 * logl + M * log(n))
}

# EM algorithm
MatManly.init <- function(Y, X = NULL, K, la = NULL, nu = NULL, Mu.type = 0, Psi.type = 0, n.start = 10, tol = 1e-05){
	
	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong dimensionality T...\n")



	i <- 0
	
	Psi.inv <- array(NA, c(T, T, K))
	if(is.null(la)){la <- matrix(0.0, K, p)}
	if(is.null(nu)){nu <- matrix(0.0, K, T)}
	
	if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
	if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")
	if(K != dim(nu)[1]) stop("Inconsistent number of mixture components K...\n")	
	if(T != dim(nu)[2]) stop("Inconsistent dimensionality T...\n")


	if(!is.null(X)){
		if(n != dim(X)[3]) stop("Inconsistent number of observations n...\n")
		if(T != dim(X)[1]) stop("Inconsistent dimensionality T...\n")

		q <- dim(X)[2]
	}


	init <- list()
	if(Psi.type == 1){	
		cat("nu's are set to be zero...\n")

	}
	repeat{
		repeat{
			s <- sample(1:n,K)
			
			mat.Y <- t(apply(Y, 3, as.matrix))
			if((p>1) || (T>1)){
				centers <- mat.Y[s,]
			}
			else{
				centers <- mat.Y[s]

			}
			D <- NULL
			if(K == 1) {D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers))}
			else{
				if((p>1) || (T>1)){
					for (k in 1:K) D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers[k,]))
				}
				else{
					for (k in 1:K) D <- cbind(D, t(mat.Y - centers[k])^2)
				}
			}
			PI <- as.matrix(D == apply(D, 1, min)) * 1
			W.result <- NULL

			iif <- 0
			for (k in 1:K){
				index <- PI[,k] == 1
				if((p>1) && (T>1)){
					var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
					var.est <- var(var.est)
					A <- try(Psi.inv[,,k] <- solve(var.est))
				}
		
				else if((p == 1) && (T > 1)){
					var.est <- var(t(Y[,,index]))
					A <- try(Psi.inv[,,k] <- solve(var.est))
				}
				else if(T == 1){

					var.est <- var(as.vector(Y[,,index]))
					A <- try(Psi.inv[,,k] <- solve(var.est))
				}
				if(class(A) == "try-error"){iif <- 1}	
			}
			if(iif == 0){
				if(is.null(X)){
					y <- as.vector(Y)
					misc_int <- c(p, T, n, K, Mu.type)
				}
				else{
					y <- as.vector(Y)
					x <- as.vector(X)
					misc_int <- c(p, T, n, q, K)
					beta1 <- rep(0, q*p*K)
				}
				gamma1 <- PI
				la1 <- as.vector(la)
				nu1 <- as.vector(nu)
				invS1 <- rep(0, p*p*K)
				tau <- rep(0, K)

				misc_double <- c(tol, 0.0, 0.0)
				Mu1 <- rep(0, p*T*K)
				invPsi1 <- Psi.inv
				detS <- rep(0, K)
				detPsi <- rep(0, K)	

			
				if(Psi.type == 0){
					if(is.null(X)){
						try0 <- try(W <- .C("run_Mstep_Manly_Full", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), nu1 = as.double(nu1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau)))	
					}
					else{
						try0 <- try(W <- .C("run_Mstep_Manly_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), nu1 = as.double(nu1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau)))					
					}
				}					
				else if(Psi.type == 1){

					if(is.null(X)){
						try0 <- try(W <- .C("run_Mstep_Manly_AR", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau)))
					}
					else{

						try0 <- try(W <- .C("run_Mstep_Manly_AR_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau)))

					}
				}
				if ((class(try0) != "try-error")){
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					if(Psi.type == 0){
						W.result$nu <- matrix(W$nu1, nrow = K)
					}
					if(is.null(X)){
						W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
					}
					else{
						W.result$beta <- array(W$beta, dim = c(q, p, K))
					}
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
					W.result$detS <- W$detS
					W.result$detPsi <- W$detPsi
					W.result$Psi.type <- Psi.type
					if(is.null(X)){
						W.result$Mu.type <- Mu.type
					}

				}
			}
			if (is.null(W.result) || is.na(W.result$detS) || any(W.result$detS < 10^(-300))|| is.na(W$misc_double[2])){
			} else {
				i <- i + 1
				init[[i]] <- W.result					
				break
			}			
		}

		if (i == n.start) break

	}
	return(init)

}





MatManly.EM <- function(Y, X = NULL, initial = NULL, id = NULL, la = NULL, nu = NULL, tau = NULL, Mu = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Mu.type = 0, Psi.type = 0, tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE){

	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong dimensionality T...\n")


	if(!is.null(initial)){
		
		if(length(initial) < 1) stop("Wrong initialization...\n")

		K <- length(initial[[1]]$tau)

		best.loglik <- -Inf
		best.BIC <- Inf

		if(!is.null(X)){

			if(n != dim(X)[3]) stop("Inconsistent number of observations n...\n")
			if(T != dim(X)[1]) stop("Inconsistent dimensionality T...\n")

			q <- dim(X)[2]
		}

		for(i in 1:length(initial)){

			result <- NULL
			M <- NA


			y <- as.vector(Y)
			gamma1 <- rep(0, n*K)
			tau <- initial[[i]]$tau
			ll <- rep(0, 3)
			if(is.null(X)){
				misc_int <- c(p, T, n, K, max.iter, initial[[i]]$Mu.type)
				Mu1 <- initial[[i]]$Mu
			}
			else{
				misc_int <- c(p, T, n, q, K, max.iter)
				beta1 <- initial[[i]]$beta
			}
			misc_double <- c(tol, 0.0, 0.0)
			conv <- rep(0, 2)
			id <- rep(0, n)
			la1 <- initial[[i]]$la
			nu1 <- initial[[i]]$nu
			invS1 <- initial[[i]]$invS
			invPsi1 <- initial[[i]]$invPsi
			detS <- initial[[i]]$detS
			detPsi <- initial[[i]]$detPsi


			if(initial[[i]]$Psi.type == 0){
				if(is.null(X)){
					if(is.null(Mu1)) stop("Provided initialization invalid...\n")
					try0 <- try(result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

				}
				else{
					if(is.null(beta1)) stop("Provided initialization invalid...\n")
					x <- as.vector(X)
					try0 <- try(result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

				}
			}

			else if(initial[[i]]$Psi.type == 1){
				if(is.null(X)){
					if(is.null(Mu1)) stop("Provided initialization invalid...\n")
					try0 <- try(result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

				}
				else{
					if(is.null(beta1)) stop("Provided initialization invalid...\n")
					x <- as.vector(X)
					try0 <- try(result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

				}
			}

			try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
			try2 <- try(invPsi <- array(result$invPsi1, dim = c(T, T, K)))
				
			try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
			try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

			if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){

				if(!is.na(result$ll[1])){

					if ((result$ll[1] > best.loglik) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(result$id) > size.control)){
						best.loglik <- result$ll[1]


						best.la <- matrix(result$la1, nrow = K)
						if(initial[[i]]$Psi.type == 0){
							best.nu <- matrix(result$nu1, nrow = K)
							quantity <- best.nu[,T]
							best.nu <- best.nu - quantity
							best.la <- best.la + quantity
						}
						best.tau <- result$tau
						best.Sigma <- Sigma
						best.Psi <- Psi

						quantity <- apply(best.Psi, 3, det)
						quantity <- quantity^(1/T)
	
						for(k in 1:K){
							best.Psi[,,k] <- best.Psi[,,k] / quantity[k]		
							best.Sigma[,,k] <- best.Sigma[,,k] * quantity[k]
						}

						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]
						ind <- as.vector(best.la)==0
						if(initial[[i]]$Psi.type == 0){
							ind2 <- as.vector(best.nu[,1:(T-1)])==0
						}
						best.id <- result$id
						best.flag <- result$conv[2]
						best.Psi.type <- initial[[i]]$Psi.type
						if(initial[[i]]$Psi.type == 0){
							if(is.null(X)){
								best.Mu <- array(result$Mu1, dim = c(p, T, K))
								best.Mu.type <- initial[[i]]$Mu.type
								M <- K - 1 + K * p + K * T + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K
							}
							else{
								best.beta <- array(result$beta1, dim = c(q, p, K))
								M <- K - 1 + K * p + K * T + K * q* p + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K

							}
						}
						else if(initial[[i]]$Psi.type == 1){
							if(is.null(X)){
								best.Mu <- array(result$Mu1, dim = c(p, T, K))
								best.Mu.type <- initial[[i]]$Mu.type
								M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2 +K - sum(ind) 							}
							else{
								best.beta <- array(result$beta1, dim = c(q, p, K))
								M <- K - 1 + K * p + K * q * p + K * p * (p + 1) / 2 +K - sum(ind)
							}

						}
						best.BIC <- Manly.bic(best.loglik, n, M)
					}
				}

			}


			if(silent == FALSE){
				if(!is.null(result)){
							
					cat("Initialization", i, ":", "ll =", result$ll[1], "(best ll =", best.loglik, "BIC =", best.BIC, ")", "\n")
					
				}
				else{

					cat("Initialization", i, ":", "ll = NA", "(best ll =", best.loglik, "BIC =", best.BIC, ")", "\n")
				}

			}

	

		}

		if(best.loglik > - Inf){
			if(best.Psi.type == 0){
				if(is.null(X)){		
					ret <- list(la = best.la, nu = best.nu, tau = best.tau, Mu = best.Mu, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, Mu.type = best.Mu.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				}
				else{
					ret <- list(la = best.la, nu = best.nu, tau = best.tau, beta = best.beta, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				}
			}
			else if(best.Psi.type == 1){
				if(is.null(X)){		
					ret <- list(la = best.la, tau = best.tau, Mu = best.Mu, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, Mu.type = best.Mu.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				}
				else{
					ret <- list(la = best.la, tau = best.tau, beta = best.beta, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				}
		
			}
			class(ret) <- "MatManlyMix"
		}
		else{
			warning("The EM algorithm does not converge...\n")
			if(best.Psi.type == 0){
				if(is.null(X)){
					ret <- list(la = NULL, nu = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}
				else{				
					ret <- list(la = NULL, nu = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}
			}

			else if(best.Psi.type == 1){
				if(is.null(X)){
					ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}
				else{				
					ret <- list(la = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}
		
			}


		}


		return(ret)


	}


	else{
		if(is.null(X)){
			if(is.null(id) && (is.null(Mu) || is.null(tau) || is.null(Sigma) || is.null(Psi))) stop("Must provide one initialization method...\n")
		}
		else{
			q <- dim(X)[2]

			if(is.null(id) && (is.null(beta) || is.null(tau) || is.null(Sigma) || is.null(Psi))) stop("Must provide one initialization method...\n")
		}	
		if(!is.null(id)){

			K <- max(id)	
			if(K < 1) stop("Wrong number of mixture components K...\n")
			if(is.null(la)){
				la <- matrix(0.0, K, p)
			}
			if(is.null(nu)){
				nu <- matrix(0.0, K, T)
			}
			if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
			if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")
			if(K != dim(nu)[1]) stop("Inconsistent number of mixture components K...\n")	
			if(T != dim(nu)[2]) stop("Inconsistent dimensionality T...\n")

			if(n != length(id)) stop("Inconsistent number of observations n...\n")


			PI <- matrix(NA, n, K)

			for(i in 1:n){
				for(k in 1:K){

					PI[i,k] <- as.numeric(id[i] == k)
				}

			}



			inv.Psi <- array(NA, c(T, T, K))
			iif <- 0
			for (k in 1:K){

				index <- PI[,k] == 1
				var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
				var.est <- var(var.est)
				A <- try(inv.Psi[,,k] <- solve(var.est))
				if(class(A) == "try-error"){
					iif <- 1
					stop("Provided initialization invalid...\n")
					
				}	

			}
			y <- as.vector(Y)
			gamma1 <- PI
			la1 <- as.vector(la)
			nu1 <- as.vector(nu)
			invPsi1 <- inv.Psi
			tau <- rep(0, K)
			if(is.null(X)){
				Mu1 <- rep(0, p*T*K)
				misc_int <- c(p, T, n, K, max.iter, Mu.type)
			}
			else{
				x <- as.vector(X)
				misc_int <- c(p, T, n, q, K, max.iter)
				beta1 <- rep(0, q*p*K)
				
			}
			misc_double <- c(tol, 0.0, 0.0)				
			invS1 <- rep(0, p*p*K)
			detS <- rep(0, K)
			detPsi <- rep(0, K)
			if(Psi.type == 0){
				if(is.null(X)){
					W <- .C("run_Mstep_Manly_Full", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), nu1 = as.double(nu1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
					Mu1 <- array(W$Mu1, dim = c(p, T, K))
				}
				else{

					W <- .C("run_Mstep_Manly_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), nu1 = as.double(nu1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
					beta1 <- array(W$beta1, dim = c(q, p, K))
				}
			}
			else if(Psi.type == 1){
				if(is.null(X)){
					W <- .C("run_Mstep_Manly_AR", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
					Mu1 <- array(W$Mu1, dim = c(p, T, K))
				}
				else{

					W <- .C("run_Mstep_Manly_AR_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
					beta1 <- array(W$beta1, dim = c(q, p, K))
				}
			}

			if (any(W$detS < 10^(-300)) || is.na(W$misc_double[2])){
				stop("Provided initialization invalid...\n")
			}


			tau <- W$tau 
			ll <- rep(0, 3)
			conv <- rep(0, 2)
			id <- rep(0, n)
			la1 <- matrix(W$la1, nrow = K)
			if(Psi.type == 0){
				nu1 <- matrix(W$nu1, nrow = K)
			}
			invS1 <- array(W$invS1, dim = c(p, p, K))
			invPsi1 <- array(W$invPsi1, dim = c(T, T, K))
			detS <- W$detS
			detPsi <- W$detPsi
			if(Psi.type == 0){	
				if(is.null(X)){
					result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
				else{
					result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
			}
			else if(Psi.type == 1){	
				if(is.null(X)){
					result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
				else{
					result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
			}	


			la <- matrix(result$la1, nrow = K)
			if(initial[[i]]$Psi.type == 0){
				nu <- matrix(result$nu1, nrow = K)
				quantity <- nu[,T]
				nu <- nu - quantity
				la <- la + quantity
			}

			ind <- as.vector(la)==0
			if(initial[[i]]$Psi.type == 0){
				ind2 <- as.vector(nu[,1:(T-1)])==0
			}
			Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 	Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

			quantity <- apply(Psi, 3, det)
			quantity <- quantity^(1/T)
			for(k in 1:K){
				Psi[,,k] <- Psi[,,k] / quantity[k]			
				Sigma[,,k] <- Sigma[,,k] * quantity[k]
			}



			if(Psi.type == 0){					
				if(is.null(X)){
					M <- K - 1 + K * p + K * T + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K
					ret <- list(la = la, nu = nu, tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}	
				else{
					M <- K - 1 + K * p + K * T + K * p * q + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K
					ret <- list(la = la, nu = nu, tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}
			}
			else if(Psi.type == 1){

				if(is.null(X)){
					M <- K - 1 + K * p + K * T + K * p * T + K * p * (p + 1) / 2 +K - sum(ind) 
					ret <- list(la = la, tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}	
				else{
					M <- K - 1 + K * p + K * T + K * p * q + K * p * (p + 1) / 2 +K  - sum(ind) 
					ret <- list(la = la, tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}

			}
			class(ret) <- "MatManlyMix"

			if(ret$flag == 1){
				warning("The EM algorithm does not converge...\n")
				if(Psi.type == 0){
					if(is.null(X)){				
						ret <- list(la = NULL, nu = NULL,tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
					else{
						ret <- list(la = NULL, nu = NULL,tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
				}
				else if(Psi.type == 1){
					if(is.null(X)){				
						ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
					else{
						ret <- list(la = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
				}				
			}

			return(ret)

		}

		else{
			K <- length(tau)
			if(is.null(la)){la <- matrix(0.0, K, p)}
			if(is.null(nu)){nu <- matrix(0.0, K, T)}
			if(is.null(X)){
				equal.K <- c(dim(la)[1], dim(Mu)[3], dim(Sigma)[3], dim(Psi)[3])
				equal.p <- c(dim(la)[2], dim(Mu)[1], dim(Sigma)[1], dim(Sigma)[2])
				equal.T <- c(dim(Mu)[2], dim(Psi)[1], dim(Psi)[2])

				if (K < 1) stop("Wrong number of mixture components K...\n")
				if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
				if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
				if ((T != equal.T[1]) || (T != equal.T[2]) || (T != equal.T[3])) stop("Inconsistent dimensionality T...\n")
			}
			else{

				q <- dim(X)[2]

				equal.K <- c(dim(la)[1], dim(beta)[3], dim(Sigma)[3], dim(Psi)[3])
				equal.p <- c(dim(la)[2], dim(beta)[2], dim(Sigma)[1], dim(Sigma)[2])
				equal.T <- c(dim(Psi)[1], dim(Psi)[2])
				equal.q <- c(dim(beta)[1])
		

				if (K < 1) stop("Wrong number of mixture components K...\n")
				if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
				if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
				if ((T != equal.T[1]) || (T != equal.T[2])) stop("Inconsistent time T...\n")
				if (q != equal.q[1]) stop("Inconsistent number of variables q...\n")

			}	


			y <- as.vector(Y)
			gamma1 <- rep(0, n*K)
			ll <- rep(0, 3)
			misc_double <- c(tol, 0.0, 0.0)
			conv <- rep(0, 2)
			id <- rep(0, n)
			la1 <- la
			nu1 <- nu

			invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
			invPsi1 <- array(apply(Psi, 3, solve), dim = c(T,T,K))
			detS <- apply(Sigma, 3, det)
			detPsi <- apply(Psi, 3, det)

			if(Psi.type == 0){			
				if(is.null(X)){
					Mu1 <- Mu
					misc_int <- c(p, T, n, K, max.iter, Mu.type)
					result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
				else{
					x <- as.vector(X)
					beta1 <- beta
					misc_int <- c(p, T, n, q, K, max.iter)
					result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), nu1 = as.double(nu1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
			}
			else if(Psi.type == 1){			
				if(is.null(X)){
					Mu1 <- Mu
					misc_int <- c(p, T, n, K, max.iter, Mu.type)
					result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
				else{
					x <- as.vector(X)
					beta1 <- beta
					misc_int <- c(p, T, n, q, K, max.iter)
					result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

				}
			}



			la <- matrix(result$la1, nrow = K)
			if(Psi.type == 0){
				nu <- matrix(result$nu1, nrow = K)
				quantity <- nu[,T]
				nu <- nu - quantity
				la <- la + quantity
			}
			ind <- as.vector(la)==0
			if(Psi.type == 0){
				ind2 <- as.vector(nu[,1:(T-1)])==0
			}
			Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 	Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

			quantity <- apply(Psi, 3, det)
			quantity <- quantity^(1/T)
	
			for(k in 1:K){
				Psi[,,k] <- Psi[,,k] / quantity[k]			
				Sigma[,,k] <- Sigma[,,k] * quantity[k]
			}


			if(Psi.type == 0){				
				if(is.null(X)){
					M <- K - 1 + K * p + K * T + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K
					ret <- list(la = la, nu = nu, tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
				}
				else{
					M <- K - 1 + K * p + K * T + K * p * q + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind) -sum(ind2) - 2*K
					ret <- list(la = la, nu = nu, tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}
			}
			else if(Psi.type == 1){				
				if(is.null(X)){
					M <- K - 1 + K * p + K * T + K * p * T + K * p * (p + 1) / 2 +K - sum(ind) 
					ret <- list(la = la, tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
				}
				else{
					M <- K - 1 + K * p + K * T + K * p * q + K * p * (p + 1) / 2 +K - sum(ind)
					ret <- list(la = la, tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, Mu.type = Mu.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])

				}
			}
			
			class(ret) <- "MatManlyMix"

			if(ret$flag == 1){
				warning("The EM algorithm does not converge...\n")
				if(Psi.type == 0){
					if(is.null(X)){				
						ret <- list(la = NULL, nu = NULL,tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
					else{
						ret <- list(la = NULL, nu = NULL,tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
				}
				else if(Psi.type == 1){
					if(is.null(X)){				
						ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}
					else{
						ret <- list(la = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, Mu.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
					}

				}
			}
			return(ret)

		}


	}


}