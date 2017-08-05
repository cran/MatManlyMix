options(show.error.messages = FALSE)
eucl.dist <- function(X, Y){

	sqrt(sum(X - Y)^2)

}


# EM algorithm

MatManly.init <- function(Y, X = NULL, K, la = NULL, Psi.type = 0, n.start = 10, tol = 1e-05){
	
	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong time sequence T...\n")


	if(is.null(X)){
		i <- 0

		Sigma.inv <- array(NA, c(p, p, K))

		Psi.inv <- array(NA, c(T, T, K))

		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}
	
		if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
		if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")


		init <- list()


		repeat{

			

			#cat("Start", i, "\n")

			repeat{

				s <- sample(1:n, K)

				mat.Y <- t(apply(Y, 3, as.matrix))
				centers <- mat.Y[s,]

				D <- NULL
				if(K == 1) {

					D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers))
				}
				else{
					for (k in 1:K) D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers[k,]))
				}
				PI <- as.matrix(D == apply(D, 1, min)) * 1



				W.result <- NULL


				if(Psi.type == 0){
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){iif <- 1}	

					}
					if(iif == 0){
			
						y <- as.vector(Y)
						gamma1 <- PI
						la1 <- as.vector(la)
						invS1 <- Sigma.inv

						tau <- rep(0, K)
						misc_int <- c(p, T, n, K)
						misc_double <- c(tol, 0.0, 0.0)
						Mu1 <- rep(0, p*T*K)
						invPsi1 <- rep(0, T*T*K)
						detS <- rep(0, K)
						detPsi <- rep(0, K)
						Psi2 <- rep(0, K)


						W <- .C("run_Mstep_Manly_Full", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))

					
						W.result$tau <- W$tau 
						W.result$la <- matrix(W$la1, nrow = K)
						W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
						W.result$invS <- array(W$invS1, dim = c(p, p, K))
						W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
						W.result$detS <- W$detS
						W.result$detPsi <- W$detPsi
						W.result$Psi.type <- Psi.type

					}
				}

				else if(Psi.type == 1){
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){iif <- 1}	

					}
					if(iif == 0){
			
						y <- as.vector(Y)
						gamma1 <- PI
						la1 <- as.vector(la)
						invS1 <- Sigma.inv

						tau <- rep(0, K)
						misc_int <- c(p, T, n, K)
						misc_double <- c(tol, 0.0, 0.0)
						Mu1 <- rep(0, p*T*K)
						invPsi1 <- rep(0, T*T*K)
						detS <- rep(0, K)
						detPsi <- rep(0, K)
						Psi2 <- rep(0, K)



						W <- .C("run_Mstep_Manly_diag", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), Psi2 = as.double(Psi2), detS = as.double(detS), tau = as.double(tau))

						
						W.result$tau <- W$tau 
						W.result$la <- matrix(W$la1, nrow = K)
						W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
						W.result$invS <- array(W$invS1, dim = c(p, p, K))
						W.result$Psi2 <- W$Psi2
						W.result$detS <- W$detS
						W.result$Psi.type <- Psi.type

					}
				}



				else if(Psi.type == 2){
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Psi.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){iif <- 1}	

					}


					if(iif == 0){
			
						y <- as.vector(Y)
						gamma1 <- PI
						la1 <- as.vector(la)
						invS1 <- rep(0, p*p*K)

						tau <- rep(0, K)
						misc_int <- c(p, T, n, K)
						misc_double <- c(tol, 0.0, 0.0)
						Mu1 <- rep(0, p*T*K)
						invPsi1 <- Psi.inv
						detS <- rep(0, K)
						detPsi <- rep(0, K)


						W <- .C("run_Mstep_Manly_AR", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))

						
						W.result$tau <- W$tau 
						W.result$la <- matrix(W$la1, nrow = K)
						W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
						W.result$invS <- array(W$invS1, dim = c(p, p, K))
						W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
						W.result$detS <- W$detS
						W.result$detPsi <- W$detPsi
						W.result$Psi.type <- Psi.type

					}
				}




				if (is.null(W.result) || any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
					#cat("warning: M-step error produced...\n")
				} else {
					i <- i + 1
					init[[i]] <- W.result					
					break
				}
			
			
			}

			if (i == n.start) break

		}


	}

	else{

		i <- 0
		Sigma.inv <- array(NA, c(p, p, K))
		Psi.inv <- array(NA, c(T, T, K))
		Psi2 <- rep(NA, K)

		if(is.null(la)){
			la <- matrix(0.0, K, p)
		}
	


		if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
		if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")
		if(n != dim(X)[3]) stop("Inconsistent number of observations n...\n")
		if(T != dim(X)[1]) stop("Inconsistent time T...\n")

		q <- dim(X)[2]

		init <- list()
		repeat{

			

			#cat("Start", i, "\n")

			repeat{
				
				s <- sample(1:n, K)

				mat.Y <- t(apply(Y, 3, as.matrix))
				centers <- mat.Y[s,]

				D <- NULL
				if(K == 1) {

					D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers))
				}
				else{
					for (k in 1:K) D <- cbind(D, apply(mat.Y, 1, eucl.dist, centers[k,]))
				}
				PI <- as.matrix(D == apply(D, 1, min)) * 1


				W.result <- NULL
				if(Psi.type == 0 || Psi.type == 2){


					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Psi.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){iif <- 1}	
					}





					if(iif == 0){
			
						y <- as.vector(Y)
						x <- as.vector(X)
						gamma1 <- PI
						la1 <- as.vector(la)
						invS1 <- rep(0, p*p*K)

						tau <- rep(0, K)
						misc_int <- c(p, T, n, q, K)
						misc_double <- c(tol, 0.0, 0.0)
						beta1 <- rep(0, q*p*K)
						invPsi1 <- Psi.inv
						detS <- rep(0, K)
						detPsi <- rep(0, K)
					

						if(Psi.type == 0){
							W <- .C("run_Mstep_Manly_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
						}
						else{
							W <- .C("run_Mstep_Manly_AR_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))

						}
					
						W.result$tau <- W$tau 
						W.result$la <- matrix(W$la1, nrow = K)
						W.result$beta <- array(W$beta, dim = c(q, p, K))
						W.result$invS <- array(W$invS1, dim = c(p, p, K))
						W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
						W.result$detS <- W$detS
						W.result$detPsi <- W$detPsi
						W.result$Psi.type <- Psi.type

						if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
							break
							#cat("warning: M-step error produced...\n")
						} else {
							i <- i + 1
							init[[i]] <- W.result
							break
						}

					} 





				}

				else if(Psi.type == 1){
		

					iif <- 0
					for (k in 1:K){
				
						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){iif <- 1}	

					}
			
					if(iif == 0){
			


						y <- as.vector(Y)
						x <- as.vector(X)
						gamma1 <- PI
						la1 <- as.vector(la)
						invS1 <- Sigma.inv

						tau <- rep(0, K)
						misc_int <- c(p, T, n, q, K)
						misc_double <- c(tol, 0.0, 0.0)
						beta1 <- rep(0, q*p*K)
						Psi2 <- rep(0, K)
						detS <- rep(0, K)
					
			

						W <- .C("run_Mstep_Manly_diag_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), Psi2 = as.double(Psi2), detS = as.double(detS), tau = as.double(tau))


			
						W.result$tau <- W$tau 
						W.result$la <- matrix(W$la1, nrow = K)
						W.result$beta <- array(W$beta, dim = c(q, p, K))
						W.result$invS <- array(W$invS1, dim = c(p, p, K))
						W.result$Psi2 <- W$Psi2
						W.result$detS <- W$detS
						W.result$Psi.type <- Psi.type
					
				
						#cat(W.result$Psi2, "\n")

						#cat(W.result$detS, "\n")
						if (is.null(W.result) || any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
							break
							#cat("warning: M-step error produced...\n")
						} else {
							i <- i + 1
							init[[i]] <- W.result
							
							break
						}
					}

				}



			
			}

			if (i == n.start) break

		}


	}

	return(init)

}


MatManly.EM <- function(Y, X = NULL, initial = NULL, id = NULL, la = NULL, tau = NULL, Mu = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = 0, tol = 1e-05, max.iter = 1000, size.control = 0, silent = TRUE){

	A <- dim(Y)
	p <- A[1]
	T <- A[2]
	n <- A[3]

	if(n < 1) stop("Wrong number of observations n...\n")
	if(p < 1) stop("Wrong dimensionality p...\n")
	if(T < 1) stop("Wrong time sequence T...\n")


	if(is.null(X)){

		if(!is.null(initial)){
		
			if(length(initial) < 1) stop("Wrong initialization...\n")

			K <- length(initial[[1]]$tau)

			best.loglik <- -Inf
			best.BIC <- Inf


			for(i in 1:length(initial)){

				result <- NULL
				M <- NA

				if(initial[[i]]$Psi.type == 0){
					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					Mu1 <- initial[[i]]$Mu
					invS1 <- initial[[i]]$invS
					invPsi1 <- initial[[i]]$invPsi
					detS <- initial[[i]]$detS
					detPsi <- initial[[i]]$detPsi

					try0 <- try(result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					try2 <- try(invPsi <- array(result$invPsi1, dim = c(T, T, K)))
					
					try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
					try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

					if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){

					if(!is.na(result$ll[1])){

					if ((result$ll[1] > best.loglik) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(result$id) > size.control)){


						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.Mu <- array(result$Mu1, dim = c(p, T, K))
						best.Sigma <- Sigma
						best.Psi <- Psi
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]

						ind <- as.vector(best.la)==0
						M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)


						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]
						best.Psi.type <- initial[[i]]$Psi.type
					}
					}

					}


				}
				else if(initial[[i]]$Psi.type == 1){

					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					Mu1 <- initial[[i]]$Mu
					invS1 <- initial[[i]]$invS
					detS <- initial[[i]]$detS
					Psi2 <- initial[[i]]$Psi2

					try0 <- try(result <- .C("run_Mat_Manly_diag", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					
					try2 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))

					if ((class(try0) != "try-error") && (class(try1) != "try-error")){

					if(!is.na(result$ll[1])){


					if ((result$ll[1] > best.loglik) && (class(try2) != "try-error") && all(table(result$id) > size.control)){

						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.Mu <- array(result$Mu1, dim = c(p, T, K))
						best.Sigma <- Sigma
						best.Psi <- result$Psi2
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]

						ind <- as.vector(best.la)==0
		
						M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K - sum(ind)

						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]
						best.Psi.type <- initial[[i]]$Psi.type
					}
					}

					}
				}







				else if(initial[[i]]$Psi.type == 2){
					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					Mu1 <- initial[[i]]$Mu
					invS1 <- initial[[i]]$invS
					invPsi1 <- initial[[i]]$invPsi
					detS <- initial[[i]]$detS
					detPsi <- initial[[i]]$detPsi

					try0 <- try(result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))

					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					try2 <- try(invPsi <- array(result$invPsi1, dim = c(T, T, K)))
					
					try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
					try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))


					if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){
					if(!is.na(result$ll[1])){


					if ((result$ll[1] > best.loglik) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(result$id) > size.control)){

						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.Mu <- array(result$Mu1, dim = c(p, T, K))
						best.Sigma <- Sigma
						best.Psi <- Psi
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]

						ind <- as.vector(best.la)==0
						M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * 2 - sum(ind)

						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]
						best.Psi.type <- initial[[i]]$Psi.type

					}

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
					
				ret <- list(la = best.la, tau = best.tau, Mu = best.Mu, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				class(ret) <- "MatManlyMix"
			}
			else{
				warning("The EM algorithm does not converge...\n")
				ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)

			}


			return(ret)


		}

		else{
			if(is.null(id) && (is.null(Mu) || is.null(tau) || is.null(Sigma) || is.null(Psi))) stop("Must provide one initialization method...\n")
			if(!is.null(id)){

				K <- max(id)	
			
				if(K < 1) stop("Wrong number of mixture components K...\n")
				if(is.null(la)){
					la <- matrix(0.0, K, p)
				}
				if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
				if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")

				if(n != length(id)) stop("Inconsistent number of observations n...\n")


				PI <- matrix(NA, n, K)

				for(i in 1:n){
					for(k in 1:K){

						PI[i,k] <- as.numeric(id[i] == k)
					}

				}

				if(Psi.type == 0){

					Sigma.inv <- array(NA, c(p, p, K))
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){
							iif <- 1
							stop("Provided initialization invalid...\n")
						
						}	

					}


					y <- as.vector(Y)
					gamma1 <- PI
					la1 <- as.vector(la)
					invS1 <- Sigma.inv
					tau <- rep(0, K)
					misc_int <- c(p, T, n, K)
					misc_double <- c(tol, 0.0, 0.0)
					Mu1 <- rep(0, p*T*K)
					invPsi1 <- rep(0, T*T*K)
					detS <- rep(0, K)
					detPsi <- rep(0, K)
				
					W <- .C("run_Mstep_Manly_Full", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))
					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
					W.result$detS <- W$detS
					W.result$detPsi <- W$detPsi
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					Mu1 <- W.result$Mu
					invS1 <- W.result$invS
					invPsi1 <- W.result$invPsi
					detS <- W.result$detS
					detPsi <- W.result$detPsi


					result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")



					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)
					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))



				}

				else if(Psi.type == 1){




					Sigma.inv <- array(NA, c(p, p, K))
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){
							iif <- 1
							stop("Provided initialization invalid...\n")
						
						}	

					}


					y <- as.vector(Y)
					gamma1 <- PI
					la1 <- as.vector(la)
					invS1 <- Sigma.inv
					tau <- rep(0, K)
					misc_int <- c(p, T, n, K)
					misc_double <- c(tol, 0.0, 0.0)
					Mu1 <- rep(0, p*T*K)
					detS <- rep(0, K)
					detPsi <- rep(0, K)
					Psi2 <- rep(0, K)


					W <- .C("run_Mstep_Manly_diag", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), Psi2 = as.double(Psi2), detS = as.double(detS), tau = as.double(tau))
	

					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$Psi2 <- W$Psi2
					W.result$detS <- W$detS
				
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					Mu1 <- W.result$Mu
					invS1 <- W.result$invS
					Psi2 <- W.result$Psi2
					detS <- W.result$detS
		


					result <- .C("run_Mat_Manly_diag", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K - sum(ind)
					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- result$Psi2

				}





				else if(Psi.type == 2){

					Psi.inv <- array(NA, c(T, T, K))
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Psi.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){
							iif <- 1
							stop("Provided initialization invalid...\n")
						
						}	

					}


					y <- as.vector(Y)
					gamma1 <- PI
					la1 <- as.vector(la)
					invS1 <- rep(0, p*p*K)
					tau <- rep(0, K)
					misc_int <- c(p, T, n, K)
					misc_double <- c(tol, 0.0, 0.0)
					Mu1 <- rep(0, p*T*K)
					invPsi1 <- Psi.inv
					detS <- rep(0, K)
					detPsi <- rep(0, K)


					W <- .C("run_Mstep_Manly_AR", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), Mu1 = as.double(Mu1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))

					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$Mu <- array(W$Mu1, dim = c(p, T, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
					W.result$detS <- W$detS
					W.result$detPsi <- W$detPsi
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					Mu1 <- W.result$Mu
					invS1 <- W.result$invS
					invPsi1 <- W.result$invPsi
					detS <- W.result$detS
					detPsi <- W.result$detPsi


					result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")


					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * 2 - sum(ind)
					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))


				}


				ret <- list(la = matrix(result$la1, nrow = K), tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
	
				class(ret) <- "MatManlyMix"

				if(ret$flag == 1){
					warning("The EM algorithm does not converge...\n")
					ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}

				return(ret)

			}
			else{



				if(Psi.type == 0){
					K <- length(tau)
					if(is.null(la)){
						la <- matrix(0.0, K, p)
					}

					equal.K <- c(dim(la)[1], dim(Mu)[3], dim(Sigma)[3], dim(Psi)[3])
					equal.p <- c(dim(la)[2], dim(Mu)[1], dim(Sigma)[1], dim(Sigma)[2])
					equal.T <- c(dim(Mu)[2], dim(Psi)[1], dim(Psi)[2])

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					if ((T != equal.T[1]) || (T != equal.T[2]) || (T != equal.T[3])) stop("Inconsistent time T...\n")


					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					Mu1 <- Mu
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					invPsi1 <- array(apply(Psi, 3, solve), dim = c(T,T,K))
					detS <- apply(Sigma, 3, det)
					detPsi <- apply(Psi, 3, det)


					result <- .C("run_Mat_Manly_Full", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}

				else if(Psi.type == 1){

					K <- length(tau)
					if(is.null(la)){
						la <- matrix(0.0, K, p)
					}

					equal.K <- c(dim(la)[1], dim(Mu)[3], dim(Sigma)[3], length(Psi))
					equal.p <- c(dim(la)[2], dim(Mu)[1], dim(Sigma)[1], dim(Sigma)[2])

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					

					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					Mu1 <- Mu
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					Psi2 <- Psi
					detS <- apply(Sigma, 3, det)

					result <- .C("run_Mat_Manly_diag", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2 + K - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- Psi2

				}



				else if(Psi.type == 2){
					K <- length(tau)
					if(is.null(la)){
						la <- matrix(0.0, K, p)
					}

					equal.K <- c(dim(la)[1], dim(Mu)[3], dim(Sigma)[3], dim(Psi)[3])
					equal.p <- c(dim(la)[2], dim(Mu)[1], dim(Sigma)[1], dim(Sigma)[2])
					equal.T <- c(dim(Mu)[2], dim(Psi)[1], dim(Psi)[2])

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					if ((T != equal.T[1]) || (T != equal.T[2]) || (T != equal.T[3])) stop("Inconsistent time T...\n")


					y <- as.vector(Y)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					Mu1 <- Mu
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					invPsi1 <- array(apply(Psi, 3, solve), dim = c(T,T,K))
					detS <- apply(Sigma, 3, det)
					detPsi <- apply(Psi, 3, det)


					result <- .C("run_Mat_Manly_AR", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), Mu1 = as.double(Mu1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * T + K * p * (p + 1) / 2  + K * 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}


				ret <- list(la = matrix(result$la1, nrow = K), tau = result$tau, Mu = array(result$Mu1, dim = c(p, T, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
				class(ret) <- "MatManlyMix"

				if(ret$flag == 1){
					warning("The EM algorithm does not converge...\n")
					ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}

				return(ret)


			}


		}


	}

	else{

		if(n != dim(X)[3]) stop("Inconsistent number of observations n...\n")
		if(T != dim(X)[1]) stop("Inconsistent time T...\n")

		q <- dim(X)[2]


		if(!is.null(initial)){
		
			if(length(initial) < 1) stop("Wrong initialization...\n")

			K <- length(initial[[1]]$tau)

			best.loglik <- -Inf
			best.BIC <- Inf

			for(i in 1:length(initial)){


				result <- NULL
				M <- NA

				if(initial[[i]]$Psi.type == 0){

					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					beta1 <- initial[[i]]$beta
					invS1 <- initial[[i]]$invS
					invPsi1 <- initial[[i]]$invPsi
					detS <- initial[[i]]$detS
					detPsi <- initial[[i]]$detPsi
		


					try0 <- try(result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))


					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					try2 <- try(invPsi <- array(result$invPsi1, dim = c(T, T, K)))
					
					try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
					try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

					
					if ((class(try0) != "try-error") && (class(try1) != "try-error") && (class(try2) != "try-error")){

					if(!is.na(result$ll[1])){	
					if ((result$ll[1] > best.loglik) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(result$id) > size.control)){

						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.beta <- array(result$beta1, dim = c(q, p, K))
						best.Sigma <- Sigma
						best.Psi <- Psi
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]
						best.Psi.type <- initial[[i]]$Psi.type
						ind <- as.vector(best.la)==0
						M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)

						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]

					}
					}
					}

				}

				else if(initial[[i]]$Psi.type == 1){

					#cat("hi this is correct")
					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					beta1 <- initial[[i]]$beta
					invS1 <- initial[[i]]$invS
					detS <- initial[[i]]$detS
					Psi2 <- initial[[i]]$Psi2
		

					try0 <- try(result <- .C("run_Mat_Manly_diag_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))


					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					
					try2 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
					
					if ((class(try0) != "try-error") && (class(try1) != "try-error")){

					if(!is.na(result$ll[1])){

	
					if ((result$ll[1] > best.loglik) && (class(try2) != "try-error") && all(table(result$id) > size.control)){

						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.beta <- array(result$beta1, dim = c(q, p, K))
						best.Sigma <- Sigma
						best.Psi <- result$Psi2
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]
						best.Psi.type <- initial[[i]]$Psi.type
						ind <- as.vector(best.la)==0
						M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K - sum(ind)

						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]

					}
					}
					}




				}



				else if(initial[[i]]$Psi.type == 2){

					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					tau <- initial[[i]]$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- initial[[i]]$la
					beta1 <- initial[[i]]$beta
					invS1 <- initial[[i]]$invS
					invPsi1 <- initial[[i]]$invPsi
					detS <- initial[[i]]$detS
					detPsi <- initial[[i]]$detPsi
		


					try0 <- try(result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix"))


					try1 <- try(invS <- array(result$invS1, dim = c(p, p, K)))
					try2 <- try(invPsi <- array(result$invPsi1, dim = c(T, T, K)))
					
					try3 <- try(Sigma <- array(apply(invS, 3, solve), dim = c(p,p,K)))
					try4 <- try(Psi <- array(apply(invPsi, 3, solve), dim = c(T,T,K)))

					if ((class(try0) != "try-error") && (class(try1) != "try-error")&& (class(try2) != "try-error")){

					if(!is.na(result$ll[1])){
	
					if ((result$ll[1] > best.loglik) && (class(try3) != "try-error") && (class(try4) != "try-error") && all(table(result$id) > size.control)){

						best.loglik <- result$ll[1]

						best.la <- matrix(result$la1, nrow = K)
						best.tau <- result$tau
						best.beta <- array(result$beta1, dim = c(q, p, K))
						best.Sigma <- Sigma
						best.Psi <- Psi
						best.gamma <- matrix(result$gamma1, nrow = n)
						best.iter <- result$conv[1]
						best.Psi.type <- initial[[i]]$Psi.type
						ind <- as.vector(best.la)==0
						M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * 2 - sum(ind)

						best.BIC <- Manly.bic(best.loglik, n, M)
						best.id <- result$id
						best.flag <- result$conv[2]

					}
					}
					}
				}

				if(silent == FALSE){
					if(!is.null(result)){
							
						cat("Initialization", i, ":", "ll =", result$ll[1], "(best ll =", best.loglik, "BIC =", best.BIC, ")", "\n")
					
					}
					else{

						cat("Initialization ", i, ":", "ll = NA", "(best ll = ", best.loglik, "BIC =", best.BIC, ")", "\n")
					}

				}



			}

			if(best.loglik > - Inf){
	
				ret <- list(la = best.la, tau = best.tau, beta = best.beta, Sigma = best.Sigma, Psi = best.Psi, Psi.type = best.Psi.type, gamma = best.gamma, id = best.id, ll = best.loglik, bic = best.BIC, iter = best.iter, flag = best.flag)
				class(ret) <- "MatManlyMix"
			}
			else{
				warning("The EM algorithm does not converge...\n")
				ret <- list(la = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)

			}
			return(ret)


		}

		else{
			if(is.null(id) && (is.null(beta) || is.null(tau) || is.null(Sigma) || is.null(Psi))) stop("Must provide one initialization method...\n")
			if(!is.null(id)){

				K <- max(id)	
			
				if(K < 1) stop("Wrong number of mixture components K...\n")
				if(is.null(la)){
					la <- matrix(0.0, K, p)
				}
				if(K != dim(la)[1]) stop("Inconsistent number of mixture components K...\n")	
				if(p != dim(la)[2]) stop("Inconsistent dimensionality p...\n")

				if(n != length(id)) stop("Inconsistent number of observations n...\n")


				PI <- matrix(NA, n, K)

				for(i in 1:n){
					for(k in 1:K){

						PI[i,k] <- as.numeric(id[i] == k)
					}

				}
				

				if(Psi.type == 0){


					Psi.inv <- array(NA, c(T, T, K))
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Psi.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){
							iif <- 1
							stop("Provided initialization invalid...\n")
						
						}	

					}


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					la1 <- as.vector(la)
					invS1 <- rep(0, p*p*K)
					tau <- rep(0, K)
					misc_int <- c(p, T, n, q, K)
					misc_double <- c(tol, 0.0, 0.0)
					beta1 <- rep(0, q*p*K)
					invPsi1 <- Psi.inv
					detS <- rep(0, K)
					detPsi <- rep(0, K)


					W <- .C("run_Mstep_Manly_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))



					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$beta <- array(W$beta1, dim = c(q, p, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
					W.result$detS <- W$detS
					W.result$detPsi <- W$detPsi
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					beta1 <- W.result$beta
					invS1 <- W.result$invS
					invPsi1 <- W.result$invPsi
					detS <- W.result$detS
					detPsi <- W.result$detPsi


					result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")


					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}

				else if(Psi.type == 1){

					Sigma.inv <- array(NA, c(p, p, K))
					iif <- 0
					for (k in 1:K){
				
						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 1, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Sigma.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){

							iif <- 1
							stop("Provided initialization invalid...\n")
	
						}	

					}
		


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					la1 <- as.vector(la)
					tau <- rep(0, K)
					misc_int <- c(p, T, n, q, K)
					misc_double <- c(tol, 0.0, 0.0)
					beta1 <- rep(0, q*p*K)
					invS1 <- Sigma.inv
					detS <- rep(0, K)
					Psi2 <- rep(0, K)
		


					W <- .C("run_Mstep_Manly_diag_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), Psi2 = as.double(Psi2), detS = as.double(detS),  tau = as.double(tau))


					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$beta <- array(W$beta1, dim = c(q, p, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$Psi2 <- W$Psi2
					W.result$detS <- W$detS
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					beta1 <- W.result$beta
					invS1 <- W.result$invS
					Psi2 <- W.result$Psi2
					detS <- W.result$detS


					result <- .C("run_Mat_Manly_diag_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- result$Psi2


				}


				else if(Psi.type == 2){


					Psi.inv <- array(NA, c(T, T, K))
					iif <- 0
					for (k in 1:K){

						index <- PI[,k] == 1

						var.est <- apply(Y[,,index], 2, as.matrix, byrow = TRUE)
						var.est <- var(var.est)
						A <- try(Psi.inv[,,k] <- solve(var.est))
						if(class(A) == "try-error"){
							iif <- 1
							stop("Provided initialization invalid...\n")
						
						}	

					}


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					la1 <- as.vector(la)
					invS1 <- rep(0, p*p*K)
					tau <- rep(0, K)
					misc_int <- c(p, T, n, q, K)
					misc_double <- c(tol, 0.0, 0.0)
					beta1 <- rep(0, q*p*K)
					invPsi1 <- Psi.inv
					detS <- rep(0, K)
					detPsi <- rep(0, K)


					W <- .C("run_Mstep_Manly_AR_Reg", y = as.double(y), misc_double = as.double(misc_double), misc_int = as.integer(misc_int), gamma1 = as.double(gamma1), la1 = as.double(la1), invS1 = as.double(invS1), x = as.double(x), beta1 = as.double(beta1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), tau = as.double(tau))



					W.result <- NULL
					W.result$tau <- W$tau 
					W.result$la <- matrix(W$la1, nrow = K)
					W.result$beta <- array(W$beta1, dim = c(q, p, K))
					W.result$invS <- array(W$invS1, dim = c(p, p, K))
					W.result$invPsi <- array(W$invPsi1, dim = c(T, T, K))
					W.result$detS <- W$detS
					W.result$detPsi <- W$detPsi
					if (any(W.result$detS < 10^(-300)) || is.na(W$misc_double[2])){
						stop("Provided initialization invalid...\n")
					}


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- PI
					tau <- W.result$tau
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- W.result$la
					beta1 <- W.result$beta
					invS1 <- W.result$invS
					invPsi1 <- W.result$invPsi
					detS <- W.result$detS
					detPsi <- W.result$detPsi


					result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")


					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}



				ret <- list(la = matrix(result$la1, nrow = K), tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
				class(ret) <- "MatManlyMix"

				if(ret$flag == 1){
					warning("The EM algorithm does not converge...\n")
					ret <- list(la = NULL, tau = NULL, Mu = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}
				return(ret)

			}
			else{
				K <- length(tau)
				if(is.null(la)){
					la <- matrix(0.0, K, p)
				}


				if(Psi.type == 0){
					equal.K <- c(dim(la)[1], dim(beta)[3], dim(Sigma)[3], dim(Psi)[3])
					equal.p <- c(dim(la)[2], dim(beta)[2], dim(Sigma)[1], dim(Sigma)[2])
					equal.T <- c(dim(Psi)[1], dim(Psi)[2])
					equal.q <- c(dim(beta)[1])
		

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					if ((T != equal.T[1]) || (T != equal.T[2])) stop("Inconsistent time T...\n")
					if (q != equal.q[1]) stop("Inconsistent number of variables q...\n")


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					beta1 <- beta
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					invPsi1 <- array(apply(Psi, 3, solve), dim = c(T,T,K))
					detS <- apply(Sigma, 3, det)
					detPsi <- apply(Psi, 3, det)


					result <- .C("run_Mat_Manly_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * T * (T + 1) / 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}

				else if(Psi.type == 1){

					equal.K <- c(dim(la)[1], dim(beta)[3], dim(Sigma)[3], length(Psi))
					equal.p <- c(dim(la)[2], dim(beta)[2], dim(Sigma)[1], dim(Sigma)[2])
					equal.q <- c(dim(beta)[1])
		

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					if (q != equal.q[1]) stop("Inconsistent number of variables q...\n")


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					beta1 <- beta
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					detS <- apply(Sigma, 3, det)
					Psi2 <- Psi
				

					result <- .C("run_Mat_Manly_diag_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), Psi2 = as.double(Psi2), detS = as.double(detS), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")


					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- result$Psi2


				}



				else if(Psi.type == 2){
					equal.K <- c(dim(la)[1], dim(beta)[3], dim(Sigma)[3], dim(Psi)[3])
					equal.p <- c(dim(la)[2], dim(beta)[2], dim(Sigma)[1], dim(Sigma)[2])
					equal.T <- c(dim(Psi)[1], dim(Psi)[2])
					equal.q <- c(dim(beta)[1])
		

					if (K < 1) stop("Wrong number of mixture components K...\n")
					if ((K != equal.K[1]) || (K != equal.K[2]) || (K != equal.K[3]) || (K != equal.K[4])) stop("Inconsistent number of mixture components K...\n")
					if ((p != equal.p[1]) || (p != equal.p[2]) || (p != equal.p[3]) || (p != equal.p[4])) stop("Inconsistent number of dimensionality p...\n")
					if ((T != equal.T[1]) || (T != equal.T[2])) stop("Inconsistent time T...\n")
					if (q != equal.q[1]) stop("Inconsistent number of variables q...\n")


					y <- as.vector(Y)
					x <- as.vector(X)
					gamma1 <- rep(0, n*K)
					ll <- rep(0, 3)
					misc_int <- c(p, T, n, q, K, max.iter)
					misc_double <- c(tol, 0.0, 0.0)
					conv <- rep(0, 2)
					id <- rep(0, n)
					la1 <- la
					beta1 <- beta
					invS1 <- array(apply(Sigma, 3, solve), dim = c(p,p,K))
					invPsi1 <- array(apply(Psi, 3, solve), dim = c(T,T,K))
					detS <- apply(Sigma, 3, det)
					detPsi <- apply(Psi, 3, det)


					result <- .C("run_Mat_Manly_AR_Reg", y = as.double(y), misc_int = as.integer(misc_int), misc_double = as.double(misc_double), tau = as.double(tau), la1 = as.double(la1), x = as.double(x), beta1 = as.double(beta1), invS1 = as.double(invS1), invPsi1 = as.double(invPsi1), detS = as.double(detS), detPsi = as.double(detPsi), gamma1 = as.double(gamma1), id = as.integer(id), ll = as.double(ll), conv = as.integer(conv), PACKAGE = "MatManlyMix")

					ind <- as.vector(matrix(result$la1, nrow = K))==0
					M <- K - 1 + K * p + K * p * q + K * p * (p + 1) / 2  + K * 2 - sum(ind)

					Sigma <- array(apply(array(result$invS1, dim = c(p, p, K)), 3, solve), dim = c(p,p,K))
		 			Psi <- array(apply(array(result$invPsi1, dim = c(T, T, K)), 3, solve), dim = c(T,T,K))

				}


				ret <- list(la = matrix(result$la1, nrow = K), tau = result$tau, beta = array(result$beta1, dim = c(q, p, K)), Sigma = Sigma, Psi = Psi, Psi.type = Psi.type, gamma = matrix(result$gamma1, nrow = n), id = result$id, ll = result$ll[1], bic = Manly.bic(result$ll[1], n, M), iter = result$conv[1], flag = result$conv[2])
				class(ret) <- "MatManlyMix"

				if(ret$flag == 1){
					warning("The EM algorithm does not converge...\n")
					ret <- list(la = NULL, tau = NULL, beta = NULL, Sigma = NULL, Psi = NULL, Psi.type = NULL, gamma = NULL, id = NULL, ll = NULL, bic = NULL, iter = NULL, flag = 1)
				}

				return(ret)


			}


		}




	}

}