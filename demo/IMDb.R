set.seed(1234)

#Application to dataset IMDb
data("IMDb")
Y <- IMDb$Y

n <- 105
p <- 2
T <- 4
BIC_comparison  <- matrix(NA, 5, 2)

colnames(BIC_comparison) <- c("mGm", "mMm")



for(K in 1:5){
	init <- MatManly.init(Y, K = K, la = matrix(0, K, p), nu = matrix(0, K, T), n.start = 50)
	G2 <- MatManly.EM(Y, initial = init, max.iter = 1000, tol = 1e-05, size.control = 10)
	BIC_comparison[K,1] <- G2$bic

	init <- MatManly.init(Y, K = K, la = matrix(0.1, K, p), nu = matrix(0.1, K, T), n.start = 50)
	M2 <- MatManly.EM(Y, initial = init, max.iter = 1000, tol = 1e-05, size.control = 10)
	BIC_comparison[K,2] <- M2$bic

}

BIC_comparison