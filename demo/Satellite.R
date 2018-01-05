set.seed(123)
library(MixSim)
#Application to dataset AIS
data("Satellite")

K <- 3
p <- 4
T <- 9
Y <- Satellite$Y
trueid <- Satellite$id

AR_comparison <- rep(NA, 2)
CP_comparison <- rep(NA, 2)
BIC_comparison <- rep(NA, 2)

names(AR_comparison) <- c("mGm", "mMm")
names(CP_comparison) <- c("mGm", "mMm")
names(BIC_comparison) <- c("mGm", "mMm")

K <- 3

init <- MatManly.init(Y, K = K, n.start = 10)
G <- MatManly.EM(Y, initial = init, max.iter = 1000)
CP_comparison[1] <- ClassProp(trueid, G$id)
AR_comparison[1] <- RandIndex(trueid, G$id)$AR
BIC_comparison[1] <- G$bic

init <- MatManly.init(Y, K = K, la = matrix(0.1, K, p), nu = matrix(0.1, K, T), n.start = 10)
M <- MatManly.EM(Y, initial = init, max.iter = 1000)
CP_comparison[2] <- ClassProp(trueid, M$id)
AR_comparison[2] <- RandIndex(trueid, M$id)$AR
BIC_comparison[2] <- M$bic



CP_comparison
AR_comparison
BIC_comparison

