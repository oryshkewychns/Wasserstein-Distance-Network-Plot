library(transport)
library(poLCA)
library(qgraph)
library(RColorBrewer)

set.seed(1234)


data(carcinoma)

# run 3 cluster model
f <- cbind(A,B,C,D)~1
m1 <- poLCA(f, nclass = 3, data = carcinoma, nrep = 10, maxiter = 100000)


# extract predicted probabilities for each cluster
c1_probs <- list(A = matrix(m1$probs$A[1,], ncol = 2),
                 B = matrix(m1$probs$B[1,], ncol = 2),
                 C = matrix(m1$probs$C[1,], ncol = 2),
                 D = matrix(m1$probs$D[1,], ncol = 2))

c2_probs <- list(A = matrix(m1$probs$A[2,], ncol = 2),
                 B = matrix(m1$probs$B[2,], ncol = 2),
                 C = matrix(m1$probs$C[2,], ncol = 2),
                 D = matrix(m1$probs$D[2,], ncol = 2))

c3_probs <- list(A = matrix(m1$probs$A[3,], ncol = 2),
                 B = matrix(m1$probs$B[3,], ncol = 2),
                 C = matrix(m1$probs$C[3,], ncol = 2),
                 D = matrix(m1$probs$D[3,], ncol = 2))


# place all probabilities in a list
all_probs <- list(c1_probs, c2_probs, c3_probs)
names(all_probs) <- c("c1", "c2", "c3")


WD_LCA <- function(S, myprobs1, myprobs2){
  c1data <- poLCA.simdata(S, probs = myprobs1, nresp = NULL, niv = 0, P = 1)$dat
  c2data <- poLCA.simdata(S, probs = myprobs2, nresp = NULL, niv = 0, P = 1)$dat
  
  x1 <- pp(c1data)
  x2 <- pp(c2data)
  
  wd <- wasserstein(x1, x2, p = 2)
  wd
}

wd_mat <- sapply(all_probs, function(x) sapply(all_probs, function(y) WD_LCA(S = 5000, x, y)))


# rescale WD to 0-1 range
rescale_01 <- function(x){
  old_min <- min(x, na.rm = T)
  old_max <- max(x, na.rm = T)
  
  rescaled_mat <- ((x - old_min)/(old_max - old_min)) * (1 - 0) + 0
  return(rescaled_mat)
}

wd_mat_lower <- wd_mat
wd_mat_lower[upper.tri(wd_mat_lower, diag = T)] <- 0
wd_mat_lower[wd_mat_lower == 0] <- NA

wd_mat_01 <- rescale_01(wd_mat_lower)


qgraph(1- wd_mat_01,
       layout = "spring", 
       directed = F,
       edge.labels = T)

