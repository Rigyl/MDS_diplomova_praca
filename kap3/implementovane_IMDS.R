##### 1. MOJE IMDS #####

cross_dist <- function(Xa, Xb) {
  m <- nrow(Xa)
  l <- nrow(Xb)
  A <- matrix(0, m, l)
  
  # vypocet vzdialenosti medzi dvoma maticami
  for (i in 1:m) {
    diffs <- sweep(Xb, 2, Xa[i, ], "-")
    A[i, ] <- sqrt(rowSums(diffs^2))
  }
  A
}

imds <- function(X, l = 200, r = 2) {
  
  X0 <- as.matrix(X)      # pre istotu prevediem na maticu, ak by to bol data.frame
  n <- nrow(X0)
  
  if (l < 2) stop("l musi byt aspon 2")
  if (l > n) stop("l nemoze byt vacsie ako n")
  if (r < 1) stop("r musi byt aspon 1")
  
  # nahodna vzorka velkosti l
  idx1 <- sample.int(n, l, replace = FALSE)
  X1_data <- X0[idx1,]
  
  # matica vzdialenosti D1
  D1 <- as.matrix(dist(X1_data))
  
  # classical MDS na D1
  X1 <- cmdscale(D1, k = r)
  
  # vypocet Q1 a q1 podla vzorca
  P <- diag(l) - (1 / l) * (matrix(1, l, 1) %*% matrix(1, 1, l))
  Q1 <- (-0.5) * P %*% (D1^2) %*% t(P)
  q1 <- diag(Q1)
  
  # kovariancna matica S1
  S1 <- cov(X1)
  S1_inv <- solve(S1)
  
  # zvysne body rozdelim na bloky (velkost max l)
  idx_rest <- (1:n)[-idx1]
  blocks <- split(idx_rest, ceiling(seq_along(idx_rest) / l))
  
  # sem budem skladat finalnu konfiguraciu
  X <- matrix(NA_real_, n, r)
  X[idx1, ] <- X1
  
  # interpolacia po blokoch
  for (b in 1:length(blocks)) {
    
    idx_j <- blocks[[b]]
    Xj_data <- X0[idx_j, ]
    m_j <- nrow(Xj_data)
    
    # A_{j1} - vzdialenosti medzi blokom a X1
    A_j1 <- cross_dist(Xj_data, X1_data)
    
    # Gowerova formula
    X_hat_j <- (1 / (2 * l)) * ((matrix(1, m_j, 1) %*% t(q1)) - (A_j1^2)) %*% X1 %*% S1_inv
    
    X[idx_j, ] <- X_hat_j
  }
  X
}


##### 2. PRIKLADY #####
### 2.0 PRIPRAVA ###
# install.packages(bigmds)
library(bigmds)
r = 2
n_cores = 1

set.seed(42)
setwd("C:/Users/furde/Downloads/Diplomovka")

### 2.1 MNIST FASHION ###
# install.packages(c("OpenML", "farff"))
library(OpenML)

d <- getOMLDataSet(data.id = 40996)
df <- d$data

y <- df$class
X <- df[, setdiff(names(df), "class")]

keep <- c("0","1","2","3","4")

df2 <- df[as.character(df$class) %in% keep, ]

y <- droplevels(df2$class)
X <- as.matrix(df2[, setdiff(names(df2), "class")])
dim(X)

idx <- sample(seq_len(nrow(X)), 10000)
X <- X[idx, ]
X <- as.matrix(X)
y <- y[idx]
dim(X)
dim(y)
class(X)
class(y)

colors_map <- c("0" = "blue", "1" = "red", "2" = "green", "3" = "orange", "4" = "purple")
colors <- colors_map[as.character(y)]

## VYPOCET ##
D <- dist(X)
l = 430
par(mfrow = c(1, 2))

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove
plot(fit_nove, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "Implementované IMDS")

s2 <- proc.time()
fit_stare <- interpolation_mds(X, l, r, n_cores)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_stare <- dist(fit_stare$points)
stress_stare <- mean((D - D_stare)^2)
t2
stress_stare
plot(fit_stare$points, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "IMDS")
par(mfrow = c(1, 1))

### 2.2 MNIST ###
M <- read.table("MNIST.txt")

# farby
library(RColorBrewer)
col.pal <- c(brewer.pal(9, "Set1"),"black")

set.seed(369)
ind <- sample.int(nrow(M), 10000, replace=TRUE)
lab <- M[ind,1]
X <- as.matrix(M[ind,-1])

colors <- col.pal[lab+1]

## VYPOCET ##
D <- dist(X)
l = 430
par(mfrow = c(1, 2))

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove
plot(fit_nove, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "Implementované IMDS")

s2 <- proc.time()
fit_stare <- interpolation_mds(X, l, r, n_cores)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_stare <- dist(fit_stare$points)
stress_stare <- mean((D - D_stare)^2)
t2
stress_stare
plot(fit_stare$points, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "IMDS")
par(mfrow = c(1, 1))

### 2.3. MEDICAL CONDITIONS ###
data <- read.csv("C:/Users/furde/Downloads/Diplomovka/Data/health_con.csv")
data[, 2:ncol(data)] <- scale(data[, 2:ncol(data)])

df <- data[sample(nrow(data), 10000),]
X <- df[,2:ncol(df)]
X <- as.matrix(X)

colors_map <- c("Arthritis" = "blue", "Asthma" = "yellow", "Cancer" = "red",
                "Diabetes" = "orange", "Obesity" = "purple",
                "Hypertension" = "black", "Healthy" = "green")
colors <- colors_map[df$Medical_Condition]

## VYPOCET ##
D <- dist(X)
l = 320
par(mfrow = c(1, 2))

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove
plot(fit_nove, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "Implementované IMDS")

s2 <- proc.time()
fit_stare <- interpolation_mds(X, l, r, n_cores)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_stare <- dist(fit_stare$points)
stress_stare <- mean((D - D_stare)^2)
t2
stress_stare
plot(fit_stare$points, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "IMDS")
par(mfrow = c(1, 1))


##### 3. POROVNANIE S CMDS (FASHION-MNIST) #####
r = 2

## Fashion-MNIST ##

# install.packages(c("OpenML", "farff"))
library(OpenML)

d <- getOMLDataSet(data.id = 40996)
df <- d$data

y <- df$class
X <- df[, setdiff(names(df), "class")]

keep <- c("0","1","2","3","4")

df2 <- df[as.character(df$class) %in% keep, ]

y <- droplevels(df2$class)
X <- as.matrix(df2[, setdiff(names(df2), "class")])
dim(X)

idx <- sample(seq_len(nrow(X)), 3000)
X <- X[idx, ]
X <- as.matrix(X)
y <- y[idx]
dim(X)
dim(y)
class(X)
class(y)

colors_map <- c("0" = "blue", "1" = "red", "2" = "green", "3" = "orange", "4" = "purple")
colors <- colors_map[as.character(y)]

## VYPOCET ##
D <- dist(X)
l = 500
par(mfrow = c(1, 2))

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove
plot(fit_nove, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "Implementované IMDS")

s2 <- proc.time()
fit_cmds <- cmdscale(D, k=2, eig=TRUE)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_cmds <- dist(fit_cmds$points)
stress_cmds <- mean((D - D_cmds)^2)
t2
stress_cmds
plot(fit_cmds$points, pch = 19, col = colors,
     xlab = "x", ylab = "y", main = "CMDS")
par(mfrow = c(1, 1))

## este zvysne dve na ukazku ##
# MNIST

M <- read.table("MNIST.txt")

set.seed(369)
ind <- sample.int(nrow(M), 3000, replace=TRUE)
lab <- M[ind,1]
X <- as.matrix(M[ind,-1])

## VYPOCET ##
D <- dist(X)
l = 500

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove

s2 <- proc.time()
fit_cmds <- cmdscale(D, k=2, eig=TRUE)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_cmds <- dist(fit_cmds$points)
stress_cmds <- mean((D - D_cmds)^2)
t2
stress_cmds

# Zdravotne stavy

data <- read.csv("C:/Users/furde/Downloads/Diplomovka/Data/health_con.csv")
data[, 2:ncol(data)] <- scale(data[, 2:ncol(data)])

df <- data[sample(nrow(data), 3000),]
X <- df[,2:ncol(df)]
X <- as.matrix(X)

## VYPOCET ##
D <- dist(X)
l = 400

s1 <- proc.time()
fit_nove <- imds(X, l, r)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_nove <- dist(fit_nove)
stress_nove <- mean((D - D_nove)^2)
t1
stress_nove

s2 <- proc.time()
fit_cmds <- cmdscale(D, k=2, eig=TRUE)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_cmds <- dist(fit_cmds$points)
stress_cmds <- mean((D - D_cmds)^2)
t2
stress_cmds