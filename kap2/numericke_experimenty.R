##### 1. DATA #####

set.seed(42)
data <- read.csv("C:/Users/furde/Downloads/Diplomovka/Data/health_con.csv")
data[, 2:ncol(data)] <- scale(data[, 2:ncol(data)])

set.seed(42)
df <- data[sample(nrow(data), 5000), ]
hod <- df[,2:ncol(df)]
hod <- as.matrix(hod)
D <- dist(hod)

# farby podla ochorenia
colors_map <- c("Arthritis" = "blue", "Asthma" = "yellow", "Cancer" = "red",
            "Diabetes" = "orange", "Obesity" = "purple",
            "Hypertension" = "black", "Healthy" = "green")
colors <- colors_map[df$Medical_Condition]

##### 2. CMDS #####
s1 <- proc.time()
dist_CMDS <- cmdscale(D, k=2, eig=TRUE)
e1 <- proc.time()
t1 <- (e1 - s1)["elapsed"]
D_cmds <- dist(dist_CMDS$points)
stress_cmds <- mean((D - D_cmds)^2)
# vykreslenie
plot(dist_CMDS$points, xlab = "X", ylab = 'Y', pch = 19, main = "Vykreslené ochorenia - CMDS", col = colors)
legend("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
                             "Obesity", "Hypertension", "Healthy"),
       col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
       pch = 19, cex = 0.6)

##### 3. Pokrocile MDS metody #####

# install.packages('bigmds')
library(bigmds)

##### 3. POROVNANIE S CMDS #####
col_big <- c("blue", "green", "red")
l <- 300
r <- 2
c <- 5*r
n_cores <- 1

# divide_conquer_mds(x, l, c_points, r, n_cores)
s2 <- proc.time()
dist_DCMDS <- divide_conquer_mds(hod, l, c, r, n_cores)
e2 <- proc.time()
t2 <- (e2 - s2)["elapsed"]
D_dcmds <- dist(dist_DCMDS$points)
stress_dcmds <- mean((D - D_dcmds)^2)
# interpolation_mds(x, l, r, n_cores)
s3 <- proc.time()
dist_IMDS <- interpolation_mds(hod, l, r, n_cores)
e3 <- proc.time()
t3 <- (e3 - s3)["elapsed"]
D_imds <- dist(dist_IMDS$points)
stress_imds <- mean((D - D_imds)^2)
# fast_mds(x, l, s_points, r, n_cores)
s4 <- proc.time()
dist_FMDS <- fast_mds(hod, l, c, r, n_cores)
e4 <- proc.time()
t4 <- (e4 - s4)["elapsed"]
D_fmds <- dist(dist_FMDS$points)
stress_fmds <- mean((D - D_fmds)^2)

par(mfrow = c(2, 2))
body <- list(dist_CMDS$points, dist_DCMDS$points,
             dist_IMDS$points,dist_FMDS$points)
mena <- c("CMDS", "Divide and Conquer MDS", "Interpolation MDS", "Fast MDS")

for (k in 1:4) {
  plot(body[[k]], xlab = "X", ylab = 'Y', pch = 19, cex = 0.6,
       main = mena[k], col = colors)
}
par(mfrow = c(1,1))

plot(c(t1,t2,t3,t4), c(stress_cmds, stress_dcmds, stress_imds, stress_fmds),
     xlab = "Čas", ylab = "Stres",
     main = paste("Graf pre: n = ", dim(hod)[1]),
     pch = 19, col = c("black", col_big),
     xlim = c(0,t1))
legend("topright",
        legend = c("CMDS", "Divide and conquer MDS", "Interpolation MDS", "Fast MDS"),
        col = c("black", col_big), pch = 19)
# tu vidime ze CMDS sa neoplati - dalej s nim nebudeme pracovat

### 3.1 3D ###
# r <- 3
# c <- 5*r
# 
# dist_CMDS_3D <- cmdscale(D, k=3, eig=TRUE)
# dist_DMDS_3D <- divide_conquer_mds(hod, l, c, r, n_cores)
# dist_IMDS_3D <- interpolation_mds(hod, l, r, n_cores)
# dist_FMDS_3D <- fast_mds(hod, l, c, r, n_cores)
# 
# library(rgl)
# plot3d(dist_CMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia CMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)
# plot3d(dist_DMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia DCMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)
# plot3d(dist_IMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia IMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)
# plot3d(dist_FMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia FMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)

##### 4. POROVNANIE PRE ROZNE l #####
r <- 2
c <- 5*r
l <- seq(20, 500, by = 20)
matica <- matrix(0, nrow = length(l), ncol = 6)

par(mfrow = c(1,2))
for (j in 1:length(l)) {
  l_now = l[j]
  
  # divide_conquer_mds(x, l, c_points, r, n_cores)
  s2 <- proc.time()
  dist_DCMDS <- divide_conquer_mds(hod, l_now, c, r, n_cores)
  e2 <- proc.time()
  t2 <- (e2 - s2)["elapsed"]
  D_dcmds <- dist(dist_DCMDS$points)
  stress_dcmds <- mean((D - D_dcmds)^2)
  # interpolation_mds(x, l, r, n_cores)
  s3 <- proc.time()
  dist_IMDS <- interpolation_mds(hod, l_now, r, n_cores)
  e3 <- proc.time()
  t3 <- (e3 - s3)["elapsed"]
  D_imds <- dist(dist_IMDS$points)
  stress_imds <- mean((D - D_imds)^2)
  # fast_mds(x, l, s_points, r, n_cores)
  s4 <- proc.time()
  dist_FMDS <- fast_mds(hod, l_now, c, r, n_cores)
  e4 <- proc.time()
  t4 <- (e4 - s4)["elapsed"]
  D_fmds <- dist(dist_FMDS$points)
  stress_fmds <- mean((D - D_fmds)^2)
  
  # plot(c(t2, t3, t4), c(stress_dcmds, stress_imds, stress_fmds),
  #      xlab = "Čas", ylab = "Stres",
  #      xlim = c(0, max(t2, t3, t4)),
  #      ylim = c(0, max(stress_dcmds, stress_imds, stress_fmds)),
  #      main = paste("Graf pre: n = ", dim(hod)[1],", l = ", l_now),
  #      pch = 19, col = col_big)
  # legend("bottomright", cex = 0.7,
  #        legend = c("Divide and Conquer MDS", "Interpolation MDS", "Fast MDS"),
  #        col = col_big, pch = 19,
  #        title = "Metódy z knižnice BIGMDS")
  
  matica[j,] <- c(stress_dcmds, stress_imds, stress_fmds, t2, t3, t4)
}
par(mfrow = c(1,2))

plot(l, matica[,1], main = "Divide and Conquer MDS - stres", ylab = "Stres",
     col = "blue", pch = 19)
plot(l, matica[,4], main = "Divide and Conquer MDS - čas", ylab = "Čas",
     col = "blue", pch = 19)

plot(l, matica[,2], main = "Interpolation MDS - stres", ylab = "Stres",
     col = "green", pch = 19)
plot(l, matica[,5], main = "Interpolation MDS - čas", ylab = "Čas",
     col = "green", pch = 19)

plot(l, matica[,3], main = "Fast MDS - stres", ylab = "Stres",
     col = "red", pch = 19)
plot(l, matica[,6], main = "Fast MDS - čas", ylab = "Čas",
     col = "red", pch = 19)

# najlepsie pre l su pre vsetky metody rozne, priblizne
l_dc <- 60
l_i  <- 420
l_f  <- 220

##### 5. POROVNAVANIE METOD PRE ROZNE n #####
library(bigmds)
r <- 2
c <- 5*r
n_cores <- 1
n <- seq(2000, 18000, by = 2000)
dlzka <- length(n)
matica <- matrix(0, nrow = dlzka, ncol = 6)

l <- seq(20, 500, by = 20)
mat_l <- matrix(0, nrow = dlzka, ncol = 4)

col_big <- c("blue", "green", "red")
colors_map <- c("Arthritis" = "blue", "Asthma" = "yellow", "Cancer" = "red",
            "Diabetes" = "orange", "Obesity" = "purple",
            "Hypertension" = "black", "Healthy" = "green")

for (i in 1:dlzka) {
  set.seed(42)
  df <- data[sample(nrow(data), n[i]), ]
  hod <- df[,2:ncol(df)]
  hod <- as.matrix(hod)
  D <- dist(hod)
  
  colors <- colors_map[df$Medical_Condition]
  
  # vypocet optimalneho l
  mat_pomocna <- matrix(0, nrow = length(l), ncol = 6)
  for (j in 1:length(l)) {
    l_now = l[j]
    
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(hod, l_now, c, r, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    D_dcmds <- dist(dist_DCMDS$points)
    stress_dcmds <- mean((D - D_dcmds)^2)
    
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(hod, l_now, r, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    D_imds <- dist(dist_IMDS$points)
    stress_imds <- mean((D - D_imds)^2)
    
    s4 <- proc.time()
    dist_FMDS <- fast_mds(hod, l_now, c, r, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    D_fmds <- dist(dist_FMDS$points)
    stress_fmds <- mean((D - D_fmds)^2)
    
    mat_pomocna[j,] <- c(stress_dcmds, stress_imds, stress_fmds, t2, t3, t4)
  
  }
  
  min_stres <- function(matica, stres, sek, limit = 5) {
    podm <- matica[, sek] < limit
    min_hod <- min(matica[podm, stres])
    which(matica[, stres] == min_hod & matica[, sek] < limit)
  }
  
  r1 <- min_stres(mat_pomocna, 1, 4)
  r2 <- min_stres(mat_pomocna, 2, 5)
  r3 <- min_stres(mat_pomocna, 3, 6)
  
  mat_l[i,] <- c(n[i], l[r1], l[r2], l[r3])
  
  # stres a cas pre optimalne l
  stress_dcmds <- mat_pomocna[r1, 1]
  stress_imds <- mat_pomocna[r2, 2]
  stress_fmds <- mat_pomocna[r3, 3]
  t2 <- mat_pomocna[r1, 4]
  t3 <- mat_pomocna[r2, 5]
  t4 <- mat_pomocna[r3, 6]
  
  matica[i,1:3] <- c(stress_dcmds,  stress_imds,  stress_fmds)
  matica[i,4:6] <- c(t2, t3, t4)

  # plot(c(t2, t3, t4), c(stress_dcmds, stress_imds, stress_fmds),
  #      xlab = "Čas", ylab = "Stres",
  #      main = paste("Graf pre: n = ", dim(hod)[1]),
  #      pch = 19, col = col_big)
  # legend("topright", cex = 0.7, pch = 19,
  #        legend = c("Divide and Conquer MDS", "Interpolation MDS","Fast MDS"),
  #        col = col_big,
  #        title = "Metódy z knižnice BIGMDS")
}

l_dc <- mat_l[dlzka, 2]
l_i <- mat_l[dlzka, 3]
l_f <- mat_l[dlzka, 4]
# divide_conquer_mds(x, l, c_points, r, n_cores)
dist_DCMDS <- divide_conquer_mds(hod, l_dc, c, r, n_cores)
# interpolation_mds(x, l, r, n_cores)
dist_IMDS <- interpolation_mds(hod, l_i, r, n_cores)
# fast_mds(x, l, s_points, r, n_cores)
dist_FMDS <- fast_mds(hod, l_f, c, r, n_cores)

# vykreslenie
par(mfrow = c(1, 3))
body <- list(dist_DCMDS$points, dist_IMDS$points, dist_FMDS$points)
mena <- c("Divide and Conquer MDS", "Interpolation MDS", "Fast MDS")
for (k in 1:3) {
  plot(body[[k]], pch = 19, xlab= "X", ylab = "Y",
       main = mena[k], col = colors)
}
par(mfrow = c(1,1))

### 5.1 3D ###
# r <- 3
# c <- 5*r
# 
# dist_DMDS_3D <- divide_conquer_mds(hod, l_dc, c, r, n_cores)
# dist_IMDS_3D <- interpolation_mds(hod, l_i, r, n_cores)
# dist_FMDS_3D <- fast_mds(hod, l_f, c, r, n_cores)
# 
# plot3d(dist_DMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia DCMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)
# plot3d(dist_IMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia IMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                                "Obesity", "Hypertension", "Healthy"),
#          col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#          pch = 19, cex = 0.6)
# plot3d(dist_FMDS_3D$points, col = colors, type = "s", size = 1,
#        xlab = "X", ylab = "Y", zlab = "Z", main = "Rotacia FMDS 3D")
# legend3d("topleft", legend = c("Arthritis", "Asthma", "Cancer", "Diabetes",
#                              "Obesity", "Hypertension", "Healthy"),
#        col = c("blue", "yellow", "red", "orange", "purple", "black", "green"),
#        pch = 19, cex = 0.6)

### 5.2 Porovnanie stresu ###
matica

plot(n, matica[,1], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[, 1:3])*1.1),
     xlab = "Veľkosť dát", ylab = "Stres", main = "Porovnanie metód - stres")
lines(n, matica[,2], type = "b", pch = 17,  lty = 3, col = "green")
lines(n, matica[,3], type = "b", pch = 18,  lty = 4, col = "red")
legend("bottomright", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

### 5.3 Porovnanie času ###
plot(n, matica[,4], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[, 4:6])*1.1),
     xlab = "Veľkosť dát", ylab = "Čas", main = "Porovnanie metód - čas")
lines(n, matica[,5], type = "b", pch = 17,  lty = 3, col = "green")
lines(n, matica[,6], type = "b", pch = 18,  lty = 4, col = "red")
legend("topleft", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

##### 6. MNIST #####

# library(dslabs)
# mnist <- read_mnist()
# X <- mnist$train$images
# lab <- mnist$train$labels
# M <- data.frame(lab = lab, X)
# write.table(M, file="MNIST.txt")

set.seed(42)
setwd("C:/Users/furde/Downloads/Diplomovka")
M <- read.table("MNIST.txt")

# farby
library(RColorBrewer)
col.pal <- c(brewer.pal(9, "Set1"),"black")

# bigmds metody
library(bigmds)
col_big <- c("blue", "green", "red")
r <- 2
c <- 5*r
n_cores <- 1

### 6.1 vypocet l's ###
n <- seq(3000, 12000, by = 3000)
dlzka <- length(n)
matica <- matrix(0, nrow = dlzka, ncol = 6)

l <- seq(20, 500, by = 20)
mat_l <- matrix(0, nrow = dlzka, ncol = 4)

for (i in 1:dlzka) {
  set.seed(42)
  ind <- sample.int(nrow(M), n[i], replace=FALSE) 
  lab <- M[ind,1]
  X <- as.matrix(M[ind,-1])
  D <- dist(X)
  
  cols <- col.pal[lab+1]
  
  # vypocet optimalneho l
  mat_pomocna <- matrix(0, nrow = length(l), ncol = 6)
  for (j in 1:length(l)) {
    l_now = l[j]
    
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(X, l_now, c, r, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    D_dcmds <- dist(dist_DCMDS$points)
    stress_dcmds <- mean((D - D_dcmds)^2)
    
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(X, l_now, r, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    D_imds <- dist(dist_IMDS$points)
    stress_imds <- mean((D - D_imds)^2)
    
    s4 <- proc.time()
    dist_FMDS <- fast_mds(X, l_now, c, r, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    D_fmds <- dist(dist_FMDS$points)
    stress_fmds <- mean((D - D_fmds)^2)
    
    mat_pomocna[j,] <- c(stress_dcmds, stress_imds, stress_fmds, t2, t3, t4)
    
  }
  
  min_stres <- function(matica, stres, sek, limit = 10) { # tu je (zatial) zmena z 5 na 10
    podm <- matica[, sek] < limit
    min_hod <- min(matica[podm, stres])
    which(matica[, stres] == min_hod & matica[, sek] < limit)
  }
  
  r1 <- min_stres(mat_pomocna, 1, 4)
  r2 <- min_stres(mat_pomocna, 2, 5)
  r3 <- min_stres(mat_pomocna, 3, 6)
  
  mat_l[i,] <- c(n[i], l[r1], l[r2], l[r3])
  
  # stres a cas pre optimalne l
  stress_dcmds <- mat_pomocna[r1, 1]
  stress_imds <- mat_pomocna[r2, 2]
  stress_fmds <- mat_pomocna[r3, 3]
  t2 <- mat_pomocna[r1, 4]
  t3 <- mat_pomocna[r2, 5]
  t4 <- mat_pomocna[r3, 6]
  
  matica[i,1:3] <- c(stress_dcmds,  stress_imds,  stress_fmds)
  matica[i,4:6] <- c(t2, t3, t4)
  
  # plot(c(t2, t3, t4), c(stress_dcmds, stress_imds, stress_fmds),
  #      xlab = "Čas", ylab = "Stres",
  #      main = paste("Graf pre: n = ", dim(X)[1]),
  #      pch = 19, col = col_big)
  # legend("topright", cex = 0.7, pch = 19,
  #        legend = c("Divide and Conquer MDS", "Interpolation MDS","Fast MDS"),
  #        col = col_big,
  #        title = "Metódy z knižnice BIGMDS")
}

l_dc <- mat_l[dlzka, 2]
l_i <- mat_l[dlzka, 3]
l_f <- mat_l[dlzka, 4]
# divide_conquer_mds(x, l, c_points, r, n_cores)
dist_DCMDS <- divide_conquer_mds(X, l_dc, c, r, n_cores)
# interpolation_mds(x, l, r, n_cores)
dist_IMDS <- interpolation_mds(X, l_i, r, n_cores)
# fast_mds(x, l, s_points, r, n_cores)
dist_FMDS <- fast_mds(X, l_f, c, r, n_cores)

## vykreslenie ##
par(mfrow = c(1, 3))
body <- list(dist_DCMDS$points, dist_IMDS$points, dist_FMDS$points)
mena <- c("Divide and Conquer MDS", "Interpolation MDS", "Fast MDS")
for (k in 1:3) {
  plot(body[[k]], pch = 19, xlab= "X", ylab = "Y",
       main = mena[k], col = cols)
}
par(mfrow = c(1,1))

### POROVNAVANIE ###

### 6.2 Porovnanie stresu ###
matica

par(mfrow = c(1,2))
plot(n, matica[,1], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[, 1:3])*1.1),
     xlab = "Veľkosť dát", ylab = "Stres", main = "Porovnanie metód - stres")
lines(n, matica[,2], type = "b", pch = 17,  lty = 3, col = "green")
lines(n, matica[,3], type = "b", pch = 18,  lty = 4, col = "red")
legend("bottomright", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

### 6.3 Porovnanie času ###
plot(n, matica[,4], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[, 4:6])*1.1),
     xlab = "Veľkosť dát", ylab = "Čas", main = "Porovnanie metód - čas")
lines(n, matica[,5], type = "b", pch = 17,  lty = 3, col = "green")
lines(n, matica[,6], type = "b", pch = 18,  lty = 4, col = "red")
legend("topleft", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))
par(mfrow = c(1,1))

##### 7. VYPOCET PRE ROZNE DIMENZIE (r) #####

set.seed(42)
setwd("C:/Users/furde/Downloads/Diplomovka")
M <- read.table("MNIST.txt")

ind <- sample.int(nrow(M), 3000, replace=FALSE)
lab <- M[ind,1]
X <- as.matrix(M[ind,-1])
D <- dist(X)

# install.packages('bigmds')
library(bigmds)

n_cores <- 1

r <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50)

min_stres <- function(matica, stres, sek, limit = 10) {
  podm <- matica[, sek] < limit
  min_hod <- min(matica[podm, stres])
  which(matica[, stres] == min_hod & matica[, sek] < limit)
}

matica <- matrix(0, nrow = length(r), ncol = 12)

for (i in 1:length(r)) {
  r_now <- r[i]
  c_now <- 5*r_now

  l <- sapply(1:21, function(i) 2*c_now + 20*(i-1))

  mat_pomocna <- matrix(0, nrow = length(l), ncol = 6)

  s1 <- proc.time()
  dist_CMDS <- cmdscale(D, k= r_now, eig=TRUE)
  e1 <- proc.time()
  t1 <- (e1 - s1)["elapsed"]
  D_cmds <- dist(dist_CMDS$points)
  stress_cmds <- mean((D - D_cmds)^2)
  
  for (j in 1:length(l)) {
    l_now <- l[j]

    # divide_conquer_mds(x, l, c_points, r, n_cores)
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(X, l_now, c_now, r_now, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    D_dcmds <- dist(dist_DCMDS$points)
    stress_dcmds <- mean((D - D_dcmds)^2)
    # interpolation_mds(x, l, r, n_cores)
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(X, l_now, r_now, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    D_imds <- dist(dist_IMDS$points)
    stress_imds <- mean((D - D_imds)^2)
    # fast_mds(x, l, s_points, r, n_cores)
    s4 <- proc.time()
    dist_FMDS <- fast_mds(X, l_now, c_now, r_now, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    D_fmds <- dist(dist_FMDS$points)
    stress_fmds <- mean((D - D_fmds)^2)
  
    mat_pomocna[j,] <- c(stress_dcmds, stress_imds, stress_fmds, t2, t3, t4)
  }
  
  r1 <- min_stres(mat_pomocna, 1, 4)
  r2 <- min_stres(mat_pomocna, 2, 5)
  r3 <- min_stres(mat_pomocna, 3, 6)
  
  l_dc <- l[r1]
  l_i  <- l[r2]
  l_f  <- l[r3]
  
  stress_dcmds <- mat_pomocna[r1, 1]
  stress_imds <- mat_pomocna[r2, 2]
  stress_fmds <- mat_pomocna[r3, 3]
  
  n_time <- 10
  time_mat <- matrix(0, nrow = n_time, ncol = 3)
  for (k in 1:n_time) {
    # divide_conquer_mds(x, l, c_points, r, n_cores)
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(X, l_dc, c_now, r_now, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    # interpolation_mds(x, l, r, n_cores)
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(X, l_i, r_now, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    # fast_mds(x, l, s_points, r, n_cores)
    s4 <- proc.time()
    dist_FMDS <- fast_mds(X, l_f, c_now, r_now, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    
    time_mat[k,] <- c(t2, t3, t4)
  }
  
  t2 <- mean(time_mat[,1])
  t3 <- mean(time_mat[,2])
  t4 <- mean(time_mat[,3])
  
  print(t2)

  matica[i,] <- c(r_now, t1, stress_cmds, l_dc, t2, stress_dcmds,
                 l_i, t3, stress_imds, l_f, t4, stress_fmds)
  
  print(r_now)
}
matica

# saveRDS(matica, "matica_pre_dimenzie.rds")
matica <- readRDS("matica_pre_dimenzie.rds")
matica

### 7.1 Porovnanie l-ka ###
plot(r, matica[,4], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[,4], matica[,7], matica[,10])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Veľkosť podvzorky", main = "Porovnanie metód - veľkosť podvzorky")
lines(r, matica[,7], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,10], type = "b", pch = 18,  lty = 4, col = "red")
legend("topleft", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

### 7.2 Porovnanie času ###
# par(mfrow = c(1, 2))
plot(r, matica[,5], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[,5], matica[,8], matica[,11])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Čas", main = "Porovnanie metód - čas")
lines(r, matica[,8], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,11], type = "b", pch = 18,  lty = 4, col = "red")
legend("bottomright", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))
## 7.2.1 Porovnanie casu CMDS ##
# plot(r, matica[,2], type = "b", pch = 15, col = "black",
#      lty = 1, ylim = c(0, max(matica[,2])*1.1),
#      xlab = "Veľkosť dimenzie", ylab = "Čas", main = "Výpočtový čas CMDS")
# legend("bottomright", legend = c("C_MDS"),
#        col = c("black"), pch = 15, lty = 1)
# par(mfrow = c(1, 1))

### 7.3 Porovnanie stresu ###
plot(r, matica[,3], type = "b", pch = 15, col = "black",
     lty = 1, ylim = c(0, max(matica[,3], matica[,6], matica[,9], matica[,12])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Stres", main = "Porovnanie metód - stres")
lines(r, matica[,6], type = "b", pch = 16,  lty = 2, col = "blue")
lines(r, matica[,9], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,12], type = "b", pch = 18,  lty = 4, col = "red")
legend("topright", legend = c("C_MDS", "DC_MDS", "I_MDS", "F_MDS"),
       col = c("black", "blue", "green", "red"), pch = c(15, 16, 17, 18), lty = c(1, 2, 3, 4))



##### OPTIMALIZACNE KRITERIUM CMDS #####

set.seed(42)
setwd("C:/Users/furde/Downloads/Diplomovka")
M <- read.table("MNIST.txt")

ind <- sample.int(nrow(M), 3000, replace=FALSE)
lab <- M[ind,1]
X <- as.matrix(M[ind,-1])
D <- dist(X)

# install.packages('bigmds')
library(bigmds)

n_cores <- 1

r <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50)

min_obj <- function(matica, obj, sek, limit = 10) {
  podm <- matica[, sek] < limit
  min_hod <- min(matica[podm, obj])
  which(matica[, obj] == min_hod & matica[, sek] < limit)
}

D_mat <- as.matrix(D)
n <- nrow(D_mat)
J <- diag(n) - matrix(1 / n, n, n)
B <- -0.5 * J %*% (D_mat^2) %*% J

objective <- function(B, X) {
  X <- as.matrix(X)
  X <- scale(X, center = TRUE, scale = FALSE)
  X <- as.matrix(X)
  B_hat <- X %*% t(X)
  sum((B - B_hat)^2)
}

matica <- matrix(0, nrow = length(r), ncol = 12)

for (i in 1:length(r)) {
  r_now <- r[i]
  c_now <- 5*r_now
  
  l <- sapply(1:21, function(i) 2*c_now + 20*(i-1))
  
  mat_pomocna <- matrix(0, nrow = length(l), ncol = 6)
  
  s1 <- proc.time()
  dist_CMDS <- cmdscale(D, k= r_now, eig=TRUE)
  e1 <- proc.time()
  t1 <- (e1 - s1)["elapsed"]
  obj_cmds <- objective(B, dist_CMDS$points)
  
  for (j in 1:length(l)) {
    l_now <- l[j]
    
    # divide_conquer_mds(x, l, c_points, r, n_cores)
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(X, l_now, c_now, r_now, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    obj_dcmds <- objective(B, dist_DCMDS$points)
    # interpolation_mds(x, l, r, n_cores)
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(X, l_now, r_now, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    obj_imds <- objective(B, dist_IMDS$points)
    # fast_mds(x, l, s_points, r, n_cores)
    s4 <- proc.time()
    dist_FMDS <- fast_mds(X, l_now, c_now, r_now, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    obj_fmds <- objective(B, dist_FMDS$points)
    
    mat_pomocna[j,] <- c(obj_dcmds, obj_imds, obj_fmds, t2, t3, t4)
  }
  
  r1 <- min_obj(mat_pomocna, 1, 4)
  r2 <- min_obj(mat_pomocna, 2, 5)
  r3 <- min_obj(mat_pomocna, 3, 6)
  
  l_dc <- l[r1]
  l_i  <- l[r2]
  l_f  <- l[r3]
  
  obj_dcmds <- mat_pomocna[r1, 1]
  obj_imds <- mat_pomocna[r2, 2]
  obj_fmds <- mat_pomocna[r3, 3]
  
  n_time <- 10
  time_mat <- matrix(0, nrow = n_time, ncol = 3)
  for (k in 1:n_time) {
    # divide_conquer_mds(x, l, c_points, r, n_cores)
    s2 <- proc.time()
    dist_DCMDS <- divide_conquer_mds(X, l_dc, c_now, r_now, n_cores)
    e2 <- proc.time()
    t2 <- (e2 - s2)["elapsed"]
    # interpolation_mds(x, l, r, n_cores)
    s3 <- proc.time()
    dist_IMDS <- interpolation_mds(X, l_i, r_now, n_cores)
    e3 <- proc.time()
    t3 <- (e3 - s3)["elapsed"]
    # fast_mds(x, l, s_points, r, n_cores)
    s4 <- proc.time()
    dist_FMDS <- fast_mds(X, l_f, c_now, r_now, n_cores)
    e4 <- proc.time()
    t4 <- (e4 - s4)["elapsed"]
    
    time_mat[k,] <- c(t2, t3, t4)
  }
  
  t2 <- mean(time_mat[,1])
  t3 <- mean(time_mat[,2])
  t4 <- mean(time_mat[,3])
  
  matica[i,] <- c(r_now, t1, obj_cmds, l_dc, t2, obj_dcmds,
                  l_i, t3, obj_imds, l_f, t4, obj_fmds)
  print(r_now)
}
matica

### 7.1 Porovnanie l-ka ###
plot(r, matica[,4], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[,4], matica[,7], matica[,10])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Veľkosť podvzorky", main = "Porovnanie metód - veľkosť podvzorky")
lines(r, matica[,7], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,10], type = "b", pch = 18,  lty = 4, col = "red")
legend("topleft", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

### 7.2 Porovnanie času ###
plot(r, matica[,5], type = "b", pch = 16, col = "blue",
     lty = 2, ylim = c(0, max(matica[,5], matica[,8], matica[,11])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Čas", main = "Porovnanie metód - čas")
lines(r, matica[,8], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,11], type = "b", pch = 18,  lty = 4, col = "red")
legend("bottomright", legend = c("DC_MDS", "I_MDS", "F_MDS"),
       col = c("blue", "green", "red"), pch = c(16, 17, 18), lty = c(2, 3, 4))

### 7.3 Porovnanie optimalizacie ###
plot(r, matica[,3], type = "b", pch = 15, col = "black",
     lty = 1, ylim = c(0, max(matica[,3], matica[,6], matica[,9], matica[,12])*1.1),
     xlab = "Veľkosť dimenzie", ylab = "Hodnota cieľovej funkcie", main = "Porovnanie metód - hodnota cieľovej funkcie")
lines(r, matica[,6], type = "b", pch = 16,  lty = 2, col = "blue")
lines(r, matica[,9], type = "b", pch = 17,  lty = 3, col = "green")
lines(r, matica[,12], type = "b", pch = 18,  lty = 4, col = "red")
legend("topright", legend = c("C_MDS", "DC_MDS", "I_MDS", "F_MDS"),
       col = c("black", "blue", "green", "red"), pch = c(15, 16, 17, 18), lty = c(1, 2, 3, 4))
# -> CMDS je skutocne najlepsie v ramci svojho kriteria