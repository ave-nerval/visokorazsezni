### Eva Lavrencic
### Januar 2022
## Permutation methods for differences between groups


library(sva)
library(permuco)

##########################################################################
# Funkcije za metode

funkcija_perm_ost <- function(data){
  pval.Gen.H0 <- NULL
  pval.Gen.H0_brezZ <- NULL
  
  if(z_value == "cat") data$z = as.factor(data$z)
  
  modelFull <- lm(data$y ~ data$x + data$z)
  modelRed <- lm(data$y ~ data$z)
  
  modelFull_brezZ <- lm(data$y ~ data$x)
  modelRed_brezZ <- lm(data$y ~ 1)
  
  # referencna F vrednost
  FstatRef <- anova(modelRed, modelFull)[2,5]
  FstatRef_brezZ <- anova(modelRed_brezZ, modelFull_brezZ)[2,5]
  
  # navadna verzija
  FstatPermSave <- NULL
  for (i in 1:m){
    # permutiramo ostanke, m permutacij
    ostanki_perm <- sample(modelRed$residuals, size=n, replace = FALSE)
    # izracunamo nov y
    y_perm <- fitted(modelRed) + ostanki_perm
    # polni in reducirani model z novim y
    modelFullPerm <- lm(y_perm ~ data$x + data$z)
    modelRedPerm <- lm(y_perm ~ data$z)
    # izracun vrednosti F
    FstatPermSave[i] <- anova(modelRedPerm, modelFullPerm)[2,5]
  }
  
  # brez z
  FstatPermSave_brezZ <- NULL
  for (i in 1:m){
    # permutiramo ostanke, m permutacij
    ostanki_perm <- sample(modelRed_brezZ$residuals, size=n, replace = FALSE)
    # izracunamo nov y
    y_perm <- fitted(modelRed_brezZ) + ostanki_perm
    # polni in reducirani model z novim y
    modelFullPerm_brezZ <- lm(y_perm ~ data$x + data$z)
    modelRedPerm_brezZ <- lm(y_perm ~ data$z)
    # izracun vrednosti F
    FstatPermSave_brezZ[i] <- anova(modelRedPerm_brezZ, modelFullPerm_brezZ)[2,5]
  }
  
  # vrednost p: delez vrednosti F, ki so vecje od referencne
  pValGen.H0 <- mean(FstatPermSave >= FstatRef)
  pValGen.H0_brezZ <- mean(FstatPermSave_brezZ >= FstatRef_brezZ)
  
  cbind(pValGen.H0, pValGen.H0_brezZ)
}




funkcijaAovperm <- function(data){
  model <- aovperm(data$y ~ data$z + data$x, np=m, method="freedman_lane")
  model_brezZ <- aovperm(data$y ~ data$x, np=m, method="freedman_lane")
  
  pValPermuco <- summary(model)[2,5]
  pValPermuco_brezZ <- summary(model_brezZ)[1,5]
  
  cbind(pValPermuco, pValPermuco_brezZ)
  
}



funkcijaPermCategorical <- function(data){
  pval_fin <- NULL
  pval_fin_brezZ <- NULL
  
  
  # permutacija oznak skupin znotraj stratumov (stratume doloca z)
  
  rez<- summary(lm(data$y~data$x+data$z))
  pvrednostStart <- rez$coefficient[2,4]
  
  dataZ0 <- data[data$z==0,]
  dataZ1 <- data[data$z==1,]
  pvrednostPerm = rep(NA, m)
  
  for(i in 1:m){
    dataZ0$x <- sample(dataZ0$x, replace=F)
    dataZ1$x <- sample(dataZ1$x, replace=F)
    dataXYZ <- rbind(dataZ0, dataZ1)
    
    
    rez<- summary(lm(dataXYZ$y~dataXYZ$x+dataXYZ$z))
    pvrednostPerm[i] <- rez$coefficient[2,4]
    
  }
  
  # permutacija oznak skupin (brez upostevanja z)
  
  rez<- summary(lm(data$y~data$x))
  pvrednostStart_brezZ <- rez$coefficient[2,4]
  
  pvrednostPerm_brezZ = rep(NA, m)
  
  for(i in 1:m){
    data$x <- sample(data$x, replace=F)
    
    rez<- summary(lm(data$y~data$x))
    pvrednostPerm_brezZ[i] <- rez$coefficient[2,4]
    
  }
  
  pval_fin <- mean(pvrednostStart >= pvrednostPerm)
  pval_fin_brezZ <- mean(pvrednostStart_brezZ >= pvrednostPerm_brezZ)
  
  cbind(pval_fin, pval_fin_brezZ)
}




funkcijaSVA <- function(k, genSVA, z_value, hypothesis, z_included){
  
  pValuesSv <- matrix(NA, nrow = k, ncol = genSVA)

  for(i in 1:k){
    
    if (z_value == "cont"){
      z <- rnorm(n,1,1) #nuisance
    } else if (z_value == "cat"){
      z <- rbinom(n,1,0.3) #nuisance; npr. spol
    }
    
    x <- rbinom(n,size=1,prob=1/(1+exp(-b1*z))) #skupina
    
    if (hypothesis == "null"){
      y<-b2*z+0*x+matrix(rnorm(n*gen,sd=sqrt(sd.eps)), nrow=n, ncol=gen)
    } else if (hypothesis == "alt"){
      y<-b2*z+b3*x+matrix(rnorm(n*gen,sd=sqrt(sd.eps)), nrow=n, ncol=gen)
    }
    
    data = data.frame(y,x,z)
    data$x <- as.factor(data$x)
    rm(x,y,z)
    
    # browseVignettes("sva")
    # transpose: variables in rows and samples in columns
    data.sva <- t(data[,1:genSVA])
    
    if (z_included == TRUE){
    # full model matrix - including both the adjustment variables and the variable of interest (group x)
    mod = model.matrix(~as.factor(x) + z, data=data)
    # The null model contains only the adjustment variables
    mod0 = model.matrix(~z, data=data)
    # apply sva function
    svobj = sva(data.sva,mod,mod0)
    # Adjusting for surrogate variables using the f.pvalue function
    # can calculate the F-test p-values for differential expression with respect to cancer status,
    # without adjusting for surrogate variables, adjust them for multiple testing,
    # and calculate the number that are significant with a Q-value less than 0.05
    modSv = cbind(mod,svobj$sv)
    mod0Sv = cbind(mod0,svobj$sv)
    pValuesSv[i,] = f.pvalue(data.sva,modSv,mod0Sv)
    }
    
    if (z_included == FALSE){
    # fitting without z
    mod = model.matrix(~as.factor(x), data=data)
    mod0 = model.matrix(~1, data=data)
    svobj = sva(data.sva,mod,mod0)
    modSv = cbind(mod,svobj$sv)
    mod0Sv = cbind(mod0,svobj$sv)
    pValuesSv[i,] = f.pvalue(data.sva,modSv,mod0Sv)
  }}
  
  pValuesSv
  
}




##########################################################################
# Funkcija za simulacijo


simulacija_k <- function(hypothesis = "null", z_value = "cont"){
  pValPermOst <- matrix(NA, nrow = k, ncol = 2)
  pValPermuco <- matrix(NA, nrow = k, ncol = 2)
  pValCat <- matrix(NA, nrow = k, ncol = 2)
  
  for(j in 1:k){
    
    if (z_value == "cont"){
      z <- rnorm(n,1,1) #nuisance
    } else if (z_value == "cat"){
      z <- rbinom(n,1,0.5) #nuisance; npr., spol
    }
    
    x <- rbinom(n,size=1,prob=1/(1+exp(-b1*z)))
      
    if (hypothesis == "null"){
      y <-b2*z+0*x+rnorm(n,sd=sd.eps)
    } else if (hypothesis == "alt"){
      y <-b2*z+b3*x+rnorm(n,sd=sd.eps)
    }  
      
      
    data = data.frame(y,x,z)
    data$x <- as.factor(data$x)
    rm(x,y,z)
      
    pValPermOst[j,] <- funkcija_perm_ost(data = data)
    pValPermuco[j,] <- funkcijaAovperm(data = data)
    
    if(z_value == "cat") {
      pValCat[j,] <- funkcijaPermCategorical(data = data)
    }
    
  }
  
  if(z_value == "cont"){
    res <- data.frame(pValPermOst, pValPermuco)
    colnames(res) = c("PermOst", "PermOst_brezZ", "Permuco", "Permuco_brezZ")}
  if(z_value == "cat"){
    res <- data.frame(pValPermOst, pValPermuco, pValCat)  
    colnames(res) = c("PermOst", "PermOst_brezZ", "Permuco", "Permuco_brezZ", "Categorical", "Categorical_brezZ")
  }
  
  res

}


##########################################################################
# Simulacija

m = 50
k = 50
n = 50
gen = 100
z_value = "cont"
hypothesis = "null"
b1 = 1
b2 = 1
b3 = 1
sd.eps = 1

res_null_cont <- simulacija_k(hypothesis = "null", z_value = "cont")
res_alt_cont <- simulacija_k(hypothesis = "alt", z_value = "cont")
res_null_cat <- simulacija_k(hypothesis = "null", z_value = "cat")
res_alt_cat <- simulacija_k(hypothesis = "alt", z_value = "cat")

res_SVA_null_cont <- funkcijaSVA(k, gen, z_value = "cont", hypothesis = "null", z_included = TRUE)
res_SVA_alt_cont <- funkcijaSVA(k, genSVA = gen, z_value = "cont", hypothesis = "alt", z_included = TRUE)
res_SVA_null_cat <- funkcijaSVA(k, genSVA = gen, z_value = "cat", hypothesis = "null", z_included = TRUE)
res_SVA_alt_cat <- funkcijaSVA(k, genSVA = gen, z_value = "cat", hypothesis = "alt", z_included = TRUE)

res_SVA_null_cont_brezZ <- funkcijaSVA(k, gen, z_value = "cont", hypothesis = "null", z_included = FALSE)
res_SVA_alt_cont_brezZ <- funkcijaSVA(k, genSVA = gen, z_value = "cont", hypothesis = "alt", z_included = FALSE)
res_SVA_null_cat_brezZ <- funkcijaSVA(k, genSVA = gen, z_value = "cat", hypothesis = "null", z_included = FALSE)
res_SVA_alt_cat_brezZ <- funkcijaSVA(k, genSVA = gen, z_value = "cat", hypothesis = "alt", z_included = FALSE)





##########################################################################
# Rezultati

# REZULTATI

# Porazdelitev vrednosti p; podatki generirani pod H0
par(mfrow=c(2,4))
hist(res_null_cat$PermOst, main="Perm. ostankov, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cat$Permuco, main="Aovperm, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_null_cat, main = "SVA, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cat$Categorical, main = "Perm. x, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cont$PermOst, main = "Perm. ostankov, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cont$Permuco, main = "Aovperm, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_null_cont, main = "SVA, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")

# Porazdelitev vrednosti p; podatki generirani pod Ha
par(mfrow=c(2,4))
hist(res_alt_cat$PermOst, main="Perm. ostankov, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cat$Permuco, main="Aovperm, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_alt_cat, main = "SVA, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cat$Categorical, main = "Perm. x, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cont$PermOst, main = "Perm. ostankov, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cont$Permuco, main = "Aovperm, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_alt_cont, main = "SVA, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")

# Porazdelitev vrednosti p; podatki generirani pod H0 BREZ Z
par(mfrow=c(2,4))
hist(res_null_cat$PermOst_brezZ, main="Perm. ostankov, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cat$Permuco_brezZ, main="Aovperm, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_null_cat_brezZ, main = "SVA, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cat$Categorical_brezZ, main = "Perm. x, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cont$PermOst_brezZ, main = "Perm. ostankov, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_null_cont$Permuco_brezZ, main = "Aovperm, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_null_cont_brezZ, main = "SVA, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")

# Porazdelitev vrednosti p; podatki generirani pod Ha BREZ Z
par(mfrow=c(2,4))
hist(res_alt_cat$PermOst_brezZ, main="Perm. ostankov, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cat$Permuco_brezZ, main="Aovperm, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_alt_cat_brezZ, main = "SVA, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cat$Categorical_brezZ, main = "Perm. x, z = kat. spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cont$PermOst_brezZ, main = "Perm. ostankov, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_alt_cont$Permuco_brezZ, main = "Aovperm, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")
hist(res_SVA_alt_cont_brezZ, main = "SVA, z = zvezna spr.", ylab = "Frekvenca", xlab = "Vrednost p")









































