## ----setup, include=FALSE--------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ---- eval=FALSE-----------------------------------------------------------------------------
## 
my.fast.t.test<-function(data.Class1, data.Class2, n1, n2){

  mean.x1<-colMeans(data.Class1)

  mean.x2<-colMeans(data.Class2)

  #izracun variance
  s2p<-(colSums( (data.Class1-rep(mean.x1, each=n1))^2) +

          colSums( (data.Class2-rep(mean.x2, each=n2))^2))/(n1+n2-2)

  t.values<-(mean.x1-mean.x2)/(sqrt(s2p*(1/n1+1/n2)))
  t.values
}

n=500 #stevilo genov
p.x=0.5 #porazdelitev X
M=500 #stevilo permutacij
B=1000 #ponovitev simulacij

k <-B #stevilo simulacij v vaji

testnaPerm <- rep(NA,500)

#prazne matrike in vektorji
tt <- rep(NA, 500)
holmV <- rep(NA,B)
bonfV <- rep(NA,B)
BHV <- rep(NA,B)
permV <- rep(NA,B)
p.value <- rep(NA,500)
pmin <- rep(NA,500)
permpval <- rep(NA,500)

for(i in 1:k){
  n=25 #velikost vzorca za 1 gen
  p.x=0.5
  y <- replicate(500,rnorm(n,0,1)) #generiraj 500 genov (po 25 ponovitev)
  x <- rbinom(n,size=1,prob=0.5) #generiraj vektor x 0/1

  n1<-sum(x==0) #stevilo 0
  n2<-n-n1 #stevilo 1
  d1 <- y[x==0,] #y-i ki imajo 0
  d2 <- y[x==1,] #y-i ki imajo 1

  tt<-my.fast.t.test(d1,d2,n1,n2) #dobis vrednost t

  p.value = 2*pt(-abs(tt), df=25-2)

  bonfp <- p.adjust(p.value, method="bonf")
  holmp <- p.adjust(p.value, method="holm")
  BHp <- p.adjust(p.value, method="BH")


for (ii in 1:M){

  razmeciX <- x[sample(1:n)]

  d1 <- y[razmeciX==0,] #y-i ki imajo 0
  d2 <- y[razmeciX==1,] #y-i ki imajo 1

  testt <- t.test(d1,d2,var.equal = T)

  tt2<-my.fast.t.test(d1,d2,n1,n2) #izracunaj p vrednost za vsak gen v permutaciji !!
  #v vsaki ponovitvi shrani samo najmanjso p vredn (M) !!
  pperm <-  2*pt(-abs(tt2), df=25-2)
  pmin[ii] <- min(pperm)
  }

for (kk in 1:500){

  permpval[kk] <-  (sum(pmin<=p.value[kk])+1)/(M+1)}


#izracun V

holmV[i] <- sum(holmp<0.05)
bonfV[i] <- sum(bonfp<0.05)
BHV[i] <- sum(BHp<0.05)
permV[i] <- sum(permpval<0.05)

}

#FWER
holmFWER <- sum(holmV>0)/B
bonfFWER <- sum(bonfV>0)/B
bhFWER <- sum(BHV>0)/B
permFWER <- sum(permV>0)/B

FWER <- c(holmFWER, bonfFWER, bhFWER, permFWER)
FWER

#FDR
holmR <- holmV
holmpFDR <- mean(holmV/holmR,na.rm=TRUE)
holmFDR <- holmpFDR * mean((holmR)>0)

funkcijaFDR <- function(x){
  xR <- x
  return(mean(x/xR,na.rm=TRUE)*mean((xR>0)))
}

funkcijaFDR(holmV)
funkcijaFDR(bonfV)
funkcijaFDR(BHV)
funkcijaFDR(permV)


## ---- eval=FALSE-----------------------------------------------------------------------------
## 
my.fast.t.test<-function(data.Class1, data.Class2, n1, n2){

  mean.x1<-colMeans(data.Class1)

  mean.x2<-colMeans(data.Class2)

  #izracun variance
  s2p<-(colSums( (data.Class1-rep(mean.x1, each=n1))^2) +

          colSums( (data.Class2-rep(mean.x2, each=n2))^2))/(n1+n2-2)

  t.values<-(mean.x1-mean.x2)/(sqrt(s2p*(1/n1+1/n2)))
  t.values
}

n=500 #stevilo genov
p.x=0.5 #porazdelitev X
M=500 #stevilo permutacij
B=1000 #ponovitev simulacij
mu <- 0.5

set.seed(1)
k <-B #stevilo simulacij v vaji
#prazne matrike, vektorji
tt <- rep(NA,500)

holmV <- rep(NA,B)
bonfV <- rep(NA,B)
BHV <- rep(NA,B)
permV <- rep(NA,B)

holmS <- rep(NA,B)
bonfS <- rep(NA,B)
BHS <- rep(NA,B)
permS <- rep(NA,B)

p.value <- rep(NA,500)
pmin <- rep(NA,500)
permpval <- rep(NA,500)

for(i in 1:k){
  n=25 #velikost vzorca za 1 gen
  p.x=0.5
  y <- replicate(500,rnorm(n,0,1)) #generiraj 500 genov (po 25 ponovitev)
  x <- rbinom(n,size=1,prob=0.5) #generiraj vektor x 0/1
  y[x==1,1:10] <- matrix(rnorm(sum(x==1)*10, mu, 1), ncol=10)

  n1<-sum(x==0) #stevilo 0
  n2<-n-n1 #stevilo 1
  d1 <- y[x==0,] #y-i ki imajo 0
  d2 <- y[x==1,] #y-i ki imajo 1

  tt<-my.fast.t.test(d1,d2,n1,n2) #dobis vrednost t!

  p.value = 2*pt(-abs(tt), df=25-2) #vrednost p

  #p.adjust(p.value, bonf) !!
  #dobimo vektor l=500

  bonfp<- p.adjust(p.value, method="bonf")
  holmp <- p.adjust(p.value, method="holm")
  BHp <- p.adjust(p.value, method="BH")

  for (ii in 1:M){

  razmeciX <- x[sample(1:n)]

  d1 <- y[razmeciX==0,] #y-i ki imajo 0
  d2 <- y[razmeciX==1,] #y-i ki imajo 1

  testt <- t.test(d1,d2,var.equal = T)

  tt2<-my.fast.t.test(d1,d2,n1,n2)#izracunaj p vrednost za vsak gen v permutaciji
  #v vsaki ponovitvi shrani samo najmanjso p vredn (M) !!
  pperm <-  2*pt(-abs(tt2), df=25-2)
  pmin[ii] <- min(pperm)
  }

for (kk in 1:500){

  permpval[kk] <-  (sum(pmin<=p.value[kk])+1)/(M+1)}


#izracun V

holmV[i] <- sum(holmp[11:500]<0.05)
bonfV[i] <- sum(bonfp[11:500]<0.05)
BHV[i] <- sum(BHp[11:500]<0.05)
permV[i] <- sum(permpval[11:500]<0.05)

holmS[i] <- sum(holmp[1:10]<0.05)
bonfS[i] <- sum(bonfp[1:10]<0.05)
BHS[i] <- sum(BHp[1:10]<0.05)
permS[i] <- sum(permpval[1:10]<0.05)

# R = V+S

holmR = holmV + holmS
bonfR = bonfV + bonfS
BHR = BHV + BHS
permR = permV + permS


#FWER
holmFWER <- sum(holmV>0)/B #FWER
bonfFWER <- sum(bonfV>0)/B #FWER
bhFWER <- sum(BHV>0)/B #FWER, tukaj je visja
permFWER <- sum(permV>0)/B

#pFDR in FDR

#pFDR IN FDR
holmpFDR <- mean(holmV/holmR,na.rm=TRUE)
holmFDR <- holmpFDR * mean((holmR)>0)

bonfpFDR <- mean(bonfV/bonfR,na.rm=TRUE)
bonfFDR <- bonfpFDR * mean((bonfR)>0)

BHpFDR <- mean(BHV/BHR,na.rm=TRUE)
BHFDR <- BHpFDR * mean(BHR>0)

permpFDR <- mean(permV/permR,na.rm=TRUE)
permFDR <- permpFDR * mean(permR>0)


rezultati <- rbind(c(holmFWER, holmpFDR, holmFDR), c(bonfFWER, bonfpFDR, bonfFDR), c(bhFWER, BHpFDR, BHFDR), c(permFWER, permpFDR, permFDR))

#moč
#mean(holmS/10)
#mean(bonfS/10)
#mean(BHS/10)
#mean(permS/10)



## ---- eval=FALSE-----------------------------------------------------------------------------
set.seed(1)

M <- 500
B <- 1000
n=25 #velikost vzorca za 1 gen

mu <- 0 # H0
#mu <- 0.5 HA
testnaPerm <- rep(NA,M)
vrednostp <- rep(NA,B)

for (k in 1:B){

y <- replicate(250,rnorm(n,0,1)) #generiraj 250 genov (po 25 ponovitev)
x <- rbinom(25,size=1,prob=0.5)
y[x==1,1:10] <- matrix(rnorm(sum(x==1)*10, mu, 1), ncol=10)

n1<-sum(x==0) #stevilo 0x
n2<-n-n1

#ponovi za razmetan X
#shrani originalno testno stat

d1 <- y[x==0,] #y-i ki imajo 0
d2 <- y[x==1,] #y-i ki imajo 1

povprx0 <- apply(d1,2,mean)
povprx1 <- apply(d2,2,mean)
razlikapovpr <- (povprx0-povprx1)
skupnaVariancakoren <- sqrt((apply(d1,2,var)*(n1-1)+apply(d2,2,var)*(n2-1))/(n1+n2-2))
stNapaka <- skupnaVariancakoren*sqrt(1/n1+1/n2)
testna <- sum((razlikapovpr/stNapaka))

#permutacije
for (i in 1:M){
  razmeciX <- x[sample(1:n)]
  d1 <- y[razmeciX==0,] #y-i ki imajo 0
  d2 <- y[razmeciX==1,] #y-i ki imajo 1
  povprx0 <- apply(d1,2,mean)
  povprx1 <- apply(d2,2,mean)
  razlikapovpr <- (povprx0-povprx1)
  skupnaVariancakoren <- sqrt((apply(d1,2,var)*(n1-1)+apply(d2,2,var)*(n2-1))/(n1+n2-2))
  stNapaka <- skupnaVariancakoren*sqrt(1/n1+1/n2)
  testnaPerm[i] <- sum(razlikapovpr/stNapaka)
}
vrednostp[k] <- (1+ (sum(abs(testnaPerm)>=abs(testna))) )/(M+1)
}


#moc testa
#mean(vrednostp <0.05)

#histogram
#hist(vrednostp)


## ---- eval=FALSE-----------------------------------------------------------------------------
set.seed(1)
M <- 500
B <- 1000
n=25 #velikost vzorca za 1 gen

mu <- 0 # H0
#mu <- 0.5 HA
vrednostp2 <- rep(NA,B)
testnaPerm <- rep(NA,M)

for (k in 1:B){
y <- replicate(250,rnorm(n,0,1))#generiraj 250 genov (po 25 ponovitev)
x <- rbinom(25,size=1,prob=0.5)
y[x==1,1:10] <- matrix(rnorm(sum(x==1)*10, mu, 1), ncol=10)

n1<-sum(x==0) #stevilo 0
n2<-n-n1

d1 <- y[x==0,] #y-i ki imajo 0
d2 <- y[x==1,] #y-i ki imajo 1

povprx0 <- apply(d1,2,mean)
povprx1 <- apply(d2,2,mean)
razlikapovpr <- (povprx0-povprx1)
stNapaka <- sqrt(apply(d1,2,var)/(n1) +apply(d2,2,var)/(n2))
testna <- sum((razlikapovpr/stNapaka))

#permutacije
for (i in 1:M){
  razmeciX <- x[sample(1:n)]
  d1 <- y[razmeciX==0,] #y-i ki imajo 0
  d2 <- y[razmeciX==1,] #y-i ki imajo 1
  povprx0 <- apply(d1,2,mean)
  povprx1 <- apply(d2,2,mean)
  razlikapovpr <- (povprx0-povprx1)
  stNapaka <- sqrt(apply(d1,2,var)/(n1) +apply(d2,2,var)/(n2))
  testnaPerm[i] <- sum((razlikapovpr/stNapaka))
}
vrednostp2[k] <-  (1+ (sum(abs(testnaPerm)>=abs(testna)) ))/(M+1)
}


#moc testa
#mean(vrednostp2<0.05)

#histogram
#hist(vrednostp2)



## ---- eval=FALSE-----------------------------------------------------------------------------
set.seed(1)
M <- 500
B <- 1000
n=25 #velikost vzorca za 1 gen

mu <- 0 # H0
#mu <- 0.5 HA
testnaPerm <- rep(NA,M)
vrednostp3 <- rep(NA,B)

for (k in 1:B){
y <- replicate(250,rnorm(n,0,1)) #generiraj 250 genov (po 25 ponovitev)
x <- rbinom(25,size=1,prob=0.5)
y[x==1,1:10] <- matrix(rnorm(sum(x==1)*10, mean=mu, sd= 1), ncol=10)
n1<-sum(x==0) #stevilo 0
n2<-n-n1

d1 <- y[x==0,] #y-i ki imajo 0
d2 <- y[x==1,] #y-i ki imajo 1

povprx0 <- apply(d1,2,mean)
povprx1 <- apply(d2,2,mean)
razlikapovpr <- (povprx0-povprx1)^2
testna = sum(razlikapovpr)

#permutacije
for (i in 1:M){
  razmeciX <- x[sample(1:n)]
  d1 <- y[razmeciX==0,] #y-i ki imajo 0
  d2 <- y[razmeciX==1,] #y-i ki imajo 1
  povprx0p <- (apply(d1,2,mean))
  povprx1p <- (apply(d2,2,mean))
  testnaPerm[i]<- sum((povprx0p-povprx1p)^2)
}
vrednostp3[k] <- (1+ (sum(abs(testnaPerm)>=abs(testna)))) /(M+1)
}

#moc testa
#mean(vrednostp3<0.05)

#histogram
#hist(vrednostp3)


## ---- eval=FALSE-----------------------------------------------------------------------------
C


## ---- echo=F---------------------------------------------------------------------------------
library(class)
library(kableExtra)
spremenljivke=1000
spremenljivkeblok=10
stBlokov = spremenljivke/spremenljivkeblok
enote = 100
mu <- c(0,0.5,1)


## ---- echo=F---------------------------------------------------------------------------------
#vse skupaj
set.seed(1)

funkcijaKNN <- function(mu,k){
PA <- rep(NA,1000)
G <- rep(NA,1000)
Fmera <- rep(NA,1000)

for (ii in 1:1000){

#simulacija podatkov
y <- replicate(1000,rnorm(100,mean=0,sd=sqrt(0.2)))
X <- sample(c(rep(0,50), rep(1,50)), replace=FALSE) #razred
XX <- sample(c(rep(0,950), rep(1,50)), replace=FALSE) #katerih 50 spremenljivk je drugace izrazenih

#popravi drugace izrazene spremenljivke
y[X==1,XX==1] <- matrix(rnorm(sum(X==1)*50, mean=mu, sd=sqrt(0.2)), ncol=50)

#pristej x3 vsaki spremenljivki znotraj bloga
#tako so spremenljivke znotraj bloka korelirane
for (i in 1:100){
  x3 <- rnorm(100,0,sqrt(0.8))
  x33 <- matrix(rep(x3,10),ncol=10)
  y[,(i*10-9):(i*10)] <- y[,(i*10-9):(i*10)]+x33
}

#izbor učne in testne množice

indeksi1 <- which(X==1)
indeksi0 <- which(X==0)

ucna0 <- sample(indeksi0, size=0.8*enote/2)
ucna1 <- sample(indeksi1, size=0.8*enote/2)
ucnaizbor <- c(ucna0, ucna1)

ucna <- y[ucnaizbor,]
testna <- y[-ucnaizbor,]

#k-nn algoritem
KNN <- knn(ucna,testna,cl=X[ucnaizbor],k)

tabela <- table(KNN,X[-ucnaizbor])

#izracun tocnosti ..

PA[ii] <- sum(diag(tabela))/(0.2*enote)
PA1 <- tabela[1,1]/sum(tabela[,1])
PA2 <- tabela[2,2]/sum(tabela[,2])

G[ii] <- sqrt(PA1*PA2)

beta <- 1
Fmera[ii] <- (1+beta^2)*tabela[1,1]/((1+beta^2)*tabela[1,1]+beta^2*tabela[2,1]+tabela[1,2])
}

return(cbind(PA,G,Fmera))
}

knn1 <- lapply(mu,k=1,funkcijaKNN)
knn1 <- matrix(unlist(knn1), ncol=9)
knn1 <- apply(knn1,2,mean)

knn3 <- lapply(mu,k=3,funkcijaKNN)
knn3 <- matrix(unlist(knn3), ncol=9)
knn3 <- apply(knn3,2,mean)

knn5 <- lapply(mu,k=5,funkcijaKNN)
knn5 <- matrix(unlist(knn5), ncol=9)
knn5 <- apply(knn5,2,mean)



## ---- echo=F---------------------------------------------------------------------------------
vrednostK <- c(1,3,5)
mu0 <- rbind(knn1[1:3],knn3[1:3],knn5[1:3])
mu05 <- rbind(knn1[4:6],knn3[4:6],knn5[4:6])
mu1 <- rbind(knn1[7:9],knn3[7:9],knn5[7:9])

#shrani podatke v tabelo
Tabela1 <- cbind(vrednostK, mu0)
colnames(Tabela1) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela1,caption = "mu=0") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")


## ---- echo=F---------------------------------------------------------------------------------
Tabela3 <- cbind(vrednostK, mu05)
colnames(Tabela3) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela3,caption = "mu=0.5") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")


## ---- echo=F---------------------------------------------------------------------------------
Tabela5 <- cbind(vrednostK, mu1)
colnames(Tabela5) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela5,caption = "mu=1") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")


## ---- eval=FALSE-----------------------------------------------------------------------------
library(class)
library(kableExtra)
spremenljivke=1000
spremenljivkeblok=10
stBlokov = spremenljivke/spremenljivkeblok
enote = 100
mu <- c(0,0.5,1)

#vse skupaj
set.seed(1)

funkcijaKNN <- function(mu,k){
PA <- rep(NA,1000)
G <- rep(NA,1000)
Fmera <- rep(NA,1000)

for (ii in 1:1000){

#simulacija podatkov
y <- replicate(1000,rnorm(100,mean=0,sd=sqrt(0.2)))
X <- sample(c(rep(0,50), rep(1,50)), replace=FALSE) #razred
XX <- sample(c(rep(0,950), rep(1,50)), replace=FALSE) #katerih 50 spremenljivk je drugace izrazenih

#popravi drugace izrazene spremenljivke
y[X==1,XX==1] <- matrix(rnorm(sum(X==1)*50, mean=mu, sd=sqrt(0.2)), ncol=50)

#pristej x3 vsaki spremenljivki znotraj bloga
#tako so spremenljivke znotraj bloka korelirane
for (i in 1:100){
  x3 <- rnorm(100,0,sqrt(0.8))
  x33 <- matrix(rep(x3,10),ncol=10)
  y[,(i*10-9):(i*10)] <- y[,(i*10-9):(i*10)]+x33
}

#izbor učne in testne množice

indeksi1 <- which(X==1)
indeksi0 <- which(X==0)

ucna0 <- sample(indeksi0, size=0.8*enote/2)
ucna1 <- sample(indeksi1, size=0.8*enote/2)
ucnaizbor <- c(ucna0, ucna1)

ucna <- y[ucnaizbor,]
testna <- y[-ucnaizbor,]

#k-nn algoritem
KNN <- knn(ucna,testna,cl=X[ucnaizbor],k)

tabela <- table(KNN,X[-ucnaizbor])

#izracun tocnosti ..

PA[ii] <- sum(diag(tabela))/(0.2*enote)
PA1 <- tabela[1,1]/sum(tabela[,1])
PA2 <- tabela[2,2]/sum(tabela[,2])

G[ii] <- sqrt(PA1*PA2)

beta <- 1
Fmera[ii] <- (1+beta^2)*tabela[1,1]/((1+beta^2)*tabela[1,1]+beta^2*tabela[2,1]+tabela[1,2])
}

return(cbind(PA,G,Fmera))
}

knn1 <- lapply(mu,k=1,funkcijaKNN)
knn1 <- matrix(unlist(knn1), ncol=9)
knn1 <- apply(knn1,2,mean)

knn3 <- lapply(mu,k=3,funkcijaKNN)
knn3 <- matrix(unlist(knn3), ncol=9)
knn3 <- apply(knn3,2,mean)

knn5 <- lapply(mu,k=5,funkcijaKNN)
knn5 <- matrix(unlist(knn5), ncol=9)
knn5 <- apply(knn5,2,mean)

vrednostK <- c(1,3,5)
mu0 <- rbind(knn1[1:3],knn3[1:3],knn5[1:3])
mu05 <- rbind(knn1[4:6],knn3[4:6],knn5[4:6])
mu1 <- rbind(knn1[7:9],knn3[7:9],knn5[7:9])

#shrani podatke v tabelo
Tabela1 <- cbind(vrednostK, mu0)
colnames(Tabela1) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela1,caption = "mu=0") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")

Tabela3 <- cbind(vrednostK, mu05)
colnames(Tabela3) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela3,caption = "mu=0.5") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")

Tabela5 <- cbind(vrednostK, mu1)
colnames(Tabela5) <- c("K", "Napovedna točnost", "G povprečje", "F1-mera")
kable(Tabela5,caption = "mu=1") %>%
  kable_styling("striped",full_width = F, font_size=10)%>%
  column_spec(1, bold = T, color = "black", background = "yellow")


