######SARIMA#######

setwd ("c:\\Users\\Ana\\Documents\\PUC\\Simulação\\Base")
  
##Bibliotecas#
library(Kendall)
library(openxlsx)
library(astsa)
library(forecast)
library(urca)
library(tseries)
##configração dos plots#
par(mfrow=c(1,1))

##Série de ena mensal#
ena<-read.xlsx("ena.xlsx")

##Definindo a série de Ena como série temporal#
x<- ts(ena[,"SUDESTE"], start=c(2000,1), end = c(2013,12), frequency=12)
x1<- ts(ena[,"PREVSD"], start=c(2014,1), end = c(2014,12), frequency=12)
plot(a, ylab= "SE/CO")

a<- ts(ena[,"SD"], start=c(2000,1), end = c(2014,12), frequency=12)
c<- matrix(ena[,"MLTSD"], nrow = 12, ncol=1)

###Estatísticas descritivas####
m = mean(a)
minx = min(a)
maxx = max(a)
medx = median(a)
DP = sd(a)
CV = DP/m *100
AS = 3*(m- medx)/DP
QT<- quantile(a)
QT

####Dados####
am<- matrix(x, 12)#matriz com os dados por mês na linha
aa <- matrix(NA, nrow = 12)
for (i in 1:12){
  aa[i,] <- mean(am[i,])
}

amx <- matrix (NA, nrow = 15, ncol = 12) #matriz com os dados por mês na coluna

for (i in 1:15){
  amx [i,] <- am[,i]  
}

MMnames= list("Jan", "Fev", "Mar", "Abr", "Mai", "Jun", "Jul", "Ago", "Set", "Out", "Nov", "Dez")

boxplot(amx, names= MMnames, ylab = "SE/CO") 

####Estatísticas descritivas dos meses####
 ED <- matrix(NA, nrow= 6, ncol = 12)

for (i in 1:12){
  ED[1,i] <- min(amx[,i])
  ED[2,i] <- mean(amx[,i])
  ED[3,i] <- median(amx[,i])
  ED[4,i] <- max(amx[,i])
  ED[5,i] <- sd(amx[,i])
}

####Testes de estacionariedade####
#Teste MannKendall - Tendência
MannKendall(x)
#Conclusão dos teste: Série estacionaria#

####Transformação de BOXCOX####
l1=forecast::BoxCox.lambda(x, method=c("loglik"),lower=-1, upper=2)
z<-forecast::BoxCox(x,l1)

#Teste MannKendall - Tendência
MannKendall(z)
hist(z)

#Conclusão dos teste: Série estacionaria#

#### Função de autocorrelação e autocorrelação parcial ts ####
## Calculando, mas não plotando, acf
acfx<-acf(x,lag.max = 24, plot = FALSE)
## Transformando lags de anos para mês
acfx$lag <- acfx$lag*12
# Plot do acf 
plot(acfx,ylab = "ACF", xlab="Lag(meses)")

## Calculando, mas não plotando, pacf
pacfx<-pacf(x,lag.max = 24, plot=FALSE)
## Transformando lags de anos para mês
pacfx$lag <- pacfx$lag*12
# Plot do pacf 
plot(pacfx, ylab = "PACF", xlab="Lag(meses)")

#### Analisando fit do modelo Sarima ####
auto.arima(x)

#### Teste de independência###
y<-sarima(z,4,0,2,1,0,1,12)
r<-y$fit$residuals
Me<- mean(r)
DPe <- sd(r)
y
###Teste de normalidade dos resíduos###
tseries::jarque.bera.test (r)
shapiro.test(r)

Box.test(r, lag = 10, type = "Ljung-Box")

density(r)
hist(r)

###Teste de homocedasticidade###
set.seed(2^31-1)
r1<-split(r,sample(r, size=2))
A1<-r1$`0.141892798144254`
A2<-r1$`-0.064956683379853`

y<- c(A1, A2)
group <- as.factor(c(rep(1, length(A1)), rep(2, length(A2))))
car::leveneTest(y,group)

####Modelo de previsão####
####Simulação####
w<- sarima.for(z,12,4,0,2,1,0,1,12)

w1<- InvBoxCox(w$pred, l1)

MAPE <- mean(abs((x1- w1)/x1))

MRSE <- sqrt(mean((x1-w1)^2))

RSQT <- (sum(x1 - w1)^2)/(sum(x1- mean(x1))^2)

####Simulação####
#Matriz de erros aleatório#
#dmn <- list(month.abb, unique(floor(time(r))))
#rc<-as.data.frame(t(matrix(r, 12, dimnames = dmn)))

rr<-matrix(r,12)

MJ<- matrix(NA, nrow = 12, ncol = 2)

for (i in 1:12) {
  MJ[i,1] <- mean(rr[i,])
  MJ[i,2] <- sd(rr[i,])
}

nrSamples = 100
SeriesSinteticas <- matrix(NA, nrow = 12, ncol = nrSamples)

set.seed(2^20-1)
for (i in 1:12) {
  SeriesSinteticas[i,]<- rnorm(n = nrSamples, mean = MJ[i,1], sd = MJ[i,2])
}

PREV<- matrix(w$pred, nrow= 12, ncol=1)

MC<- matrix(NA, nrow = 12, ncol = nrSamples)
for (i in 1:nrSamples) {
  MC[,i]<- SeriesSinteticas[,i] + PREV
}
INVMC <- InvBoxCox(MC, l1)

####Plot#####
MonteC <- matrix(NA, nrow = 12, ncol = 1)
for (i in 1:12) {
  MonteC[i,]<- mean(MC[i,])
}

MCw1<- InvBoxCox(MonteC, l1)

MAPE4 <- mean(abs((x1- MCw1)/x1))

MRSE1 <- sqrt(mean((x1-MCw1)^2))

RSQT1 <- (sum(x1 - MCw1)^2)/(sum(x1- mean(x1))^2)

#### Gráficos dos resultados####
###Intervalo de confiança da previsão####
IC <- matrix(NA, nrow=12, ncol=2)
ICs=w$pred + 1.96 * w$se 
ICi= w$pred - 1.96*w$se
IC[,1] =InvBoxCox(ICs,lambda = l1)
IC[,2] =InvBoxCox(ICi,lambda = l1)

#gráfico com previsão e série histórica
MCmt <- matrix(x1, ncol=1)

plot(MCmt, main= "",type= "l", ylab= "SE/CO", xlab= "Período", col= "black", ylim = c(0,120000), lwd=2.5)
for (i in 1:ncol(INVMC)){
  lines(INVMC[,i], col = "gray",lwd=2.5)
  lines(MCmt, col="black",lwd=2.5)
  lines(IC[,1], col='dimgray', lwd=2.5, lty=2)
  lines(IC[,2], col='dimgray',lwd=2.5, lty=2)
  }
legend(6, 120000, legend=c("Cenários","Dados observados","Intervalo de confiança - 95%"),lty=c(1,1,2), lwd=c(2.5,2.5,2.5),col=c("gray","black","dimgray"), pt.cex=10, cex=0.8)


##gráfico com média de cenários de monte carlo e intervalo de confiança
plot(MCmt,type= "l", ylab= "SE/CO", xlab= "Período", col ="black", ylim = c(0,80000),lwd=2.5, lty=5)
lines(InvBoxCox(MonteC, lambda = l1),col='gray59', lwd=2.5, lty=1)
lines(InvBoxCox(PREV,l1), col='gray30',lwd=2.5, lty=3)
legend(6.5, 75000, legend=c("Dados observados","Média dos cenários","previsão"),lty=c(5,1,3), lwd=c(2.5,2.5,2.5), col=c("black","gray", "gray30"), pt.cex=10, cex=0.7)

###GRÁFICO COLORIDO
plot(MCmt, main= "NE",type= "l", ylab= "NE", xlab= "Período", col= "black", ylim = c(0,40000), lwd=2.5)
for (i in 1:ncol(INVMC)){
  lines(INVMC[,i], col = "gray",lwd=2.5)
  lines(MCmt, col="black",lwd=2.5)
}
legend(6, 30000, legend=c("série histórica","Cenários"),lty=c(1,1), lwd=c(2.5,2.5), col=c("gray","black"), pt.cex=10, cex=0.8)


##gráfico com média de cenários de monte carlo e intervalo de confiança
plot(MCmt,type= "l", ylab= "NE", xlab= "Período", col ="black", ylim = c(0,25000),lwd=2.5, lty=5)
lines(InvBoxCox(MonteC, lambda = l1),col='red', lwd=2.5, lty=1)
lines(IC[,1], col='dimgray', lwd=2.5, lty=3)
lines(IC[,2], col='dimgray',lwd=2.5, lty=3)
lines(InvBoxCox(PREV,l1), col='green',lwd=2.5, lty=1)
legend(6.5, 25000, legend=c("série histórica"," Média dos Cenários","Intervalo de confiança - 95%", "previsão"),lty=c(5,1,1,2), lwd=c(2.5,2.5,2.5,2.5), col=c("black","red","dimgray", "blue4"), pt.cex=10, cex=0.6)


#### Análise dos resultados####

ks.test(aa, w1, alternative = c("two.sided"))  # previsão e série histórica
ks.test(aa, MCw1, alternative = c("two.sided"))  # previsão e MC

t.test(aa, w1, alternative = c("two.sided")) # previsão e série histórica
t.test(aa, MCw1, alternative = c("two.sided")) # previsão e MC

LEV1<- c(aa, w1)
group2 <- as.factor(c (rep(1, length(aa)), rep(2, length(w1)))) 
car::leveneTest(LEV1,group2)   # previsão e série histórica

LEV2<- c(aa, MCw1)
group3 <- as.factor(c(rep(1, length(aa)), rep(2, length(MCw1))))
car::leveneTest(LEV2,group3) # previsão e MC


####MLT aderencia
ks.test(c, w1, alternative = c("two.sided"))  # previsão e série histórica
ks.test(c, MCw1, alternative = c("two.sided"))  # previsão e MC

t.test(c, w1, alternative = c("two.sided")) # previsão e série histórica
t.test(c, MCw1, alternative = c("two.sided")) # previsão e MC

LEV1<- c(c, w1)
group2 <- as.factor(c (rep(1, length(c)), rep(2, length(w1)))) 
car::leveneTest(LEV1,group2)   # previsão e série histórica

LEV2<- c(c, MCw1)
group3 <- as.factor(c(rep(1, length(c)), rep(2, length(MCw1))))
car::leveneTest(LEV2,group3) # previsão e MC

