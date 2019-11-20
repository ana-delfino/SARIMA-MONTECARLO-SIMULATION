 ######SARIMA#######

setwd ("c:\\Users\\Ana\\Documents\\PUC\\Simulação\\Base")
  
##Bibliotecas#
library(Kendall)
library(openxlsx)
library(astsa)
library(forecast)
library(urca)
library(tseries)
library(ggplot2)
  
##configração dos plots#
par(mfrow=c(1,1))

##Série de ena mensal#
ena<-read.xlsx("ena.xlsx")

##Definindo a série de Ena como série temporal#
x<- ts(ena[,"NORDESTE"], start=c(2000,1), end = c(2013,12), frequency=12)
x1<- ts(ena[,"PREVND"], start=c(2014,1), end = c(2014,12), frequency=12)
b<- ts(ena[,"ND"], start=c(2000,1), end = c(2014,12), frequency=12)
c<- matrix(ena[,"MLTND"], nrow = 12, ncol=1)

###Estatísticas descritivas####
m = mean(b)
minx = min(b)
maxx = max(b)
medx = median(b)
DP = sd(b)
CV = DP/m *100
AS = 3*(m- medx)/DP
QT<- quantile(b)
AS1 = (9404.75+ 2748 - 2*medx)/(9404.75- 2728)

plot(b, ylab="NE")

bm<- matrix(x, 12)#matriz com os dados por mês na linha

bmx <- matrix (NA, nrow = 15, ncol = 12) #matriz com os dados por mês na coluna

for (i in 1:15){
bmx [i,] <- bm[,i]  
}

MMnames= list("Jan", "Fev", "Mar", "Abr", "Mai", "Jun", "Jul", "Ago", "Set", "Out", "Nov", "Dez")


boxplot(bmx, names= MMnames, ylab = "NE")     
am<- matrix(x, 12)#matriz com os dados por mês na linha
aa <- matrix(NA, nrow = 12)
for (i in 1:12){
  aa[i,] <- mean(am[i,])
}

####Estatísticas descritivas dos meses####
EDb <- matrix(NA, nrow= 6, ncol = 12)


for (i in 1:12){
  EDb[1,i] <- min(bmx[,i])
  EDb[2,i] <- mean(bmx[,i])
  EDb[3,i] <- median(bmx[,i])
  EDb[4,i] <- max(bmx[,i])
  EDb[5,i] <- sd(bmx[,i])
  EDb[6,i] <- EDb[5,i]/EDb[2,i]
}

format(EDb, decimal.mark = ',')

####Testes de estacionariedade####
#Teste MannKendall - Tendência
MannKendall(x)
#Conclusão dos teste: Série estacionaria#

####Transformação de BOXCOX####
l1=forecast::BoxCox.lambda(x, method=c("loglik"),lower=0, upper=2)
z<-forecast::BoxCox(x,l1)

#Teste MannKendall - Tendência
MannKendall(z)
hist(z)

#Conclusão dos teste: Série estacionaria#

#### Função de autocorrelação e autocorrelação parcial ts ####
## Calculando, mas não plotando, acf
acfx<-acf(X12,lag.max = 24, plot = FALSE)
## Transformando lags de anos para mês
acfx$lag <- acfx$lag*12
# Plot do acf 
plot(acfx,ylab = "ACF", xlab="Lag(meses)")

## Calculando, mas não plotando, pacf
pacfx<-pacf(X12,lag.max = 24, plot=FALSE)
## Transformando lags de anos para mês
pacfx$lag <- pacfx$lag*12
# Plot do pacf 
plot(pacfx, ylab = "PACF", xlab="Lag(meses)")

#### Analisando fit do modelo Sarima ####
auto.arima(z, trace=TRUE, ic=c("bic"))

#### Teste de independência###
fitM<-sarima(z,2,0,0,2,0,0,12)
r<-fitM$fit$residuals
Me<- mean(r)
DPe <- sd(r)

Box.test(r, lag = 1, type = "Ljung-Box")

fitM
###Teste de normalidade dos resíduos####
tseries::jarque.bera.test (r)
shapiro.test(r)

###Teste de homocedasticidade###
set.seed(2^31-1)
r1<-split(r,sample(r, size=2))
A1<-r1$`-0.206021116179237`
A2<-r1$`-0.133163471154843`
y<- c(A1, A2)
group <- as.factor(c(rep(1, length(A1)), rep(2, length(A2))))
car::leveneTest(y,group)

###forecast####
w<- sarima.for(z,12,2,0,0,2,0,0,12)
w1
w1<- InvBoxCox(w$pred, l1)

MAPE <- mean(abs((x1- w1)/x1))

MRSE <- sqrt(mean((x1-w1)^2))

MCmt <- matrix(x1, ncol=1)

MADD <- mean(abs(x1- w1))

####Simulação####
#Matriz de erros aleatório#

rr<-matrix(r,12) # matrix com os resíduos dos modelos ajustado, tem 12lx14c

MJ<- matrix(NA, nrow = 12, ncol = 2) # matriz para receber a média e dp de cada mês

##loop para inserir média e dp na matriz MJ
for (i in 1:12) {
  MJ[i,1] <- mean(rr[i,])
  MJ[i,2] <- sd(rr[i,])
}
nrSamples = 100 # números de séries sintáticas

SeriesSinteticas <- matrix(NA, nrow = 12, ncol = nrSamples) #matriz de séries sintáticas

set.seed(2^20-1) # semente
#loop para conseguir 100 séries para casa mÊs com distribuição de acordo com MJ
for (i in 1:12) {
SeriesSinteticas[i,]<- rnorm(n = nrSamples, mean = MJ[i,1], sd = MJ[i,2])
}

PREV<- matrix(w$pred, nrow= 12, ncol=1) ## matriz com previsão do modelo 12 mês a frente
MCW2 <- InvBoxCox(PREV,l1) ### inversa da matriz de previsão

#matriz com as séries sintéticas - Simulação de Monte carlo
MC<- matrix(NA, nrow = 12, ncol = nrSamples)
for (i in 1:nrSamples) {
  MC[,i]<- PREV + SeriesSinteticas[,i] 
}

##Matriz anterior, pós inversão da transformação box.cox
INVMC <- InvBoxCox(MC, l1)

MonteC <- matrix(NA, nrow = 12, ncol = 1) # média das simulações de monte carlo
for (i in 1:12) {
  MonteC[i,]<- mean(MC[i,])
}

MCw1<- InvBoxCox(MonteC, l1) # média das simulação 

MAPE1 <- mean(abs((x1- MCw1)/x1)) #MAPE usando a média das simulações de MC

MRSE1 <- sqrt(mean((x1-MCw1)^2)) # RMSE usando a média das simualações

MADD1 <- mean(abs(x1- MCw1))


#### Gráficos dos resultados####
###Intervalo de confiança da previsão####
IC <- matrix(NA, nrow=12, ncol=2)
ICs=w$pred + 1.96 * w$se 
ICi= w$pred - 1.96*w$se
IC[,1] =InvBoxCox(ICs,lambda = l1)
IC[,2] =InvBoxCox(ICi,lambda = l1)

#gráfico com previsão e série histórica
plot(MCmt, main= "NE",type= "l", ylab= "NE", xlab= "Período", col= "black", ylim = c(0,40000), lwd=2.5)
for (i in 1:ncol(INVMC)){
  lines(INVMC[,i], col = "gray",lwd=2.5)
  lines(MCmt, col="black",lwd=2.5)
  lines(IC[,1], col='dimgray', lwd=2.5, lty=2)
  lines(IC[,2], col='dimgray',lwd=2.5, lty=2)
}
legend(6, 40000, legend=c("série histórica","Cenários","Intervalo de confiança - 95%"),lty=c(1,1,2), lwd=c(2.5,2.5,2.5),col=c("gray","black","dimgray"), pt.cex=10, cex=0.8)


##gráfico com média de cenários de monte carlo e intervalo de confiança
plot(MCmt,type= "l", ylab= "NE", xlab= "Período", col ="black", ylim = c(0,15000),lwd=2.5, lty=5)
lines(InvBoxCox(MonteC, lambda = l1),col='gray59', lwd=2.5, lty=1)
lines(InvBoxCox(PREV,l1), col='gray30',lwd=2.5, lty=3)
legend(6.5, 20000, legend=c("série histórica"," Média dos Cenários", "previsão"),lty=c(5,1,3), lwd=c(2.5,2.5,2.5), col=c("black","gray", "gray30"), pt.cex=10, cex=0.7)

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


###Aderencia a série historica
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
