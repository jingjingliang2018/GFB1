# Load packages 
library(nlme)

# Load plot-level data
# Download GFB1_data_figshare.xlxs from Figshare and convert to a csv file

data<- read.csv("GFB1_data_figshare.csv")
data <- subset(data, P>0)
data <- subset(data, S>0)

quantile(data$S,0.99996)
quantile(data$P,0.99996)

data1 <- subset(data,data$S<=270 & data$P<=533 & data$S >0 & data$P>0)   # removed 894 plots with 0 or extreme S and P values

###################################################################
############## Derive Global GeoRF Estimation #####################
###################################################################


logP <- log(data1$P)
# jig coordinates to avoid duplicated values
Lon1 <- data1$Lon+ runif(length(data1$Lon),-0.0001,0.0001)
Lat1 <- data1$Lat+ runif(length(data1$Lat),-0.0001,0.0001)
data1 <- cbind.data.frame(data1, logP, Lat1, Lon1)


############ Loop ##################
coef <- matrix(0, nrow=10000, ncol=20)	# Coef Matrix	


for(i in 1: 10000) {
  tryCatch({

	training <- data1[sample(1:nrow(data1), 500, replace=FALSE),]   # turn 'replace' off to maximize inclusion of new plots
	
	logS <- log(training$S)
	training <- cbind.data.frame(training, logS)

	gls1 <- gls(logP~ logS + G + T3 + C1 + C3 + PET + IAA + E, data=training, method="ML", corr= corSpher(form = ~ Lon1 + Lat1, nugget = TRUE), control=glsControl(singular.ok=TRUE))

	coef[i,3] <- i
	coef[i,4] <- logLik (gls1)
	coef[i,5] <- AIC (gls1)
	coef[i,6]<- BIC (gls1)
		#Generalized coefficient of determination
		gls0 <- gls(logP~ 1, data=training, method="ML")  
		R2   <- 1-exp(logLik(gls0)-logLik(gls1))^(2/500)
		coef[i,7]<- R2

	coef[i,8]  <- coef(gls1)[1]  	
	coef[i,9]  <- coef(gls1)[2]
	coef[i,10] <- coef(gls1)[3]
	coef[i,11] <- coef(gls1)[4]
	coef[i,12] <- coef(gls1)[5]
	coef[i,13] <- coef(gls1)[6]
	coef[i,14] <- coef(gls1)[7]
	coef[i,15] <- coef(gls1)[8]
	coef[i,16] <- coef(gls1)[9]
	coef[i,17] <- 0

	# Baseline (S=1) productivity
	# logS + B1 + T3 + C1 + C3 + PET + IAA + E
	newdata <- data.frame(logS=0, G=mean(training$G), T3=mean(training$T3), C1=mean(training$C1), C3=mean(training$C3),PET=mean(training$PET), IAA=mean(training$IAA), E=mean(training$E))
	coef[i,20]  <- exp(predict(gls1,newdata))

	#counter
	cat(i, " of ", 1000, date(),"Theta=",coef(gls1)[2], "R2=", R2, "\n" )

	#remove files
	rm(training, newdata, gls1, R2)

  }, error=function(e){})
	}

coef_df <- as.data.frame(coef)

names(coef_df) <- c("0", "0", "i", "Loglik", "AIC", "BIC", "R2","const","theta", "B", "T3", "C1", "C3", "PET", "IAA", "E", "0", "0", "0", "P_1")

write.csv(coef_df, "global_estimates.csv")



#######################################################################################################
############## Draw estimated Biodiversity-Productivity Relationship (BPR) curves #####################
#######################################################################################################

data<- read.csv("global_estimates.csv")

theta <- data$theta
mean(theta)
P_base <- mean(data$P_1)

# Predict P over an increased S from 1 to global max (271), which corresponds to S_hat from 100/271 to 100
S     <- seq(1,271,1)
S_hat <- S*100/271

P_est <- data.frame(matrix(0, 10000, ncol =273))
P_est[,1] <- P_base
P_est[,2] <- theta

for (i in 1:10000){
  P_est[i,3:273] <- P_est[i,1] * S ^ P_est[i,2]
  }

# demosntration plot only shows the first 18 iterations
plot(S_hat,colMeans(P_est[,3:273]), ylim=c(0,20), type="l",col = "blue", ylab="P")
for (i in 1:16){
  P_est[i,3:273] <- P_est[i,1] * S ^ P_est[i,2]
  lines(S_hat,P_est[i,3:273],col = "green")
}
# Confidence intervals
lines(S_hat,colMeans(P_est[,3:273])+1.96*apply(P_est[,3:273], 2, sd)/sqrt(10000), ylim=c(0,20), type="l",col = "red")
lines(S_hat,colMeans(P_est[,3:273])-1.96*apply(P_est[,3:273], 2, sd)/sqrt(10000), ylim=c(0,20), type="l",col = "red")



# End of the code
