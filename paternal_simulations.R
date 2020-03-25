library(scales)

N <- 2000 # This needs to be an even number for clean pairing
precis <- 100000
mutate <- 0.0001
V0 <- 0.0001

phis <- 0.5
Y_D <- 2.5
Y_Cboost <- 2.5

shiftsize <- 1/4
t1 <- -2
k1 <- 12.5
t2 <- -0.4
k2 <- 25

delta_D <- 1
X <- 5
A <- 20

#########

lwd <- 3
repeet <- 1

ts <- seq(-3,0,length.out=precis)

dev.new(width=7, height=8)
par(mfcol=c(2,1))
par(mar=c(4,3,0,0.5)+0.5)

######################################################################
######################################################################

Y_C <- Y_D

fracDs <- rep(NA,precis)
meanV_D <- rep(NA,precis)
meanV_C <- rep(NA,precis)
runDs <- matrix(ncol=length(phis) * repeet,nrow=length(ts))
runphis <- rep(NA,length(phis) * repeet)
runcolors <- rep(NA,length(phis) * repeet)
kappas <- rep(NA,length(ts))

l <- 1
counter <- 0	
for (l in 1:repeet) {
	for (k in 1:length(phis)) {
		counter <- 1 + counter
		phi <- phis[k]
		Ds <- rep(0,N)
		cat(paste(counter,"/",repeet*length(phis)," phi =",phi," ###########\n"))
		for (i in 1:length(ts)) {
			
			kappa <- 0 + shiftsize / (1 + exp(-k1*(ts[i]-t1))) + shiftsize / (1 + exp(-k2*(ts[i]-t2)))
			mu <- kappa
			kappas[i] <- kappa
						
			####
			
			fracDs[i] <- mean(Ds)
			
			partners <- sample(1:N,N,replace=F)
			part1 <- partners[1:floor(N/2)]
			part2 <- partners[floor(N/2 + 1):N]
			
			D1 <- Ds[part1]
			D2 <- Ds[part2]
			Ds <- c(D1,D2)
			
			V1 <- rep(NA,length(D1))
			V2 <- rep(NA,length(D1))
			j <- 1
			for (j in 1:length(D1)) {
			
				V1[j] <-  Pi_function4(phi=phi, kappa=kappa, mu=mu, delta_D=delta_D,Y_D=Y_D, Y_C=Y_C, X=X, D_i=D1[j], D_j=D2[j],A=A)
				V2[j] <-  Pi_function4(phi=phi, kappa=kappa, mu=mu, delta_D=delta_D,Y_D=Y_D, Y_C=Y_C, X=X, D_i=D2[j], D_j=D1[j],A=A)	
			}
			Vraw <- c(V1,V2)
		
			V <- Vraw + V0
			V[V<0] <- 0
			V <- V / mean(V)
			
			wh <- which(Ds == 1)
			wh <- which(Ds == 0)
			
			Ds <- sample(Ds,N,replace=T,prob=V)
			Nmutants <- rpois(1,lambda=mutate * N)
			mutants <- sample(1:N,Nmutants,replace=FALSE)
			Ds[mutants] <- 1 - Ds[mutants]
			
			if(i %% 1000 == 0) {cat(paste("     ", round(i,2),"; kappa=",kappa,"; mu=",mu,"; D=",round(fracDs[i],3),"; ", Nmutants,"mutants\n"))}
		}
		runDs[,counter] <- fracDs
		runphis[counter] <- phi
		
	}
}

runDsBOT <- runDs
runphisBOT <- runphis
runcolorsBOT <- runcolors

runDs[1:10]

######################################################################
######################################################################

plot(c(-3,0),c(0,0),bty="n",col="white",ylab="",xlab="Time",yaxt="n",xaxt="n",type="l",lty=1,ylim=c(-0.25,1),lwd=2)
axis(1,c(-3,-2,-1,0),c("3mya","2mya","1mya","present"))
mtext("Dads",2,at=0.5,line=2)
axis(2,c(0,1))

alph<-0.25
polygon(c(-2.5,-1.5,-1.5,-2.5),c(-0.25,-0.25,-0.05,-0.05),col=alpha("black", alph),border=NA)
text(-2.0,-0.15,"Increasing\ncomplementarities",cex=1)
polygon(c(-0.6,-0.2,-0.2,-0.6),c(-0.25,-0.25,-0.05,-0.05),col=alpha("black",alph),border=NA)
text(-0.39,-0.15,"Increasing\n comple-\nmentarities",cex=0.7)

#lines(ts,kappas,pch=16,lty=2,lwd=1,col="black")

runcolorsBOT[1] <- "red"
for (i in 1:ncol(runDsBOT)) {
	lines(ts,runDsBOT[,i],pch=16,lwd=lwd,col=runcolorsBOT[i])
}

###

write.csv(data.frame(ts,runDsBOT[,1]),"run_on-off.csv")

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################


fracDs <- rep(NA,precis)
meanV_D <- rep(NA,precis)
meanV_C <- rep(NA,precis)

runDs <- matrix(ncol=length(phis) * repeet,nrow=length(ts))
runphis <- rep(NA,length(phis) * repeet)
runcolors <- rep(NA,length(phis) * repeet)
kappas <- rep(NA,length(ts))


Y_C <- Y_D + Y_Cboost

l <- 1
counter <- 0	
for (l in 1:repeet) {
	for (k in 1:length(phis)) {
		counter <- 1 + counter
		phi <- phis[k]
		Ds <- rep(0,N)
		cat(paste(counter,"/",repeet*length(phis)," phi =",phi," ###########\n"))
		for (i in 1:length(ts)) {
			
			kappa <- 0 + shiftsize / (1 + exp(-k1*(ts[i]-t1))) + shiftsize / (1 + exp(-k2*(ts[i]-t2)))
			mu <- kappa
			kappas[i] <- kappa
						
			####
			
			fracDs[i] <- mean(Ds)
			
			partners <- sample(1:N,N,replace=F)
			part1 <- partners[1:floor(N/2)]
			part2 <- partners[floor(N/2 + 1):N]
			
			D1 <- Ds[part1]
			D2 <- Ds[part2]
			Ds <- c(D1,D2)
			
			V1 <- rep(NA,length(D1))
			V2 <- rep(NA,length(D1))
			j <- 1
			for (j in 1:length(D1)) {
				V1[j] <-  Pi_function4(phi=phi, kappa=kappa, mu=mu, delta_D=delta_D,Y_D=Y_D, Y_C=Y_C, X=X, D_i=D1[j], D_j=D2[j],A=A)
				V2[j] <-  Pi_function4(phi=phi, kappa=kappa, mu=mu, delta_D=delta_D,Y_D=Y_D, Y_C=Y_C, X=X, D_i=D2[j], D_j=D1[j],A=A)	
			}							
			Vraw <- c(V1,V2)
		
			V <- Vraw + V0
			V[V<0] <- 0
			V <- V / mean(V)
			
			wh <- which(Ds == 1)
			wh <- which(Ds == 0)
			
			
			Ds <- sample(Ds,N,replace=T,prob=V)
			Nmutants <- rpois(1,lambda=mutate * N)
			mutants <- sample(1:N,Nmutants,replace=FALSE)
			Ds[mutants] <- 1 - Ds[mutants]
			
			if(i %% 1000 == 0) {cat(paste("     ", round(i,2),"; kappa=",kappa,"; mu=",mu,"; D=",round(fracDs[i],3),"; ", Nmutants,"mutants\n"))}
		}
		
		runDs[,counter] <- fracDs
		runphis[counter] <- phi
	}
}

runDsBOT <- runDs
runphisBOT <- runphis
runcolorsBOT <- runcolors

runDs

######################################################################
######################################################################
######################################################################

plot(c(-3,0),c(0,0),bty="n",col="white",ylab="",xlab="Time",yaxt="n",xaxt="n",type="l",lty=1,ylim=c(-0.25,1),lwd=2)
axis(1,c(-3,-2,-1,0),c("3mya","2mya","1mya","present"))
mtext("Dads",2,at=0.5,line=2)
axis(2,c(0,1))

alph<-0.25
polygon(c(-2.5,-1.5,-1.5,-2.5),c(-0.25,-0.25,-0.05,-0.05),col=alpha("black", alph),border=NA)
text(-2.0,-0.15,"Increasing\ncomplementarities",cex=1)
polygon(c(-0.6,-0.2,-0.2,-0.6),c(-0.25,-0.25,-0.05,-0.05),col=alpha("black",alph),border=NA)
text(-0.39,-0.15,"Increasing\n comple-\nmentarities",cex=0.7)

runcolorsBOT[1] <- "red"
for (i in 1:ncol(runDsBOT)) {
	lines(ts,runDsBOT[,i],pch=16,lwd=lwd,col=runcolorsBOT[i])
}

###

write.csv(data.frame(ts,runDsBOT[,1]),"run_polymorphic.csv")

