par(mfrow=c(1,2))
plot(ts(log(out_h[,2]), start=0, deltat=0.1), col="blue", lwd="1.4",
     xlab="Day", main="Human", ylab="Population (log scale)", ylim=c(0,20))
lines(ts(log(out_h[,3]), start=0, deltat=0.1), col="orange", lwd="1.4")
lines(ts(log(out_h[,4]), start=0, deltat=0.1), col="red", lwd="1.4")
lines(ts(log(out_h[,5]), start=0, deltat=0.1), col="cyan", lwd="1.4")
lines(ts(log(out_h[,6]), start=0, deltat=0.1), col="green3", lwd="1.4")
legend(x=47,y=16.15, lwd=1.5, cex=0.65, box.lty=0.6,
       legend = c("Susceptible","Exposed","Infectious","Quarantine","Recovered"),
       col = c("blue","orange","red","cyan","green3")) 

plot(ts(log(out_r[,2]), start=0, deltat=0.1), col="blue", lwd="1.4",
     xlab="Day", main="Rodent", ylab="Population (log scale)", ylim=c(13,15))
lines(ts(log(out_r[,3]), start=0, deltat=0.1), col="orange", lwd="1.4")
lines(ts(log(out_r[,4]), start=0, deltat=0.1), col="red", lwd="1.4")
legend(x=47,y=14.8, lwd=1.5, cex=0.65, box.lty=0.6,
       legend = c("Susceptible","Exposed","Infectious"),
       col = c("blue","orange","red")) 


par(mfrow=c(1,1))
plot(ts(changeI, start=0, deltat=0.1), lwd="1.4", xlab="Day", main="New infection cases", ylab="Cases")
