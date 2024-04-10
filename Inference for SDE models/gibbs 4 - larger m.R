# m=20
set.seed(1109)
out20 <- gibbs(N=5000, m=20, Vtune=0.12^2*diag(2), data, par, alphaSIR, betaSIR)
post20 <- out20[[1]][-(1:200),]
plot(ts(post20, start=0))
effectiveSize(post20) #96.82968 114.56149 


###Susceptible trace plot###
par(mfrow=c(1,3))
plot(ts(st2[,8], start=0), ylim=c(80,145), xlab="Iteration",
     ylab="Susceptible", main="m=2")
s2 <- st2[,seq(10,14,2)]
for(i in 1:3) {
  lines(ts(s2[,i], start=0))}

plot(ts(st5[,19], start=0), ylim=c(80,145), xlab="Iteration",
     ylab="Susceptible", main="m=5")
s5 <- st5[,c(24,29,34)]
for(i in 1:3) {
  lines(ts(s5[,i], start=0))}

plot(ts(out20[[2]][,71], start=0), ylim=c(80,145), xlab="Iteration",
     ylab="Susceptible", main="m=20")
s20 <- out20[[2]][,c(91,111,131)]
for(i in 1:3) {
  lines(ts(s20[,i], start=0))}

###Acceptance rates###
#m=2
#0.3696739
out2[[4]]
# [1] 0.4514903 0.4614923 0.4060812 0.3626725 
# [5] 0.1524305 0.2258452 0.2448490

#m=5
#0.3626725
out5[[4]]
# [1] 0.1694339 0.1858372 0.1974395 0.1942388
# [5] 0.1058212 0.1064213 0.0520104

#m=20
#0.3034607
out20[[4]]
# [1] 0.03880776 0.04780956 0.06461292 0.06581316
# [5] 0.04060812 0.02440488 0.00320064
