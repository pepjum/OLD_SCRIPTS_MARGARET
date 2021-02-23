#plot espectro


pdf("/home/margaret/data/pepe/espectros.pdf")
plot(x=500:1500, y= espectro_real$Z, type="l")
lines(x=500:1500, y=espectro_simulado$normalized, type="l", col="red")
legend(1200, 0.6, legend=c("real", "simulado"),
       col=c("black", "red"), lty=1:2, cex=0.8)

dev.off()
