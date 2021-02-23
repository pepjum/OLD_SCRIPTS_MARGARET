
pdf('wavelet.pdf', width=20, height=14, onefile=FALSE)
m<-matrix(terrain.colors(20), ncol=1)
m<-matrix(m[nrow(m):1, ], ncol=1)
legend_image <- as.raster(m)
layout(matrix(1:2,ncol=2), width = c(3.5,1),height = c(1,2))
plot(CWT_signal, col=terrain.colors(length(CWT_signal)))
plot(c(0,3),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=1.5, y = seq(0,1,l=5), labels = seq(round(min(as.matrix(CWT_signal))),round(max(as.matrix(CWT_signal))),l=5))
rasterImage(legend_image, 0, 0, 1,1)

dev.off()


pdf("/home/nostromo/data/pepe/01_CHIPSEQ/EGR1/wavelet_series.pdf",width=20, height=14, onefile=FALSE)
plot(CWT_signal_plotting, col=terrain.colors(length(CWT_signal_plotting)), series=TRUE)
dev.off()


pdf("/home/nostromo/data/pepe/01_CHIPSEQ/EGR1/zero_crossing_lines.pdf",width=20, height=14, onefile=FALSE)

for (i in 1:length(ZCL_signal_plotting$zerolineTrace)){
	if (unlist(ZCL_signal_plotting$zeroline)[i]>0) {
		print(i)
		x<-c(unlist(ZCL_signal_plotting$zerolineTrace[[i]]))
		y<-c(seq(1:length(x)))
		if (i == which(unlist(ZCL_signal_plotting$zeroline) != 0)[1]) {

            plot(x,y, type="l", xlim = c(1, 14000), ylim = c(1, 30))
        }
         else {
				lines(x,y, type="l")
			}
		}
}

dev.off()
