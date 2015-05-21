#source("http://bioconductor.org/biocLite.R")
#biocLite("flowCore")
#library("flowCore")
#biocLite("flowViz")
#library("flowViz")
library(prada)
#help.start()
file.name <- "multiColor.fcs"
#file.name <- "one_well.fcs"
#file.name <- "Specimen_001_Tube_001_(1).fcs"

x <- read.fcs(file.name)#FCS(file.name), transformation=FALSE)
#summary(x)
times = exprs(x)[,"Time"]
s = 100; 
a = hist(times,breaks=seq(from=min(times)-s,to=max(times)+s,by=s)) # per second
counts = a$counts
#counts[counts<50]=0
library(quantmod)
library(caTools)
library(signal)


# define a plate region
counts[counts<=15]=0
counts.smooth = runmean(counts,20,align='center')
possible.breakpoints = which(counts.smooth==0)


distance = diff(possible.breakpoints)
reasonable.distances = which(distance>1000&distance<2000)
actual.breakpoints.start = possible.breakpoints[reasonable.distances]
actual.breakpoints.end =   possible.breakpoints[reasonable.distances + 1] 

pdf('plate_identification.pdf')
plot(counts)
lines(counts.smooth,col='orange')
points(which(counts.smooth==0),counts.smooth[counts.smooth==0],pch=16,col='blue')
points(actual.breakpoints.start,counts.smooth[actual.breakpoints.start],pch=16,col='green')
points(actual.breakpoints.end,counts.smooth[actual.breakpoints.start],pch=16,col='red')
dev.off()

for( i in 1:length(actual.breakpoints.start)){
    plate.counts = counts[actual.breakpoints.start[i]:actual.breakpoints.end[i]]
    plate.counts.smooth = runmean(plate.counts,1,align='center')
    possible.plate.breakpoints = which(plate.counts.smooth==0)
    plate.distance = diff(possible.plate.breakpoints)
    reasonable.plate.distances = which(plate.distance>4&plate.distance<15)
    actual.plate.breakpoints.start = possible.plate.breakpoints[reasonable.plate.distances]
    actual.plate.breakpoints.end =   possible.plate.breakpoints[reasonable.plate.distances + 1] 
    
    plot(plate.counts)
    lines(plate.counts.smooth,col='orange')
    points(possible.plate.breakpoints,plate.counts.smooth[possible.plate.breakpoints],pch=16,col='blue')
    points(actual.plate.breakpoints.start,plate.counts.smooth[actual.plate.breakpoints.start],pch=16,col='green')
    points(actual.plate.breakpoints.end,plate.counts.smooth[actual.plate.breakpoints.start],pch=16,col='red')
    
    
    

counts.smooth = runmean(counts,600) # 20 min
counts.smooth.filtered = filter(butter(4, c(1/5500), type="low"),counts.smooth)
plot(counts.smooth,type='l')
lines(counts.smooth.filtered,col='blue')
counts.smooth.filtered.d = runmean(diff(counts.smooth.filtered),2000)
lines(1000*counts.smooth.filtered.d+3,col='red')
peaks = findValleys(counts.smooth.filtered.d)
points(peaks,counts.smooth.filtered[peaks],col='green',pch=16)

sep = runmean(peaks,2)
sep[1]=0
a= matrix(nrow=1,ncol=9)
a[]=25
points(sep,a,pch=16)


plot(counts)
lines(10*runmean(counts,60,align='center'),col='green') # seems to work the best for seperating plates
lines(10*runmean(counts,120,align='center'),col='red')
lines(10*runmean(counts,30,align='center'),col='blue')
lines(10*runmean(counts,180,align='center'),col='brown')
lines(10*runmean(counts,300,align='center'),col='orange')


peaks = findPeaks(counts.smooth)
points(peaks,counts.smooth[peaks],col='red3',pch=16)







counts.smooth.filtered = filter(butter(4, c(1/250,1/50), type="pass"),counts.smooth)
peaks = findPeaks(counts.smooth.filtered)
lines(counts.smooth.filtered,col='blue')
points(peaks,counts.smooth.filtered[peaks],col='green',pch=16)
