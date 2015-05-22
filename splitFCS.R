file.name <- "one_well.fcs"
#file.name <- "multiColor.fcs"

split.fcs = function(counts, smooth, reasonable.time.min, reasonable.time.max, squash, plot=F, title='default title'){
    counts[counts<=squash]=0
    counts.smooth = runmean(counts,smooth,align='center')
    possible.breakpoints = which(counts.smooth==0)
    distance = diff(possible.breakpoints)
    reasonable.distances = which(distance>reasonable.time.min&distance<reasonable.time.max)
    actual.breakpoints.start = possible.breakpoints[reasonable.distances]
    actual.breakpoints.end =   possible.breakpoints[reasonable.distances + 1] 

    if(plot){
    plot(counts,main=paste(title,': ', toString(length(actual.breakpoints.start)), ' detected'))
    lines(counts.smooth,col='orange')
    points(which(counts.smooth==0),counts.smooth[counts.smooth==0],pch=16,col='blue')
    points(actual.breakpoints.start,counts.smooth[actual.breakpoints.start],pch=16,col='green')
    points(actual.breakpoints.end,counts.smooth[actual.breakpoints.start],pch=16,col='red')
    }    
    output = list()
    output$start = actual.breakpoints.start
    output$end = actual.breakpoints.end
    return(output)
}
library(prada)
library(quantmod)
library(caTools)
library(signal)

fcs <- exprs(read.fcs(file.name))#FCS(file.name), transformation=FALSE)
times = fcs[,"Time"]
s = 100; # Convert to seconds
limits = seq(from=min(times)-s,to=max(times)+s,by=s)
time.buckets = hist(times,breaks=limits)
counts = time.buckets$counts # per second

# define a plate region
pdf('test.pdf')
 
plate.splits = split.fcs(counts, 20, 1000, 2000, squash=15, plot=T, title='test')

plate = list()
for( i in 1:length(plate.splits$start)){
    well.splits = split.fcs(counts[plate.splits$start[i]:plate.splits$end[i]], 1, 5, 15, 10, plot=T, title=paste('plate',toString(i)))
    plate[[i]] = list()
    for (j in 1:length(well.splits$start)){
        plate[[i]][[j]]=fcs[which(times>limits[well.splits$start[j]]&times<limits[well.splits$start[j+1]]),]
    }
}



dev.off()
