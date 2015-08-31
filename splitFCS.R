rm(list=ls())
## This is where you specifiy the fcs file
#file.name <- "msg5ste4_Tube_001.fcs"
file.name <- commandArgs(TRUE)[1]
# print(test)
# q('yes')
#file.name <- "dig1ste4_Tube_001.fcs"
#file.name <- "one_well.fcs"
#file.name <- "multiColor.fcs"
#file.name <- "ste4gpa1_di_robot_t2.fcs"
#file.name <- "Andres_FACS_Timecourse.fcs"
#file.name <- "patrick.fcs"
#file.name <- "20150618robot.fcs"

## This is the function to split up plates. I tried using it to split up wells as well, but it didnt works so well.

split.fcs = function(first.index, counts, smooth, reasonable.time.min, reasonable.time.max, squash, plot=F, title='default title'){
    # first.index: starting point to do the splitting, must be before the split actually takes place
	# counts: this is the event rate (usually per second but not nessisarily)
	# smooth: how many points will get averaged in the smoothing
	# reasonable.time.min, reasonable.time.max: what is considered a reasonable time between the events: in the case of a plate a reasonable time is between 1000 seconds and 2000 seconds
	# squash: low event rates need to be zero for this to work, this number should change based on the dencity of cells you are using, I found that 15 works quite well with event rates of about 200/sec this is going to take all the event rates <15 and make them zero so for higher ODs you may need a higher number
	# plot: do you want the function to plot out what it did, very useful for first time users
	# title: the title on the plot you create.
	counts[counts<=squash]=0
    counts.smooth = runmean(counts,smooth,align='center')
    possible.breakpoints = which(counts.smooth==0)
    distance = diff(possible.breakpoints)
    reasonable.distances = which(distance>=reasonable.time.min&distance<=reasonable.time.max)
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
    output$start = actual.breakpoints.start + first.index - 1
    output$end = actual.breakpoints.end + first.index - 1
    return(output)
}

## These are the libraries that you need to have loaded. you can install them with install.packages('quantmod') etc. for the flowcore you need to install that from bioconductor. make sure you have version 1.35.5 or later (bioconductor 3.2 will get you there)
library(flowCore)
library(quantmod)
library(caTools)
library(signal)
library(ggplot2)

# read fcs file
fcs <- exprs(read.FCS(file.name, transformation=FALSE, truncate_max_range = FALSE))
# extract time of events
times = fcs[,"Time"]
s = 100; # Conversion into seconds 
limits = seq(from=min(times)-s,to=max(times)+s,by=s)
pdf('plots.pdf')
time.buckets = hist(times,breaks=limits)
counts = time.buckets$counts # per second

# define a plate region
#pdf('test.pdf') only if you want to save the image, remember to add dev.off() at the end
 
plate.splits = split.fcs(0, counts, 20, 1000, 2000, squash=15, plot=T, title='test')


## this is the old code that used split.fcs on wells as well as plates. it needs tweakeing and doing a simple breakdown works well enough
# plate = list()
# for( i in 1:length(plate.splits$start)){
#     well.splits = split.fcs(first.index = plate.splits$start[i], counts[plate.splits$start[i]:plate.splits$end[i]], smooth=1, reasonable.time.min=5, reasonable.time.max=14, 25, plot=T, title=paste('plate',toString(i)))
#     plate[[i]] = list()
#     for (j in 1:length(well.splits$start))
#     {
#         plate[[i]][[j]]=fcs[which(times>limits[well.splits$start[j]]&times<limits[well.splits$start[j+1]]),]
#         fname = paste('timepoint_',as.character(i),'_well_', as.character(j), '_', file.name, sep='')
#         write.FCS(flowFrame(plate[[i]][[j]]),file = fname)
#     }
# }
# this is the dataframe that I will exort at the end of the forloop
boundries = data.frame(plate_num= integer(0), well_num= integer(0), start_time = numeric(0), end_time = numeric(0))
counter = 1
plate.events = numeric()
for (plate.num in 1:length(plate.splits$start)) {
	# break down a single plate
	plate = fcs[fcs[,'Time']>limits[plate.splits$start[plate.num]] & fcs[,'Time']<limits[plate.splits$end[plate.num]],]
	plate = data.frame(plate[,c("SSC-H","mCherry-H","FITC-H","Time")])
	plate.events[plate.num] =  dim(plate)[1]
	# define well boundries simply on breaking the time down into 96 pieces
	wells.bound = seq(from=min(plate[,"Time"]),to=max(plate[,"Time"]),length.out=96*1+1)
	for (well.num in 1:96) {
		# now break down the well into 10 time segments, so that we can get rid of contamination from previous and next well
		tenths = seq(from=wells.bound[well.num],to=wells.bound[well.num+1],length.out=11)
		# take the middle 6 tenths (leave out 2 on each side)
		boundries[counter,] = c(plate.num, well.num, tenths[3],tenths[9])
		counter = counter + 1
	}
	## Print out the wells for sanity sake, only really makes sense in the validation script
	# d = ggplot(plate,aes(mCherry.H,FITC.H))
	# d + geom_point() + facet_wrap(~Time,ncol=12)
	
}
y = diff(log(plate.events))#+log(530/30)
plate.time = (1:length(plate.events))*20
plot(y, type='l')
dev.off()

outfile = paste(strsplit(file.name, "\\.")[[1]][1],'well_splits.csv',sep='_')

write.table(boundries, outfile,row.names=FALSE,sep=',')


#dev.off()








