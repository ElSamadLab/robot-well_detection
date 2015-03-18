import os, FlowCytometryTools
from FlowCytometryTools import FCMeasurement
import matplotlib.pyplot as plt
import numpy as np
import scipy
#from numpy import convolve
import math
import pandas as pd


datadir='/Users/alain/Documents' #directory to facs file(s)
datafile=r'/Users/alain/Documents/Specimen_001_Tube_001_(1).fcs' #the location of the actual FACS file of interest
#datafile1 = os.path.join(datadir, 'Specimen_001_Tube_001_(1).fcs')
sample=FCMeasurement(ID='Test Sample',datafile=datafile)
#print sample.channel_names
#print "The number of events is equal to %s" % sample.shape[0]
Event_Rate=sample.data['Time']
SSC1 =[]
SSC = sample.data['SSC-H']
# SSC1 = SSC.values.tolist()
FITC=sample.data['FITC-H']
xaxis = range(114275)
Diff=[]
def Slope(a,b):
	return float((b-a)/25)
print Event_Rate
Well_Index= []
Wells=[]
Pre_ImageSC = []
ImageSC_Norm=[]
ImageSC = []
#plt.figure(1)
#plt.subplot(1,2,1)
#plt.gca().set_aspect('equal',adjustable='box')
#sample.plot(['FSC-H','SSC-H'], kind = 'scatter',color = 'red')


#plt.subplot(1,2,2)
#plt.gca().set_aspect('equal',adjustable='box')
#sample.plot('FITC-H',bins=1000)
#plt.axis([0.0, 2000.0, 0.0, 100000.0])
#print len(Event_Rate)

Diff = np.diff(Event_Rate)
Sec_Diff = np.diff(Diff)


for x in range(1900,len(Diff)-1):
	if Diff[x] > 60 and Diff[x] > 400:
		if Diff[x+1] < 22 and Diff[x+2] < 22 and Diff[x+3] < 22 and Diff[x+4] < 22 and Diff[x+5] < 22 and Diff[x+6] < 22 and Diff[x+7] < 22 and Diff[x+8] < 22:
			Wells.append(range(Event_Rate[x+1675],Event_Rate[x+4080]))
			Well_Index.append(range(x+1675,x+4080))
print Well_Index
#Histogram= np.array([[23*np.array([[len(SSC[Well_Index[0]]),1]])

for i in range(len(Well_Index)):
	for k in Well_Index[i]:
		ImageSC.append(SSC[k])
	#ImageSC_Norm.append(ImageSC[i]/sum(ImageSC[i]))

#ImageSC_Norm = np.histogram(ImageSC)/sum(ImageSC)
#print "the number of wells is %d " %len(Well_Index)
#print ImageSC[22]
#print ImageSC_Norm[0]
#print SSC[Well_Index[0][0]]
#print SSC[Well_Index[0][1]]
#print len(Well_Index[0])
#for j in range(len(Well_Index[0])):
#	print SSC[range(Well_Index[0][j])]
#print Histogram

plt.figure()
plt.subplot(1,2,1)
for n in range(len(Wells)):
	plt.plot(Event_Rate[Well_Index[n]])
#plt.subplot(1,2,2)
#plt.plot(Event_Rate)

#plt.figure()
#plt.plot(Diff)
#for j in range(len(Well_Index)):
#	for n in range(len(Well_Index[j])):
#		Histogram.append(SSC[Well_Index[j][n]])
#		ImageSC.append([Histogram])
Histogram= np.zeros([23,1])
#for z in range(len(ImageSC[0])):
#print Histogram

#Histogram = np.histogram([ImageSC[0]], bins = 100)

heat = []
for i in range(23):
	dat = ImageSC[i*2405 : (i+1)*2405]
	hist = np.histogram(dat, bins = 100)
	if len(dat) != 2405:
		print "nooo"
	col = hist[0]/float(len(dat))
	heat.append(col)
heat = np.array(heat)
heat = np.transpose(heat)


# plt.figure()
# for k in range(2):#
# 	plt.subplot(3,8,k)
# 	plt.hist(ImageSC[k],bins= 100,normed = True)
# 	#ImageSC_Norm=np.array([[ImageSC[k-1]]])
# 	ImageSC_Norm = np.append(ImageSC_Norm,ImageSC[k])
# 	ImageSC_Norm_TP = np.transpose(ImageSC_Norm)


#ImageSC_Norm1 = np.histogram(ImageSC)

#ImageSC_Norm = np.append(ImageSC_Norm, [[ImageSC[k]]])
plt.figure()
#plt.hist(ImageSC[0],bins= 100)

plt.figure()
plt.pcolor(heat)
#print SSC
#print Event_Rate
#Image_NP = np.array(ImageSC_Norm)
#print ImageSC_Norm[0]
#plt.figure()
#plt.imshow(float(ImageSC[0]))

#print len(ImageSC[4])
#plt.figure()

#for b in range(len(ImageSC[0])):
#	plt.pcolor(np.array([[ImageSC[b]]]))


plt.show()








