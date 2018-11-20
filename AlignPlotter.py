#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Alignment Measurement Plotter
#3/21/18

#############################
#							#
#		  Imports			#
#							#
#############################

import ROOT as R
import os.path
import sys
from array import array
#import matplotlib.pyplot as plt
#import matplotlib.ticker as mtick

#############################
#							#
#		  Globals			#
#							#
#############################


def my_index(L, obj):
     try:
             return L.index(obj)
     except ValueError:
             return -1

def readvar(L, obj):
     if obj == -1:
             return -1
     else:
             return float(L[obj])


alignDataFiles = sys.argv[1]
PSS_fid = sys.argv[2]
PSP_fid = sys.argv[3]

nameForSave = alignDataFiles.replace('.txt','')
print nameForSave
if PSS_fid == "0":
	pss_xoff = 420
	pss_yoff = -185
else:
	pss_xoff = 525
	pss_yoff = -110
if PSP_fid == "0":
	psp_xoff = 420
	psp_yoff = -185
else:
	psp_xoff = 525
	psp_yoff = -110

print "Offsets", pss_xoff, pss_yoff, psp_xoff, psp_yoff

x = array('f',[ 0. ])
y = array('f',[ 0. ])
z = array('f',[ 0. ])
ImgX = array('f',[ 0. ])
ImgY = array('f',[ 0. ])
score = array('f',[ 0. ])
edge = array('i',[ 0 ])
fidX = array('f',[ 0. ])
fidY = array('f',[ 0. ])
goodpt = array('f',[ 0. ])

RTree = R.TTree('align','align')
RTree.Branch('x',x,'x/F')
RTree.Branch('y',y,'y/F')
RTree.Branch('z',z,'z/F')
RTree.Branch('ImgX',ImgX,'ImgX/F')
RTree.Branch('ImgY',ImgY,'ImgY/F')
RTree.Branch('fidX',fidX,'fidX/F')
RTree.Branch('fidY',fidY,'fidY/F')
RTree.Branch('score',score,'score/F')
RTree.Branch('edge',edge,'edge/I')
RTree.Branch('goodpt',goodpt,'goodpt/F')
	
with open(alignDataFiles, "r") as iv:
            txtLines = [line for line in iv]
	    print "Reading headers..."
            idx = [i for i, line in enumerate(txtLines) if "edge" in line][0]
            headers = txtLines[idx].replace('\n','').replace('\r','').split('\t')
	    #print headers

            x_idx = my_index(headers,"x")
            y_idx = my_index(headers,"y")
            z_idx = my_index(headers,"zFocus")
            ImgX_idx = my_index(headers,"ImgX")
            ImgY_idx = my_index(headers,"ImgY")
            score_idx = my_index(headers,"score")
            edge_idx = my_index(headers,"edge")
            goodpt_idx = my_index(headers,"GoodPt")

            xlist = [[],[],[],[]]
            ylist = [[],[],[],[]]
            zlist = [[],[],[],[]]
	    crnrs = [[],[],[],[]]
	    crnrxy = [[],[],[],[]]
            #ImgX = []
            #ImgY = []
            #score = []
            #edge = []

	    mmX = 0.98
	    mmY = 0.99
	    dzdx = -.0009
	    dzdy = .0005

            data = txtLines[idx+1:]
	    print "Reading data..."

	    lastcorner = [0,0,0,0]
	    lastcornerxy = [0,0]
	    linenum = 0
	    lastedge = -1
	    lastxy = -1
            for line in data:
              words = line.replace('\n','').replace('\r','').split('\t')
	      #print words
	      if len(words)>4:
                edge[0]=int(readvar(words,edge_idx))
                x[0]=readvar(words,x_idx)
                y[0]=readvar(words,y_idx)
		if linenum==0:
			x0 = x[0]
			y0 = y[0]
                z[0]=10160-readvar(words,z_idx)+(x[0]-x0)*dzdx + (y[0]-y0)*dzdy
		#print "diff: ", (x[0]-x0)*dzdx, (y[0]-y0)*dzdy
                ImgX[0]=readvar(words,ImgX_idx)
                ImgY[0]=readvar(words,ImgY_idx)
                score[0]=readvar(words,score_idx)
                goodpt[0]=readvar(words,goodpt_idx)

		fidX[0] = x[0] + ImgX[0]*mmX
		fidY[0] = y[0] + ImgY[0]*mmY
		if (edge[0] == lastedge): 
		    lastcornerxy = [x[0],y[0]]
		    if (goodpt[0]>0):
			lastcorner = [x[0],y[0],z[0],fidX[0],fidY[0]]

		#print edge[0], lastedge, x[0], y[0], lastcorner
		if (edge[0] != lastxy):
			crnrxy[edge[0]].append([x[0],y[0]])
			if (lastxy>=0): crnrxy[lastxy].append(lastcornerxy)
		if (edge[0] != lastedge) and (goodpt[0]>0): 
			crnrs[edge[0]].append([x[0],y[0],z[0],fidX[0],fidY[0]])
			if lastedge>=0: crnrs[lastedge].append(lastcorner)
		if (goodpt[0]>0):
			lastedge = edge[0]
		lastxy = edge[0]

		xlist[edge[0]].append(fidX[0])
		ylist[edge[0]].append(fidY[0])
		zlist[edge[0]].append(z[0])
		RTree.Fill()
		linenum += 1

	    crnrs[lastedge].append(lastcorner)
	    crnrxy[lastxy].append(lastcornerxy)
	    print "File Reading complete"
	    #print "Corners Fiducials", crnrs
	    #print "Corners XY", crnrxy



if os.path.isfile(nameForSave+".root"): RFile = R.TFile(nameForSave+".root",'UPDATE')
else: RFile = R.TFile(nameForSave+".root",'RECREATE')
RTree.Write("",R.TObject.kOverwrite)

Z1 = R.TCanvas( 'Plots', 'Plots', 0, 0, 400, 400 )
#c1.SetTitle(prefix+"_qvt")
Z1.Divide(2,2)
Z1.Update()

RFit = R.TF1()
flatZ = []
flatZerr = []
edges = []
edgeserr = []

for i in xrange(4):    
	RTree.Draw("z:x","edge==" + str(i) + " && score>0.4")
	ZGr = R.TGraph(RTree.GetSelectedRows(), RTree.GetV2(), RTree.GetV1())
	ZGr.SetName("Z Edge " + str(i))
	if (i > 0) and (i < 3): sensor = "PSP"
	else: sensor = "PSS"
	if i < 2: color = 1
	else: color = 2
	ZGr.SetMarkerColor(color)
	ZGr.SetLineColor(color)
	ZGr.SetTitle(sensor + " Edge " + str(i))
	ZGr.GetXaxis().SetTitle("x (um)")
	ZGr.GetYaxis().SetTitle("z (um)")
	#Z1.cd(i+1)
	ZGr.Draw("*L")
	#Z1.Update()
	Z1.SaveAs(nameForSave+"_edge"+str(i)+".jpg")
	#RFit = ZGr.Fit("pol1")
	#flatZ.append(R.gROOT.GetFunction("pol1").GetParameter(1))
	#flatZerr.append(R.gROOT.GetFunction("pol1").GetParError(1))
	ZGr.Write("",R.TObject.kOverwrite)

	RTree.Draw("fidY:fidX","edge==" + str(i) + " && goodpt>0")
	FGr = R.TGraph(RTree.GetSelectedRows(), RTree.GetV2(), RTree.GetV1())
	FGr.SetName("Fid Edge " + str(i))
	FGr.SetTitle("Fiducials")
	FGr.GetXaxis().SetTitle("x (um)")
	FGr.GetYaxis().SetTitle("y (um)")
	FGr.Draw("*L")
	FGr.Fit("pol1")
	edges.append(R.gROOT.GetFunction("pol1").GetParameter(1))
	edgeserr.append(R.gROOT.GetFunction("pol1").GetParError(1))
	FGr.Write("",R.TObject.kOverwrite)

#c1.Print(nameForSave+"_zProfile.pdf")
Z1.Update()
Z1.SaveAs(nameForSave+"_zProfile.jpg")

ChuckMap = R.TH2F("ChuckMap","Chuck Z height",10,121500,211500,10,58500,148500)
RTree.Draw("x:y>>ChuckMap","z","colz")
ChuckMap.Write("",R.TObject.kOverwrite)

#print "Flatness"
#print flatZ
#print flatZerr

print "Corner Fiducials" 
for i in xrange(4): 
	print crnrs[i]

corners = [[[0,0],[0,0]],[[0,0],[0,0]],[[0,0],[0,0]],[[0,0],[0,0]]]

print "Corners (edges)"
corners[0][0] = [crnrs[0][0][3]+crnrxy[0][0][0]-crnrs[0][0][0]-pss_xoff, crnrs[0][0][4]+crnrxy[0][0][1]-crnrs[0][0][1]+pss_yoff]
corners[1][0] = [crnrs[1][0][3]+crnrxy[1][0][0]-crnrs[1][0][0]-psp_xoff, crnrs[1][0][4]+crnrxy[1][0][1]-crnrs[1][0][1]+psp_yoff]
corners[2][0] = [crnrs[2][0][3]+crnrxy[2][0][0]-crnrs[2][0][0]-psp_xoff, crnrs[2][0][4]+crnrxy[2][0][1]-crnrs[2][0][1]-psp_yoff]
corners[3][0] = [crnrs[3][0][3]+crnrxy[3][0][0]-crnrs[3][0][0]-pss_xoff, crnrs[3][0][4]+crnrxy[3][0][1]-crnrs[3][0][1]-pss_yoff]
corners[0][1] = [crnrs[0][1][3]+crnrxy[0][1][0]-crnrs[0][1][0]+pss_xoff, crnrs[0][1][4]+crnrxy[0][1][1]-crnrs[0][1][1]+pss_yoff]
corners[1][1] = [crnrs[1][1][3]+crnrxy[1][1][0]-crnrs[1][1][0]+psp_xoff, crnrs[1][1][4]+crnrxy[1][1][1]-crnrs[1][1][1]+psp_yoff]
corners[2][1] = [crnrs[2][1][3]+crnrxy[2][1][0]-crnrs[2][1][0]+psp_xoff, crnrs[2][1][4]+crnrxy[2][1][1]-crnrs[2][1][1]-psp_yoff]
corners[3][1] = [crnrs[3][1][3]+crnrxy[3][1][0]-crnrs[3][1][0]+pss_xoff, crnrs[3][1][4]+crnrxy[3][1][1]-crnrs[3][1][1]-pss_yoff]

distS = []
distL = []
for i in xrange(4): 
	print corners[i]
	if i<3: j = i+1
	else: j = 0
	distS.append(((corners[i][0][0]-corners[i][1][0])**2+(corners[i][0][1]-corners[i][1][1])**2)**0.5)

for i in xrange(2): 
	if i == 0:  j = 3
	else: j = 1
	distL.append(((corners[2*i][0][0]-corners[j][0][0])**2+(corners[2*i][0][1]-corners[j][0][1])**2)**0.5)
	distL.append(((corners[2*i][1][0]-corners[j][1][0])**2+(corners[2*i][1][1]-corners[j][1][1])**2)**0.5)

deltaZ = [crnrs[0][0][2]-crnrs[1][0][2], crnrs[0][1][2]-crnrs[1][1][2], crnrs[3][0][2]-crnrs[2][0][2], crnrs[3][1][2]-crnrs[2][1][2]]
overhang = [corners[0][0][1]-corners[1][0][1], corners[0][1][1]-corners[1][1][1], corners[3][0][1]-corners[2][0][1], corners[3][1][1]-corners[2][1][1]]
shift = [corners[0][0][0]-corners[1][0][0], corners[0][1][0]-corners[1][1][0], corners[3][0][0]-corners[2][0][0], corners[3][1][0]-corners[2][1][0]]
print "Overhangs: ", overhang
print "Shifts: ", shift
print "Delta Z:", deltaZ

print "Dimensions:"
print "\t\t",distS[1]
print "\t\t",distS[0]
print "PSS: ",distL[0],"\t",distL[1]
print "PSP: ",distL[2],"\t",distL[3]
print "\t\t",distS[3]
print "\t\t",distS[2]

print "\n","Offsets:"
print "\t", overhang[0], overhang[1], "\t"
print shift[0], "\t\t", shift[1] 
print shift[2], "\t\t", shift[3] 
print "\t", overhang[2], overhang[3], "\t"

print "\nSlope of Edges"
print edges[0], edges[1], edges[2], edges[3]
print "Error in Fit:"
print edgeserr

zPSS = R.TMultiGraph()
zPSS.SetName("zPSS")
zPSS.SetTitle("PSS Z Profile")
#zPSS.GetXaxis().SetTitle("x Pos (um)")
#zPSS.GetYaxis().SetTitle("z Pos (um)")
zPSP = R.TMultiGraph()
zPSP.SetName("zPSP")
zPSP.SetTitle("PSP Z Profile")
#zPSP.GetXaxis().SetTitle("x Pos (um)")
#zPSP.GetYaxis().SetTitle("z Pos (um)")

zPSS.Add(RFile.Get("Z Edge 0"))
zPSS.Add(RFile.Get("Z Edge 3"))

zPSP.Add(RFile.Get("Z Edge 1"))
zPSP.Add(RFile.Get("Z Edge 2"))



TopEdge = (edges[0]-edges[1])*1000000
BotEdge = (edges[3]-edges[2])*1000000
TopTitle = "Top Edge Parallelism " + str(int(TopEdge)) + " urad"
BotTitle = "Bottom Edge Parallelism " + str(int(BotEdge)) + " urad"
print TopTitle, BotTitle
GrTop = R.TMultiGraph()
GrTop.SetName("Top Edges")
GrTop.SetTitle("Top Edge Parallelism " + str(int(TopEdge)) + " urad")
#GrTop.GetXaxis().SetTitle("x Pos (um)")
#GrTop.GetYaxis().SetTitle("y Pos (um)")
GrBot = R.TMultiGraph()
GrBot.SetName("Bottom Edges")
GrBot.SetTitle("Bottom Edge Parallelism " + str(int(BotEdge)) + " urad")
#GrBot.GetXaxis().SetTitle("x Pos (um)")
#GrBot.GetYaxis().SetTitle("y Pos (um)")

c3 = R.TCanvas( 'Fiducials', 'Fiducials', 0, 0, 1200, 1200 )
c3.Divide(1,2)
c3.cd(1)
GrTop.Add(RFile.Get("Fid Edge 0"))
GrTop.Add(RFile.Get("Fid Edge 1"))
GrTop.Draw()
c3.SaveAs(nameForSave+"_Edges.jpg")

c3.cd(2)
GrBot.Add(RFile.Get("Fid Edge 2"))
GrBot.Add(RFile.Get("Fid Edge 3"))
GrBot.Draw()
#c3.Update()

zPSS.Write("",R.TObject.kOverwrite)
zPSP.Write("",R.TObject.kOverwrite)
GrTop.Write("",R.TObject.kOverwrite)
GrBot.Write("",R.TObject.kOverwrite)
Z1.Write("",R.TObject.kOverwrite)
#c2.Write("zProfile",R.TObject.kOverwrite)
#c3.Write("",R.TObject.kOverwrite)

#c2 = R.TCanvas( 'Z Profile', 'Z Profile', 0, 0, 400, 400 )
#c2 = R.TCanvas()
#R.gStyle.SetOptStat(0)
#c2.Clear()
#c2.Divide(1,2)
#c2.Update()
#Z1.Update()
Z1.cd(1) 
zPSS.Draw()
#RFile.Get("zPSS").Draw()
#c2.cd(2)
#zPSP.Draw()
#Z1.Update()
Z1.SaveAs(nameForSave+"_zProfile.jpg")

RFile.Close()

#for i in xrange(4):
#	end = len(xlist[i])
#	print "Edge " + str(i) + ": (Num points: " + str(end) + ")"
#	print xlist[i][0],ylist[i][0],zlist[i][0]
#	print xlist[i][end-1],ylist[i][end-1],zlist[i][end-1]

rootFile = R.TFile()

textfile = open(nameForSave+"_summary.txt","w")
textfile.write ("Parallelism:\n")
textfile.write ("Top Edge: " + str(TopEdge) + " urad\n")
textfile.write ("Bottom Edge: " + str(BotEdge) + " urad\n")
textfile.write ("\n")
textfile.write ("Slopes:\n")
for i in edges: textfile.write(str(i)+"\n")
textfile.write ("\n")
textfile.write("Dimensions:\n")
textfile.write("\t\t"+str(distS[1])+"\n")
textfile.write("\t\t"+str(distS[0])+"\n")
textfile.write("PSS: "+str(distL[0])+"\t"+str(distL[1])+"\n")
textfile.write("PSP: "+str(distL[2])+"\t"+str(distL[3])+"\n")
textfile.write("\t\t"+str(distS[3])+"\n")
textfile.write("\t\t"+str(distS[2])+"\n")

textfile.write("\nOffsets:\n")
textfile.write("\t"+str(overhang[0])+"\t"+str(overhang[1])+ "\t"+"\n")
textfile.write(str(shift[0])+ "\t\t\t" + str(shift[1])+"\n") 
textfile.write(str(shift[2])+ "\t\t\t" + str(shift[3])+"\n")
textfile.write("\t"+str(overhang[2]) +"\t"+ str(overhang[3])+"\t\n\n")

textfile.write ("Corner Fiducials:\n")
for i in crnrs: textfile.write(str(i)+"\n")
textfile.write ("\n")
textfile.write ("PSS Fiducial Offset: (" + str(pss_xoff) + ", " + str(pss_yoff) + ")\n")
textfile.write ("PSP Fiducial Offset: (" + str(psp_xoff) + ", " + str(psp_yoff) + ")\n")
textfile.write ("\n")
textfile.write ("Corners Silicon:\n")
for i in corners: textfile.write(str(i)+"\n")
textfile.write ("\n")
#textfile.write ("Overhangs:\n")
#for i in overhang: textfile.write(str(i)+"\n")
#textfile.write ("\n")
#textfile.write ("Shift:\n")
#for i in shift: textfile.write(str(i)+"\n")
#textfile.write ("\n")
textfile.write ("Delta Z:\n")
for i in deltaZ: textfile.write(str(i)+"\n")
textfile.close()
