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

nameForSave = alignDataFiles.replace('.txt','')
print nameForSave

x = array('f',[ 0. ])
y = array('f',[ 0. ])
z = array('f',[ 0. ])
score = array('f',[ 0. ])
goodpt = array('f',[ 0. ])

RTree = R.TTree('align','align')
RTree.Branch('x',x,'x/F')
RTree.Branch('y',y,'y/F')
RTree.Branch('z',z,'z/F')
RTree.Branch('score',score,'score/F')
RTree.Branch('goodpt',goodpt,'goodpt/F')

dxdzPlot = R.TProfile('dxdz','dxdz',10,-45000,45001)
dydzPlot = R.TProfile('dydz','dydz',10,-45000,45001)
	
with open(alignDataFiles, "r") as iv:
            txtLines = [line for line in iv]
	    print "Reading headers..."
            idx = [i for i, line in enumerate(txtLines) if "edge" in line][0]
            headers = txtLines[idx].replace('\n','').replace('\r','').split('\t')
	    #print headers

            x_idx = my_index(headers,"x")
            y_idx = my_index(headers,"y")
            z_idx = my_index(headers,"zFocus")
            score_idx = my_index(headers,"score")
            goodpt_idx = my_index(headers,"GoodPt")

	    mmX = -0.98
	    mmY = 0.99

            data = txtLines[idx+1:]
	    print "Reading data..."

	    linenum = 0
            for line in data:
              words = line.replace('\n','').replace('\r','').split('\t')
	      #print words
	      if len(words)>4:
                x[0]=readvar(words,x_idx)
                y[0]=readvar(words,y_idx)
		if linenum==0:
			x0 = x[0]
			y0 = y[0]
                z[0]=readvar(words,z_idx)
                score[0]=readvar(words,score_idx)
                goodpt[0]=readvar(words,goodpt_idx)
		RTree.Fill()
		if score[0] > 2:
			dxdzPlot.Fill(x[0],z[0])
			dydzPlot.Fill(y[0],z[0])
		linenum += 1

	    print "File Reading complete"

if os.path.isfile(nameForSave+".root"): RFile = R.TFile(nameForSave+".root",'UPDATE')
else: RFile = R.TFile(nameForSave+".root",'RECREATE')
RTree.Write("",R.TObject.kOverwrite)

dxdzPlot.Fit("pol1")
dydzPlot.Fit("pol1")
dxdzPlot.Write("",R.TObject.kOverwrite)
dydzPlot.Write("",R.TObject.kOverwrite)

Z1 = R.TCanvas( 'Plots', 'Plots', 0, 0, 400, 400 )
#c1.SetTitle(prefix+"_qvt")
Z1.Divide(2,2)
Z1.Update()

RFit = R.TF1()
flatZ = []
flatZerr = []


ChuckMap = R.TH2F("ChuckMap","Chuck Z height",10,-45000,45001,10,-45000,45001)
ChuckMap.SetMinimum(9600)
ChuckMap.SetMaximum(9710)
ChuckMap.SetOption("colz")
ChuckMap.SetStats(0)
RTree.Draw("x:y>>ChuckMap","z*(score>2)","colz")
ChuckMap.Write("",R.TObject.kOverwrite)

#print "Flatness"
#print flatZ
#print flatZerr

RFile.Close()

rootFile = R.TFile()

#textfile = open(nameForSave+"_summary.txt","w")
#textfile.close()
