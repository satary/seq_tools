# -*- coding: utf-8 -*-
"""
This library plots profiles on top of sequences.
Using combination of R and python.
It can generate visual representation of sequences by itself.

"""

from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys
from Bio import AlignIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import csv
import collections
from Bio import Entrez
import  pickle
from Bio import SeqIO

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

import numpy as np
import pandas as pd
from io import StringIO

#from hist_ss import get_hist_ss
#from hist_ss import get_hist_ss_in_aln, get_hist_ss_in_aln_for_shade
from Bio.Align.AlignInfo import SummaryInfo
from seq_tools import aln_tools
from seq_tools.shade_aln import shade_aln2png


from seq_tools import CONFIG


def plot_mat4seq(filename='default',data=[],seq1=[],seq1lab='Chain 1',seq1offset=0,features1=[],seq2=[],seq2lab='Chain 2',seq2offset=0,features2=[],title=''):
	"""
	will plot a 2D matrix for intersection of two sequences
	data is a dataframe with three columns Resid1, Resid2, Value.
	Resids - 0 based numbering, or offset specified
	you have to have zeros for non interacting residue(!)
	"""
	tempdf=os.path.join(CONFIG.TEMP_DIR,'temp.csv')
	temppng1=os.path.join(CONFIG.TEMP_DIR,'tempprofseq1.png')
	temppng2=os.path.join(CONFIG.TEMP_DIR,'tempprofseq2.png')

	lenseq1=len(seq1[0])
	lenseq2=len(seq2[0])
	# print profile
	#convert matrix to a dataframe
	# df=pd.DataFrame(np.array([profile,range(len(profile))]).T,columns=[axis,'Resid'])
	data.to_csv(tempdf,index=False)
	shade_aln2png(seq1,filename=temppng1,shading_modes=['charge_functional'], legend=False, features=features1,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
	shade_aln2png(seq2,filename=temppng2,shading_modes=['charge_functional'], legend=False, features=features2,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False,rotate=True)


	#let's write an R-script
	a=open(os.path.join(CONFIG.TEMP_DIR,'mat4seq.r'),'w')

	a.write(r"""
	library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img1 <- readPNG("%s")
img2 <- readPNG("%s")


"""%(tempdf,temppng1,temppng2))

	a.write("""
	seqimg1 <- rasterGrob(img1, interpolate=TRUE,width=1)
	seqimg2 <- rasterGrob(img2, interpolate=TRUE,height=1)


theme_set(theme_bw(base_size=24)+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=df,aes(x=Resid1+1-%d,y=Resid2+1-%d,color=Value))+

xlab("%s")+
ylab("%s")+
xlim(%f,%f)+
ylim(%f,%f)+
ggtitle("%s")+
# geom_hline(yintercept = c_c$PROT2_resid, colour="green", linetype = "longdash",size=0.5)+
# geom_vline(xintercept = h2azimp, colour="green", linetype = "longdash",size=0.5)+
geom_point(size=5)+scale_colour_gradient(low="blue", high="red",name="Value")+

annotation_custom(seqimg1, ymin=%f, ymax=0, xmin=0.5,xmax=%f)+
annotation_custom(seqimg2, xmin=%f,xmax=0,ymin=0.5, ymax=%f)
"""%(seq1offset,seq2offset,seq1lab,seq2lab,\
	-lenseq1*0.05,lenseq1*1.01,-lenseq2*0.05,lenseq2*1.01,\
	title[-lenseq1:len(title)],\
	-lenseq2*0.05,lenseq1+0.5,\
	-lenseq1*0.05,lenseq2+0.5))

	a.write("""
	ggsave("%s",plot=a,height=%f,width=%f)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),12./60.*lenseq2+3,12./60.*lenseq1+3))

	a.close()
	os.system('R --vanilla --slave < '+os.path.join(CONFIG.TEMP_DIR,'mat4seq.r'))
	os.system('rm Rplots.pdf')
	os.remove(tempdf)
	os.remove(temppng1)
	os.remove(temppng2)


def plot_manyprof4seq(filename='default',profiledf=[[]],seqmsa=[],features=[],axis='X',title='',offset=0.,yoffset=0.,funcgroups=None,seqontop=False,ruler=False,htune=1.0,ltune=1.0,wtune=1.0,oy=0.0,type='bar',bwidth=0.2,xbreaksby=2,base_size=24,psize=1.0,lsize=1.0,seqmargin=None,spacing=1,dropnan=False,vline=None,vlsize=0.0,xlab=None,scale_color=None,scale_x=None,scale_y=None,dpi=300):
	"""
	will plot several profile for every position in a sequence or small msa
	profiledf now is a pandas dataframe with columns to be plotted.
	oy - y lower lim
	type='bar' or 'point'
	profiledf - is a dataframe
	spacing just divides all indexes (1-based) in array by spacing
	"""
	tempdf=TEMP_DIR+'/temp.csv'
	temppng=TEMP_DIR+'/tempprofseq.png'
	# print profile
	#convert profile to a dataframe
	ldf=len(profiledf)
	mdf=max(list(profiledf.max()))
	df=pd.concat([profiledf,pd.DataFrame(np.arange(1.0,float(ldf)+1.0).T/spacing,columns=['Resid'])],axis=1)
	df.to_csv(tempdf,index=False)
	shade_aln2png(seqmsa,filename=temppng,shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=ruler,show_seq_names=False,show_seq_length=False,funcgroups=funcgroups,margins=seqmargin)

	if('expression' not in axis):
		axis='"'+axis+'"'
	#let's write an R-script
	a=open(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r','w'))

	a.write(r"""
	library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img <- readPNG("%s")

"""%(tempdf,temppng))

	a.write("""
	seqimg <- rasterGrob(img, interpolate=TRUE,width=1)

theme_set(theme_bw(base_size=%d)+theme(plot.margin=unit(c(0,0,0,0),"mm"),panel.border =element_rect(linetype = "dashed", colour = "white")))

ndf <- melt(df, id="Resid")
if(%s){
ndf=ndf[!is.na(ndf$value),]}

a<-ggplot(data=ndf,aes(x=Resid+%f,y=value+%f,fill=variable,color=variable))+
%s geom_vline(xintercept = %s,size=%f)+
geom_%s(stat='identity'%s)+#scale_y_continuous(limits=c(-5,6),breaks=seq(0,6),labels=seq(0,6,by=1))+
%s geom_line(stat='identity',position='identity',size=%f)+
# scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_color_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
xlab('%s')+
%s+
ylim(%f,%f)+
%s+
# geom_point(data=h3data_amore,aes(color=color),fill='red',size=2)+
ylab(%s)+
ggtitle("%s")+
scale_x_continuous(limits=c(%d,%d),breaks = seq(%d, %d, by = %d))+
%s+
annotation_custom(seqimg, ymin=%f, ymax=%f, xmin=%f,xmax=%f)"""%(base_size,'TRUE' if dropnan else 'FALSE',\
	offset,yoffset, '' if vline!=None else '#', vline if vline!=None else '0.0',vlsize,type,',position=\'dodge\',width=%f'%bwidth if type == 'bar' else ',\
	size=%f'%psize,' ' if (type == 'point') or (type=='line') else '#',lsize,\
	xlab if xlab else 'Resid',\
	scale_color if scale_color else '#',\
	(-mdf*1.01*ltune+oy)+yoffset,mdf*1.01+yoffset,\
	scale_y if scale_y else '#',\
	axis,title,offset+0.5,offset+ldf/spacing+0.5,\
	offset-(offset%xbreaksby),offset+ldf+xbreaksby,xbreaksby,\
	scale_x if scale_x else '#',\
	-mdf*ltune+oy+yoffset,oy+yoffset,offset+0.5,ldf/spacing+0.5+offset))

	if(seqontop):
		a.write("""+geom_text(aes(x=Resid+%f-0.25,label=seq),hjust=0, vjust=%f)"""%(offset,-mdf*0.5))
	else:
		a.write("""
			""")
	a.write("""
	ggsave("%s",plot=a,height=%f,width=%f,dpi=%d)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),4.0*htune,12./60.*ldf*wtune,dpi))

	a.close()
	os.system('R --vanilla --slave < '+os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.remove(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.remove(tempdf)
	os.remove(temppng)





def plot_prof4seq(filename='default',profile=[],seqmsa=[],features=[],axis='X',title='',offset=0.,funcgroups=None,seqontop=False,ruler=False,htune=1.0,ltune=1.0,oy=0.0,type='bar',bwidth=0.7,fontsize=18):
	"""
	will plot a profile for every position in a sequence or small msa
	profile is a list of values.
	oy - y lower lim
	type='bar' or 'point'
	"""
	tempdf=os.path.join(CONFIG.TEMP_DIR,'temp.csv')
	temppng=os.path.join(CONFIG.TEMP_DIR,'tempprofseq.png')
	# print profile
	#convert profile to a dataframe
	if(seqontop):
		df=pd.DataFrame(np.array([profile,range(len(profile)),list(str(seqmsa[0].seq))]).T,columns=['axis','Resid','seq'])
	else:
		df=pd.DataFrame(np.array([profile,range(len(profile))]).T,columns=['axis','Resid'])
	df.to_csv(tempdf,index=False)
	shade_aln2png(seqmsa,filename=temppng,shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=ruler,show_seq_names=False,show_seq_length=False,funcgroups=funcgroups)

	#let's write an R-script
	a=open(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'),'w')

	a.write(r"""
	library(ggplot2)
#library(fitdistrplus)
#library(gtools)
library(png)
library(grid)
#library(reshape2)
#library(gridExtra)
#library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img <- readPNG("%s")

"""%(tempdf,temppng))

	a.write("""
	seqimg <- rasterGrob(img, interpolate=TRUE,width=1)

theme_set(theme_bw(base_size=%d)+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))


a<-ggplot(data=df,aes(x=Resid+1+%f,y=axis))+
geom_%s(stat='identity',position='identity',width=%f)+#scale_y_continuous(limits=c(-5,6),breaks=seq(0,6),labels=seq(0,6,by=1))+
%s geom_line(stat='identity',position='identity')+
# scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_color_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
xlab('Resid')+
ylim(%f,%f)+
# geom_point(data=h3data_amore,aes(color=color),fill='red',size=2)+
ylab("%s")+
ggtitle("%s")+
annotation_custom(seqimg, ymin=%f, ymax=%f, xmin=0.5,xmax=%f)"""%(fontsize,offset,type,bwidth,' ' if type == 'point' else '#',(-max(profile)*1.01*ltune+oy),max(profile)*1.01,axis,title[-len(profile):len(title)],-max(profile)*ltune+oy,oy,len(profile)+0.5))
	if(seqontop):
		a.write("""+geom_text(aes(x=Resid+1+%f-0.25,label=seq),hjust=0, vjust=%f)"""%(offset,-max(profile)*0.5))
	else:
		a.write("""
			""")
	a.write("""
	ggsave("%s",plot=a,height=%f,width=%f)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),4.0*htune,12./60.*len(profile)))

	a.close()
	os.system('R --vanilla --slave < %s'%os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.system('rm Rplots.pdf')
	os.remove(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.remove(tempdf)
	os.remove(temppng)


def plot_2prof4seq(filename='default',profile=[],profile2=[],seqmsa=[],features=[],axis='X',axis2='X2',title='',offset=0.,funcgroups=None,ruler=False,htune=1.0,ltune=1.0):
	"""
	thes same as plot_prof4seq but two profiles
	profiles should be of the same length!
	"""
	tempdf=os.path.join(CONFIG.TEMP_DIR,'temp.csv')
	temppng=os.path.join(CONFIG.TEMP_DIR,'tempprofseq.png')
	# print profile
	#convert profile to a dataframe
	df=pd.DataFrame(np.array([profile,profile2,range(len(profile))]).T,columns=[axis,axis2,'Resid'])
	df.to_csv(tempdf,index=False)
	shade_aln2png(seqmsa,filename=temppng,shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=ruler,show_seq_names=False,show_seq_length=False,funcgroups=funcgroups)

	#let's write an R-script
	a=open(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'),'w')

	a.write(r"""
	library(ggplot2)
library(fitdistrplus)
library(gtools)
library(png)
library(grid)
library(reshape2)
library(gridExtra)
library(plyr)""")

	a.write("""
	df<-read.csv("%s",skip=0,header=TRUE,check.name=FALSE)
img <- readPNG("%s")

"""%(tempdf,temppng))

	a.write("""
	seqimg <- rasterGrob(img, interpolate=TRUE,width=1)

theme_set(theme_bw()+theme(panel.border =element_rect(linetype = "dashed", colour = "white")))

ndf <- melt(df, id="Resid")

a<-ggplot(data=ndf,aes(x=Resid+1+%f,y=value,fill=variable))+
geom_bar(stat='identity',position='dodge')+#scale_y_continuous(limits=c(-5,6),breaks=seq(0,6),labels=seq(0,6,by=1))+
# scale_x_continuous(limits=c(0,136),labels=c(),breaks=c(),expand=c(0,0))+
# scale_fill_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
# scale_color_manual(breaks=c('CA','bSCH','amore'),values=c('amore'='red','CA'='blue','bSCH'='green'),labels=c(expression(paste('C',alpha,'-atoms')),'Side chain','> 6Å      '),name='')+
xlab('Resid')+
ylim(%f,%f)+
# geom_point(data=h3data_amore,aes(color=color),fill='red',size=2)+
#ylab("RMSD, Å")+
ggtitle("%s")+
annotation_custom(seqimg, ymin=%f, ymax=0, xmin=0.5,xmax=%f)
"""%(offset,-max(profile)*1.01*ltune,max(profile)*1.01,title[-len(profile):len(title)],-max(profile)*ltune,len(profile)+0.5))

	a.write("""
	ggsave("%s",plot=a,height=%f,width=%f)

"""%(filename if filename[-3:]=='png' else (filename+'.png'),4.0*htune,12./60.*len(profile)))

	a.close()
	os.system('R --vanilla --slave < '+os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.remove(os.path.join(CONFIG.TEMP_DIR,'prof4seq.r'))
	os.remove(tempdf)
	os.remove(temppng)



if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A')])
    features=get_hist_ss_in_aln_for_shade(msa,below=True)
    plot_prof4seq('default',map(np.abs,map(np.sin,range(len(msa[0])))),msa,features,axis='conservation')
    # print features
    # shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    


