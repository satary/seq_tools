# -*- coding: utf-8 -*-
"""
This is a library that makes good images of shaded alignments 
through TeXShade.

Input:
1) alignment (might be a single sequence).
2) Shading options
3) Features list
Output:
pdf or image (needed for automaitc plotting)


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
import cPickle as pickle
from Bio import SeqIO
import math
from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

from StringIO import StringIO

from hist_ss import get_hist_ss
from hist_ss import get_hist_ss_in_aln, get_hist_ss_in_aln_for_shade
from Bio.Align.AlignInfo import SummaryInfo
import L_aln_tools

import CONFIG
TEMP_DIR=CONFIG.TEMP_DIR

def seqfeat2shadefeat(msa,seqref=None,idseqref=True):
    """
    converts SeqFeature records for every sequence in msa to our style feature list
    """
    features=[]
    #In texshade [1,1] will be colored as one residue. this is different from biopython where [1,1] selects nothing!
    for m,i in zip(msa,range(len(msa))):
        for f in m.features:
            if seqref:
                sr=seqref
            elif idseqref:
                sr=m.id
            else:
                sr=i+1

            if f.type=='motif':
                features.append({'style':'shaderegion','seqref':sr,'sel':[f.location.start,f.location.end-1],'shadcol':f.qualifiers['color']})
            if f.type=='domain':
                features.append({'style':'shaderegion','seqref':sr,'sel':[f.location.start,f.location.end-1],'shadcol':f.qualifiers['color']})
            if f.type=='frameblock':
                features.append({'style':'frameblock','seqref':sr,'sel':[f.location.start,f.location.end-1],'color':f.qualifiers['color']})
            if f.type=='---':
                features.append({'style':'---','seqref':sr,'sel':[f.location.start,f.location.end-1],'color':f.qualifiers['color'],'text':f.qualifiers['text']})

    return features

def shade_aln2png(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,rotate=False,threshold=[80,50],resperline=0,margins=None, density=150):
    intf=TEMP_DIR+'/tempshade.pdf'
    shade_aln2pdf(msa,intf,shading_modes,features,title,legend, logo,hideseqs,splitN,setends,ruler,show_seq_names,show_seq_length,funcgroups,threshold,resperline=resperline)
    #let's use imagemagic
    #margins - add on each side margins %
    if margins:
        m=margins
    else:
        m=0
    print("Converting PDF to PNG")
    if rotate:
        cmd='convert -density %d '%density+intf+' -trim -bordercolor White -border %.3f%%x0%% -rotate -90 %s'%(m,filename if filename[-3:]=='png' else filename+'.png')
        print(cmd)
        os.system(cmd)
    else:
        cmd='convert -density %d '%density+intf+' -trim -bordercolor White -border %.3f%%x0%% %s'%(m,filename if filename[-3:]=='png' else filename+'.png')
        print(cmd)
        os.system(cmd)


def shade_aln2pdf(msa,filename='default',shading_modes=['similar'],features=[],title='',legend=True, logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=True,show_seq_length=True,funcgroups=None,threshold=[80,50],resperline=0):
    """
will convert msa to a shaded pdf.
shading_modes: similar, ... see write_texshade code
features - a list of dictionaries:
{'style':'shaderegion','shadeblock','frameblock','helix','loop','-->','---','<--',',-,' - all texshade features types of features
+ frameblock + shaderegion + shadeblock
'position':'top','bottom','ttop','bbottom', etc. if no - automatic
'seqref':number of sequence for selection - default consensus
'sel':[begin,end] - region selection, this is in 0 based numbering (as in msa - we override here the TexShade 1-based numbering)
'text':'text'
'color':'Red' this is for frame block
'thickness':1.5 -1.5pt also for frameblock
'rescol':'Black' this is for \shaderegion{ seqref  }{ selection }{ res.col. }{ shad.col. }
'shadcol':'Green' - the same
funcgroup example fg="\\funcgroup{xxx}{CT}{White}{Green}{upper}{up} \\funcgroup{xxx}{GA}{White}{Blue}{upper}{up}"
}

    """
    ns='consensus'
    hideseqs_by_name=[]
    if(len(msa)==1):
        msa=msa[:]
        msa.extend([SeqRecord(msa[0].seq,id='dum',name='dum')])
        hideseqs_by_name.append('dum')
    #####if we are splitting the alignment into blocks - get number of blocks

    a_len=len(msa)
    num=int(a_len/splitN)+1
    while ((a_len-(num-1)*splitN)<2):
        splitN=splitN+1
        num=int(a_len/splitN)+1
    print "Chosen splitting parameters"    
    print a_len, splitN

    ####iterate over blocks and create alignment fasta files
    for i in range(num):
        t_aln=msa[(i*splitN):((i+1)*splitN)]
        AlignIO.write(t_aln,open(TEMP_DIR+'/alignment%d.fasta'%i,'w'), 'fasta')
    if resperline==0:
        res_per_line=len(msa[0])
    else:
        res_per_line=resperline

    #prepare feature section

    #alias dict
    aliasf={'alpha':'helix','beta':'-->','domain':'loop'}


    features_dict={}
    for i in features:
        sr=str(i.get('seqref','consensus'))
        if(i['style']=='block' or i['style']=='frameblock'):
            features_dict[sr]=features_dict.get(sr,'')+"\\frameblock{%s}{%d..%d}{%s[%.1fpt]}"%(sr,i['sel'][0]+1,i['sel'][1]+1,i.get('color','Red'),i.get('thickness',1.5))
        elif(i['style']=='shaderegion' or i['style']=='shadeblock'):
            features_dict[sr]=features_dict.get(sr,'')+"\\%s{%s}{%d..%d}{%s}{%s}"%(i['style'],sr,i['sel'][0]+1,i['sel'][1]+1,i.get('rescol','Black'),i.get('shadcol','Green'))
        else:
            features_dict[sr]=features_dict.get(sr,'')+"\\feature{%s}{%s}{%d..%d}{%s}{%s}"%(i.get('position','top'),sr,i['sel'][0]+1,i['sel'][1]+1,aliasf.get(i.get('style','loop'),i.get('style','loop')),i.get('text','').replace('_','-'))
            
    a=open(TEMP_DIR+'/align.tex','w')

    a.write(r"""\documentclass[11pt,landscape]{article}
%\documentclass{standalone}
%\usepackage[a0paper]{geometry}
\usepackage{hyperref}
""")

    lmult=math.ceil(float(len(msa[0]))/float(res_per_line))
    ca_len=a_len*lmult
    h=((ca_len/30.*18 + (2.5 if legend else 0.0)) if (ca_len/30.*18 + (2.5 if legend else 0.0) <18.0) else 18)
    w=(22/200.*res_per_line+2.5)

    if title:
        h+=1
        if w<len(title)*0.4:
            w=len(title)*0.4
    a.write("""
\\usepackage[paperwidth=%fin, paperheight=%fin]{geometry}
        """%(w,h))

    a.write(r"""
\usepackage{texshade}

\begin{document}""")
    a.write("""
\\pagenumbering{gobble}
\\centering
\\Huge{%s}
\\vspace{-0.5in}"""%title)

    for i in range(num-1):
        features_code=features_dict.get('consensus','')
        for s,ns in zip(msa[(i*splitN):((i+1)*splitN)],range((i*splitN),((i+1)*splitN))):
            features_code+=features_dict.get(s.id,'')
            features_code+=features_dict.get(str(ns+1),'')
        write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%i , features_code, res_per_line,False,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name,funcgroups=funcgroups,threshold=threshold)
    
    i=num-1
    features_code=features_dict.get('consensus','')
    for s,ns in zip(msa[(i*splitN):((i+1)*splitN)],range((i*splitN),((i+1)*splitN))):
        features_code+=features_dict.get(s.id,'')
        features_code+=features_dict.get(str(ns+1),'')
    write_texshade(a,TEMP_DIR+'/alignment%d.fasta'%(num-1) , features_code, res_per_line,legend,shading_modes,logo,hideseqs,setends,ruler,numbering_seq='consensus',hide_ns=False,show_seq_names=show_seq_names,show_seq_length=show_seq_length,hideseqs_by_name=hideseqs_by_name,funcgroups=funcgroups,threshold=threshold)

    a.write(r"""
\end{document} """)
    a.close()

    command='/usr/bin/env pdflatex --file-line-error --synctex=1 -output-directory=%s --save-size=10000  %s/align.tex > /dev/null'%(TEMP_DIR,TEMP_DIR)

    print('Launcning Latex:')
    print(command)
    os.system(command)
    os.system('mv '+TEMP_DIR+'/align.pdf %s'%(filename if filename[-3:]=='pdf' else (filename+'.pdf')))




def write_texshade(file_handle,aln_fname,features,res_per_line=120,showlegend=True,shading_modes=['similar'],logo=False,hideseqs=False,setends=[],ruler=False,numbering_seq='consensus',hide_ns=False,show_seq_names=True,show_seq_length=True,hideseqs_by_name=[],funcgroups=None,threshold=[80,50]):

    for shading in shading_modes:
        shading=str(shading)
        file_handle.write("""
    \\begin{texshade}{%s}
    \\residuesperline*{%d}
    """%(aln_fname,res_per_line))
       


        file_handle.write(r"""
    \seqtype{P}
    \defconsensus{{}}{*}{upper}
    """)
        #a very dirty hack
        if(setends):
            if(numbering_seq=='consensus'):
                file_handle.write("""
    \\setends[%d]{%s}{%d..%d}
    """%(setends[0]+1,numbering_seq,setends[0]+setends[0],setends[1]+setends[0]+setends[0]-1))
            else:
                    file_handle.write("""
    \\setends{%s}{%d..%d}
    """%(numbering_seq,setends[0],setends[1]))
           
        if(ruler):
            if(ruler=='bottom'):
                file_handle.write("""
    \\showruler{bottom}{%s}
    """%numbering_seq)
            else:
                file_handle.write("""
    \\showruler{top}{%s}
    """%numbering_seq)
        if(hide_ns):
            file_handle.write("""
    \\hideseq{%s}
    """%numbering_seq)

        file_handle.write(r"""
    \seqtype{P}
    """)
        if(not show_seq_names):
            file_handle.write(r"""
    \hidenames
    """) 
        if(not show_seq_length):
            file_handle.write(r"""
    \hidenumbering
    """) 
        if((shading=='similar')|(shading=='0')):
            file_handle.write(r"""
    \shadingmode{similar}
    \threshold[%d]{%d}
    """%(threshold[0],threshold[1]))
            if(logo):
                file_handle.write(r"""
    \showsequencelogo{top} \showlogoscale{leftright}
    \namesequencelogo{logo}
    """)

        if((shading=='hydropathy_functional')|(shading=='1')):
            file_handle.write(r"""
    \shadingmode[hydropathy]{functional}
    \shadeallresidues
    \threshold[%d]{%d}
    
    """%(threshold[0],threshold[1]))     
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[hydropathy]{top} \showlogoscale{leftright}

    """)
        if((shading=='chemical_functional')|(shading=='2')):
            file_handle.write(r"""

    \shadingmode[chemical]{functional}
    \shadeallresidues
    """)
            if(logo):
                file_handle.write(r"""
    \showsequencelogo[chemical]{top} \showlogoscale{leftright}

    """)
        if((shading=='structure_functional')|(shading=='3')):
            file_handle.write(r"""
    \shadingmode[structure]{functional}
    \shadeallresidues
    """)

        if((shading=='charge_functional')|(shading=='4')):
            file_handle.write(r"""
    \shadingmode[charge]{functional}
    \shadeallresidues
    """)

        if((shading=='diverse')|(shading=='5')):
            file_handle.write(r"""
    \shadingmode{diverse}

    """)
        if(funcgroups):
            file_handle.write(funcgroups)
            
        if(hideseqs):
            file_handle.write(r"""
    \hideseqs

    """)
        for s in hideseqs_by_name:
            file_handle.write("""
    \\hideseq{%s}
    """%s)

        file_handle.write(r"""

    %\setends{consensus}{1..160}
    %\setends{consensus}{1..160}


    %\feature{ttop}{1}{1..160}{bar:conservation}{}
    %\showfeaturestylename{ttop}{conserv.}
    \ttopspace{-\baselineskip}

    %\feature{top}{1}{1..160}{color:charge}{}
    %\showfeaturestylename{top}{charge}

    %\feature{bottom}{1}{1..160}{color:molweight[ColdHot]}{}
    %\showfeaturestylename{bottom}{weight}

    %\bbottomspace{-\baselineskip}
    %\feature{bbottom}{2}{1..160}{bar:hydrophobicity[Red,Gray10]}{}
    %\showfeaturestylename{bbottom}{hydrophob.}

    %\bargraphstretch{3}
    %\featurestylenamescolor{Red}
    %\featurestylenamesrm  \featurestylenamesit

    %\showsequencelogo{top}


    %\showconsensus[ColdHot]{bottom}
    \showconsensus[black]{top}

    %\defconsensus{.}{lower}{upper}
    %\defconsensus{{}}{lower}{upper}
    %\defconsensus{{}}{*}{upper}

    """)
        file_handle.write(features)
        file_handle.write(r"""
    \featurerule{1mm}""")
        if(showlegend):
            file_handle.write(r"""
    \showlegend""")

        align = AlignIO.read(open(aln_fname,'r'), "fasta")
        for a,i in zip(align,range(len(align))):
            # print a.id.replace('|',' | ')
            file_handle.write("""
    \\nameseq{%d}{%s}"""%(i+1,a.id.replace('|',' | ')))

        file_handle.write(r"""

    \end{texshade}
    """)

def feature_str2dict(featurestring,position='top'):
    """converts string of secondary structure annotation (like in VMD) to our type of dict"""
    #numbering should be 0 based
    #HHHHHHEEEEEBBBBBBCCCCCbTGI
    features=[]
    style=''
    begin=-1
    for i,s in enumerate(featurestring):
        if(s in ['H','G','I']): #helices
            if style!='helix':
                style='helix'
                begin=i
        else:
            if style=='helix':
                style=''
                end=i-1
                features.append({'style':'helix','sel':[begin,end],'position':position})

    for i,s in enumerate(featurestring):
        if(s in ['E','B','b']): #helices
            if style!='-->':
                style='-->'
                begin=i
        else:
            if style=='-->':
                style=''
                end=i-1
                features.append({'style':'-->','sel':[begin,end],'position':position})
    return features




    #prof=cons_prof(alignment)
    #pylab.plot(prof)
if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A')])
    features=get_hist_ss_in_aln_for_shade(msa,below=True)
    # features=[{'style':'fill:$\uparrow$','sel':[5,10],'text':'test'}]
    print features
    shade_aln2png(msa,filename='default',shading_modes=['charge_functional'], legend=False, features=features,title='',logo=False,hideseqs=False,splitN=20,setends=[],ruler=False,show_seq_names=False,show_seq_length=False)
    



    
            