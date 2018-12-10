# -*- coding: utf-8 -*-
"""
Tools to work with sequence sets and alignments.

Retrive from NCBI by gis, align with MUSCLE,
profiles with AL2CO.
Trim alignemnt to show only postitions found in another sequence.

Now extended to do taxo_msa
and also just to display a set of annotated sequences on taxonomy 

Requirements: ete2==2.3.9 and PyQt
The taxo_seq_architecture to fully work requires the modified version of ete2 v.2.3.
on top of the installed ete2 with PyQt
https://github.com/molsim/ete/tree/2.3

"""
import uuid
from Bio import ExPASy
from Bio import SwissProt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
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
from Bio.SeqUtils.CheckSum import seguid
import numpy as np

from Bio.Align import MultipleSeqAlignment
import re
from Bio import AlignIO
from Bio.Align.Applications import MuscleCommandline
import subprocess

#from ete2 import NCBITaxa
# import pylab
from StringIO import StringIO
import L_shade_hist_aln
from Bio.Align.AlignInfo import SummaryInfo
from Bio.Emboss.Applications import NeedleCommandline
#from pylab import *
# sys.path.append('../sec_str')

from hist_ss import get_hist_ss_in_aln
 
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist


import CONFIG
Entrez.email = CONFIG.EMAIL

TEMP_DIR=CONFIG.TEMP_DIR
MUSCLE_BIN=CONFIG.MUSCLE_BIN_PATH
BLOSSUM_PATH=CONFIG.BLOSSUM_PATH

PATH_TO_AL2CO=CONFIG.PATH_TO_AL2CO
os.environ['PATH']='/Users/alexeyshaytan/soft/hmmer3.0/bin:/Users/alexeyshaytan/soft/al2co:/Users/alexeyshaytan/soft/x3dna-v2.1/bin:/Users/alexeyshaytan/soft/amber12/bin:/Users/alexeyshaytan/soft/sratoolkit/bin:/Users/alexeyshaytan/soft/bins/gromacs-4.6.3/bin:/opt/local/bin:/opt/local/sbin:/Users/alexeyshaytan/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/X11/bin:/usr/local/ncbi/blast/bin:/usr/texbin'+os.path.sep+os.environ['PATH']


def get_prot_seq_by_gis(gi_list):
    """
    Download a dictionary of Seqs (not fasta seqrecs - no identifiers) from NCBI given a list of GIs.
    """

    print("Downloading FASTA SeqRecords by GIs from NCBI")
    num=len(gi_list)
    fasta_seq=dict()
    for i in range(int(num/1000)+1):
        while i*1000 < num:
            try:
                print("Fetching %d th thousands from %d"%(i,num))
                strn = ",".join(map(str,gi_list[i*1000:(i+1)*1000]))
                request=Entrez.epost(db="protein",id=strn)
                result=Entrez.read(request)
                # print result
                webEnv=result["WebEnv"]
                queryKey=result["QueryKey"]
                handle=Entrez.efetch(db="protein",rettype='fasta',retmode='text',webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle,'fasta'):
                    fasta_seq[r.id.split('|')[1]]=r.seq
            except Exception as e:
                    print e
                    continue
            break
    print("Sequences downloaded:")
    print(len(fasta_seq))
    return(fasta_seq)




def get_prot_seqrec_by_gis(gi_list):
    """
    Download a dictionary of fasta SeqsRec from NCBI given a list of GIs.
    """

    print("Downloading FASTA SeqRecords by GIs from NCBI")
    num=len(gi_list)
    fasta_seqrec=dict()
    for i in range(int(num/1000)+1):
        while i*1000 < num:
            try:
                print("Fetching %d th thousands from %d"%(i,num))
                strn = ",".join(gi_list[i*1000:(i+1)*1000])
                request=Entrez.epost(db="protein",id=strn)
                result=Entrez.read(request)
                webEnv=result["WebEnv"]
                queryKey=result["QueryKey"]
                handle=Entrez.efetch(db="protein",rettype='fasta',retmode='text',webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle,'fasta'):
                    fasta_seqrec[r.id.split('|')[1]]=r
            except:
                    continue
            break
    print("FASTA Records downloaded:")
    print(len(fasta_seqrec))
    return(fasta_seqrec)


def get_prot_gbrec_by_gis(gi_list):
    """
    Download a dictionary of fasta GenBank Rec from NCBI given a list of GIs.
    """

    print("Downloading FASTA SeqRecords by GIs from NCBI")
    num=len(gi_list)
    fasta_seqrec=dict()
    for i in range(int(num/1000)+1):
        while i*1000 < num:
            try:
                print("Fetching %d th thousands from %d"%(i,num))
                strn = ",".join(gi_list[i*1000:(i+1)*1000])
                request=Entrez.epost(db="protein",id=strn)
                result=Entrez.read(request)
                webEnv=result["WebEnv"]
                queryKey=result["QueryKey"]
                handle=Entrez.efetch(db="protein",rettype='gb',retmode='text',webenv=webEnv, query_key=queryKey)
                for r in SeqIO.parse(handle,'gb'):
                    # print r.features[0].qualifiers['db_xref']
                    # fasta_seqrec[r.annotations['gi']]=r
                    fasta_seqrec[r.annotations['accessions'][0]]=r
            except Exception as e:
                    print "Exception detected: %s"%e
                    print r.annotations
                    continue
            break
    print("FASTA Records downloaded:")
    print(len(fasta_seqrec))
    return(fasta_seqrec)



# def del_gaps(alignment,key_sequence):
#     """Deletes columns in alignment that do not correspond to key sequence"""

#     #Let's get the key sequence id in the alignmnent
#     for i in range(len(alignment)):
#         if(str(alignment[i].seq).replace('-','')==str(key_sequence)):
#             key_index=i
#             print("Key index %d found"% i)
#             break

#     new_aln=alignment[:,0:1]
#     for i in range(len(alignment[key_index].seq)):
#         if(alignment[key_index][i] != '-'):
#             new_aln=new_aln+alignment[:,i:i+1]

#     return(new_aln[:,1:])



def aln_undup(alignment):
    """Removes duplicate keys"""
    aln=MultipleSeqAlignment([])
    checksums = set()
    for record in alignment:
        checksum = seguid(record.seq)
        if checksum in checksums:
            print "Ignoring %s" % record.id
            continue
        checksums.add(checksum)
        aln.append(record)

    return aln

def cons_prof(alignment,f=2,c=2,m=0,norm='F'):
    """Uses al2co to build conservation profiles"""
# do Loading    
 # -f    Weighting scheme for amino acid frequency estimation [Integer] Optional
 #        Options:
 #        0=unweighted,
 #        1=weighted by the modified method of Henikoff & Henikoff (2)(3),
 #        2=independent-count based (1)(4)
 #        Default = 2

   # -c    Conservation calculation method [Integer] Optional
   #      Options:
   #      0=entropy-based    C(i)=sum_{a=1}^{20}f_a(i)*ln[f_a(i)], where f_a(i)
   #        is the frequency of amino acid a at position i,
   #      1=variance-based   C(i)=sqrt[sum_{a=1}^{20}(f_a(i)-f_a)^2], where f_a
   #        is the overall frequency of amino acid a,
   #      2=sum-of-pairs measure   C(i)=sum_{a=1}^{20}sum_{b=1}^{20}f_a(i)*f_b(i)*S_{ab},
   #        where S_{ab} is the element of a scoring matrix for amino acids a and b
   #      Default = 0
    
  # -m    Scoring matrix transformation [Integer] Optional
  #       Options:
  #       0=no transformation,
  #       1=normalization S'(a,b)=S(a,b)/sqrt[S(a,a)*S(b,b)],
  #       2=adjustment S"(a,b)=2*S(a,b)-(S(a,a)+S(b,b))/2
  #       Default = 0

 # -n    Normalization option [T/F] Optional
 #        Subtract the mean from each conservation index and divide by the
 #        standard deviation.
 #        Default = T

 #          -s    Input file with the scoring matrix [File in] Optional
 #        Format: NCBI
 #        Notice: Scoring matrix is only used for sum-of-pairs measure
 #        with option -c  2.
 #        Default = identity matrix
    os.environ['PATH']+=":"+PATH_TO_AL2CO

    AlignIO.write(alignment,TEMP_DIR+"/tmp.aln","clustal")
    proc=subprocess.Popen(["al2co","-i",TEMP_DIR+"/tmp.aln","-f",str(f),'-n',norm,'-c',str(c),'-m',str(m),'-s',BLOSSUM_PATH],shell=False, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return_code=proc.wait()
    # print proc
    # print proc.stdout
    os.remove(TEMP_DIR+'/tmp.aln')
    p=re.compile("\d+\s+\S\s+\S+")
    cons_prof=[]
    for line in proc.stdout:
        print line
        if(p.match(line)):
            cons_prof.append(line.split()[2])
    return cons_prof    
     

def trim_aln_gaps(alignment,threshold=0.8):
    """Removes positions with more than threshold gaps in alignment"""
    a=SummaryInfo(alignment)
    cons=a.gap_consensus(threshold=threshold, ambiguous='X')
    new_aln=alignment[:,0:0]
    for c,i in zip(cons,range(len(cons))):
        if(c=='-'):
            continue
        else:
            new_aln+=alignment[:,i:i+1]

    return new_aln


def trim_aln_to_key_seq(alignment,key_sequence):
    """Deletes columns in alignment that does not correspond to key sequence, sequence should be present in the alignment already"""

    #Let's get the key sequence id in the alignmnent
    for i in range(len(alignment)):
        if(str(alignment[i].seq).replace('-','')==str(key_sequence)):
            key_index=i
            # print("Key index %d found"% i)
            break

    new_aln=alignment[:,0:1]
    for i in range(len(alignment[key_index].seq)):
        if(alignment[key_index][i] != '-'):
            new_aln=new_aln+alignment[:,i:i+1]

    return(new_aln[:,1:])


def trim_aln_to_seq(alignment,sequence):
    """Trim alignment to a sequence, i.e. leave only postions that correspond to this sequence
    Note that seqeuence should be incorportatable into alignment without additional gaps in alignment.
    """

    n1=str(uuid.uuid4())
    n2=str(uuid.uuid4())
    #Get consensus
    a=SummaryInfo(alignment)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')
    #Needle it
    SeqIO.write([SeqRecord(cons,id='CONS',name='CONS')],n1+'.fasta','fasta')
    SeqIO.write([SeqRecord(sequence,id='KEY',name='KEY')],n2+'.fasta','fasta')

#Now we will redo it with Needlman Wunsh - the global alignment
    needle_cline = NeedleCommandline(asequence=n1+".fasta", bsequence=n2+".fasta",gapopen=10, gapextend=0.5, outfile=n1+".txt")
    stdout, stderr = needle_cline()
# print('Needle alignment')

    align = AlignIO.read(n1+".txt", "emboss")
    os.system('rm %s.fasta %s.fasta %s.txt'%(n1,n2,n1))
    # print align
    # print alignment
    align.extend(alignment)
    a=align[1:,:]

    return trim_aln_to_key_seq(a,sequence)[1:,:]

def trim_aln_to_seq_length(alignment,sequence):
    """Trim alignment to a sequence, i.e. leave only postions that correspond to this sequence span"""

    n1=str(uuid.uuid4())
    n2=str(uuid.uuid4())
    #Get consensus
    a=SummaryInfo(alignment)
    cons=a.dumb_consensus(threshold=0.1, ambiguous='X')
    #Needle it
    SeqIO.write([SeqRecord(cons,id='CONS',name='CONS')],n1+'.fasta','fasta')
    SeqIO.write([SeqRecord(sequence,id='KEY',name='KEY')],n2+'.fasta','fasta')

#Now we will redo it with Needlman Wunsh - the global alignment
    needle_cline = NeedleCommandline(asequence=n1+".fasta", bsequence=n2+".fasta",gapopen=10, gapextend=0.5, outfile=n1+".txt")
    stdout, stderr = needle_cline()
# print('Needle alignment')

    align = AlignIO.read(n1+".txt", "emboss")
    os.system('rm %s.fasta %s.fasta %s.txt'%(n1,n2,n1))
    # print align
    # print alignment
    #first seq is consensus, we need to get borders useing second one.
    seq=str(align[1,:].seq)
    # print seq
    begin=seq.index(str(sequence[0]))
    end=len(seq)-seq[::-1].index(str(sequence[-1]))
    print begin
    print end
    return alignment[:,begin:end]

    
def muscle_aln(seqreclist,**kwargs):
    """Align with muscle
    muscle_aln(f_fasta_dict.values(),gapopen=float(-20))
    """
            #let's write to file
    s=str(uuid.uuid4())

    output_handle = open(TEMP_DIR+"/%s.fasta"%s, "w")
    SeqIO.write(seqreclist, output_handle, "fasta")
    output_handle.close()
        
    muscle_cline = MuscleCommandline(MUSCLE_BIN,input=TEMP_DIR+"/%s.fasta"%s,**kwargs)
    # print muscle_cline
    stdout, stderr = muscle_cline()
      #  # print stderr
    # print stdout
    msa = AlignIO.read(StringIO(stdout), "fasta")
    os.system("rm "+TEMP_DIR+"/%s.fasta"%s)
    return msa
        

def add_consensus(alignment,threshold=0.9, ambiguous='-',name='consensus'):
    """Add a consensus line"""
    a=SummaryInfo(alignment)
    # cons=a.dumb_consensus(threshold, ambiguous)
    cons=a.gap_consensus(threshold, ambiguous)
    alignment.extend([SeqRecord(cons,id=name,name=name)])
    return alignment


def cluster_seq_support_nw(seq_dict,ident_thresh=0.90):
    matrix = matlist.blosum62
    items=seq_dict.items()
    ident_matrix=np.identity(len(items))

    for ind1 in range(len(items)):
        (gi1,sr1)=items[ind1]
        # print ind1,' from ',len(items)
        for ind2 in range(ind1):
            (gi2,sr2)=items[ind2]
            # pairwise2.align.globalds(p53_human, p53_mouse, matrix, gap_open, gap_extend)
            # alns = pairwise2.align.globalds(sr1.seq, sr2.seq, matrix, -10, -0.5)
            # alns = pairwise2.align.globalxx(sr1.seq, sr2.seq)
            needle_cline = NeedleCommandline(asequence="asis::"+sr1.seq, bsequence="asis::"+sr2.seq,gapopen=10, gapextend=0.5, outfile=TEMP_DIR+"/needle.txt")
            stdout, stderr = needle_cline()
            align = AlignIO.read(TEMP_DIR+"/needle.txt", "emboss")
            # print align
            # l1,l2=alns[0][0:2]
            l1=align[0].seq
            l2=align[1].seq

            matches = sum(aa1 == aa2 for aa1, aa2 in zip(l1, l2))
            identity = matches / float(len(l1))
            # print identity
            ident_matrix[ind1,ind2]=identity
            ident_matrix[ind2,ind1]=identity

    #crude clustering
    # print ident_matrix
    support=dict()
    # print ident_matrix
    for i in range(len(items)):
        support[items[i][0]]=0
        for k in range(len(items)):
            if(ident_matrix[i,k]>ident_thresh):
                support[items[i][0]]+=1

    return support


def cluster_seq_support(seq_dict,ident_thresh=0.90):
    """
    This function for a set of sequences - returns a dictionary,
    where for every sequence there is a number - corresponding to the number
    of similar sequences in this sequence set.
    """
    matrix = matlist.blosum62
    items=seq_dict.items()
    ident_matrix=np.identity(len(items))
    # if(513031220 in seq_dict.keys()):
        # print "kuku"
    msa=muscle_aln([SeqRecord(v.seq,id=str(k)) for k,v in seq_dict.iteritems()])
    seq_dict={int(sr.id):sr.seq for sr in msa}
    items=seq_dict.items()

    for ind1 in range(len(items)):
        (gi1,s1)=items[ind1]
        # print ind1,' from ',len(items)
        for ind2 in range(ind1):
            (gi2,s2)=items[ind2]
            # print s1
            # print s2
            matches = sum(1 if ((aa1==aa2)&(aa1!='-')) else 0 for aa1, aa2 in zip(s1, s2))
            length = float(sum(0 if ((aa1=='-') and (aa2=='-')) else 1 for aa1, aa2 in zip(s1, s2)))
            identity = matches/length 
            # print matches,length
            ident_matrix[ind1,ind2]=identity
            ident_matrix[ind2,ind1]=identity

    #crude clustering
    # print ident_matrix
    support=dict()
    # print ident_matrix
    for i in range(len(items)):
        support[items[i][0]]=0
        for k in range(len(items)):
            if(ident_matrix[i,k]>ident_thresh):
                support[items[i][0]]+=1

    return support




def taxo_msa(outfile='taxo_msa.svg',taxids=[],annotation='',msa=[],title='',width=2000):
    """
    Visualize MSA together with a taxonomy tree
    taxids - list of taxids in the same order as seqs in msa
    """
    from ete2 import NCBITaxa, Tree, SeqMotifFace, TreeStyle, add_face_to_node,AttrFace,TextFace

    # taxid2gi={f_df.loc[f_df.gi==int(gi),'taxid'].values[0]:gi for gi in list(f_df['gi'])}
    # gi2variant={gi:f_df.loc[f_df.gi==int(gi),'hist_var'].values[0] for gi in list(f_df['gi'])}

    # msa_dict={i.id:i.seq for i in msa_tr}
    ncbi = NCBITaxa()
    taxids=map(int,taxids)

    t = ncbi.get_topology(taxids,intermediate_nodes=False)
    a=t.add_child(name='annotation')
    a.add_feature('sci_name','annotation')
    t.sort_descendants(attr='sci_name')
    ts = TreeStyle()
    def layout(node):
        # print node.rank
        # print node.sci_name
        if getattr(node, "rank", None):
            if(node.rank in ['order','class','phylum','kingdom']):   
                rank_face = AttrFace("sci_name", fsize=7, fgcolor="indianred")
                node.add_face(rank_face, column=0, position="branch-top")
        if node.is_leaf():
            sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
            node.add_face(sciname_face, column=0, position="branch-right")
        if node.is_leaf() and not node.name=='annotation':
            s=str(msa[taxids.index(int(node.name))].seq)
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            # gi=taxid2gi[int(node.name)]
            add_face_to_node(TextFace(' '+msa[taxids.index(int(node.name))].id),node,column=1, position = "aligned")
            # add_face_to_node(TextFace('      '+str(int(node.name))+' '),node,column=2, position = "aligned")
            # add_face_to_node(TextFace('      '+str(gi2variant[gi])+' '),node,column=3, position = "aligned")

        if node.is_leaf() and node.name=='annotation':
            if(annotation):
                s=annotation
                # get_hist_ss_in_aln_as_string(msa_tr)
            else:
                s=' '*len(msa[0].seq)
            seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            add_face_to_node(seqFace, node, 0, position="aligned")
            add_face_to_node(TextFace(' '+'SEQ_ID'),node,column=1, position = "aligned")
            # add_face_to_node(TextFace('       '+'NCBI_TAXID'+' '),node,column=2, position = "aligned")
            # add_face_to_node(TextFace('       '+'Variant'+'       '),node,column=3, position = "aligned")



    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(title, fsize=20), column=0)
    t.render(outfile, w=width, dpi=300, tree_style=ts)

def gen_fake_msa(seqreclist):
    """
    Fake an MSA by adding dashes to the end of sequences up to maximum length and 
    return it as an MSA
    """
    mlength=max(map(len,seqreclist))
    newseqreclist=list()
    for s in seqreclist:
        l=mlength-len(s.seq)
        newseqreclist.append(SeqRecord(id=s.id,name=s.name,seq=Seq(str(s.seq)+'-'*l)))
    msa=MultipleSeqAlignment(newseqreclist)
    return msa


def features_via_hmm(seq,hmmdb,eval_thresh=1.0,ftype='domain',conv_name_dict={}):
    
    """
    This function takes a Seq, runs hmmscan against a compressed hmmdb (prepare with hmmpress)
    and output a list of biobython SeqFeature.
    #Needs strictly HMMER 3.0!!!!
    """
    def get_color(str):
        colorlist=['Red','Green','Yellow','LightBlue','Cyan','Magenta','Orange','Pink','LightGreen']
        return colorlist[hash(str)%9]

    features=list()
    ufn=str(uuid.uuid4())
    SeqIO.write([SeqRecord(seq,id='QUERY',name='QUERY',description="QUERY")],ufn+'.fasta','fasta')
    subprocess.call(["hmmscan","-o",ufn+".out","--tblout",ufn+".tbl","--domtblout",ufn+".dtbl",hmmdb,ufn+'.fasta'])
    #Now let's read it

    for v in SearchIO.parse(ufn+".dtbl", "hmmscan3-domtab"):
        for hit in v:
            for h in hit.hsps:
                # print h
                if h.evalue<eval_thresh:
                    features.append(SeqFeature(FeatureLocation(h.query_start,h.query_end), type=ftype,qualifiers={'name':conv_name_dict.get(h.hit_id,h.hit_id),'evalue':h.evalue,'color':get_color(h.hit_id),'text':conv_name_dict.get(h.hit_id,h.hit_id)}))


    os.system("rm %s %s %s %s"%(ufn+'.fasta',ufn+'.out',ufn+'.tbl',ufn+'.dtbl'))
    return features


def taxo_seq_architecture(seqreclist=[],outfile='taxo_arch.svg',taxids=[],annotation='',title='',width=2000):
    """
    Visualize sequence architecture together with a taxonomy tree
    seqreclist - contains a list of seqres.
    each seqrec should have a list of features in biobython SeqFeature format.

    features of type "domain" will be plotted as boxes 
    features of type "xxxx" will be plotted as ...

    taxids - list of taxids in the same order as seqs in msa, if now provided will assume that seqrecs
    are in genbank format and attempt to get taxids from there.
    """
    # print sys.path
    # import ete2
    # print ete2.__version__
    from ete2 import NCBITaxa, Tree, SeqMotifFace, TreeStyle, add_face_to_node,AttrFace,TextFace

    aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','.','-']
    
    def get_color(str):
        colorlist=['red','green','yellow','lightblue','cyan','magenta','orange','pink','lightgreen']
        return colorlist[hash(str)%9]

    if len(taxids)==0:
        taxids=map(get_taxid_from_gbrec,seqreclist)

    ncbi = NCBITaxa()
    taxids=map(int,taxids)

    t = ncbi.get_topology(taxids,intermediate_nodes=False)
    # a=t.add_child(name='annotation')
    # a.add_feature('sci_name','annotation')
    t.sort_descendants(attr='sci_name')
    ts = TreeStyle()
    def layout(node):
        # print node.rank
        # print node.sci_name
        if getattr(node, "rank", None):
            if(node.rank in ['order','class','phylum','kingdom']):   
                rank_face = AttrFace("sci_name", fsize=7, fgcolor="indianred")
                node.add_face(rank_face, column=0, position="branch-top")
        if node.is_leaf():
            sciname_face = AttrFace("sci_name", fsize=9, fgcolor="steelblue")
            node.add_face(sciname_face, column=0, position="branch-right")
        if node.is_leaf() and not node.name=='annotation':
            #here we are adding faces and we need to play with seqmotif face
            seq=str(seqreclist[taxids.index(int(node.name))].seq)
            motifs=[]#[[0,len(seq), "seq", 10, 10, None, None, None]]
            for f in seqreclist[taxids.index(int(node.name))].features:
                if f.type=='domain':
                    motifs.append([f.location.start,f.location.end,"[]",None,10,"blue", f.qualifiers.get('color',get_color(f.qualifiers['name'])).lower(), "arial|8|black|%s"%f.qualifiers['name']])
                if f.type=='motif':
                    #It turns out that we need to solve overlap problem here, here it is solved only in case of one overlap
                    s=f.location.start
                    e=f.location.end
                    flag=True
                    overlappedm=[]
                    for m in motifs:
                        if m[2]=='seq' and m[0]<e and m[1]>s: #we have an overlap, four cases, preceding motife always is on top
                            flag=False
                            overlappedm.append(m)
                    if not flag: #we have to solve multiple overlap problem
                    #let's do it by scanning
                        sflag=False
                        eflag=False
                        for x in range(s,e+1):
                            if not sflag: #check if we can start
                                overlap=False
                                for m in overlappedm:
                                    if x>=m[0] and x<m[1]:
                                        overlap=True
                                if not overlap:
                                    ts=x
                                    sflag=True

                            #check if is time to end
                            if sflag and not eflag:
                                overlap=False
                                for m in overlappedm:
                                    if x==m[0]:
                                        overlap=True
                                if overlap or x==e:
                                    te=x
                                    eflag=True


                            if sflag and eflag:
                                motifs.append([ts,te,"seq",10,10,"black",f.qualifiers.get('color',get_color(f.qualifiers['name'])).lower(),None])
                                sflag=False
                                eflag=False
                    if flag:
                        motifs.append([f.location.start,f.location.end,"seq",10,10,"black",f.qualifiers.get('color',get_color(f.qualifiers['name'])).lower(),None])
            seqFace = SeqMotifFace(seq,motifs,scale_factor=1,seq_format="[]")
            seqFace.overlaping_motif_opacity = 1.0
            # seqFace.fg=aafgcolors
            # seqFace.bg=aabgcolors_gray

            add_face_to_node(seqFace, node, 0, position="aligned")
            # gi=taxid2gi[int(node.name)]
            add_face_to_node(TextFace(' '+seqreclist[taxids.index(int(node.name))].id+'         '),node,column=1, position = "aligned")
            # add_face_to_node(TextFace('      '+str(int(node.name))+' '),node,column=2, position = "aligned")
            # add_face_to_node(TextFace('      '+str(gi2variant[gi])+' '),node,column=3, position = "aligned")


        #We currently disable annotation
        if node.is_leaf() and node.name=='annotation':
            if(annotation):
                s=annotation
                # get_hist_ss_in_aln_as_string(msa_tr)
            else:
                s=' '*max(map(lambda x: len(x.seq),seqreclist))
            # seqFace = SeqMotifFace(s,[[0,len(s), "seq", 10, 10, None, None, None]],scale_factor=1)
            # add_face_to_node(seqFace, node, 0, position="aligned")
            # add_face_to_node(TextFace(' '+'SEQ_ID'),node,column=1, position = "aligned")


            # add_face_to_node(TextFace('       '+'NCBI_TAXID'+' '),node,column=2, position = "aligned")
            # add_face_to_node(TextFace('       '+'Variant'+'       '),node,column=3, position = "aligned")



    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.title.add_face(TextFace(title, fsize=20), column=0)
    t.render(outfile, w=width, dpi=300, tree_style=ts)


def debump_features(seqrec):
    """
    if features ovelap this is not good for graphics programs.
    here we debump them and return a new set of features
    numbering conventions here should as that of biopython (0 based, and [1,1] - nothing)
    """
    l=len(seqrec.seq)*[-1]
    for f,n in zip(seqrec.features,range(len(seqrec.features))):
        # print f
        for i in range(f.location.start,f.location.end):
            l[i]=n
    # newseqrec=SeqRecord(id=seqrec.id,description=seqrec.description,name=seqrec.name, seq=seqrec.seq,annotations=seqrec.annotations,dbxref=seqrec.dbxref)
    newfeat=[]
    flag=l[0]
    start=0
    for s,i in zip(l,range(len(l))):
        if(s!=flag):
            n=l[i-1]
            if n>=0:
                feat=seqrec.features[n]
                newfeat.append(SeqFeature(FeatureLocation(start,i), type=feat.type,qualifiers=feat.qualifiers))
            start=i
            flag=s
    return newfeat

if __name__ == '__main__':
    human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
    xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    test_h2a_core=Seq('TRSTRAHLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
    ts=Seq('GAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLG')
    
    # human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLIKATIAGGGVIPHIHKSLIG')
    msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2A.Z',name='H2A.Z'),SeqRecord(test_h2a_core,id='test',name='test')])
    # print trim_aln_gaps(msa,threshold=0.1)
    # print trim_aln_to_seq(msa,ts)
    print cons_prof(msa)

    # get_pdf('H2A',MultipleSeqAlignment([SeqRecord(human_h2a_z_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='1H2A.Z',name='H2A.Z')]),'H2AvsH2A.Z',[0,5,1])

    
            