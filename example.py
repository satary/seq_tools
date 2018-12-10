import seq_tools.plot4seq as p
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seq_tools.hist_ss as st
import numpy as np
#import get_hist_ss_in_aln_for_shade

human_h2a_z_core=Seq('SRSQRAGLQFPVGRIHRHLKSRTTSHGRVGATAAVYSAAILEYLTAEVLELAGNASKDLKVKRITPRHLQLAIRGDEELDSLI-KATIAGGGVIPHIHKSLIG')
xenopus_h2a_core=Seq('TRSSRAGLQFPVGRVHRLLRKGNYAE-RVGAGAPVYLAAVLEYLTAEILELAGNAARDNKKTRIIPRHLQLAVRNDEELNKLLGRVTIAQGGVLPNIQSVLLP')
  
msa=MultipleSeqAlignment([SeqRecord(xenopus_h2a_core,id='H2A',name='H2A'),SeqRecord(human_h2a_z_core,id='H2Ah',name='H2Ah')])
features=st.get_hist_ss_in_aln_for_shade(msa,below=True)
p.plot_prof4seq('default',list(map(np.abs,map(np.sin,range(len(msa[0]))))),msa,features,axis='conservation')
