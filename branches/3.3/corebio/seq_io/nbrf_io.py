
#  Copyright (c) 2006, The Regents of the University of California, through 
#  Lawrence Berkeley National Laboratory (subject to receipt of any required
#  approvals from the U.S. Dept. of Energy).  All rights reserved.

#  This software is distributed under the new BSD Open Source License.
#  <http://www.opensource.org/licenses/bsd-license.html>
#
#  Redistribution and use in source and binary forms, with or without 
#  modification, are permitted provided that the following conditions are met: 
#
#  (1) Redistributions of source code must retain the above copyright notice, 
#  this list of conditions and the following disclaimer. 
#
#  (2) Redistributions in binary form must reproduce the above copyright 
#  notice, this list of conditions and the following disclaimer in the 
#  documentation and or other materials provided with the distribution. 
#
#  (3) Neither the name of the University of California, Lawrence Berkeley 
#  National Laboratory, U.S. Dept. of Energy nor the names of its contributors 
#  may be used to endorse or promote products derived from this software 
#  without specific prior written permission. 
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
#  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
#  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
#  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
#  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
#  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
#  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
#  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
#  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
#  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
#  POSSIBILITY OF SUCH DAMAGE. 

"""Sequence IO for NBRF/PIR format.

The format is similar to fasta. The header line consistins of '>', a two-
letter sequence type (P1, F1, DL, DC, RL, RC, or XX), a semicolon, and a
sequence ID. The next line is a textual description of the sequence, 
followed by one or more lines containing the sequence data. The end of 
the sequence is marked by a "*" (asterisk) character.

type_code -- A map between NBRF two letter type codes and Alphabets.


see:  http://www.cmbi.kun.nl/bioinf/tools/crab_pir.html

--- Example NBRF File ---

>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR 
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH 
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ 
  SDVPERSIPI TREEKPAIAG AQRK*

>P1;CRAB_BOVIN
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDIAIHHPWI RRPFFPFHSP SRLFDQFFGE HLLESDLFPA STSLSPFYLR 
  PPSFLRAPSW IDTGLSEMRL EKDRFSVNLD VKHFSPEELK VKVLGDVIEV 
  HGKHEERQDE HGFISREFHR KYRIPADVDP LAITSSLSSD GVLTVNGPRK 
  QASGPERTIP ITREEKPAVT AAPKK*

"""

from corebio.utils import *
from corebio.seq import *
from corebio.seq_io import *

names = ("nbrf", "pir",)
extensions = ('nbrf', 'pir', 'ali')




type_code = {
    'P1' : protein_alphabet,   # Protein (complete)
    'F1' : protein_alphabet,   # Protein (fragment)
    'DL' : dna_alphabet,       # DNA (linear)
    'DC' : dna_alphabet,       # DNA (circular)
    'RC' : rna_alphabet,       # RNA (linear)
    'RL' : rna_alphabet,       # RNA (circular)
    'N3' : rna_alphabet,       # tRNA
    'N1' : rna_alphabet,       # other functional RNA
    'XX' : generic_alphabet
    }

def read(fin, alphabet=None):  
    """Read and parse a NBRF seqquence file. 

    Args:
        fin -- A stream or file to read
        alphabet -- The expected alphabet of the data. If not supplied, then
                an appropriate alphabet will be inferred from the data.
    Returns: 
        SeqList -- A list of sequences
    Raises: 
        ValueError -- If the file is unparsable        
    """
    seqs = [ s for s in iterseq(fin, alphabet)]
    return SeqList(seqs)


        
def iterseq(fin, alphabet=None):
    """ Generate sequences from an NBRF file.
    
    arguments:
        fin -- A stream or file to read
        alphabet --    
    yields :
        Seq
    raises :
        ValueError -- On a parse error.
    """
        
    body, header,sequence = range(3) # Internal states
    
    state = body
    seq_id = None
    seq_desc = None
    seq_alpha = None
    seqs = []

    for lineno, line in enumerate(fin) :
        if state == body :
            if line == "" or line.isspace() :
                continue
            if line[0] == '>':
                seq_type, seq_id = line[1:].split(';')
                if alphabet :
                    seq_alpha = alphabet
                else :
                    seq_alpha = type_code[seq_type]
                state = header
                continue
            raise ValueError("Parse error on line: %d" % lineno)

        elif state == header :
            seq_desc = line.strip()
            state = sequence
            continue
                          
        elif state == sequence :
            data = "".join(line.split()) # Strip out white space
            if data[-1] =='*' :
                # End of sequence data
                seqs.append(data[:-1])   

                seq = Seq( "".join(seqs), name = seq_id.strip(), 
                    description = seq_desc, alphabet = seq_alpha)

                yield seq
                state= body
                seq_id = None
                seq_desc = None
                seqs = []
                continue
            else :
                seqs.append(data)
                continue
        else :       
            # If we ever get here something has gone terrible wrong
            assert(False)

    # end for


