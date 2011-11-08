 
#  Copyright (c) 2005 Gavin E. Crooks <gec@threeplusone.com>
#
#  This software is distributed under the MIT Open Source License.
#  <http://www.opensource.org/licenses/mit-license.html>
#
#  Permission is hereby granted, free of charge, to any person obtaining a 
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
#  THE SOFTWARE.
#

""" Read and write the CLUSTAL sequence file format.
 
See :
- http://www.cmpharm.ucsf.edu/~goh/Treecorr/sampleAlignment.html
- http://www.bioperl.org/wiki/ClustalW_multiple_alignment_format
  
Ref :
-   Higgins D., Thompson J., Gibson T., Thompson J.D., Higgins D.G., Gibson 
    T.J. (1994). CLUSTAL W: improving the sensitivity of progressive multiple
    sequence alignment through sequence weighting, position-specific gap
    penalties and weight matrix choice. Nucleic Acids Res. 22:4673-4680.
"""

# TODO: What happens if CLUSTAL is not the first line of the file?


import re

from corebio.utils import *
from corebio.seq import *
from corebio.seq_io import *

__all__ = ('example', 'names', 'extensions', 'read')

example = """
CLUSTAL W (1.81) multiple sequence alignment


CXCR3_MOUSE       --------------------------LENSTSPYDYGENESD-------FSDSPPCPQDF
BLR_HUMAN         --------------------------LENLEDLF-WELDRLD------NYNDTSLVENH-
CXCR1_HUMAN       --------------------------MSNITDPQMWDFDDLN-------FTGMPPADEDY
CXCR4_MURINE      -----------------------------------YTSDN---------YSGSGDYDSNK
                                                     :  :          :..     ..
 
CXCR3_MOUSE       -SL-------NFDRTFLPALYSLLFLLGLLGNGAVAAVLLSQRTALSSTDTFLLHLAVAD
BLR_HUMAN         --LC-PATMASFKAVFVPVAYSLIFLLGVIGNVLVLVILERHRQTRSSTETFLFHLAVAD
CXCR1_HUMAN       -SPC-MLETETLNKYVVIIAYALVFLLSLLGNSLVMLVILYSRVGRSVTDVYLLNLALAD
CXCR4_MURINE      -EPC-RDENVHFNRIFLPTIYFIIFLTGIVGNGLVILVMGYQKKLRSMTDKYRLHLSVAD
                             :.  .:   * ::** .::**  *  ::   :   * *: : ::*::**

CXCR3_MOUSE       VLLVLTLPLWAVDAA-VQWVFGPGLCKVAGALFNINFYAGAFLLACISFDRYLSIVHATQ
BLR_HUMAN         LLLVFILPFAVAEGS-VGWVLGTFLCKTVIALHKVNFYCSSLLLACIAVDRYLAIVHAVH
CXCR1_HUMAN       LLFALTLPIWAASKV-NGWIFGTFLCKVVSLLKEVNFYSGILLLACISVDRYLAIVHATR
CXCR4_MURINE      LLFVITLPFWAVDAM-ADWYFGKFLCKAVHIIYTVNLYSSVLILAFISLDRYLAIVHATN
                  :*:.: **: ...     * :*  ***..  :  :*:*.. ::** *:.****:****..
"""



names = ("clustal", "clustalw",)
extensions = ('aln',)


header_line = re.compile(r'(CLUSTAL.*)$')

# (sequence_id) (Sequence) (Optional sequence number)
seq_line   = re.compile(r'(\s*\S+\s+)(\S+)\s*(\d*)\s*$')

# Saved group includes variable length leading space.
# Must consult a seq_line to figure out how long the leading space is since
# the maximum CLUSTAL ids length (normally 10 characters) can be changed.
match_line = re.compile(r'([\s:\.\*]*)$') 


def iterseq(fin, alphabet=None):
    """Iterate over the sequences in the file."""
    # Default implementation
    return iter(read(fin, alphabet) )


def read(fin, alphabet=None) :  
    alphabet = Alphabet(alphabet)      
    seq_ids = []
    seqs = []
    block_count = 0
    data_len =0

    for token in _scan(fin):
        if token.typeof== "begin_block":
            block_count = 0
        elif token.typeof == "seq_id":
            if len(seqs) <= block_count :
                seq_ids.append(token.data)
                seqs.append([])
        elif token.typeof == "seq":
            if not alphabet.alphabetic(token.data) :
                raise ValueError(
                    "Character on line: %d not in alphabet: %s : %s" % (
                    token.lineno, alphabet, token.data) )
            seqs[block_count].append(token.data)
            if block_count==0 :
                data_len = len(token.data) 
            elif data_len != len(token.data) :
                raise ValueError("Inconsistent line lengths")
                
            block_count +=1

          
    seqs = [ Seq("".join(s), alphabet, name= i) for s,i in zip(seqs,seq_ids)]
    return SeqList(seqs)


# 1) The word "CLUSTAL" should be the first word on the first line of the file.
# (But sometimes isn't.)
# 2) The alignment is displayed in blocks of fixed length.
# 3) Each line in the block corresponds to one sequence.
# 4) Each sequence line starts with a sequence name followed by at least one
#     space and then the sequence.

def _scan( fin ):
    """Scan a clustal format MSA file and yield tokens.
        The basic file structure is
            begin_document
                header?     
               (begin_block
                   (seq_id seq seq_index?)+
                   match_line?
               end_block)*
            end_document     
    
        Usage:
        for token in scan(clustal_file):
            do_something(token)
    """
    header, body, block = range(3)
    
    yield Token("begin")
    leader_width = -1
    state = header
    for L, line in enumerate(fin):
        if state==header :
            if line.isspace() : continue
            m = header_line.match(line)
            state = body
            if m is not None :
                yield Token("header", m.group() )
                continue
            # Just keep going and hope for the best.
            #else :
                #raise ValueError("Cannot find required header")
                
                
        
        if state == body :
            if line.isspace() : continue
            yield Token("begin_block")
            state = block
            # fall through to block
        
        if state ==  block:
            if line.isspace() :
                yield Token("end_block")
                state = body
                continue
            
            m = match_line.match(line)
            if m is not None :
                yield Token("match_line", line[leader_width:-1])
                continue
     
            m = seq_line.match(line) 
            if m is None: 
                raise ValueError("Parse error on line: %d (%s)" % (L,line))
            leader_width = len(m.group(1))
            yield Token("seq_id", m.group(1).strip() )
            yield Token("seq", m.group(2).strip() )
            if m.group(3)  :
                yield Token("seq_num", m.group(3)) 
            continue

        # END state blocks. If I ever get here something has gone terrible wrong
        raise RuntimeError()
    
    if state==block:
         yield Token("end_block")
    yield Token("end")     
    return

def write(fout, seqs) :
    """Write 'seqs' to 'fout' as text in clustal format"""
    header = "CLUSTAL W (1.81) multiple sequence alignment"
    name_width = 17
    seq_width = 60
    
    print >>fout, header
    print >>fout
    print >>fout 
    
    L = 0
    for s in seqs: L = max(L, len(s))
    
    for block in range(0, L, seq_width):
        for s  in seqs :
            start = min(block, len(s))
            end = min( start+seq_width, len(s))
            print  >>fout, s.name.ljust(name_width),
            print  >>fout, s[start:end]
        print  >>fout 
        
        
        



