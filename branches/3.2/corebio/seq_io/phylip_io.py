#!/usr/bin/env python

#  Copyright (c) 2005 David D. Ding <dding@berkeley.edu>
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

"""Reads Sequences in interleaved Phylip format (not sequential) and returns a
list of sequences. Phylip is a very common phylogeny generating sequence type
that has the following traits:
1) First line contains number of species and number of characters in a species'
sequence. Options may follow, and they can be spaced or unspaced. Options are
simply letters such as A and W after the number of characters.
2) Options don't have to contain U in order for a usertree to appear.
3) If there are options then options appear first, then the sequences. For the
first iteration of sequences the first ten spaces are reserved for names of
options and species, the rest is for sequences.
4) For the second and following iterations the names are removed, only
sequence appears
4) At end of file a usertree may appear. First there is a number that indicts
the number of lines the usertree will take, and then the usertrees follow.

Examples:
  6   50   W
W         0101001111 0101110101 01011	
dmras1    GTCGTCGTTG GACCTGGAGG CGTGG	
hschras   GTGGTGGTGG GCGCCGGCCG TGTGG
ddrasa    GTTATTGTTG GTGGTGGTGG TGTCG
spras     GTAGTTGTAG GAGATGGTGG TGTTG
scras1    GTAGTTGTCG GTGGAGGTGG CGTTG
scras2    GTCGTCGTTG GTGGTGGTGG TGTTG

0101001111 0101110101 01011	
GTCGTCGTTG GACCTGGAGG CGTGG	
GTGGTGGTGG GCGCCGGCCG TGTGG
GTTATTGTTG GTGGTGGTGG TGTCG
GTAGTTGTAG GAGATGGTGG TGTTG
GTAGTTGTCG GTGGAGGTGG CGTTG
GTCGTCGTTG GTGGTGGTGG TGTTG

1					
((dmras1,ddrasa),((hschras,spras),(scras1,scras2)));


"""

from corebio.seq import *

names = ( 'phylip',)
extensions = ('phy',)

def iterseq(fin, alphabet=None):
    """Iterate over the sequences in the file."""
    # Default implementation
    return iter(read(fin, alphabet) )


#Read takes in a phylip file name, reads it, processes it, and returns a SeqList
def read(fin, alphabet=None):
    
   
    sequence=[] #where sequences are stored
    idents=[]
    num_seq=0
    num_total_seq=0 #length of sequence of 1 species
    tracker=0 #track what sequence the line is on
    usertree_tracker=0 #track usertree lines
    options='' #options
    num_options=0 #number/lens of options - U

    line=fin.readline()
    while line:
        s_line=line.split() #for ease of use, not used in all scenarios, but easier on the eye

        if s_line == []: #see nothing do nothing
            pass

        elif (s_line[0].isdigit() and len(s_line) == 1 and len(sequence)==num_seq and len(sequence[0])==num_total_seq):    #identifies usertree
            usertree_tracker = int(s_line[0])
            pass
        
        elif num_options > 0:
            if len(sequence) < num_seq:
                if s_line[0][0] in options:
                    num_options -= 1
                    pass
                else:
                    raise ValueError('Not an option, but it should be one')
            else:
                num_options -= 1
                pass
    
        elif usertree_tracker > 0:                    #bascally skip usertree
            if len(sequence[num_seq-1]) == num_total_seq:
                usertree_tracker -=1
                pass
            else:
                raise ValueError('User Tree in Wrong Place')
    
        #####problems parse error unexpected
        elif s_line[0].isdigit():
            if len(s_line) >= 2 and len(sequence) == 0: #identifies first line of file
                num_seq = int(s_line[0])               #get number of sequences
                num_total_seq = int(s_line[1])         #get length of sequences
                if len(s_line) > 2:                   #takes care of the options
                    options= (''.join(s_line[2:]))
                    num_options=len(options) - options.count('U')
            else:
                raise ValueError('parse error') 


    #when options end, this takes care of the sequence
        elif num_options == 0:
            if (num_seq==0):
                raise ValueError("Empty File, or possibly wrong file")
            elif tracker < num_seq:
                if num_seq > len(sequence):
                    sequence.append(''.join(line[10:].split())) #removes species name
                    idents.append(line[0:10].strip())
                    tracker +=1
                    
                else:
                    sequence[tracker] += (''.join(s_line))      
                    tracker +=1
                    
                if tracker == num_seq:
                    tracker = 0
                    num_options = len(options)-options.count('U')
        
        line=fin.readline()
   
    if len(sequence) != len(idents) or len(sequence)!=num_seq:
        raise ValueError("Number of different sequences wrong") 
    
    seqs = []
    for i in range (0, len(idents)):
        if len(sequence[i])==num_total_seq:
            seqs.append(Seq(sequence[i], alphabet, idents[i]))
        else:
            raise ValueError("extra sequence in list")
    
    return SeqList(seqs)


        
    
    
   
