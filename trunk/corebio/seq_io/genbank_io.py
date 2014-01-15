#!/usr/bin/env python
 

"""Read GenBank flat files. 

Currently only reads sequence data and not annotations.

"""
from __future__ import absolute_import

from ..utils import *
from ..seq import *

  
names = ( 'genbank',)
extensions = ('gb','genbank', 'gbk')



def read(fin, alphabet=None): 
    """Read and parse a file of genbank records. 

    Args:
    fin -- A stream or file to read
    alphabet -- The expected alphabet of the data, if given
    
    Returns: 
    SeqList -- A list of sequences
    
    Raises: 
    ValueError -- If the file is unparsable
    """         
    seqs = [ s for s in iterseq(fin, alphabet)]
    return SeqList(seqs)

    
def iterseq(fin, alphabet=None):
    """ Iterate over genbank records
    
    Args:
    fin -- A stream or file to read
    alphabet -- The expected alphabet of the data, if given    
    
    Yeilds: 
    Seq -- One alphabetic sequence at a time.
    
    Raises: 
    ValueError -- If the file is unparsable
    """
    alphabet = Alphabet(alphabet)

    seq = []
    
    def notblank(string) :
        return not isblank(string)

    lines = Reiterate(iter(fin))
    
    
    while True :
        line = lines.filter( notblank )
        if not line.startswith('LOCUS') :
            raise ValueError(
                "Cannot find start of record at line %d"% lines.index() )

        line = lines.filter(lambda s : s.startswith('ORIGIN') 
                                            or  s.startswith('//') )

        if line.startswith('//') :
            # No sequence data    
            yield Seq( '', alphabet)
        else:
            for line in lines :
                if line.startswith('//') :
                    yield Seq( ''.join(seq), alphabet)
                    seq = []
                    break    
                seq.extend( line.split()[1:] )
       
    
        
    


     
     
     