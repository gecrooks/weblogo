
#  Copyright (c) 2006 John Gilman
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

""" Parse the output of BLAST and similar sequence search analysis reports. 

The result of a sequence database search is represented by the Report class.
    o Each Report contains one or more results, one for each database query.
    o Each Result contains one or more hits
    o Each Hit may contain one or more Alignments (High scoring Sequence pairs)

CoreBio is often capable of guessing the correct format:
>>> from corebio import ssearch_io
>>> afile = open("test_corebio/data/ssearch/ssearch_out.txt")
>>> report = ssearch_io.read(afile)
>>> print report

Alternatively, each report type has a separate module. Each module defines a
read(fin) method that can parse that report format.

>>> from corebio.ssearch_io import fasta
>>> report = fasta.read( open("test_corebio/data/ssearch/ssearch_out.txt") )
>>> print report

Module          Application            Comments  
---------------------------------------------------------------------------
fasta           FASTA / SSEARCH     Default (-m 1) or compact (-m 9 -d 0)
blastxml        NCBI Blast          NCBI XML format                 

Status: Beta 
"""
# Dev. References :
#   Inspired by Bioperl's searchIO system
#       http://www.bioperl.org/wiki/HOWTO:SearchIO
from __future__ import absolute_import

__all__ = ['read', 'Report', 'Result', 
            'Hit','Annotation', 'Alignment']


from ..utils import stdrepr

def read(fin) :
    """ Read and parse an analysis report. 
    
    returns :
        A database search Report.
    raises :
        ValueError - If the file cannot be parsed
    """

    from . import fasta, blastxml

    parsers = (fasta, blastxml)
    for p in parsers:
        try:    
            return p.read(fin)
        except ValueError as e:
            pass
        fin.seek(0)             # FIXME. Non seakable stdin?
            
    raise ValueError("Cannot parse sequence file: Tried fasta and blastxml")
      


class Report(object) :
    """The results of a database search. The Report contains a list of 1 or more
    Results, one for each query. Each query result contains a list of hits. 
    Each Hit contains a list of HSP's (High scoring segment pairs).
    
    The structure of the report will vary somewhat depending on the source.
    
	 algorithm	          -- e.g. 'BLASTX'
	 algorithm_version	  -- e.g. '2.2.4 [Aug-26-2002]'
     algorithm_reference	  -- 	 
	 database_name	      -- e.g. 'test.fa'
	 database_letters	  -- number of residues in database e.g. 1291	 
	 database_entries	  -- number of database entries

     parameters           -- Dictionary of parameters used in search
	 
	 results              -- A list of lists of Results, one per query
	 """
    __slots__ = ['algorithm', 'algorithm_version', 'algorithm_reference','database_name', 
                'database_letters', 'database_entries', 'parameters', 'results']

    def __init__(self) :
        for name in self.__slots__ : setattr(self, name, None)
        self.parameters = {}
        self.results = []

    def __repr__(self):
        return stdrepr(self)


class Result(object) :
    """ The result from searching a database with a single query sequence.
    
    query        -- Information about the query sequence
    statistics     -- A dictionary of search statistics
    hits         -- A list of Hits
    """
    __slots__ = ['query', 'statistics', 'hits']

    def __init__(self) :
        for name in self.__slots__ : setattr(self, name, None)
        self.query = Annotation() 
        self.statistics = {}
        self.hits = []

    def __repr__(self):
        return stdrepr(self)        

        
class Hit(object) :
    """ A search hit between a query sequence and a subject sequence.
    Each hit may have one or more Alignments
    
    target       -- Information about the target sequence. 
    raw_score	 -- Typically the significance of the hit in bits, e.g. 92.0
    significance -- Typically evalue. e.g '2e-022' 
    alignments   -- A list of alignments between subject and target
    """
    __slots__ =['target', 'raw_score', 'bit_score', 'significance', 
                'alignments']
    def __init__(self) :
        for name in self.__slots__ : setattr(self, name, None)
        self.target      = Annotation()
        self.alignments	 = []
        
    def __repr__(self):
        return stdrepr(self) 

class Annotation(object) :
    """ Information about a subject or query sequence.
    
    name	     -- subject sequence name, e.g. '443893|124775'
    description	 -- e.g.  'LaForas sequence'
    length	     -- subject sequence length, e.g. 331
    locus	     -- e.g. '124775'
    accession	 -- e.g. '443893'
    """
    # Fixme: change into generic sequence annotation class?
    __slots__ = ['name', 'description', 'length', 'locus', 'accession', ]

    def __init__(self):
        for name in self.__slots__ :
            setattr(self, name, None)
             
    def __repr__(self):
        return stdrepr(self) 

class Alignment(object):
    """An alignment between query and subject sequences. 
    For BLAST, these are High scoring Segment pairs (HSPs)
  
    raw_score	     -- Typically significance of the hit in bits, e.g. 92.0
    significance     -- Typically evalue. e.g '2e-022' 

    similar	          -- number of conserved residues #FIXME either frac or num
    identical	      -- number of identical residues
    gaps              -- number of gaps    
    length            -- length of the alignment
    
    query_seq	      -- query string from alignment
    target_seq	      -- hit string from alignment
    mid_seq	          --
    
    query_start       --
    query_frame       --

    target_start      --
    target_frame      --
    
    """
    __slots__ = ['raw_score', 'bit_score', 'significance', 'similar',
     'identical', 'gaps', 'length', 'query_seq', 'target_seq', 'mid_seq',
      'query_start', 'query_frame', 'target_start', 
      'target_frame']
      
    def __init__(self):
        for name in self.__slots__ :
            setattr(self, name, None)
    
    def __repr__(self):
        return stdrepr(self)
        
                
