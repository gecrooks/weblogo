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




"""Read the output of a fasta sequence similarity search. 

FASTA is a DNA and Protein sequence alignment software package first described
by David J. Lipman and William R. Pearson in 1985. In addition to rapid 
heuristic search methods, the FASTA package provides SSEARCH, an implementation
of the optimal Smith Waterman algorithm.

The module can parse the output from fasta, ssearch and other search programs
in the fasta collection. It will parse both default ('-m 1') and compact 
('-m 9 -d 0') output. 
    
Refs:
    ftp.virginia.edu/pub/fasta
    http://en.wikipedia.org/wiki/FASTA
"""



from corebio.utils import Reiterate, Token, isblank
from corebio.ssearch_io import Report, Result, Hit, Annotation, Alignment
from math import floor
import re

__all__ = 'read'

_rangere = re.compile(r"\((\d+)-\d+:(\d+)-\d+\)")


def read(fin) :
    """Read and parse a fasta search output file.
    
    returns: a list of Results
    """    
    scanner = _scan(fin) 

    report = None
    result = None
    hit = None
    #query_seq = None
    #target_seq = None
    alignment_num = 0
    
    
    for token in scanner :
        #print token
        typeof = token.typeof
        value = token.data
        
        if typeof == 'begin_report' :
            report = Report()
        elif   typeof == 'algorithm' : 
            report.algorithm = value
        elif typeof == 'algorithm_version' : 
            report.algorithm_version = value
        elif typeof == 'algorithm_reference' : 
            report.algorithm_reference = value
        elif typeof == 'database_name' : 
            report.database_name = value
        elif typeof == 'database_letters' : 
            report.database_letters = value
        elif typeof == 'database_entries' : 
            report.database_entries = value
        elif typeof == 'end_report' :
            # Final sanity checking
            break
        elif typeof == 'parameter' : 
            key = value[0]
            value = value[1]
            report.parameters[key] = value           
        
        elif typeof == 'begin_result' :
            result = Result()
            report.results.append(result)            
            
        elif typeof == 'query_name' :
            result.query.name = value
        elif typeof == 'query_description' :
            result.query.description = value
        elif typeof == 'end_result' :
            pass
            
        elif typeof == 'begin_hit' :
            hit = Hit()
        elif typeof == 'target_name' :
            hit.target.name = value
        elif typeof == 'target_description' :
            hit.target.description = value
        elif typeof == 'target_length' :
            hit.target.length = value
        elif typeof == 'raw_score' :
            hit.raw_score = value
        elif typeof == 'bit_score' :
            hit.bit_score = value
        elif typeof == 'significance' :
            hit.significance = value               
        elif typeof == 'end_hit' :
            result.hits.append(hit)
            hit = None
            
        elif typeof == 'begin_alignment' :
            alignment = Alignment()
            tseq = []
            qseq = []
        elif typeof == 'end_alignment' :
            tseq = ''.join(tseq)
            qseq = ''.join(qseq)
            L = max (len(tseq), len(qseq) )
            tseq = tseq.ljust(L).replace(' ', '.')
            qseq = qseq.ljust(L).replace(' ', '.')
            alignment.query_seq = tseq
            alignment.target_seq = qseq
            result.hits[alignment_num].alignments.append(alignment)    
            alignment_num+=1                   
            tseq = None
            qseq = None
        elif typeof == 'target_seq' :
            tseq += value
        elif typeof == 'query_seq' :
            qseq += value
        elif typeof == 'alignment_raw_score' :
            alignment.raw_score = value

        elif typeof == 'alignment_bit_score' :
            alignment.bit_score = value
        elif typeof == 'alignment_significance' :
            alignment.significance = value
        elif typeof == 'alignment_length' :
            alignment.length = value
        elif typeof == 'alignment_similar' :
            alignment.similar = value
        elif typeof == 'alignment_identical' :
            alignment.identical = value
        elif typeof == 'alignment_query_start' :
            alignment.query_start = value
        elif typeof == 'alignment_target_start' :
            alignment.target_start = value

        else: 
            # Should never get here.
            raise RuntimeError("Unrecoverable internal parse error (SPE)")
            pass


    return report
# End method read()


def _scan(fin) :

    def next_nonempty(i) :
        L = i.next()
        while L.strip() == '':  L = i.next()
        return L
    

    lines = Reiterate(iter(fin))
    try :    
    
        yield Token("begin_report", lineno= lines.index())
        
        # find header line : "SSEARCH searches a sequence data bank"
        L = lines.next()
        
        if L[0] == '#' :
            yield Token("parameter", ("command", L[1:].strip()), lines.index())
            L = lines.next()
            
        while not L : L= lines.next()
        algorithm = L.split()[0]
        expected = [ "SSEARCH", "FASTA","TFASTA","FASTX",
                        "FASTY","TFASTX","TFASTY"]          
        if algorithm not in expected:
            raise ValueError("Parse failed: line %d" % lines.index() ) 
        yield Token ("algorithm", algorithm, lines.index() )
              
        # Next line should be the version
        L = lines.next()
        if not L.startswith(" version") : 
            raise ValueError("Parse failed: Cannot find version.")
        yield Token( "algorithm_version", L[8:].split()[0].strip(), lines.index())
        
        # Algorithm reference        
        L = lines.next()
        if not L.startswith("Please cite:") : 
            raise ValueError("Parse failed: Expecting citation" + L)
        cite = lines.next().strip() + ' ' + lines.next().strip()            
        yield Token( "algorithm_reference", cite)

        # Find line "searching testset.fa library"
        L = lines.next()
        while not L.startswith("searching") : L = lines.next()
        yield Token("database_name", L[10:-8], lines.index() )
        
        # Results 
        L = lines.next()
        while isblank(L) : L = lines.next()
        if ">>>" not in L :
            raise ValueError("Parse failed on line %d: " % lines.index())

        while ">>>" in L :
            yield Token("begin_result", lineno= lines.index())
            index = L.find('>>>')
            (name, description) = L[index+3:].split(' ',1) 
            yield Token("query_name", name, lines.index())
            yield Token("query_description", description, lines.index())

            while not L.startswith("The best scores are:") :
                L = lines.next()
            L = lines.next()
            # hits
            while not isblank(L) :
                lineno = lines.index()
                desc = L[0:49]
                yield Token("begin_hit", lineno= lineno)
                yield Token("target_description", desc, lineno, 0)
                yield Token("target_name",  desc.split(' ',1)[0], lineno, 0)
                yield Token("target_length", int(L[52:56]), lineno, 52)
                fields = L[57:].split()
                raw, bit, sig  = fields[0], fields[1], fields[2]
                #print raw, bit, sig 
                yield Token("raw_score",    float(raw), lineno, 57)
                yield Token("bit_score",    float(bit), lineno)
                yield Token("significance", float(sig), lineno)
                yield Token("end_hit", lineno=lineno)
                L = lines.next()
    
            # Optimal alignment information
            L = next_nonempty(lines)
            #print ">>>", L, L.startswith('>>')
            while L.startswith('>>'):
                if  L.startswith('>>>') : break
                
                yield Token("begin_alignment", lineno=lines.index() )

                #          1         2         3         4
                #01234567890123456789012345678901234567890123456789
                # s-w opt:  46  Z-score: 70.7  bits: 18.5 E():  3.6
                L = lines.next()
                fields = L.split()
                raw, bit, sig = fields[2], fields[6], fields[8]
                yield Token("alignment_raw_score",    float(raw), lineno)
                yield Token("alignment_bit_score",    float(bit), lineno)
                yield Token("alignment_significance", float(sig), lineno)

                #Smith-Waterman score: 46;  38.095% identity (71.429% similar) in 21 aa overlap (2-22:36-56)
                L = lines.next()
                lineno = lines.index()
                fields = L.split()
                assert( len(fields) ==12)
                alen = int(fields[8])
                identical = int( floor(0.5+alen* float(fields[3][:-1])/100.))
                similar = int( floor(0.5+alen* float(fields[3][:-1])/100.))
                yield Token("alignment_length", alen, lineno)
                yield Token("alignment_similar", similar, lineno)
                yield Token("alignment_identical", identical, lineno)                
                
                m = _rangere.match( fields[11])
                assert (m is not None)
                yield Token("alignment_query_start", int(m.group(1))-1, lineno)
                yield Token("alignment_target_start", int(m.group(2))-1, lineno) 
                
                
                count = 1                
                while True:
                    L = lines.next()
                    count += 1


                    
                    if L.startswith('>>'): break
                    if '>>>' in L:
                        lines.push(L)
                        break
                    if 'residues' in L and 'sequences' in L :
                        lines.push(L)
                        break
                    if not L or L[0].isspace() : continue
  
  
                    # there are 2 lines before the first query sequence (but
                    # we started the count at 1). There is 1 line between query
                    # and target, 3 lines between target and query, unless the
                    # query ends before the ends and the target wraps onto another
                    # Then there are two lines between target and target.
                    
# Smith-Waterman score: 34;  35.294% identity ...
# 
#           30        40        50        60         70              
# d1pfsa EGFLHLEDKPHPLQCQFFVESVIPAGSYQVPYRINVNNG-RPELAFDFKAMKRA      
#                                      : . . .:: .:  .::             
# d8rxna                 MKKYVCTVCGYEYDPAEGDPDNGVKPGTSFDDLPADWVCPVCGA
#                                10        20        30        40    
# 
# d8rxna PKSEFEAA
#            50                      
                                        
                    lineno=lines.index()
                    if count==4 :
                        yield Token("query_seq", L[7:].rstrip(),  lineno)
                    else : 
                        yield Token("target_seq", L[7:].rstrip(),lineno)
                    count = 0
                    
                yield Token("end_alignment", lineno=lines.index() )
            yield Token("end_result", lineno= lines.index())
            L = next_nonempty(lines)
        # End results
        
        # "13355 residues in 93 query   sequences"
        # "13355 residues in 93 library sequences"
        #print '>>', L
        LL = L.split()
        yield Token("database_letters",int(LL[0]), lines.index() )
        yield Token("database_entries", int(LL[3]), lines.index() )
    
        yield Token("end_report", lineno= lines.index())
    except StopIteration :
        raise ValueError("Premature end of file ")


    


