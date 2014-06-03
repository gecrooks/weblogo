
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


"""Read BLAST XML output.

The DTD is available at
http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.mod.dtd

"""
# See also
# 
# http://bugzilla.open-bio.org/show_bug.cgi?id=1933
#http://portal.open-bio.org/pipermail/biojava-dev/2004-December/002513.html
from __future__ import absolute_import

from . import Report, Result, Hit, Alignment

import xml.sax
from xml.sax.handler import ContentHandler

__all__ = 'read'

def read(fin):
    """Read BLAST xml output and return a list of Result objects.
    """
    parser = xml.sax.make_parser()
    handler = _BlastHandler() 
    parser.setContentHandler(handler)
        
    #To avoid ValueError: unknown url type: NCBI_BlastOutput.dtd
    parser.setFeature(xml.sax.handler.feature_validation, 0)
    parser.setFeature(xml.sax.handler.feature_namespaces, 0)
    parser.setFeature(xml.sax.handler.feature_external_pes, 0)
    parser.setFeature(xml.sax.handler.feature_external_ges, 0)

    try :
        parser.parse(fin)
    except xml.sax.SAXParseException as e :
        raise ValueError("Cannot parse file; " + str(e))
    return handler.report

class _BlastHandler( ContentHandler) :
    def __init__(self):
        """
        """
        ContentHandler.__init__(self)
        self._content = []
        self.report = None
        self._result = None
        self._hit = None
        self._hsp = None

        
    def characters(self, ch):
        self._content.append(ch) 
   
    def startDocument(self):
        self.report = Report()
        
    def endDocument(self) :
        pass
        
    def startElement(self, name, attr):
        if name == 'BlastOutput' :
            pass
        elif name == 'Iteration' :
            result = Result()
            self._result = result
            self.report.results.append(result)
        elif name == 'Parameters' :
            pass
        elif name == 'Statistics' :
            pass
        elif name == 'Hit' :
            self._hit = Hit()
            self._result.hits.append(self._hit)
        elif name == 'Hsp' :
            self._hsp = Alignment()
            self._hit.alignments.append(self._hsp)
        else :
            pass


    def endElement(self, name):
        content = ''.join(self._content).strip()
        self._content = []

        report = self.report
        result = self._result
        hsp = self._hsp
        hit = self._hit
        
        if name == 'BlastOutput' : 
            pass
        elif name == 'BlastOutput_program' :
            report.algorithm = content
        elif name == 'BlastOutput_version' :
            report.algorithm_version = content.split()[1]
        elif name == 'BlastOutput_reference' :
            report.algorithm_reference = content
        elif name == 'BlastOutput_db' :
            report.database_name = content
        elif name == 'BlastOutput_query-ID' : pass
        elif name == 'BlastOutput_query-def' : pass
        elif name == 'BlastOutput_query-len' : pass
        elif name == 'BlastOutput_query-seq' : pass            
        elif name == 'BlastOutput_param' : pass
        elif name == 'BlastOutput_iterations' : pass
        elif name == 'BlastOutput_mbstat' : pass
            
        elif name == 'Iteration' : pass
        elif name == 'Iteration_iter-num' : pass            
        elif name == 'Iteration_query-ID' :  
            result.query.name = content
        elif name == 'Iteration_query-def' :             
            result.query.description = content
        elif name == 'Iteration_query-len' : 
            result.query.length = int(content)            
        elif name == 'Iteration_hits' : pass            
        elif name == 'Iteration_stat' : pass            
        elif name == 'Iteration_message' : pass  
                      
        elif name == 'Parameters' : 
            pass        
        elif name == 'Parameters_matrix' :
            report.parameters['matrix'] = content            
        elif name == 'Parameters_expect' :
            report.parameters['expect'] = content              
        elif name == 'Parameters_include' :
            report.parameters['include'] = content              
        elif name == 'Parameters_sc-match' :
            report.parameters['sc-match'] = content              
        elif name == 'Parameters_sc-mismatch' :
            report.parameters['sc-mismatch'] = content              
        elif name == 'Parameters_gap-open' :
            report.parameters['gap-open'] = content              
        elif name == 'Parameters_gap-extend' :
            report.parameters['gap-extend'] = content              
        elif name == 'Parameters_filter' :
            report.parameters['filter'] = content              
        elif name == 'Parameters_pattern' :
            report.parameters['pattern'] = content              
        elif name == 'Parameters_entrez-query' :
            report.parameters['entrez-query'] = content  

        elif name == 'Statistics' :
            pass              
        elif name == 'Statistics_db-num' :
            result.statistics['db-num'] = int(content)            
        elif name == 'Statistics_db-len' :
            result.statistics['db-len'] = int(content)              
        elif name == 'Statistics_hsp-len' :
            result.statistics['hsp-len'] = int(content)            
        elif name == 'Statistics_eff-space' :
            result.statistics['eff-space'] = float(content)            
        elif name == 'Statistics_kappa' :
            result.statistics['kappa'] = float(content)            
        elif name == 'Statistics_lambda' :
            result.statistics['lambda'] = float(content)          
        elif name == 'Statistics_entropy' :
            result.statistics['entropy'] = float(content)            

        elif name == 'Hit' :
            self._hit = None
        elif name == 'Hit_num' :
            pass            
        elif name == 'Hit_id' :
            hit.target.name = content            
        elif name == 'Hit_def' :
            hit.target.description = content
        elif name == 'Hit_accession' :
            hit.target.accession = content            
        elif name == 'Hit_len' :
            hit.target.length = int(content)             
        elif name == 'Hit_hsps' :
            pass            

        elif name == 'Hsp' :
            self._hsp = None                
        elif name == 'Hsp_num' :
            pass            
        elif name == 'Hsp_bit-score' :
            hsp.bit_score = float(content)            
        elif name == 'Hsp_score' :
            hsp.raw_score = float(content)             
        elif name == 'Hsp_evalue' :
            hsp.significance = float(content)             
        elif name == 'Hsp_query-from' :
            hsp.query_start = int(content) -1           
        elif name == 'Hsp_query-to' :
            #hsp.query_end= int(content)              
            pass
        elif name == 'Hsp_hit-from' :
            hsp.target_start = int(content) -1              
        elif name == 'Hsp_hit-to' :
            #hsp.target_end = int(content)             
            pass
        elif name == 'Hsp_pattern-from' :
            pass            
        elif name == 'Hsp_pattern-to' :
            pass            
        elif name == 'Hsp_query-frame' :
            hsp.query_frame = int(content)              
        elif name == 'Hsp_hit-frame' :
            hsp.target_frame = int(content)            
        elif name == 'Hsp_identity' :
            hsp.identical = int(content)            
        elif name == 'Hsp_positive' :
            hsp.similar = int(content)             
        elif name == 'Hsp_gaps' :
            hsp.gaps = int(content)            
        elif name == 'Hsp_align-len' :
            hsp.length = int(content)             
        elif name == 'Hsp_density' :
            pass            
        elif name == 'Hsp_qseq' :
            hsp.query_seq = content               
        elif name == 'Hsp_hseq' :
            hsp.target_seq = content            
        elif name == 'Hsp_midline' :
            hsp.mid_seq = content     
        else :
            pass       
                


