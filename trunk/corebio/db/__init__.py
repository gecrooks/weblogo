#  Copyright (c) 2006, John X. Gilman

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

"""Database access. In biological parlance, a 'database' is an organized  
collection of biological data. 

    -- Dbxref - A database cross-reference, a resource identifier.
    -- DataSource - A source of biological data.
    -- databases - A list of standard databases.

Usage:
    ref = Dbxref("Swiss-Prot:P12345")
    data = ref.

"""

import re
import urllib

from corebio.data import data_stream
from corebio.utils import update, stdrepr




__all__ = ['Dbxref','DataSource', 'default_registry', 
            'astral', 'scop' ]

_dbxref = re.compile(r'\s*(?:/db_xref\s*=\s*\")?(\w+):(\w+)\"?\s*')

class Dbxref(object) :

    """A database cross-reference, a kind of Universal Resource Name (URN) used 
    to locate information in biological databases.

    -- database  - the name of the database 
    -- identifier - the internal identifier of the related information
                    according to the naming conventions of the cross-referenced
                    database.

    Examples:        

    cross reference to GDB identifier:            /db_xref="GDB:39999"   
    cross reference to Swiss-Prot entry:          /db_xref="Swiss-Prot:P12345" 

    Refs:
        http://www.ncbi.nlm.nih.gov/collab/db_xref.html
        http://www.ebi.ac.uk/cgi-bin/emblfetch
        http://www.expasy.org/cgi-bin/lists?dbxref.txt
    """
    
    __slots__ = ["database", "identifier"]
    
    def __init__(self, ref, identifier=None) :
        """ Dbxref('/db_xref="GDB:39999"')
            Dbxref("GDB:39999")
            Dbxref('GDB','39999')
        """
        if identifier:
            self.database = ref
            self.identifier = identifier
            return
        m = _dbxref.match(ref)
        if m is None : 
            raise ValueError, "Cannot parse database cross-reference."
        self.database = m.group(1)
        self.identifier = m.group(2)
    
    def __repr__(self) :
        return "Dbxref(%s, %s)" %(self.database, self.identifier)
    
    def __str__(self) :
        return "%s:%s" %(self.database, self.identifier)
        
    def datasource(self, registry = None) :
        if registry is None: registry = default_registry
        database = self.database.lower()
        if not (database in registry) :
            raise ValueError, "Unknown database: %s" % self.database
        return registry[database]
        
    def data_url(self) :
        return self.datasource().data_url(self.identifier)    
        
    def data_stream(self) :
        return self.datasource().data_stream(self.identifier)            
        
    def data_string(self) :
        return self.datasource().data_string(self.identifier)    
     
     
class DataSource(object) :
    """A collection of indexed biological data. 
    
    Attributes:
        -- name - The name of the database.
        -- abbrev - Standard abrreviation
        -- alt_abbrev - A list of alternative abrreviations or names.
        -- url - Main URL of the database.
        -- resource_url - The URL of individual database resources, e.g.    
                    'biocyc.org/getid?id=%s', where '%s' is replaced by the
                    resource identifier.
        -- parser - A function that can convert raw data into an object.
        -- description - Description of the type of biological data in this
             database. E.g. ''Organism-specific gene databases'
    """
    __slots__=['abbrev', 'alt_abbrev', 'name', 
                'url', 'resource_url', 'parser', 'description']
    
    def __init__(self, **kwarg) :
        for s in self.__slots__ : setattr(self, s, None)
        update(self, **kwarg)

    def __repr__(self) : 
        return stdrepr(self)
      
    def data_url(self, identifier) :    
        url = self.resource_url
        if url is None : 
            raise ValueError, "Cannot resolve resource."
        url = url % identifier # Fixme
        return url
        
    def data_stream(self, identifier) :
        url = self.data_url(identifier) 
        return urllib.urlopen(url)
 
    def data_string(self, identifier) :
        return self.data_stream(identifier).read()
    
    def data_object(self, identifier) :
        parser =  self.parser
        if parser == None :    
            raise ValueError, "No parser for resource."
        stream = self.data_stream(identifier)
        return parser(stream)
  
  

def _build_default_registry() :
    registry = dict()
    
    # Parse database specifications from standard 'dbxref' file.
    ds = _read_dbxref( data_stream('dbxref.txt') )

    ds.append(DataSource(
        name= 'Swiss-Prot',
        abbrev = 'swissprot',
        alt_abbrev = ('uniprotkb', 'uniprot'),
        url = 'http://ca.expasy.org/sprot/',
        resource_url = 'http://www.expasy.ch/cgi-bin/get-sprot-raw.pl?%s',
        description='Protein knowledgebase',
        ) )
            
    for d in ds :
        registry[d.abbrev.lower() ] = d
        if d.alt_abbrev :
            for alt in d.alt_abbrev : 
                registry[alt.lower()] = ds                 
    return registry


def _read_dbxref(stream) :
    """ Read the dbxref file and return a list of datasources.""" 
    namemap = dict(Abbrev = 'abbrev', 
                    Name = 'name', 
                    Server = 'url',
                    Db_URL = 'resource_url',
                    Cat    = 'description'
                    )
                    
    ds_list = []
    attr = None
    ds = None

    for line in stream :
        if len(line)<8 : continue
        typeof = line[0:6].strip()
        content = line[8:].strip()

        if typeof=='Abbrev' : # New record
            ds = DataSource()
            ds_list.append(ds)
        
        if content == 'None' : continue # No data
        
        if ds is None : continue # No record started
        
        if typeof == '' : # Continuation line
            if not attr is None : 
                content = getattr(ds, attr) +' '+ content
                setattr(ds, attr, content)
        else :
            attr = None
            if typeof in namemap:
                attr = namemap[typeof]
                if attr == 'resource_url' or attr=='url' : 
                    content = ''.join(('http://',content))
                    
                setattr(ds, attr, content)
       
    return ds_list


default_registry = _build_default_registry()
