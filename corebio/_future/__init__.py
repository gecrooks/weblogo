
"""
Private compatability module for running under python version 2.3.
    
Replacement for 
  o  string.Template    -- introduced in python 2.4
  o  subprocess         -- introduced in python 2.4
  
  o  resource_string    -- introduced in pkg_resource of setuptools
  o  resource_stream
  o  resource_filename

from string import Template -> from corebio._future import Template

"""


try :
    import pkg_resources
except ImportError :
    pkg_resources = None


try :
    from string import Template
except ImportError :
    from _string import Template




def resource_string( modulename, resource, basefilename = None):
    """Locate and return a resource as a string.
    >>> f = resource_string( __name__, 'somedatafile', __file__)
    """
    if pkg_resources : 
        return pkg_resources.resource_string(modulename, resource)    
    
    f = resource_stream( modulename, resource, basefilename)
    return f.read()
    
def resource_stream( modulename, resource, basefilename = None):
    """Locate and return a resource as a stream.
    >>> f = resource_stream( __name__, 'somedatafile', __file__)
    """    
    if pkg_resources :  
        return pkg_resources.resource_stream(modulename, resource)    
    
    return open( resource_filename( modulename, resource, basefilename) )

def resource_filename( modulename, resource, basefilename = None): 
    """Locate and return a resource filename.
    >>> f = resource_filename( __name__, 'somedatafile', __file__)

    A resource is a data file stored with the python code in a package.
    All three resource methods (resource_string,  resource_stream,
    resource_filename) call the corresponding methods in the 'pkg_resources'
    module, if installed. Otherwise, we resort to locating the resource
    in the local filesystem. However, this does not work if the package
    is located inside a zip file. 
    """ 
    if pkg_resources : 
        return pkg_resources.resource_filename(modulename, resource)            
    
    if basefilename is None :
        raise NotImplementedError(
            "Require either basefilename or pkg_resources")
    
    import os
    return os.path.join(os.path.dirname(basefilename), resource)







