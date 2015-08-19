import os

import urllib #for urlretrieve 
import urllib2
from urlparse import urlparse, urlunparse


#will need to install pyopenssl ndg-httpsclient pyasn1 to local machine to avoid SSL error

#weblogo license.txt file share url
target_url = 'https://www.dropbox.com/s/qgfto1ipfnxh0pq/LICENSE.txt?dl=0'

#google target_url weblogo_changelog.txt share url
#target_url = "https://drive.google.com/file/d/0B_QbDSPaPwJzVVh2ZW1YZmVTOUU/view?usp=sharing"

def _from_URL_fileopen(target_url):
    """allows WebLogo to open files from URL locations"""
   

    #parsing url in component parts
    (scheme, net_location, path, param, query, frag) = urlparse(target_url)

    # checks if string is URL link
    if scheme != "http" and scheme !="https" and scheme != "ftp":
        raise ValueError("Cannot open url: %s" , target_url)

    # checks for dropbox link
    if net_location == 'www.dropbox.com':
        #changes dropbox http link into download link
        if query == "dl=0":    query2 = "dl=1"

        #rebuild download URL, with new query2 variable
            target_url = urlunparse((scheme, net_location, path, param, query2,""))

    # checks for google drive link
    if net_location == 'drive.google.com':
        import shutil, tempfile

        #link configuration for direct download instead of html frame
        google_directdl_frag = "https://docs.google.com/uc?export=download&id="
        
        #pull file id
        (scheme, net_location, path_raw, param, query, frag) = urlparse(target_url)
        path = path_raw.split('/')
        id_file = path[3]

        #rebuild URL for direct download
        target_url = google_directdl_frag + id_file

    # save url to temporary file
    req = urllib2.Request(target_url)
    res = urllib2.urlopen(req)
    temp = tempfile.TemporaryFile()
    shutil.copyfileobj(res, temp)
    temp.seek(0)
    return temp


_from_URL_fileopen(target_url)