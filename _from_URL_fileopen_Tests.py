import unittest

class _from_URL_fileopen_Tests(unittest.TestCase):


    def test_URLscheme(self):
        """test for http, https, or ftp scheme"""
        broken_url = "file://foo.txt"
        self.assertRaises(ValueError, _from_URL_fileopen, (broken_url))
                         
    def test_URLfile_extension(self):
        """tests if file extension from URL is accepted format"""
        broken_text_ext = "https://www.fakebox.com/x/xxxxxxxxxx/LICENSE.pdf"
        self.assertRaises(ValueError, _from_URL_fileopen, (broken_text_ext))



