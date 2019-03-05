===========
WebLogo API
===========

.. contents:: :local:
.. currentmodule:: weblogo


To create a logo in python code:

.. code-block:: python

    >>> from weblogo import *
    >>> fin = open('cap.fa')
    >>> seqs = read_seq_data(fin)
    >>> logodata = LogoData.from_seqs(seqs)
    >>> logooptions = LogoOptions()
    >>> logooptions.title = "A Logo Title"
    >>> logoformat = LogoFormat(logodata, logooptions)
    >>> eps = eps_formatter(logodata, logoformat)

