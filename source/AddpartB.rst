

Welcome to Add part B
=====================

Definition
----------

term 3 : classifier
    Definition 3.

term 4 : classifier one : classifier two
    Definition 4.


Usage
-----

*Usage:* 

``python *.py -i INPUT -o OUTDIR``

- aaaaaa

  + a111
  + a222

- bbbbbb

  * b1
  * b2

- cccccc

  # ccccc
  # ccc

usage: 

``aviris_tcld_py5.py [-h] -i INPUT -o OUTDIR``

This code is for $cloud masking

optional arguments:

  -h, --help
                   show this help message and exit

  -i INPUT, --input INPUT   Input file name

  -o OUTDIR, --outdir OUTDIR  Output dir

cloud masking
--------------

*drmvgbmbhmhnmmnjnmjmn*


Example
-------

``python *.py -i INPUT -o OUTDIR``

Library
-------

.. code-block:: python
  :linenos:

  import numpy as np
  def get_envi_header_dict(hdr):
    #Get all "key = {val}" type matches
    regex=re.compile(r"^(.+?)\s*=\s*{\n((?:.|\n)*?)}$",re.M|re.I)#

    matches = regex.findall(hdr)

    #Remove them from the header
    subhdr2=regex.sub('',hdr)
    #Get all "key = val" type matches
    regex=re.compile(r'^(.+?)\s*=\s*(.*?)$',re.M|re.I)

    matches.extend(regex.findall(subhdr2))

    return dict(matches)
