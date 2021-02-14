

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

  -d OUTDIR, --outdir OUTDIR  Output dir

  -ss  cc aa bb

  simple table 

  -aa   bb

  -cc   dd
  
  simple table2   

    -sw   cc   aa sb 

    -sw   cc aa   sb inline Eq :math:`a^2 + b^2 = c^2`

    -ss   cc aa bb   nn

    -ss   cc aa bb nn

  grid table

  :math:`a^2 + b^2 = c^2`
  
  grid table

  =====  =====  =======
    A      B    A and B
  =====  =====  =======
  False  False  False
  True   False  False
  False  True   False
  True   True   True
  =====  =====  =======


cloud masking
--------------

  ``dedef``

*drmvgbmbhmhnmmnjnmjmn*

.. image::./ownimg/f130612_flightline.png
    :scale: 50%
    :align: center
    :alt: alternate text

.. math:: 
  :nowrap:

  \begin{align}
    y    & =  ax^2+bx+c \\
    f(x) & =  x^2+2xy+y^2
  \end{align}

  \begin{align*}
    y    & = ax^2 +bx + c \\
    f(x) & = x^2 +2xy + y^2
  \end{align*}

Example
-------

``python *.py -i INPUT -o OUTDIR``

Library
-------

.. code-block:: python
  :linenos:
  :dedent: 0

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
