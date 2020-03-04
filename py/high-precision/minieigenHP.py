##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

"""
When yade uses high-precision number as ``Real`` type the usual (old) ``from minieigen import *``
has to be replaced with:

.. code-block:: python

	from yade.minieigenHP import *

This command ensures backward compatibility between both. It is then guaranteed that python uses the same number of decimal places as yade is using everywhere else.

Please note that used precision can be very arbitrary, because ``cpp_bin_float``
or ``mpfr`` take it as a :ref:`compile-time argument <highPrecisionReal>`. Hence such ``yade.minieigenHP`` cannot be separately precompiled as a package.
Though it could be precompiled for some special types such as ``boost::multiprecision::float128``.
"""

# all C++ functions are accessible now:
import yade
from yade._minieigenHP import *

