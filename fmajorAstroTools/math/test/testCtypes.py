from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np
import ctypes
c_int = ctypes.c_int
c_int_p = ctypes.POINTER(c_int)
c_int_pp = ctypes.POINTER(c_int_p)
c_int_ppp = ctypes.POINTER(c_int_pp)

a = np.arange(24).reshape((2,3,4))
a = a.astype(np.int, order="C")
#a_c = a.ctypes.data_as(c_int_ppp)
a_c = np.ctypeslib.as_ctypes(a)

print(a_c)
for i in range(2):
    print(a_c[i])
    for j in range(3):
        print(a_c[i][j])
        for k in range(4):
            print(a_c[i][j][k])
