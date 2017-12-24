from __future__ import (absolute_import, division, print_function, unicode_literals)
import numpy as np

a = np.arange(120).reshape((4,5,6))
b = a.T

z, m, n = a.shape
for eachM in range(m):
    for eachN in range(n):
        print(a[:, eachM, eachN])
        assert np.all(a[:, eachM, eachN]==b[eachN, eachM, :])
