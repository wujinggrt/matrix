#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import time

a0 = np.mat(np.random.randint(0, 10, (100, 500)))
a1 = np.mat(np.random.randint(0, 10, (500, 10)))
a2 = np.mat(np.random.randint(0, 10, (10, 75)))
a3 = np.mat(np.random.randint(0, 10, (75, 55)))
a4 = np.mat(np.random.randint(0, 10, (55, 10)))
a5 = np.mat(np.random.randint(0, 10, (10, 100)))
a6 = np.mat(np.random.randint(0, 10, (100, 74)))
a7 = np.mat(np.random.randint(0, 10, (74, 45)))
a8 = np.mat(np.random.randint(0, 10, (45, 200)))
a9 = np.mat(np.random.randint(0, 10, (200, 100)))
a10 = np.mat(np.random.randint(0, 10, (100, 200)))
a11 = np.mat(np.random.randint(0, 10, (200, 700)))

start = time.process_time()
for i in range(1000):
    tmp = a0 * a1 * a2 * a3 * a4 * a5 * a6 * a7 * a8 * a9 * a10 * a11
end = time.process_time()

print((end - start) / 1000.)