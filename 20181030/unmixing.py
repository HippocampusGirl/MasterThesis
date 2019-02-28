#!/usr/bin/env python3

from multiprocessing import Pool
from multiprocessing.sharedctypes import RawArray
import numpy as np
from scipy.optimize import nnls
from skimage.external.tifffile import (
    TiffFile,
    TiffWriter
)
# import matplotlib.pyplot as plt
import os
from os import path as op
import sys
from glob import glob

# A = np.loadtxt("A10_14_e.txt")
A = np.loadtxt("A11HB_16_e.txt")
np.set_printoptions(suppress = True, linewidth = 200)
print(A)

fname = "Images/Alexa350Alexa405GFPAlexa555Alexa647FM4-64FX_0.tif"

# basedir = "E:\masterThesis\data\Tiling"
# fnames = glob(op.join(basedir, "**", "img.tif"), recursive=True)

fnames = [fname]

for fname in fnames:
    with TiffFile(fname) as im:
        arr = im.asarray()
        tags = im.pages[0].tags

    sizeC = arr.shape[0]

    sizeXY = arr.shape[1:]
    sizeP = int(np.prod(sizeXY))

    y = RawArray("f", sizeP*sizeC)
    yarr = np.frombuffer(y, dtype = np.float32).reshape((sizeP, sizeC))

    arrt = arr.T.reshape((sizeP, sizeC))
    np.copyto(yarr, arrt, casting = "safe")
    del arr, arrt

    sizeD = A.shape[0]

    x = RawArray("f", sizeP*sizeD)
    xarr = np.frombuffer(x, dtype = np.float32).reshape((sizeP, sizeD))

    chunksize = 1 << 18

    def init_worker(x_, y_, A_, chunksize_, sizeP_):
        global x, y, A, chunksize, sizeC, sizeD, sizeP
        x = x_
        y = y_
        A = A_
        sizeC = A.shape[1]
        sizeD = A.shape[0]
        sizeP = sizeP_

    def unmix(i):
        global x, y, A, chunksize, sizeC, sizeD, sizeP
        print("Chunk %d" % i)
        yarr = np.frombuffer(y, dtype = np.float32).reshape((sizeP, sizeC))
        xarr = np.frombuffer(x, dtype = np.float32).reshape((sizeP, sizeD))
        for j in range(chunksize * i, min(chunksize * (i + 1), sizeP)):
            xarr[j, :] = nnls(A.T, yarr[j, :])[0]

    with Pool(initializer = init_worker, initargs = (x, y, A, chunksize, sizeP)) as pool:
        chunks = (sizeP + chunksize - 1) // chunksize
        pool.map(unmix, range(chunks))

    xarrt = xarr.reshape(sizeXY[::-1]+(sizeD,)).T * 4096

    with TiffWriter(op.join(op.dirname(fname), "unmix.tif"), imagej = True) as im:
        im.save(xarrt, photometric = "minisblack")
        # im.save(xarr, compress = 4, photometric = "minisblack", extratags = tags)




