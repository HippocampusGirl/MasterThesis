#!/usr/bin/env python

print("mixing11")

import numpy as np
from skimage.external.tifffile import TiffFile
# import matplotlib.pyplot as plt
import os

import uuid

import sys

import tensorflow as tf
tf.logging.set_verbosity(4)
# import tensorflow_probability as tfp
# tf.enable_eager_execution()

##  _  _  __  __  __  __    __   
## ( \( )(  )(  )(  \/  )  /__\  
##  )  (  )(__)(  )    (  /(__)\ 
## (_)\_)(______)(_/\/\_)(__)(__)

from multiprocessing import cpu_count

from ctypes import CDLL

libnuma = CDLL("libnuma.so.1")

n_numanodes = libnuma.numa_num_task_nodes()
n_cpuspernumanode = cpu_count() // n_numanodes

config = tf.ConfigProto()
config.intra_op_parallelism_threads = n_cpuspernumanode 
config.inter_op_parallelism_threads = n_numanodes 
config.log_device_placement = True

# os.environ["MPI_OPTIMAL_PATH"] = "1"

# print(n_cpuspernumanode)

os.environ["OMP_NUM_THREADS"] = str(n_cpuspernumanode)
os.environ["MKL_NUM_THREADS"] = str(n_cpuspernumanode)

os.environ["OMP_DYNAMIC"] = "FALSE"
os.environ["MKL_DYNAMIC"] = "FALSE"

os.environ["OMP_NESTED"] = "FALSE"
os.environ["OMP_PLACES"] = "threads"
os.environ["OMP_PROC_BIND"] = "close"

os.environ["KMP_BLOCKTIME"] = "0"

for k, v in os.environ.items():
    if "MP" in k:
        print("%s = %s" % (k, v))

##  ____    __   ____   __   
## (  _ \  /__\ (_  _) /__\  
##  )(_) )/(__)\  )(  /(__)\ 
## (____/(__)(__)(__)(__)(__)

# args = ["0b0000101", "small_Alexa405.tif",
#         "0b0001001", "small_Alexa488.tif",
#         "0b0010001", "small_Alexa555.tif",
#         "0b0100001", "small_Alexa647.tif",
#         "0b1000001", "small_FM4-64FX.tif",
#         "0b1111111", "small_Alexa350Alexa405GFPAlexa555Alexa647FM4-64FX.tif"]

args = ["0b0000001", "Images/Background_0.tif",
        "0b0000101", "Images/Alexa405_0.tif",
        "0b0001001", "Images/Alexa488_0.tif",
        "0b0010001", "Images/Alexa555_0.tif",
        "0b0100001", "Images/Alexa647_0.tif",
        "0b0111011", "Images/Alexa350Alexa488Alexa555Alexa647_0.tif",
        "0b0111011", "Images/Alexa350Alexa488Alexa555Alexa647_1.tif",
        "0b0111101", "Images/Alexa405GFPAlexa555Alexa647_0.tif",
        "0b1000101", "Images/Alexa405FM4-64FX_0.tif",
        "0b1000101", "Images/Alexa405FM4-64FX_1.tif",
        "0b1000101", "Images/Alexa405FM4-64FX_2.tif",
        "0b1000001", "Images/FM4-64FX_0.tif",
        "0b1000001", "Images/FM4-64FX_1.tif",
        "0b1000001", "Images/FM4-64FX_2.tif",
        "0b1111111", "Images/Alexa350Alexa405GFPAlexa555Alexa647FM4-64FX_0.tif",
        "0b1111111", "Images/Alexa350Alexa405GFPAlexa555Alexa647FM4-64FX_1.tif"
]
        
sizeD = len(args[0]) - 2

argpairs = zip(args[::2], args[1::2])    

P = None
Y = None

sizeC = 0

skip = 1.75 #int(np.random.randint(2, 6)) 
# 
# print(skip)

for i, (b, fname) in enumerate(argpairs):
    print(fname)
    sys.stdout.flush()
    im = TiffFile(fname)
    arr = im.asarray()
    sizeC = arr.shape[0]
    y = arr.T.reshape((-1, sizeC)).astype(np.float64)
    nchoice = int(y.shape[0] / skip) #* 2 
    y = y[np.argsort(np.sum(np.square(y), axis = 1))[-nchoice:], :]
    # y = y[np.random.choice(y.shape[0], nchoice // 2, replace = False), :]
    sizeP = y.shape[0]
    
    p = np.zeros((sizeP, sizeD), dtype = bool)
    
    b = int(b, 0)
    for i in range(sizeD):
        if b & (1 << i) != 0:
            p[:, i] = 1
            
    if P is None:
        P = p
    else :
        P = np.vstack((P, p))
    
    if Y is None:
        Y = y
    else :
        Y = np.vstack((Y, y))
        
    im.close()

Y = Y.astype(np.float64)

# normalize ##
# M = Y.mean(axis = 0) 
# Y /= M

sizeP = Y.shape[0]
sizeX = P.sum()
sizeA = sizeD*sizeC

print("sizeC %d" % sizeC)
print("sizeD %d" % sizeD)
print("sizeP %d" % sizeP)

##  __  __  _____  ____  ____  __   
## (  \/  )(  _  )(  _ \( ___)(  )  
##  )    (  )(_)(  )(_) ))__)  )(__ 
## (_/\/\_)(_____)(____/(____)(____)
# 

pind = tf.constant(np.vstack(np.where(P)).T, dtype = tf.int64, name = "indP")
    
# def lsq(x):
# np.loadtxt("Ainv_theoretical.txt").astype(np.float64) / 4096
# tf.abs(tf.random_normal(shape = (sizeC, sizeD,), stddev = 1 / Y.std(), dtype = np.float64))
A = tf.Variable(
    tf.abs(tf.random_normal(shape = (sizeD, sizeC,), stddev = Y.std(), dtype = np.float64)), 
    dtype = np.float64, name = "A")
x = tf.Variable(
    tf.random_uniform(shape = (sizeX,), maxval = 1.0, dtype = np.float64), 
    dtype = np.float64, name = "x")

# A = tf.reshape(x[:sizeA], (sizeD, sizeC,), name = "reshapeA")

clA = tf.clip_by_value(A, 1, np.inf, name = "clipA")

clx = tf.clip_by_value(x, 0, np.inf, name = "clipX")
X = tf.scatter_nd(pind, clx, (sizeP, sizeD), name = "scatterX")

Yhat = tf.matmul(X, clA, name = "calc_Yhat")

loss = tf.sqrt(tf.reduce_mean(tf.squared_difference(Y, Yhat)))
    
    # return loss, tf.gradients(loss, x)[0]

# optim_results = tfp.optimizer.bfgs_minimize(
#     lsq, initial_position = start, tolerance = 1e-8,
#     parallel_iterations = 4)

# # 
optimizer = tf.contrib.opt.ScipyOptimizerInterface(loss, 
    var_to_bounds = {A: (1, np.infty), x: (0, 1.0)},
    method = "L-BFGS-B", options = {"maxiter": 500, "disp": True})

# optimizer = tf.contrib.opt.ScipyOptimizerInterface(loss, 
#     var_to_bounds = {A: (-np.infty, np.infty), x: (0, np.infty)},
#     method = "TNC", options = {"maxiter": 10000, "disp": True})
    
# train_step = tf.train.AdamOptimizer(learning_rate = 0.001, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-08).minimize(loss)
# train_step = tf.train.MomentumOptimizer(learning_rate = 0.005, momentum = 0.5).minimize(loss)

# iterations = 4000

def loss_callback(A = None):
    print(A)
    sys.stdout.flush()

np.set_printoptions(suppress = True, linewidth = 200)
with tf.Session(config = config) as session:
    # results = session.run(optim_results)
    # print (sess.run(normal_rv))
    session.run(tf.global_variables_initializer())
    A_value = session.run(clA)
    print(A_value)
    # Ainv_value = np.linalg.pinv(A_value)
    # print(Ainv_value)
    sys.stdout.flush()
    
    # min_loss_value = np.inf
    # for i in range(iterations):
    #     _, loss_value, iteration_A_value = session.run([train_step, loss, clA])  
    #     print("Iteration %4d/%4d %f" % (i + 1, iterations, loss_value))  
    #     if i == 10 and loss_value > 400:
    #         sys.exit(0)
    #     if loss_value < min_loss_value:
    #         min_loss_value = loss_value
    #         A_value = iteration_A_value
    #         print(A_value)
    # 
    #     sys.stdout.flush()
    optimizer.minimize(session, fetches = [A], loss_callback = loss_callback)

    A_value = session.run(clA)
    loss_value = session.run(loss)
    ofname = "A11HB_%d_%s.txt" % (int(loss_value), str(uuid.uuid4()).split("-")[0])
    print(ofname)
    np.savetxt(ofname, A_value)
    print(A_value)
    sys.stdout.flush()








