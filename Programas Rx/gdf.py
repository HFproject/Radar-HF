import glob
import numpy
import math
import datetime
import time
import pickle
import os
import astropy.io.ascii

debug = False

def new_gdf(dirn,dtype="<i2",itemsize=4):
    result = {}
    files =  glob.glob("%s/*.gdf"%(dirn))
    files.sort()
    files = files[1:] #In order to avoid the first elememt, which can be incompleted 
    result['file_size'] = os.path.getsize(files[0])/itemsize
    result['file_list'] = files
    result['max_n'] = result['file_size']*len(files)
    result['dtype'] = dtype
    result['cache'] = numpy.zeros(result['file_size'],dtype=numpy.complex64)
    result['cache_idx'] = -1
    result['scale'] = 1.0
    result['re_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2
    result['im_idx'] = numpy.arange(result['file_size'],dtype=numpy.int64)*2+1
    tstamp_file = open("%s/timestamps.log"%(dirn),'r')
    result['t0'] = float(tstamp_file.readline().split(" ")[1])
    tstamp_file.close()
    if debug:
        print("Read gdf dir. Nfiles=%d Fsize=%d t0 %1.2f"%(len(result['file_list']),result['file_size'],result['t0']))
    return(result)

def read_vec(gdf, idx, length):
    idx = int(idx)
    length = int(length)
    files = gdf['file_list']
    f0_idx = int(math.floor(idx/gdf['file_size']))
    fn_idx = int(math.floor((idx+length-1)/gdf['file_size']))+1
    res_vec = numpy.zeros([length],dtype=numpy.complex64)
    if debug:
        print "f0 ",f0_idx," ",fn_idx
    c0_idx = idx % gdf['file_size']
    c1_idx = (idx+length-1) % gdf['file_size']+1
    if debug:
        print c0_idx," ",c1_idx
    n_read = 0
    for f_idx in range(f0_idx,fn_idx):
        c0 = 0
        c1 = gdf['file_size']
        if f_idx == f0_idx:
            c0 = c0_idx
        if f_idx+1 == fn_idx:
            c1 = c1_idx
        if debug:
            print "open file ",f_idx,"\n"

        if gdf['cache_idx'] != f_idx:
            a = numpy.fromfile(files[f_idx],dtype=gdf['dtype'])
            gdf['cache'] = numpy.array(a[gdf['re_idx']]*gdf['scale'] + 1j*a[gdf['im_idx']]*gdf['scale'],dtype=numpy.complex64)
            gdf['cache_idx'] = f_idx
        if debug:
            print len(res_vec)," c0 ",c0," c1 ",c1
        res_vec[numpy.arange(c1-c0)+n_read] = gdf['cache'][c0:c1]
        n_read = n_read + (c1-c0)
        if debug:
            print "indices ",c0,"-",c1," n_read ",n_read
    return(res_vec)