#!/usr/bin/python
#This process analyze is from RX HFA
import matplotlib
from matplotlib.pyplot import cohere
#from scipy.linalg.fblas import snrm2
matplotlib.use('Agg')
import stuffr
import matplotlib.pyplot as plt
import numpy
import gdf
import os
import math
import glob
from matplotlib import dates
import datetime
import time
import cPickle as pickle

import h5py
import traceback
from optparse import OptionParser

import sys



## scan directories and analyze
def analyze_dirs_ric(dirn="/data0", proc_folder='', freq="", delete_old=False,old_threshold=4.0,phase_cal=0.0,reanalyze=False, max_N_analyze=30,inc_int=0,profiles_number=600,stationtx_codes='0,2'):
    print "NEW ANALYZE DIRS"

    d = dirn
    print d
    t_now = time.time()
    t_created = os.path.getctime(d)
    data_age = (t_now-t_created)/3600.0/24.0
    print "created: %s age %1.2f days" % (time.ctime(t_created), data_age)

    an_len= (profiles_number/100)*1000000

    try:
        print "This try"
        ch0 = freq.channel_dirs[0]
        ch1 = freq.channel_dirs[1]

        analyze_all2(dirn0="%s/%s"%(d,ch0), dirn1="%s/%s"%(d,ch1), rawdata_doy_folder=dirn, proc_folder=proc_folder, freq=freq,an_len=an_len,phase_cal=phase_cal, reanalyze=reanalyze, max_N_analyze=max_N_analyze,inc_int= inc_int,stationtx_codes=stationtx_codes)

    except:
	print "ESTOY AQUI"
        print "This except"
        print "Error processing %s"%d

    """ FIXME: Check this option """
    if delete_old:
        if data_age > old_threshold:
            # delete *.gdf files and mv directory to old subdirectory to avoid further processing.
            #
            cmd = "find %s/0 -name \*.gdf |sed -e 's/.*/rm \\0/ge'" % (d)
            os.system(cmd)
            cmd = "find %s/1 -name \*.gdf |sed -e 's/.*/rm \\0/ge'" % (d)
            os.system(cmd)
            os.system("mv %s %s/old/"%(d,dirn))




## dual channel
def analyze_all2(dirn0,
                 dirn1,
                 rawdata_doy_folder,
                 proc_folder,
                 freq,
                 idx0=0,
                 i0=0,
                 i1=None,
                 an_len=6000000,
                 clen=10000,
                 station=0,
                 reanalyze=False,
                 thresh_phase=0.5,
                 threshv=0.1,
                 rfi_rem=True,
                 dcblock=False,
                 phase_cal=0.0,
                 Nranges=1000,
                 max_N_analyze=0,
		 inc_int=0,  #Max number of hdf5 files to create each bucle
		 stationtx_codes='0,2'):

    #codigos="0,1,2"
    # codigos_test="0,2"
    #print "codigos_test",codigos_test
    codigos=stationtx_codes###################### SELECCION DE CODIGO COMO PARAMETROS 0,1,2 Y  EN EL ORDEN QUE ME PROVOQUE
    print "codigos",codigos
    splitcodes=codigos.split(',')
    codeList=[]
    for eachsplitcode in splitcodes:
        codeList.append(eachsplitcode)
    #####HASTA AQUI TENGO LA LISTA DE CODIGO A DECODIFICAR :D
    os.system("mkdir -p %s/cspec"%(dirn0))
    #os.system("mkdir -p %s/%s"%(proc_folder, freq.procdata_folder))
    #Create just the neccesary folders

#    print "test 1"
    for code in enumerate(codeList):
        os.system("mkdir -p %s/%s"%(proc_folder, freq.procdata_folder.replace('sp01_','sp%s1_'%(code[1]))))

    if reanalyze:
        os.system("rm %s/cspec/res*.hdf5"%(dirn0))
    g0 = gdf.new_gdf(dirn0,dtype=numpy.float32,itemsize=8)
    g1 = gdf.new_gdf(dirn1,dtype=numpy.float32,itemsize=8)

    #print g0['file_list']
    print type(g0['file_list'])
    g0_list = g0['file_list']
    #print g1['file_list']

    if i1 == None:
        i1 = math.floor(g0['max_n']/an_len)-1
        print i1
    N = i1

    print "==> Total gdf files founded: %d"%(N)
    Ns = an_len/clen
    if i0 == 0:
        if not reanalyze:
            try:
                #i0 = numpy.max(stuffr.load_object("%s/cspec/index.bin"%(dirn0)),0) #6/9/17 cambie index
		#print dirn0 , "DIRECTORIO"
                i0 = numpy.max(stuffr.load_object("%s/cspec/index.bin"%(dirn0))+1,0)
            except:
                print "No index file."

    #Check if rawdata_end_flag exists and if I've already reached the total of hdf5 processed files
    print "VALOR INC_INT: ",inc_int ## INTEGRACIONES INCOHERENTES 6 (60seg) o 0 (10seg)##aqui el parametro debe venir de arriba
    thrshold = inc_int-1
    if inc_int==0:
	thrshold=1
    if os.path.isfile("%s/rawdata_end_flag"%(rawdata_doy_folder)) and i0+thrshold >= N:
        print "DOY FOLDER: %s ALREADY COMPLETE"%(rawdata_doy_folder)
        print "Index bin mark is %d and N: %d"%(i0, N)
        print "Creating file %s/%s"%(rawdata_doy_folder,  freq.end_flag)
        os.system("touch %s/%s"%(rawdata_doy_folder,  freq.end_flag))
        return 0
    print i0,'test2'
    print N,'test3'
    if i0+1 >= N:
        print "Nothing to process yet !!!"
	time.sleep(60)#300
        return 5


    """ En caso no se indique se analizaran todos los archivos disponibles """
    if max_N_analyze == 0:
        None
    else:
        # The program process no more than 5 min of data per cycle (N-i0 equal to 30 is 5m=300s)
        if N-i0 > max_N_analyze:
            N = i0 + max_N_analyze

    if inc_int != 0 :
        print "Int Inc: %d "%(inc_int)
        save_int_data = [1,1,1]
        dic={} # temporal dictionary to save spectralmoments.
    print "i0: %d"%(i0)
    print "N: %d"%(N)
    start = 0
    plot_deco_power = 0
    for i in numpy.arange(i0,N):
        print "%d/%d"%(i,N-1)
	print "i","VALOR DE i: ......",i
        sttime=time.time()
        print g0_list[10*int(i)*(an_len/1000000)],""
        temp_a = g0_list[10*int(i)*(an_len/1000000)]
        mark = int((temp_a.split('/')[-1]).split('.')[-2][-12:])
        print "mark: %d"%(mark)
        for code in enumerate(codeList):

            if inc_int != 0:
                a0 = analyze_prc(dirn=g0,idx0=an_len*i,an_len=an_len,Nranges=Nranges,station=int(code[1]),rfi_rem=rfi_rem) #
                a1 = analyze_prc(dirn=g1,idx0=an_len*i,an_len=an_len,Nranges=Nranges,station=int(code[1]),rfi_rem=rfi_rem)
                s0= numpy.abs(a0['spec'])**2.0
                s1= numpy.abs(a1['spec'])**2.0
                cspec01=a0['spec']*numpy.conjugate(a1['spec'])

                if save_int_data[code[0]]==(inc_int):
                    os.system("rm -f %s/%s/spec-%012d.hdf5"%(proc_folder, freq.procdata_folder.replace('sp01_','sp%s1_'%(code[1])),mark))
                    res = h5py.File("%s/%s/spec-%012d.hdf5"%(proc_folder, freq.procdata_folder.replace('sp01_','sp%s1_'%(code[1])),mark))
                    #res['ch0_C%s'%(code[1])] = dic["ch0_spec_C%s"%(code[1])]/inc_int
                    #res['ch1_C%s'%(code[1])] = dic["ch1_spec_C%s"%(code[1])]/inc_int
                    res['pw0_C%s'%(code[1])] = dic["s0_C%s"%(code[1])]/inc_int # power spectra density ch0
                    res['pw1_C%s'%(code[1])] = dic["s1_C%s"%(code[1])]/inc_int # power spectra density ch1
                    res['cspec01_C%s'%(code[1])] = dic["cspec01_C%s"%(code[1])]/inc_int # Cross spectra ch0.ch1*
                    res['dc0_C%s'%(code[1])] =dic["dc0_C%s"%(code[1])]/inc_int
                    res['dc1_C%s'%(code[1])] =dic["dc1_C%s"%(code[1])]/inc_int
                    res['i'] = i
                    res['t'] = g0['t0']+float(i)*an_len/100e3
                    res['image0_C%s'%(code[1])] = spec2imgcolors(dic["s0_C%s"%(code[1])]/inc_int,threshv=threshv)
                    res['image1_C%s'%(code[1])] = spec2imgcolors(dic["s1_C%s"%(code[1])]/inc_int,threshv=threshv)
                    res['version'] = '3'
                    res.close()
                    print '------Saving data for code %s...%d/%d'%(code[1],save_int_data[code[0]],inc_int)
                    if code[0]==(len(codeList)-1): #Just save at the last code
                        stuffr.save_object(i,"%s/cspec/index.bin"%(dirn0))

                    save_int_data[code[0]]=1
                else:
                    if save_int_data[code[0]]==1:
                       # dic["ch0_spec_C%s"%(code[1])]=a0['spec']
                       # dic["ch1_spec_C%s"%(code[1])]=a1['spec']
                        dic["s0_C%s"%(code[1])]=s0
                        dic["s1_C%s"%(code[1])]=s1
                        dic["dc0_C%s"%(code[1])]=a0['spec'][0,:]
                        dic["dc1_C%s"%(code[1])]=a1['spec'][0,:]
                        dic["cspec01_C%s"%(code[1])]=cspec01

                    else:
                       # dic["ch0_spec_C%s"%(code[1])]+=a0['spec']
                       # dic["ch1_spec_C%s"%(code[1])]+=a1['spec']
                        dic["s0_C%s"%(code[1])]+=s0
                        dic["s1_C%s"%(code[1])]+=s1
                        dic["dc0_C%s"%(code[1])]+=a0['spec'][0,:]
                        dic["dc1_C%s"%(code[1])]+=a1['spec'][0,:]
                        dic["cspec01_C%s"%(code[1])]+=cspec01


                    print 'Integrating data for code %s... %d/%d'%(code[1],save_int_data[code[0]],inc_int)
                    save_int_data[code[0]]+=1
                    #time.sleep(1)
            else:
			print 'Code{0}:', code[0] #extract the string format
		        print 'Code{1}:', code[1] #extract the string format
		        a0 = analyze_prc(dirn=g0,idx0=an_len*i,an_len=an_len,Nranges=Nranges,station=int(code[1]),rfi_rem=rfi_rem)
		        a1 = analyze_prc(dirn=g1,idx0=an_len*i,an_len=an_len,Nranges=Nranges,station=int(code[1]),rfi_rem=rfi_rem)
		        s0= numpy.abs(a0['spec'])**2.0
		        s1= numpy.abs(a1['spec'])**2.0
		        cspec01=a0['spec']*numpy.conjugate(a1['spec'])
			#print "VOY A GUARDAR"
		        os.system("rm -f %s/%s/spec-%012d.hdf5"%(proc_folder, freq.procdata_folder.replace('sp01_','sp%s1_'%(code[1])),mark))
		        res = h5py.File("%s/%s/spec-%012d.hdf5"%(proc_folder, freq.procdata_folder.replace('sp01_','sp%s1_'%(code[1])),mark))
		        res['pw0_C%s'%(code[1])] = s0 # power spectra density ch0
		        res['pw1_C%s'%(code[1])] = s1 # power spectra density ch1
		        res['cspec01_C%s'%(code[1])] = cspec01 # Cross spectra ch0.ch1*
		        res['dc0_C%s'%(code[1])] = a0['spec'][0,:]
		        res['dc1_C%s'%(code[1])] = a1['spec'][0,:]

		        res['i'] = i
		        res['t'] = g0['t0']+float(i)*an_len/100e3
		        res['image0_C%s'%(code[1])] = spec2imgcolors(s0,threshv=threshv)#AnADI image0
			res['image1_C%s'%(code[1])] = spec2imgcolors(s1,threshv=threshv)#AnADI ESTA LINEA
		        res['version'] = '3'
		        res.close()
		        #print '------Saving data for code %s...%d/%d'%(code[1],save_int_data[code[0]],inc_int)
		        if code[0]==(len(codeList)-1): #Just save at the last code
		            stuffr.save_object(i,"%s/cspec/index.bin"%(dirn0))

            #if start == 0 and plot_deco_power == 1:
            #    coded=a0['res'][50,:]*numpy.conjugate(a0['res'][50,:])
            #    spect=a0['spec']
            #    for ii in numpy.arange(1000):
            #        spect[:,ii] = numpy.fft.ifft(spect[:,ii])
            #    #volts=numpy.fft.ifft(spect)
            #    power=spect[50,:]*numpy.conjugate(spect[50,:])
            #    #print 'plotting some files for debugg.'
            #    plt.clf()
            #    plt.plot(numpy.log10(numpy.real(coded))*10.0,'k')
            #    plt.plot(numpy.log10(numpy.real(power))*10.0,'--r')
            #    plt.ylabel('Red is from spectra')
            #    plt.savefig("/home/igp-114/Documents/DecoPowerPlots/HFTtest/%d.png"%(mark))
            #    start=start+1
            #    print 'plotting some files for debugg.'

        print 'time taked> ', time.time()-sttime



#
# Interferometry:
# colorsys.hsv_to_rgb(hue,sat,val)
# colorsys.hls_to_rgb
#
#
def analyze_prc(dirn="",idx0=0,an_len=6000000,clen=10000,station=0,Nranges=1000,rfi_rem=True,cache=True):
    g = []
    if type(dirn) is str:
        g = gdf.new_gdf(dirn,dtype=numpy.float32,itemsize=8)
    else:
        g = dirn
    ####################################################################3
    '''

    AQUI ES DONDE UTILIZA EL CODIGO Y SE OBTIENE EL HDF5 A PARTIR DEL  GDF

    '''
    ##########################################################################
#    print "STATION:", station
    code = stuffr.create_pseudo_random_code(len=clen,seed=station)
    #print "CODE",code
    N = an_len/clen
    res = numpy.zeros([N,Nranges],dtype=numpy.complex64)
    Z = numpy.zeros([clen,N],dtype=numpy.complex64)
    r = stuffr.create_estimation_matrixNEWH(code=code,cache=cache,rmax=Nranges)
    #raw_input('2 gardenias')
    if station == 0:
        B = r['H'][:,0:10000]
    if station == 1:
        B = r['H'][:,10000:20000]
    if station == 2:
        B = r['H'][:,20000:30000]
    #print "STATION:", station,B
    spec = numpy.zeros([N,Nranges],dtype=numpy.complex64)
    # update 17/02/2020
    clip = 2.0e-5 # interference threshold
    #make spectral window
    window = numpy.hanning(N)# have N samples in each range gate
    freq = numpy.fft.fftfreq(N)

    for i in numpy.arange(N):
        z = gdf.read_vec(g,idx0+i*clen,clen)
        Z[:,i] = z

    Z = numpy.clip(Z,-clip,clip) # remove interference
    res = numpy.transpose(numpy.dot(B,Z))

    for i in numpy.arange(Nranges):
        spec[:,i] = numpy.fft.fft(window*res[:,i]) #update spec

    #print station

    #It was commented by M Milla
    #if rfi_rem:
    #    median_spec = numpy.zeros(N,dtype=numpy.float32)
    #    for i in numpy.arange(N):
    #        median_spec[i] = numpy.median(numpy.abs(spec[i,:]))
    #    for i in numpy.arange(Nranges):
    #            spec[:,i] = spec[:,i]/median_spec[:]
    ret = {}
    ret['res'] = res
    ret['spec'] = spec
    return(ret)


def spec2imgcolors(s,threshv=0.1):
    L = s.shape[0]
    N = s.shape[1]
    # | g |  b  |  r  | g |
    i0l = math.floor(L*threshv)
    i0h = L-math.floor(L*threshv)
    im = math.floor(L/2)
    colm = numpy.zeros([N,3],dtype=numpy.float32)
    for ri in numpy.arange(N):
        colm[ri,1] = numpy.sum(s[0:i0l,ri]) + numpy.sum(s[i0h:L,ri])
        colm[ri,2] = numpy.sum(s[i0l:im,ri])
        colm[ri,0] = numpy.sum(s[im:i0h,ri])
    return(colm)

def phase2imgcolors(s,thresh_phase=0.2,phase_cal=0.0):
    L = s.shape[0]
    N = s.shape[1]
    # | g |  b  |  r  | g |
    ph = numpy.angle(s*numpy.exp(1.0j*phase_cal))
    mag = numpy.abs(s)
    colm = numpy.zeros([N,3],dtype=numpy.float32)
    for ri in numpy.arange(N):
        ph_row = ph[:,ri]
        mag_row = mag[:,ri]
        colm[ri,0] = numpy.sum(mag_row[numpy.where(ph_row < -1.0*thresh_phase)])
        colm[ri,2] = numpy.sum(mag_row[numpy.where(ph_row > thresh_phase)])
        colm[ri,1] = numpy.sum(mag_row)-colm[ri,0]-colm[ri,2]

    return(colm)


def max_phase(s):
    L = s.shape[0]
    N = s.shape[1]
    # | g |  b  |  r  | g |
    ph = numpy.angle(s)
    mag = numpy.abs(s)
    res = numpy.zeros([N],dtype=numpy.float32)
    for ri in numpy.arange(N):
        res[ri] = ph[numpy.argmax(mag[:,ri]),ri]
    return(res)
