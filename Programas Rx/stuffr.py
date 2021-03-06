#
# An attempt to translate the main functionality my main
# R radio signal packages gursipr and stuffr to python.
# Nothing extremely complicated, just conveniece functions
#
#
import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import datetime
import time
import pickle

# fit_velocity
#import scipy.constants
#import scipy.optimize

# seed is a way of reproducing the random code without
# having to store all actual codes. the seed can then
# act as a sort of station_id.
def create_pseudo_random_code(len=10000,seed=0):
    numpy.random.seed(seed)
    phases = numpy.array(numpy.exp(1.0j*2.0*math.pi*numpy.random.random(len)),
                         dtype=numpy.complex64)
    return(phases)

def periodic_convolution_matrix(envelope,rmin=0,rmax=100):
    # we imply that the number of measurements is equal to the number of elements in code
    L = len(envelope)
    ridx = numpy.arange(rmin,rmax)
    A = numpy.zeros([L,rmax-rmin],dtype=numpy.complex64)
    for i in numpy.arange(L):
        A[i,:] = envelope[(i-ridx)%L]
    result = {}
    result['A'] = A
    result['ridx'] = ridx
    return(result)

def analyze_prc_file(fname="data-000001.gdf",clen=10000,station=0,Nranges=1000):
    z = numpy.fromfile(fname,dtype=numpy.complex64)
    code = create_pseudo_random_code(len=clen,seed=station)
    N = len(z)/clen
    res = numpy.zeros([N,Nranges],dtype=numpy.complex64)
    idx = numpy.arange(clen)
    r = create_estimation_matrix(code=code,cache=True)
    B = r['B']
    spec = numpy.zeros([N,Nranges],dtype=numpy.float32)

    for i in numpy.arange(N):
        res[i,:] = numpy.dot(B,z[idx + i*clen])
    for i in numpy.arange(Nranges):
        spec[:,i] = numpy.abs(numpy.fft.fft(res[:,i]))
    r['res'] = res
    r['spec'] = spec
    return(r)

B_cache = 0
C_cache = 0
D_cache = 0
H_cache = 0
r_cache = 0
r_cache1= 0
r_cache2= 0
B_cached = False
def create_estimation_matrix(code,rmin=0,rmax=1000,cache=True):
    '''
    global B_cache	
    global C_cache
    global D_cache 
    global r_cache
    global r_cache1
    global r_cache2
    global B_cached
    '''
    if cache == False or B_cached == False:
        r_cache = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
        A = r_cache['A']
        Ah = numpy.transpose(numpy.conjugate(A))
        B_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Ah,A)),Ah)
        r_cache['B'] = B_cache

	code = create_pseudo_random_code(len=10000,seed=1)
	
	r_cache1 = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
	C = r_cache1['A']
	Ch = numpy.transpose(numpy.conjugate(C))
        C_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Ch,C)),Ch)
	r_cache['C'] = C_cache

	code = create_pseudo_random_code(len=10000,seed=2)
	r_cache2 = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
	D = r_cache2['A']
	Dh = numpy.transpose(numpy.conjugate(D))
        D_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Dh,D)),Dh)
	r_cache['D'] = D_cache

        B_cached = True
        return(r_cache)
    else:
        return(r_cache)

def create_estimation_matrixNEWH(code,rmin=0,rmax=1000,cache=True):
    global r_cache
    global r_cache1
    global r_cache2
    global H_cache
    global B_cached
    print code.shape,"1"

    if cache == False or B_cached == False:
        r_cache = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
        A = r_cache['A']
	print A.shape
        #Ah = numpy.transpose(numpy.conjugate(A))
        #B_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Ah,A)),Ah)
        #r_cache['B'] = B_cache

	code = create_pseudo_random_code(len=10000,seed=1)
	
	r_cache1 = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
	C = r_cache1['A']

	#Ch = numpy.transpose(numpy.conjugate(C))
        #C_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Ch,C)),Ch)
	#r_cache['C'] = C_cache

	code = create_pseudo_random_code(len=10000,seed=2)
	r_cache2 = periodic_convolution_matrix(envelope=code,rmin=rmin,rmax=rmax)
	D = r_cache2['A']
	#Dh = numpy.transpose(numpy.conjugate(D))
        #D_cache = numpy.dot(numpy.linalg.inv(numpy.dot(Dh,D)),Dh)
	#r_cache['D'] = D_cache

	H=numpy.concatenate((A ,C,D))
	print H.shape,"1.3"
        Hh = numpy.transpose(numpy.conjugate(H))
        #print Hh.shape,"1.5"
        #temp=numpy.linalg.inv(numpy.dot(Hh,H))
	#print temp.shape,"1.7"
	#print Hh.shape , "1.9"
	#H_cache= numpy.dot(temp,Hh)
        H_cache= numpy.dot(numpy.linalg.inv(numpy.dot(Hh,H)),Hh)
#       print H_cache.shape,"2"
	r_cache['H']=H_cache
        B_cached = True
        return(r_cache)
    else:
        return(r_cache)

def grid_search1d(fun,xmin,xmax,nstep=100):
    vals = numpy.linspace(xmin,xmax,num=nstep)
    min_val=fun(vals[0])
    best_idx = 0
    for i in range(nstep):
        try_val = fun(vals[i])
        if try_val < min_val:
            min_val = try_val
            best_idx = i
    return(vals[best_idx])

def fit_velocity(z,t,var,frad=440.2e6):
    zz = numpy.exp(1.0j*numpy.angle(z))
    def ssfun(x):
        freq = 2.0*frad*x/scipy.constants.c
        model = numpy.exp(1.0j*2.0*scipy.constants.pi*freq*t)
        ss = numpy.sum((1.0/var)*numpy.abs(model-zz)**2.0)
        #        plt.plot( numpy.real(model))
        #plt.plot( numpy.real(zz), 'red')
        #plt.show()
        return(ss)
    v0 = grid_search1d(ssfun,-800.0,800.0,nstep=50)
    #    print v0
    #    v = scipy.optimize.fmin(ssfun,numpy.array([v0]),full_output=False,disp=False,retall=False)
    return(v0)

def fit_velocity_and_power(z,t,var,frad=440.2e6):
    zz = numpy.exp(1.0j*numpy.angle(z))
    def ssfun(x):
        freq = 2.0*frad*x/scipy.constants.c
        model = numpy.exp(1.0j*2.0*scipy.constants.pi*freq*t)
        ss = numpy.sum((1.0/var)*numpy.abs(model-zz)**2.0)
        return(ss)
    v0 = grid_search1d(ssfun,-800.0,800.0,nstep=50)
    v0 = scipy.optimize.fmin(ssfun,numpy.array([v0]),full_output=False,disp=False,retall=False)
    freq = 2.0*frad*v0/scipy.constants.c
    dc = numpy.real(numpy.exp(-1.0j*2.0*scipy.constants.pi*freq*t)*z)
    p0 = (1.0/numpy.sum(1.0/var))*numpy.sum((1.0/var)*dc)

    #plt.plot( dc)
    #    plt.show()

    #    print v0
    #
    return([v0,p0])

def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    with open(filename, 'rb') as input:
        return(pickle.load(input))

def unix2date(x):
    return datetime.datetime.utcfromtimestamp(x)

def unix2datestr(x):
    return(unix2date(x).strftime('%Y-%m-%d %H:%M:%S'))

def compr(x,fr=0.001):
    sh = x.shape
    x = x.reshape(-1)
    xs = numpy.sort(x)
    mini = xs[int(fr*len(x))]
    maxi = xs[int((1.0-fr)*len(x))]
    mx = numpy.ones_like(x)*maxi
    mn = numpy.ones_like(x)*mini
    print int(fr*len(x))," ",int((1.0-fr)*len(x))
    x = numpy.where(x < maxi, x, mx)
    x = numpy.where(x > mini, x, mn)
    x = x.reshape(sh)
    return(x)

def comprz(x):
    """ Compress signal in such a way that elements less than zero are set to zero. """
    zv = x*0.0
    return(numpy.where(x>0,x,zv))

def comprz_dB(xx,fr=0.05):
    """ Compress signal in such a way that is logarithmic but also avoids negative values """
    x = numpy.copy(xx)
    sh = xx.shape
    x = x.reshape(-1)
    x = comprz(x)
    x = numpy.setdiff1d(x,numpy.array([0.0]))
    xs = numpy.sort(x)
    mini = xs[int(fr*len(x))]
    mn = numpy.ones_like(xx)*mini
    xx = numpy.where(xx > mini, xx, mn)
    xx = xx.reshape(sh)
    return(10.0*numpy.log10(xx))

def decimate(x,dec=2):
    Nout = int(math.floor(len(x)/dec))
    idx = numpy.arange(Nout)*dec
    #    print idx
    res = x[idx]*0.0
    #print res
    for i in numpy.arange(dec):
        res = res + x[idx+i]
    return(res/float(dec))

def plot_cts(x,plot_abs=False,plot_show=True):
    time_vec = numpy.linspace(0,len(x)-1,num=len(x))
    plt.clf()
    plt.plot(time_vec,numpy.real(x),"blue")
    plt.plot(time_vec,numpy.imag(x),"red")
    if plot_abs:
        plt.plot(time_vec,numpy.abs(x),"black")
    if plot_show:
        plt.show()

def hanning(L=1000):
    n = numpy.linspace(0.0,L-1,num=L)
    return(0.5*(1.0-numpy.cos(2.0*scipy.constants.pi*n/L)))

def spectrogram(x,window=1024,wf=hanning):
    wfv = wf(L=window)
    Nwindow = int(math.floor(len(x)/window))
    res = numpy.zeros([Nwindow,window])
    for i in range(Nwindow):
        res[i,] = numpy.abs(numpy.fft.fftshift(numpy.fft.fft(wfv*x[i*window + numpy.arange(window)])))**2
    return(res)

