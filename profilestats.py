#!/usr/bin/env python
import argparse
import numpy as np
import GPy
import matplotlib.pyplot as plt
import scipy.spatial as sp

def KD(X1, X2, theta):

    X1, X2 = np.matrix(X1), np.matrix(X2) # ensure both sets of inputs are matrices
    D2 = sp.distance.cdist(X1, X2, 'sqeuclidean') # calculate squared Euclidean distance
    D1 = np.zeros((X1.shape[0],X2.shape[0]))
    for i in range(X1.shape[0]):
        for j in range(X2.shape[0]):
            D1[i,j] = X1[i] - X2[j]
    K = - theta[0] * D1 * np.exp(- D2 / (2*(theta[1]**2))) / theta[1]**2
    return np.matrix(K)

print "Start"
# Read command line arguments
#------------------------------
parser = argparse.ArgumentParser(description='Pulsar nudot variability studies using GPs')
parser.add_argument('-f','--filename', help='File containing ascii profile', required=True)
parser.add_argument('-o','--outdir', help='output folder', required=True)
parser.add_argument('-i','--interpulse',action='store_true',help='run interpulses')
parser.add_argument('-i2','--interpulse2',action='store_true',help='run interpulse 2')
#------------------------------
args = parser.parse_args()
filename = args.filename
outdir = args.outdir
interpulse = args.interpulse
interpulse2 = args.interpulse2
print filename
with open(filename, 'r') as f:
    psrname = f.readline().split()[3]
if len(psrname) < 10:
    psrname = 'J'+filename[0:9]
print psrname
data = np.loadtxt(filename, comments='F')
print data.shape
stokesI = data[:,3]
peaknorm = max(stokesI)
peakbin = np.argmax(stokesI)
bins = len(stokesI)
stokesI = np.roll(stokesI, -peakbin + bins/2)/peaknorm
xaxis=np.arange(bins)
# I have the data, so train a GP
k1 = GPy.kern.RBF(1)
k2 = GPy.kern.RBF(1)
kernel = k1 
#kernel = k1 + k2
if interpulse:
    xaxis=np.arange(int(bins/2))
    xtraining = xaxis.reshape(xaxis.shape[0],1)
    newStokesI = stokesI[int(0.25*bins):int(0.75*bins)] 
    stokesI = newStokesI
    ytraining = stokesI.reshape(stokesI.shape[0],1)
elif interpulse2:
    stokesI = np.roll(stokesI, bins/2)
    newStokesI = stokesI[int(0.25*bins):int(0.75*bins)] 
    stokesI = newStokesI
    peaknorm = max(stokesI)
    peakbin = np.argmax(stokesI)
    bins = len(stokesI)
    stokesI = np.roll(stokesI, -peakbin + bins/2)/peaknorm
    xaxis=np.arange(bins)
    xtraining = xaxis.reshape(xaxis.shape[0],1)
    ytraining = stokesI.reshape(stokesI.shape[0],1)
else:
    xtraining = xaxis.reshape(xaxis.shape[0],1)
    ytraining = stokesI.reshape(stokesI.shape[0],1)
model = GPy.models.GPRegression(xtraining, ytraining ,kernel)
model.optimize()
#model.optimize_restarts(num_restarts = 5)
ypredict, yvariance = model.predict(xtraining)
print model
#model.plot()
#plt.show()
# Now use the Gaussian Process to obtain the derivative
par = np.zeros(2)
K1 = kernel.K(xtraining, xtraining)
K1invOut = GPy.util.linalg.pdinv(np.matrix(K1))
K1inv = K1invOut[1]
xtraining_d = np.matrix([xaxis]).T
xpredict_d = xtraining_d
ytraining_d = np.matrix(np.array(stokesI.flatten())).T
# First lengthscale kernel
par[0] = model['rbf.variance'] # use the optimized amplitude
par[1] = model['rbf.lengthscale'] # use the optimized lengthscale
K_prime = KD(xpredict_d, xtraining_d, par) # training points
# Second lengthscale kernel
#par[0] = model['sum.rbf_1.variance'] # use the optimized amplitude
#par[1] = model['sum.rbf_1.lengthscale'] # use the optimized lengthscale
#K_prime += KD(xpredict_d, xtraining_d, par) # training points

#K_prime_p = par[0]/par[1]**2 # These are the diagonal elements of the
#variance for the error in the derivative

KiKx, _ = GPy.util.linalg.dpotrs(K1inv, np.asfortranarray(K_prime.T), lower = 1)
stokes_dot = np.array(np.dot(KiKx.T, ytraining_d))
peak = max(ypredict)
limit = 0.1 * peak[0]
limit1 = 0.01 * peak[0]
w10 = sum(i > 0.1 * peak for i in ypredict)
w1 = sum(i > 0.01 * peak for i in ypredict)
masked_derivative = np.ma.masked_where(ypredict<limit, stokes_dot)
dmax = masked_derivative.max(fill_value=-1.0)
dbmax = masked_derivative.argmax(fill_value=-1.0)
dmin = masked_derivative.min(fill_value=1.0)
dbmin = masked_derivative.argmin(fill_value=1.0)
if np.sqrt(model['Gaussian_noise.variance']) < 0.01:
    print psrname, " w1: ", float(w1)/float(len(ypredict))*360.0
    with open("w1.dat", "a") as myfile:
        myfile.write(psrname+" w1: "+str(float(w1)/float(len(ypredict))*360.0)+"\n")

if np.sqrt(model['Gaussian_noise.variance']) < 0.1:
    print psrname, " w10: ", float(w10)/float(len(ypredict))*360.0
    with open("w10.dat", "a") as myfile:
        myfile.write(psrname+" w10: "+str(float(w10)/float(len(ypredict))*360.0)+"\n")

print " Steepest + gradient: ", dmax, " at ", dbmax
print " Steepest - gradient: ", dmin, " at ", dbmin
print "peak: ", peak
fig, ax = plt.subplots()
ax.plot(xaxis, stokesI, 'g.')
ax.fill_between(xaxis, np.add(ypredict.flatten(),1.0*np.sqrt(yvariance.flatten())),  np.add(ypredict.flatten(),-1.0*np.sqrt(yvariance.flatten())), facecolor='blue', alpha=0.5)
#ax.plot(xaxis, ypredict, 'b-')
ax.plot(xaxis, stokes_dot, 'r-')
ax.axvline(dbmax)
ax.axvline(dbmin)
ax.fill_between(xaxis, ypredict.flatten(), limit , where=ypredict.flatten()>limit,facecolor='orange', alpha=0.5)
ax.fill_between(xaxis, ypredict.flatten(), limit1 , where=ypredict.flatten()>limit1,facecolor='yellow', alpha=0.5)
plt.savefig('./'+outdir+'/'+filename+'.png')

