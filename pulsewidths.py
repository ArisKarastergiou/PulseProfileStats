#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants

def rho(P, height):
    """Function to determine the opening angle rho given the rotational
       period and emission height.

       Args:
       -----
       P      : rotational period (seconds)
       height   : emission height (in km).

       Returns:
       --------
       rho    : the opening angle (rad)

    """
    rho = np.sqrt((9.0 * np.pi * height * 1000.0) / (2.0 * constants.c * P))

    return rho

def width(rho, alpha, beta, zeta):
    """Function to return the pulse width as a function of inclination
    angle and impact parameter
    
    Args:
    rho: beam opening angle (rad)
    alpha: inclination angle (rad)
    beta: impact parameter (rad)

    Returns:
    width: pulse width (degrees)
    """
    width = np.degrees(2.0 * np.arccos((np.cos(rho) - np.cos(alpha) * np.cos(zeta))/(np.sin(alpha)*np.sin(zeta))))
    return width

# Main
# Read in arguments
# Read command line arguments
#------------------------------
parser = argparse.ArgumentParser(description='Pulsar width simulation')
parser.add_argument('-P','--period', help='constant period', type = float)
parser.add_argument('-z','--height', help='constant emission height', type = float)
parser.add_argument('-n','--number', help='number to simulate', type = int, required=True)
parser.add_argument('-f','--filename', help='list of periods')
parser.add_argument('-i','--interpulse',action='store_true',help='run interpulses')
#------------------------------
args = parser.parse_args()
period = args.period
height = args.height
population = args.number
file = args.filename
interpulse = args.interpulse
# Compute rho
if period==None:
    data=np.loadtxt(file)
    period = data[:,1]
    filewidths = data[:,0]
    print period.shape
    population = len(period)

# Check if height is given
if height==None:
    height = 0.02 * constants.c * period / 2.0 * np.pi / 1000.0

my_rho = rho(period, height)

#
if interpulse:
    print 'Orthogonal alphas'
    alpha = np.zeros(population)
    beta = np.zeros(population)
    found = 0
    while (found < population):
        b = np.random.uniform(-my_rho[found], my_rho[found], 1)
        a = np.arccos(np.random.uniform(0.0, np.cos(my_rho[found]), 1))
        if np.abs(np.pi-2*a-b) < my_rho[found]:
            alpha[found] = a
            beta[found] = b
            found += 1
        
    n, bins, patches = plt.hist(alpha*180.0/np.pi, 50, facecolor='green', alpha=0.75)
    plt.show()
else:
    alpha = np.arccos(np.random.uniform(0.0, np.cos(my_rho), population))
    
beta = np.random.uniform(-my_rho, my_rho, population)
zeta = np.add(alpha, beta)
widths = width(my_rho, alpha, beta, zeta)
logwidths = np.log10(widths)
expectedlogwidth = np.log10(6.0 / np.sqrt(period))
logfilewidths = np.log10(filewidths)
diff = np.subtract(logwidths, logfilewidths)
residuals = logfilewidths - expectedlogwidth
residuals2 = logwidths - expectedlogwidth
print "difference mean and std: ", np.mean(diff), np.std(diff)
print np.std(logwidths), np.mean(logwidths)
n, bins, patches = plt.hist(residuals, 30, facecolor='green', alpha=0.5)
n, bins, patches = plt.hist(residuals2, 30, facecolor='yellow', alpha=0.5)
plt.xlabel('log10 widths - model')
plt.ylabel('number')
plt.title('width simulation')
plt.grid(True)
plt.show()
plt.scatter(np.log10(period), logwidths, c='y', s=100, edgecolors='none',label='Simulated')
plt.scatter(np.log10(period), logfilewidths, c='b', s=100, edgecolors='none', label='Observed')
plt.legend()
#plt.plot(np.log10(period), diff, 'g.')
plt.xlabel('log10 period (s)')
plt.ylabel('log10 width (degrees)')
plt.title('Period width: Simulated versus observed')
plt.grid(True)
plt.show()
