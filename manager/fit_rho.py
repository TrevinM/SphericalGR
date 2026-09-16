from scipy.optimize import curve_fit
import numpy as np
import sys

import matplotlib.pyplot as plt

lower = float(sys.argv[1])
lower_decimals = len(str(lower).split('.')[1])

upper = float(sys.argv[2])
upper_decimals = len(str(lower).split('.')[1])

def power_law(x, delta_c, gamma, a):
    x_c = lower + delta_c
    return a * np.greater(x_c - x, 0) * np.abs(x_c - x)**(-2*gamma)

def main():
    eta = []
    data = []


    precision = max(lower_decimals, upper_decimals)

    with open("data_rho", 'r') as f:
        for line in f.readlines():
            vals = line.split()
            eta.append(vals[0])
            data.append(vals[1])
        
    p0=[(upper - lower)/2.1, 0.3, 1]
    lows=[1.e-12, 0.02, 0]
    upps=[(upper - lower), 1, 10]
    scale=[10**(-precision), 1, 1]
    popt, pcov = curve_fit(power_law, eta, data, p0=p0, bounds=(lows, upps), method='trf', x_scale=scale, nan_policy='omit')
    print(popt[0] + lower)

    print(f"{popt[0]} {popt[1]} {popt[2]}")

if __name__ == "__main__":
    main()