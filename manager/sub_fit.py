from scipy.optimize import curve_fit
import numpy as np
import sys

import matplotlib.pyplot as plt

def power_law(x, x_c, gamma):
    return np.greater(x_c - x, 0) * np.abs(x_c - x)**(-2*gamma)

def main():
    eta = []
    data = []

    lower = float(sys.argv[1])
    upper = float(sys.argv[2])

    with open("data_sub", 'r') as f:
        for line in f.readlines():
            vals = line.split()
            eta.append(vals[0])
            data.append(vals[1])

    popt, pcov = curve_fit(power_law, eta, data, p0=[lower, 0.377], bounds=([lower, 0], [upper, 1]))

    print(f"{popt[0]} {popt[1]}")

if __name__ == "__main__":
    main()