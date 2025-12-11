# alpha (y), L/Reff (x)
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np
plt.style.use('../mplstyles/standard.mplstyle')

def make_alpha_correlation_plots(fit_filename, luminosity_filename, legend_text, output_filename):
    df_lum = pd.read_csv(luminosity_filename, sep="\t", skiprows=1)
    df_fit = pd.read_csv(fit_filename)
    
    merged = df_lum.merge(df_fit, on="Galaxy")

    # take logs for the fit
    logx = np.log(merged["L/R sof"].values)
    logy = np.log(merged["alpha"].values)

    # linear regression: log y = m log x + b
    slope, intercept, r_value, p_value, stderr = linregress(logx, logy)

    k = slope
    A = np.exp(intercept)

    print("k (power-law exponent) =", k)
    print("A (prefactor) =", A)
    print("R^2 =", r_value**2)

    xx = np.linspace(min(merged["L/R sof"].values), max(merged["L/R sof"].values), 200)
    yy = A * xx**k

    midblue = '#1f77b4'
    redorange = '#ff7f0e'

    # plot the fit
    plt.loglog(xx, yy, '-', color=redorange,
               label=rf"fit: $\alpha = {A:.3g}\,(L/R)^{{{k:.3g}}}$")

    # plot the data
    plt.loglog(merged["L/R sof"].values, merged["alpha"].values, 'o', color='darkblue',
               label=legend_text)

    plt.xlabel("L/R")
    plt.ylabel(r"$\alpha$")
    plt.legend()
    plt.grid(True, which="both")
    #plt.xlim(10E-2,10E2)
    plt.ylim(0.5,10E5)

    plt.savefig(output_filename, dpi=300, bbox_inches="tight")