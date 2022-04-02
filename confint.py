#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Program to calculate estimates and confidence intervals of portions 
Arguments: number_of_hits sample_size conf_int_percentage
           the last argument is optional, 95 is assumed if missing

Example: python3 confint.py 32 48 99
          32 positive cases in 48 data points, 99% conf. interval requested

@author: laszlo
"""

import sys
import numpy as np
from scipy import stats as st


def z_from_percent(pct, onesided=False):
    """
    Z value for one/two-sided normal distribution from % of conf. interval
    """

    alpha = 1.0 - pct / 100.0

    if onesided:
        z = st.norm.ppf(1.0 - alpha)
    else:
        z = st.norm.ppf(1.0 - alpha/2.0)

    return z


def normal_approx(x, n, pct):
    """
    Normal distribution approximation of confidence interval
    (sometimes called Wald)
    Not useful for n < 20 (or x < 10 or n - x < 10)
    """

    p = x/n
    sigma = np.sqrt(p*(1.0-p)/n)
    hw = z_from_percent(pct) * sigma
    low = max(0.0, p - hw)
    high = min(1.0, p + hw)
    return p, low, high


def agresti_coull(x, n, pct):
    """
    Agresti-Coull interval (gives the same p estimate as Wilson score)
    Recommended for samples < 150; gives *on average* correct coverage
    (meaning on average 95% of values are really in the 95% conf. interval)
    """

    z = z_from_percent(pct)
    return(normal_approx(x + z*z/2.0, n + z*z, pct))


def clopper_pearson(x, n, pct):
    """
    Clopper-Pearson interval (often called 'exact')
    Sometimes excessively conservative, but guarantees *at least* the
    specified coverage
    Gives the best c.i. values for x = 0 or x = n
    """

    halpha = 0.5 - pct/200.0

    # use analytical solution for extreme cases for numerical stability
    if x == 0:
        q = 1.0 - halpha**(1.0/n)
        return q/2.0, 0.0, q

    if x == n:
        q = halpha**(1.0/n)
        return (q + 1.0) / 2.0, q, 1.0

    low = st.beta.ppf(halpha, x, n-x+1)
    high = st.beta.ppf(1.0-halpha, x+1, n-x)
    p = (high + low) / 2.0
    return p, low, high


def wilson_score(x, n, pct):
    """
    Wilson score interval
    Based on the routines in EBCIC (https://github.com/KazKobara/ebcic)
    Similar to Agresti-Coull in usage, except for weak coverage at certain
    probabilities (e.g. around 0.98 and 0.02). Better than A-C at close
    to 0 and 100%
    """
    z = z_from_percent(pct)
    p = x / n
    z2 = z * z
    mu = 2.0 * x + z2
    hw = z * np.sqrt(z2 + 4.0 * x * (1.0 - p))
    denom = 2.0 * (n + z2)
    low = max(0, (mu - hw) / denom)
    high = min(1, (mu + hw) / denom)
    return mu/denom, low, high


def wilson_score_cc(x, n, pct):
    """
    Wilson score interval with continuity correction.
    Based on the routines in EBCIC (https://github.com/KazKobara/ebcic)
    """

    z = z_from_percent(pct)
    p = x / n
    z2 = z * z
    mu = 2.0 * x + z2
    hw = 1.0 + z * np.sqrt(z2 - 1 / n + 4 * x * (1 - p) + (4 * p - 2))
    denom = 2.0 * (n + z2)
    low = max(0, (mu - hw) / denom)
    high = min(1, (mu + hw) / denom)
    return mu/denom, low, high


if len(sys.argv) < 3:
    print("Calculates confidence intervals of portions\n")
    print("Usage: {} nr_positive nr_total".format(sys.argv[0]), end="")
    print(" [confidence_percent]\n  95% confidence interval is", end="")
    print(" calculated if the last argument is omitted")
    quit()

x = int(sys.argv[1])
n = int(sys.argv[2])

if len(sys.argv) > 3:
    pct = float(sys.argv[3])
else:
    pct = 95.0

print("\nEstimates of probability\n")
print("Max. likelihood (x/n) = {}".format(x/n))

# Laplace useful if there are no positive or negative observations
print("Laplace   (x+1)/(n+2) = {}".format((x+1.0)/(n+2.0)))

# Wilson pulls towards 0.5; also helpful with extreme proportions
# Unlike Laplace, this depends on the chosen probability, with the
# higher probabity scores pulling more strongly towards 0.5
z2 = z_from_percent(pct)**2
print("Wilson (x+z2/2)/(n+z2)= {}".format((x+z2/2.0)/(n+z2)))

print("\nEstimates of {}% confidence interval (mid, min, max)\n".format(pct))
print("Normal approx:    ", end="")
print(normal_approx(x, n, pct))
print("Clopper-Pearson:  ", end="")
print(clopper_pearson(x, n, pct))
print("Agresti-Coull:    ", end="")
print(agresti_coull(x, n, pct))
print("Wilson score:     ", end="")
print(wilson_score(x, n, pct))
print("Wilson cont corr: ", end="")
print(wilson_score_cc(x, n, pct))
