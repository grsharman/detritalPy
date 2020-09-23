"""
Created on Sat Feb 18 07:43:18 2017

This package implements adaptive kernel density estimation algorithms for 1-dimensional 
signals developed by Hideaki Shimazaki. This enables the generation of smoothed histograms
that preserve important density features at multiple scales, as opposed to naive
single-bandwidth kernel density methods that can either over or under smooth density
estimates. These methods are described in Shimazaki's paper:

 H. Shimazaki and S. Shinomoto, "Kernel Bandwidth Optimization in Spike Rate Estimation," 
 in Journal of Computational Neuroscience 29(1-2): 171–182, 2010 
 http://dx.doi.org/10.1007/s10827-009-0180-4.
 
License:
All software in this package is licensed under the Apache License 2.0.
See LICENSE.txt for more details.
 
Authors:
Hideaki Shimazaki (shimazaki@jhu.edu) shimazaki on Github
Lee A.D. Cooper (cooperle@gmail.com) cooperlab on GitHub
Subhasis Ray (ray.subhasis@gmail.com)

Note: This code was modified slightly by Glenn R. Sharman (gsharman@uark.edu) grsharman on GitHub
for the purpose of use with detritalPy
 
Three methods are implemented in this package:
1. sshist - can be used to determine the optimal number of histogram bins for independent 
identically distributed samples from an underlying one-dimensional distribution. The
principal here is to minimize the L2 norm of the difference between the histogram and the
underlying distribution.

2. sskernel - implements kernel density estimation with a single globally-optimized 
bandwidth.

3. ssvkernel - implements kernel density estimation with a locally variable bandwidth.
 
Dependencies: These functions in this package depend on NumPy for various operations 
including fast-fourier transforms and histogram generation.
"""

import numpy as np

def sshist(x, N=range(2, 501), SN=30):
    """
    Returns the optimal number of bins in a histogram used for density
    estimation.

    Optimization principle is to minimize expected L2 loss function between
    the histogram and an unknown underlying density function.
    An assumption made is merely that samples are drawn from the density
    independently each other.

    The optimal binwidth D* is obtained as a minimizer of the formula,
    (2K-V) / D^2,
    where K and V are mean and variance of sample counts across bins with width
    D. Optimal number of bins is given as (max(x) - min(x)) / D.

    Parameters
    ----------
    x : array_like
        One-dimensional data to fit histogram to.
    N : array_like, optional
        Array containing number of histogram bins to evaluate for fit.
        Default value = 500.
    SN : double, optional
        Scalar natural number defining number of bins for shift-averaging.

    Returns
    -------
    optN : int
        Optimal number of bins to represent the data in X
    N : double
        Maximum number of bins to be evaluated. Default value = 500.
    C : array_like
        Cost function C[i] of evaluating histogram fit with N[i] bins

    See Also
    --------
    sskernel, ssvkernel

    References
    ----------
    .. [1] H. Shimazaki and S. Shinomoto, "A method for selecting the bin size
           of a time histogram," in  Neural Computation 19(6), 1503-1527, 2007
           http://dx.doi.org/10.1162/neco.2007.19.6.1503
    """

    # determine range of input 'x'
    x_min = np.min(x)
    x_max = np.max(x)

    # get smallest difference 'dx' between all pairwise samples
    buf = np.abs(np.diff(np.sort(x)))
    dx = min(buf[buf > 0])

    # setup bins to evaluate
    N_MIN = 2
    N_MAX = min(np.floor((x_max - x_min) / (2*dx)), max(N))
    N = range(N_MIN, N_MAX+1)
    D = (x_max - x_min) / N

    # compute cost function over each possible number of bins
    Cs = np.zeros((len(N), SN))
    for i, n in enumerate(N):  # loop over number of bins
        shift = np.linspace(0, D[i], np.int(SN))
        for p, sh in enumerate(shift):  # loop over shift window positions

            # define bin edges
            edges = np.linspace(x_min + sh - D[i]/2,
                                x_max + sh - D[i]/2, np.int(N[i]+1))

            # count number of events in these bins
            ki = np.histogram(x, edges)

            # get mean and variance of events
            k = ki[0].mean()
            v = np.sum((ki[0] - k)**2) / N[i]

            Cs[i, p] = (2*k - v) / D[i]**2

    # average over shift window
    C = Cs.mean(axis=1)

    # get bin count that minimizes cost C
    idx = np.argmin(C)
    optN = N[idx]
    optD = D[idx]
    edges = np.linspace(x_min, x_max, np.int(optN))

    return optN, optD, edges, C, N

def sskernel(x, tin=None, W=None, nbs=1000):
    """
    Generates a kernel density estimate with globally-optimized bandwidth.

    The optimal bandwidth is obtained as a minimizer of the formula, sum_{i,j}
    \int k(x - x_i) k(x - x_j) dx  -  2 sum_{i~=j} k(x_i - x_j), where k(x) is
    the kernel function.


    Parameters
    ----------
    x : array_like
        The one-dimensional samples drawn from the underlying density
    tin : array_like, optional
        The values where the density estimate is to be evaluated in generating
        the output 'y'.
    W : array_like, optional
        The kernel bandwidths to use in optimization. Should not be chosen
        smaller than the sampling resolution of 'x'.
    nbs : int, optional
        The number of bootstrap samples to use in estimating the [0.05, 0.95]
        confidence interval of the output 'y'

    Returns
    -------
    y : array_like
        The estimated density, evaluated at points t / tin.
    t : array_like
        The points where the density estimate 'y' is evaluated.
    optw : double
        The optimal global kernel bandwidth.
    W : array_like
        The kernel bandwidths evaluated during optimization.
    C : array_like
        The cost functions associated with the bandwidths 'W'.
    confb95 : array_like
        The 5% and 95% confidence interval of the kernel density estimate 'y'.
        Has dimensions 2 x len(y). confb95[0,:] corresponds to the 5% interval,
        and confb95[1,:] corresponds to the 95% interval.
    yb : array_like
        The bootstrap samples used in estimating confb95. Each row corresponds
        to one bootstrap sample.

    See Also
    --------
    sshist, ssvkernel

    References
    ----------
    .. [1] H. Shimazaki and S. Shinomoto, "Kernel Bandwidth Optimization in 
           Spike Rate Estimation," in Journal of Computational Neuroscience 
           29(1-2): 171–182, 2010 http://dx.doi.org/10.1007/s10827-009-0180-4
    """

    # set argument 't' if not provided
    if tin is None:
        T = np.max(x) - np.min(x)
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        tin = np.linspace(np.min(x), np.max(x), np.int(min(np.ceil(T / dt_samp), 1e3)))
        t = tin
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
    else:
        T = np.max(x) - np.min(x)
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        if dt_samp > min(np.diff(tin)):
            t = np.linspace(min(tin), max(tin), np.int(min(np.ceil(T / dt_samp), 1e3)))
        else:
            t = tin

    # calculate delta t
    dt = min(np.diff(t))

    # create the finest histogram
    thist = np.concatenate((t, (t[-1]+dt)[np.newaxis]))
    y_hist = np.histogram(x_ab, thist-dt/2)[0]
    N = sum(y_hist).astype(np.float)
    y_hist = y_hist / N / dt

    # global search if input 'W' is defined
    if W is not None:
        C = np.zeros((1, len(W)))
        C_min = np.Inf
        for k, w in enumerate(W):
            C[k], yh = CostFunctionSS(y_hist, N, w, dt)
            if(C[k] < C_min):
                C_min = C[k]
                optw = w
                y = yh
    else:  # optimized search using golden section
        k = 0
        C = np.zeros((20, 1))
        W = np.zeros((20, 1))
        Wmin = 2*dt
        Wmax = (np.max(x) - np.min(x))
        tol = 10e-5
        phi = (5**0.5 + 1) / 2
        a = ilogexpSS(Wmin)
        b = ilogexpSS(Wmax)
        c1 = (phi - 1) * a + (2 - phi) * b
        c2 = (2 - phi) * a + (phi - 1) * b
        f1, dummy = CostFunctionSS(y_hist, N, logexpSS(c1), dt)
        f2, dummy = CostFunctionSS(y_hist, N, logexpSS(c2), dt)
        while (np.abs(b-a) > tol * (np.abs(c1) + np.abs(c2))) & (k < 20):
            if f1 < f2:
                b = c2
                c2 = c1
                c1 = (phi - 1) * a + (2 - phi) * b
                f2 = f1
                f1, yh1 = CostFunctionSS(y_hist, N, logexpSS(c1), dt)
                W[k] = logexpSS(c1)
                C[k] = f1
                optw = logexpSS(c1)
                y = yh1 / np.sum(yh1 * dt)
            else:
                a = c1
                c1 = c2
                c2 = (2 - phi) * a + (phi - 1) * b
                f1 = f2
                f2, yh2 = CostFunctionSS(y_hist, N, logexpSS(c2), dt)
                W[k] = logexpSS(c2)
                C[k] = f2
                optw = logexpSS(c2)
                y = yh2 / np.sum(yh2 * dt)

            # increment iteration counter
            k = k + 1

        # discard unused entries in gs, C
        C = C[0:k]
        W = W[0:k]

    # estimate confidence intervals by bootstrapping
    nbs = np.asarray(nbs)
    yb = np.zeros((int(nbs),len(tin)))
    for i in range(nbs):
        idx = np.random.randint(0, len(x_ab)-1, len(x_ab))
        xb = x_ab[idx]
        thist = np.concatenate((t, (t[-1]+dt)[np.newaxis]))
        y_histb = np.histogram(xb, thist - dt / 2)[0] / dt / N
        yb_buf = fftkernel(y_histb, optw / dt)
        yb_buf = yb_buf / np.sum(yb_buf * dt)
        yb[i, ] = np.interp(tin, t, yb_buf)
    ybsort = np.sort(yb, axis=0)
    y95b = ybsort[np.int(np.floor(0.05 * nbs)), :]
    y95u = ybsort[np.int(np.floor(0.95 * nbs)), :]
    confb95 = np.concatenate((y95b[np.newaxis], y95u[np.newaxis]), axis=0)

    # return outputs
    y = np.interp(tin, t, y)
    t = tin

    return y, t, optw, W, C, confb95, yb

def CostFunctionSS(y_hist, N, w, dt):

    # build normal smoothing kernel
    yh = fftkernel(y_hist, w / dt)

    # formula for density
    C = np.sum(yh**2) * dt - 2 * np.sum(yh * y_hist) * dt + 2 \
        / (2 * np.pi)**0.5 / w / N
    C = C * N**2

    return C, yh

def logexpSS(x):
    if x < 1e2:
        y = np.log(1 + np.exp(x))
    else:
        y = x
    return y


def ilogexpSS(x):
    # ilogexp = log(exp(x)-1);
    if x < 1e2:
        y = np.log(np.exp(x) - 1)
    else:
        y = x
    return y

def ssvkernel(x, tin=None, M=80, nbs=100, WinFunc='Boxcar'):
    """
    Generates a locally adaptive kernel-density estimate for one-dimensional
    data.

    The user provides a one-dimensional vector of samples drawn from some
    underlying unknown distribution, and optionally the values where they want
    to estimate the probability density of that distribution. The algorithm
    solves an optimization problem to identify variable bandwidths across the
    domain where the data is provided.

    The optimization is based on a principle of minimizing expected L2 loss
    function between the kernel estimate and an unknown underlying density
    function. An assumption is merely that samples are drawn from the density
    independently of each other.

    The locally adaptive bandwidth is obtained by iteratively computing optimal
    fixed-size bandwidths wihtihn local intervals. The optimal bandwidths are
    selected such that they are selected in the intervals that are gamma times
    larger than the optimal bandwidths themselves. The paramter gamma is
    optimized by minimizing the L2 risk estimate.

    Parameters
    ----------
    x : array_like
        The one-dimensional samples drawn from the underlying density
    tin : array_like, optional
        The values where the density estimate is to be evaluated in generating
        the output 'y'. Default value = None.
    M : int, optional
        The number of window sizes to evaluate. Default value = 80.
    nbs : int, optional
        The number of bootstrap samples to use in estimating the [0.05, 0.95]
        confidence interval of the output 'y'.
    WinFunc : string, optional
        The type of window function to use in estimating local bandwidth.
        Choose from one of 'Boxcar', 'Laplace', 'Cauchy' and 'Gauss'. Default
        value = 'Gauss'.

    Returns
    -------
    y : array_like
        The estimated density, evaluated at points t / tin.
    t : array_like
        The points where the density estimate 'y' is evaluated.
    optw : array_like
        The optimal local kernel bandwidths at 't'.
    gs : array_like
        The stiffness constants of the variables bandwidths evaluated.
    C : array_like
        Cost functions associated with stiffness constraints.
    confb95 : array_like
        The 5% and 95% confidence interval of the kernel density estimate 'y'.
        Has dimensions 2 x len(y). confb95[0,:] corresponds to the 5% interval,
        and confb95[1,:] corresponds to the 95% interval.
    yb : array_like
        The bootstrap samples used in estimating confb95. Each row corresponds
        to one bootstrap sample.

    See Also
    --------
    sshist, sskernel

    References
    ----------
    .. [1] H. Shimazaki and S. Shinomoto, "Kernel Bandwidth Optimization in 
           Spike Rate Estimation," in Journal of Computational Neuroscience 
           29(1-2): 171–182, 2010 http://dx.doi.org/10.1007/s10827-009-0180-4
    """

    # set argument 't' if not provided
    if tin is None:
        T = np.max(x) - np.min(x)
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        tin = np.linspace(np.min(x), np.max(x), np.int(min(np.ceil(T / dt_samp), 1e3)))
        t = tin
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
    else:
        T = np.max(x) - np.min(x)
        x_ab = x[(x >= min(tin)) & (x <= max(tin))]
        dx = np.sort(np.diff(np.sort(x)))
        dt_samp = dx[np.nonzero(dx)][0]
        if dt_samp > min(np.diff(tin)):
            t = np.linspace(min(tin), max(tin), np.int(min(np.ceil(T / dt_samp), 1e3)))
        else:
            t = tin

    # calculate delta t
    dt = min(np.diff(t))

    # create the finest histogram
    thist = np.concatenate((t, (t[-1]+dt)[np.newaxis]))
    y_hist = np.histogram(x_ab, thist-dt/2)[0] / dt
    L = y_hist.size
    N = sum(y_hist * dt).astype(np.float)

    # initialize window sizes
    W = logexpSSV(np.linspace(ilogexpSSV(5 * dt), ilogexpSSV(T), np.int(M)))

    # compute local cost functions
    c = np.zeros((M, L))
    for j in range(M):
        w = W[j]
        yh = fftkernel(y_hist, w / dt)
        c[j, :] = yh**2 - 2 * yh * y_hist + 2 / (2 * np.pi)**0.5 / w * y_hist

    # initialize optimal ws
    optws = np.zeros((M, L))
    for i in range(M):
        Win = W[i]
        C_local = np.zeros((M, L))
        for j in range(M):
            C_local[j, :] = fftkernelWin(c[j, :], Win / dt, WinFunc)
        n = np.argmin(C_local, axis=0)
        optws[i, :] = W[n]

    # golden section search for stiffness parameter of variable bandwidths
    k = 0
    gs = np.zeros((30, 1))
    C = np.zeros((30, 1))
    tol = 1e-5
    a = 1e-12
    b = 1
    phi = (5**0.5 + 1) / 2
    c1 = (phi - 1) * a + (2 - phi) * b
    c2 = (2 - phi) * a + (phi - 1) * b
    f1 = CostFunctionSSV(y_hist, N, t, dt, optws, W, WinFunc, c1)[0]
    f2 = CostFunctionSSV(y_hist, N, t, dt, optws, W, WinFunc, c2)[0]
    while (np.abs(b-a) > tol * (abs(c1) + abs(c2))) & (k < 30):
        if f1 < f2:
            b = c2
            c2 = c1
            c1 = (phi - 1) * a + (2 - phi) * b
            f2 = f1
            f1, yv1, optwp1 = CostFunctionSSV(y_hist, N, t, dt, optws, W,
                                           WinFunc, c1)
            yopt = yv1 / np.sum(yv1 * dt)
            optw = optwp1
        else:
            a = c1
            c1 = c2
            c2 = (2 - phi) * a + (phi - 1) * b
            f1 = f2
            f2, yv2, optwp2 = CostFunctionSSV(y_hist, N, t, dt, optws, W,
                                           WinFunc, c2)
            yopt = yv2 / np.sum(yv2 * dt)
            optw = optwp2

        # capture estimates and increment iteration counter
        gs[k] = c1
        C[k] = f1
        k = k + 1

    # discard unused entries in gs, C
    gs = gs[0:k]
    C = C[0:k]

    # estimate confidence intervals by bootstrapping
    nbs = np.asarray(nbs)
    yb = np.zeros((nbs, tin.size))
    for i in range(nbs):
        Nb = np.random.poisson(lam=N)
        idx = np.random.randint(0, N, Nb)
        xb = x_ab[idx]
        thist = np.concatenate((t, (t[-1]+dt)[np.newaxis]))
        y_histb = np.histogram(xb, thist - dt / 2)[0]
        idx = y_histb.nonzero()
        y_histb_nz = y_histb[idx]
        t_nz = t[idx]
        yb_buf = np.zeros((L, ))
        for k in range(L):
            yb_buf[k] = np.sum(y_histb_nz * Gauss(t[k] - t_nz, optw[k])) / Nb
        yb_buf = yb_buf / np.sum(yb_buf * dt)
        yb[i, :] = np.interp(tin, t, yb_buf)
    ybsort = np.sort(yb, axis=0)
    y95b = ybsort[np.int(np.floor(0.05 * nbs)), :]
    y95u = ybsort[np.int(np.floor(0.95 * nbs)), :]
    confb95 = np.concatenate((y95b[np.newaxis], y95u[np.newaxis]), axis=0)

    # return outputs
    y = np.interp(tin, t, yopt)
    optw = np.interp(tin, t, optw)
    t = tin

    return y, t, optw, gs, C, confb95, yb


def CostFunctionSSV(y_hist, N, t, dt, optws, WIN, WinFunc, g):

    L = y_hist.size
    optwv = np.zeros((L, ))
    for k in range(L):
        gs = optws[:, k] / WIN
        if g > np.max(gs):
            optwv[k] = np.min(WIN)
        else:
            if g < min(gs):
                optwv[k] = np.max(WIN)
            else:
                idx = np.max(np.nonzero(gs >= g))
                optwv[k] = g * WIN[idx]

    # Nadaraya-Watson kernel regression
    optwp = np.zeros((L, ))
    for k in range(L):
        if WinFunc == 'Boxcar':
            Z = Boxcar(t[k]-t, optwv / g)
        elif WinFunc == 'Laplace':
            Z = Laplace(t[k]-t, optwv / g)
        elif WinFunc == 'Cauchy':
            Z = Cauchy(t[k]-t, optwv / g)
        else:  # WinFunc == 'Gauss'
            Z = Gauss(t[k]-t, optwv / g)
        optwp[k] = np.sum(optwv * Z) / np.sum(Z)

    # speed-optimized baloon estimator
    idx = y_hist.nonzero()
    y_hist_nz = y_hist[idx]
    t_nz = t[idx]
    yv = np.zeros((L, ))
    for k in range(L):
        yv[k] = np.sum(y_hist_nz * dt * Gauss(t[k]-t_nz, optwp[k]))
    yv = yv * N / np.sum(yv * dt)

    # cost function of estimated kernel
    cg = yv**2 - 2 * yv * y_hist + 2 / (2 * np.pi)**0.5 / optwp * y_hist
    Cg = np.sum(cg * dt)

    return Cg, yv, optwp


def fftkernel(x, w):
    # forward padded transform
    L = x.size
    Lmax = L + 3 * w
    n = 2 ** np.ceil(np.log2(Lmax))
    X = np.fft.fft(x, n.astype(np.int))

    # generate kernel domain
    f = np.linspace(0, n-1, np.int(n)) / n
    f = np.concatenate((-f[0: np.int(n / 2 + 1)],
                        f[1: np.int(n / 2 - 1 + 1)][::-1]))

    # evaluate kernel
    K = np.exp(-0.5 * (w * 2 * np.pi * f) ** 2)

    # convolve and transform back from frequency domain
    y = np.real(np.fft.ifft(X * K, n))
    y = y[0:L]

    return y


def fftkernelWin(x, w, WinFunc):
    # forward padded transform
    L = x.size
    Lmax = L + 3 * w
    n = 2 ** np.ceil(np.log2(Lmax))
    X = np.fft.fft(x, n.astype(np.int))

    # generate kernel domain ###
    f = np.linspace(0, n-1, np.int(n)) / n
    f = np.concatenate((-f[0: np.int(n / 2 + 1)],
                        f[1: np.int(n / 2 - 1 + 1)][::-1]))
    t = 2 * np.pi * f

    # determine window function - evaluate kernel
    if WinFunc == 'Boxcar':
        a = 12**0.5 * w
        K = 2 * np.sin(a * t / 2) / (a * t)
        K[0] = 1
    elif WinFunc == 'Laplace':
        K = 1 / (1 + (w * 2 * np.pi * f)**2 / 2)
    elif WinFunc == 'Cauchy':
        K = np.exp(-w * np.abs(2 * np.pi * f))
    else:  # WinFunc == 'Gauss'
        K = np.exp(-0.5 * (w * 2 * np.pi * f)**2)

    # convolve and transform back from frequency domain
    y = np.real(np.fft.ifft(X * K, n))
    y = y[0:L]

    return y


def Gauss(x, w):
    y = 1 / (2 * np.pi)**2 / w * np.exp(-x**2 / 2 / w**2)
    return y


def Laplace(x, w):
    y = 1 / 2**0.5 / w * np.exp(-(2**0.5) / w / np.abs(x))
    return y


def Cauchy(x, w):
    y = 1 / (np.pi * w * (1 + (x / w)**2))
    return y


def Boxcar(x, w):
    a = 12**0.5 * w
    y = 1 / a
    y[np.abs(x) > a / 2] = 0
    return y


def logexpSSV(x):
    y = np.zeros(x.shape)
    y[x < 1e2] = np.log(1+np.exp(x[x < 1e2]))
    y[x >= 1e2] = x[x >= 1e2]
    return y


def ilogexpSSV(x):
    y = np.zeros(x.shape)
    y[x < 1e2] = np.log(np.exp(x[x < 1e2]) - 1)
    y[x >= 1e2] = x[x >= 1e2]
    return y
