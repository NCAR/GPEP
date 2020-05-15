import numpy as np
import sys
from scipy import io


def bma(D, obs, w0 = None):
    # the matlab version is used for Chen, Yang, et al. "Using Bayesian model averaging to estimate terrestrial
    # evapotranspiration in China." Journal of Hydrology 528 (2015): 537-549.
    # Ma, Yingzhao, et al. "Comprehensive evaluation of ensemble multi-satellite precipitation dataset using the
    # dynamic bayesian model averaging scheme over the Tibetan Plateau." Journal of hydrology 556 (2018): 634-644.

    # D: model data [samples, models]
    # obs: observation vector
    # w0: initial weight

    if np.ndim(D) == 1:
        sys.exit('Insufficient input models')
    ind = (~np.isnan(obs)) & (~np.isnan(np.sum(D, axis=1)))
    D = D[ind, :]
    obs = obs[ind]
    K = np.shape(D)[1]

    if len(obs)<2:
        print('Not enough sample inputs for bma')
        w = np.nan * np.ones(K)
        sigma = np.nan * np.zeros(K)
        sigma_s = np.nan
        return w, sigma, sigma_s

    N = np.shape(D)[0]

    # Start EM algorithm

    # Step 1 - Initial guess of sigma and weights
    sigma = np.zeros(K)
    sigma2 = np.zeros(K)
    for k in range(K):
        sigma[k] = np.std(D[:, k] - obs)
        sigma2[k] = np.var(D[:, k] - obs)

    if w0 == None:
        w = np.ones(K)/K
    else:
        w = w0
    sigma_s = np.mean(sigma)  # A single sigma

    # Step 2 - Compute probability density g(:,k) given initial guess of sigma
    g = np.zeros([N, K])
    for k in range(K):
        ## commented out:  g(:,k) = normpdf(obs,D(:,k),sigma(k));
        ## for n=1:N;
        ##    g(n,k) = normpdf(obs(n),D(n,k),sigma_s);
        ## end;
        ## Use least square as a proxy to likelihood:
        sls = (D[:, k] - obs) ** 2
        ## g(:,k) = exp(-sls/(2*sigma(k)^2))/(sqrt(2*pi)*sigma(k));
        g[:, k] = np.exp(-sls / (2 * (sigma_s ** 2))) / (np.sqrt(2 * np.pi) * sigma_s)
    g[g==0] = 1e-20
    gsum = np.sum(g, axis=1)
    gsum2 = np.tile(gsum, (K, 1)).T

    # Step 3 - Compute the initial log-likelihood;
    wtemp = np.tile(w, [N, 1])
    l = np.sum(wtemp * g, axis=1)
    L = np.sum(np.log(l))

    # Step4 - Iterate until maximum likelihood estimates are obtained

    iter = 0
    L_diff = 1e+5
    while L_diff >= 0.0001:
        # print('iter', iter)
        iter = iter + 1

        L_old = L
        w_old = w
        sigma_old = sigma
        sigma_s_old = sigma_s

        # Step 5 - Compute z(:, k) and w(k);

        z = g / gsum2
        w1 = np.mean(z, axis=0)

        zm = np.tile(np.max(z, axis=1), (K, 1)).T
        z2 = np.zeros([N, K])
        z2[z == zm] = 1
        z2sum = np.tile(np.sum(z2, axis=1), (K, 1)).T
        z2 = z2 / z2sum

        w2 = np.mean(z2, axis=0)
        w = w1  # w = w2; z = z2;

        # Step 6 - Compute updated sigma2;

        for k in range(K):
            sigma2[k] = np.sum(z[:, k] * ((obs - D[:, k]) ** 2))
        sigma = np.sqrt(sigma2)
        sigma2_s = np.sum(sigma2) / N
        sigma_s = np.sqrt(sigma2_s)

        # Step 7 - Update g(:, k);

        for k in range(K):
            # Compute probability density g(:, k)
            # # for n=1:N;
            # # g(n, k) = normpdf(obs(n), D(n, k), sigma_s);
            # # end;
            # # Use least square as a proxy to likelihood:
            sls = (D[:, k] - obs) ** 2
            g[:, k] = np.exp(-sls / (2 * (sigma_s ** 2))) / (np.sqrt(2 * np.pi) * sigma_s)
        g[g == 0] = 1e-20
        gsum = np.sum(g, axis=1)
        gsum2 = np.tile(gsum, (K, 1)).T

        # Step 8 - Compute the updated log-likelihood;

        wtemp = np.tile(w, [N, 1])
        l = np.sum(wtemp * g, axis=1)
        L = np.sum(np.log(l))

        # Step 9 - Compute the improvement in log - likelihood
        L_diff = L - L_old

        if L < L_old:
            w = w_old
            L = L_old
            sigma = sigma_old
            sigma_s = sigma_s_old
            break

    # Computes the BMA time series and variance
    # BMA_TS = np.zeros(N)
    # BMA_VAR = np.zeros(N)
    # for n in range(N):
    #     BMA_TS[n] = np.sum(w * D[n,:])
    #     # BMA_VAR[n] = sum(W * (D[n,:] - BMA_TS[n]) ** 2) + sum(W * SIGMA);
    #     BMA_VAR[n] = np.sum(w * (D[n,:] - BMA_TS[n])** 2) + sigma_s

    return w, sigma, sigma_s


if __name__ == "__main__":
    # datatemp = io.loadmat('bma_data.mat')
    # D = datatemp['drea']  # [:,K]: Kth model
    # obs = datatemp['dobs']  # observation: vector
    # obs = np.squeeze(obs)
    # w0 = datatemp['w0']  # initial weight
    # del datatemp
    # w, sigma, sigma_s = bma(D, obs)
    datatemp = np.load('bmadata.npz')
    dobs = datatemp['dobs']
    drea = datatemp['drea']
    drea=drea[dobs>0,:]
    dobs=dobs[dobs>0]
    del datatemp
    w, sigma, sigma_s = bma(drea, dobs)
    print(w)