import numpy as np
import sys

D = []  # [:,K]: Kth model
obs = []  # observation: vector
w0 = []  # initial weight

if np.ndim(D) == 1:
    sys.exit('Insufficient input models')

N, K = np.shape(D)

# Start EM algorithm

# Step 1 - Initial guess of sigma and weights
sigma = np.zeros(K)
sigma2 = np.zeros(K)
for k in range(K):
    sigma[k] = np.std(D[:, k] - obs)
    sigma2[k] = np.var(D[:, k] - obs)

# w=np.ones(K)/K;
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
gsum = np.sum(g,axis=1)
gsum2 = np.tile(gsum,(K,1)).T

# Step 3 - Compute the initial log-likelihood;
wtemp = np.tile(w,[N,1])
l = np.sum(w * g, axis=1)
L = np.sum(np.log(l))

# Step4 - Iterate until maximum likelihood estimates are obtained

iter = 0
L_diff = 1e+5
while L_diff >= 0.0001:

    iter = iter + 1

    L_old = L
    w_old = w
    sigma_old = sigma
    sigma_s_old = sigma_s

    # Step 5 - Compute z(:, k) and w(k);

    z = g / gsum2
    w1 = np.mean(z, axis=0)

    z2 = z.copy()
    for j in range(N):
        a = [z2(j,:) >= max(z2(j,:))]
        z2(j,:)=a
        z2(j,:)= z2(j,:) / sum(z2(j,:))

    sum(z2);
    w2 = sum(z2) / N;
    w = w1; # # w = w2;
    z = z2;

    # Step
    6 - Compute
    updated
    sigma2;

    for k=1:K;
    sigma2(k) = sum(z(:, k).*(obs - D(:, k)).^ 2);
    end;
    sigma = sqrt(sigma2);
    sigma2_s = sum(sigma2) / N;
    sigma_s = sqrt(sigma2_s);

    # Step
    7 - Update
    g(:, k);

    for k=1:K;
    # Compute
    probability
    density
    g(:, k)
    # # for n=1:N;
    # # g(n, k) = normpdf(obs(n), D(n, k), sigma_s);
    # # end;
    # # Use
    least
    square as a
    proxy
    to
    likelihood:
    sls = (D(:, k) - obs ).^ 2;
    g(:, k) = exp(-sls / (2 * sigma_s ^ 2)) / (sqrt(2 * pi) * sigma_s);
    end;
    gsum = (sum(g'))';

    # Step 8 - Compute the updated log-likelihood;

    for n=1:N;
    l(n) = sum(w. * g(n,:));
    end;
    L = sum(log(l));

    # Step
    9 - Compute
    the
    improvement in log - likelihood
    L_diff = L - L_old;

    if L < L_old;
    w = w_old;
    L = L_old;
    sigma = sigma_old;
    sigma_s = sigma_s_old;
    break;
    end;

# End iteration
