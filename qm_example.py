import numpy as np

def empirical_cdf(data, prob):
    # estimate emperical cdf according to a user-defined probability
    data2 = data[~np.isnan(data)]
    if len(data2) > 0:
        ds = np.sort(data2)
        probreal = np.arange(len(data2)) / (len(data2) + 1)
        ecdf_out = np.interp(prob, probreal, ds)
    else:
        ecdf_out = np.nan * np.zeros(len(prob))
    return ecdf_out

def cdf_correction(prob_ref, value_ref, prob_raw, value_raw, value_tar):
    # _ref: probability (prob) and its corresponding cdf value (value) of the reference dataset
    # _raw: probability (prob) and its corresponding cdf value (value) of the raw/target dataset
    # value_tar: original time series of the target dataset
    prob_tar = np.interp(value_tar, value_raw, prob_raw)
    value_out = np.interp(prob_tar, prob_ref, value_ref)
    return value_out

# usage
data_ref = np.arange(10) # just an example
data_tar = np.arange(10) # just an example
prob = np.arange(0, 1, 0.1)

ecdf_ref = empirical_cdf(data_ref, prob)
ecdf_tar = empirical_cdf(data_tar, prob)

# Quantile-mapping or CDF-matching
qm_tar = cdf_correction(prob, ecdf_ref, prob, ecdf_tar, data_tar)
