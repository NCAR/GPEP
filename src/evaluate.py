# evaluate estimates using some common metrics

import numpy as np

def evaluate_allpoint(obs, est, threshold=0.1):
    # obs/est: [number of points, sample]
    nstn = obs.shape[0]
    metric_values = np.nan * np.zeros([nstn, 16])
    for i in range(nstn):
        metric_values[i, :], metric_names = evaluate(obs[i, :], est[i, :], threshold)
    return metric_values, metric_names


def evaluate(obs, est, threshold=0.1):
    # threshold: used for probability evaluation

    obs, est = data_preprocess(obs, est)

    metric_names = ['CC', 'ME', 'RB', 'MAE', 'ARB', 'RMSE', 'NSE', 'KGE_G2009', 'KGE_K2012', 'KGE_P2018', 'KGE_T2021', 'POD', 'FOH', 'FAR', 'CSI', 'HSS']
    metric_values = np.nan * np.zeros(16)
    if len(obs) > 1:
        metric_values[0] = cal_CC(obs, est, preprocess=False)
        metric_values[1] = cal_ME(obs, est, preprocess=False)
        metric_values[2] = cal_RB(obs, est, preprocess=False)
        metric_values[3] = cal_MAE(obs, est, preprocess=False)
        metric_values[4] = cal_ARB(obs, est, preprocess=False)
        metric_values[5] = cal_RMSE(obs, est, preprocess=False)
        metric_values[6] = cal_NSE(obs, est, preprocess=False)
        metric_values[7] = cal_KGE_G2009(obs, est, preprocess=False)['KGE']
        metric_values[8] = cal_KGE_K2012(obs, est, preprocess=False)['KGE']
        metric_values[9] = cal_KGE_P2018(obs, est, preprocess=False)['KGE']
        metric_values[10] = cal_KGE_T2021(obs, est, preprocess=False)['KGE']

        if ~np.isnan(threshold):
            icont = cal_contingency(obs, est, threshold)
            metric_values[11] = icont['POD']
            metric_values[12] = icont['FOH']
            metric_values[13] = icont['FAR']
            metric_values[14] = icont['CSI']
            metric_values[15] = icont['HSS']

    return metric_values, metric_names


def data_preprocess(obs, est):
    ind_invalid = np.isnan(obs+est) | np.isinf(obs+est)
    obs = obs[~ind_invalid]
    est = est[~ind_invalid]
    return obs, est


def cal_CC(obs, est, preprocess=True):
    # correlation coefficient
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        temp = np.corrcoef(obs, est)
        CC = temp[0][1]
    else:
        CC = np.nan
    return CC


def cal_ME(obs, est, preprocess=False):
    # mean error
    if not preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        ME = np.nanmean(est - obs)
    else:
        ME = np.nan
    return ME


def cal_MAE(obs, est, preprocess=False):
    # mean absolute error
    if not preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        MAE = np.nanmean(np.abs(est - obs))
    else:
        MAE = np.nan
    return MAE


def cal_RB(obs, est, preprocess=False):
    # relative bias
    if not preprocess:
        obs, est = data_preprocess(obs, est)
    temp1 = np.nansum(est)
    temp2 = np.nansum(obs)
    if temp2 != 0 and len(obs) > 1:
        RB = (temp1 - temp2) / temp2
    else:
        RB = np.nan
    return RB


def cal_ARB(obs, est, preprocess=False):
    # absolute relative error
    if not preprocess:
        obs, est = data_preprocess(obs, est)
    obss = np.nansum(obs)
    if obss != 0 and len(obs) > 1:
        ARB = np.nansum(np.abs(est - obs)) / obss
    else:
        ARB = np.nan
    return ARB


def cal_RMSE(obs, est, preprocess=True):
    # root mean square error
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        RMSE = np.sqrt(np.square(obs - est).sum() / len(obs))
    else:
        RMSE = np.nan
    return RMSE


def cal_NSE(obs, est, preprocess=True):
    # Nash-Sutcliffe Efficiency (Nash and Sutcliffe 1970 -
    # https://doi.org/10.1016/0022-1694(70)90255-6)
    # a common index for hydrological obs
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        NSE = 1 - (np.sum((obs - est) ** 2, axis=0, dtype=np.float64) /
                   np.sum((obs - np.mean(obs)) ** 2, dtype=np.float64))
    else:
        NSE = np.nan
    return NSE


def cal_KGE_G2009(obs, est, preprocess=True):
    # Original Kling-Gupta Efficiency (Gupta et al. 2009 -
    # https://doi.org/10.1016/j.jhydrol.2009.08.003)
    # calculate error in timing and dynamics r (Pearson's correlation
    # coefficient) which can also be calculated using np.corrcoef
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs):
        est_mean = np.mean(est, axis=0, dtype=np.float64)
        obs_mean = np.mean(obs, dtype=np.float64)
        r = np.sum((est - est_mean) * (obs - obs_mean), axis=0, dtype=np.float64) / \
            np.sqrt(np.sum((est - est_mean) ** 2, axis=0, dtype=np.float64) *
                    np.sum((obs - obs_mean) ** 2, dtype=np.float64))
        # calculate error in sestad of flow alpha
        alpha = np.std(est, axis=0) / np.std(obs, dtype=np.float64)
        # calculate error in volume beta (bias of mean discharge)
        beta = np.sum(est, axis=0, dtype=np.float64) / \
            np.sum(obs, dtype=np.float64)
        # calculate the Kling-Gupta Efficiency KGE
        KGE = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        KGEgroup = {'KGE': KGE, 'r': r, 'alpha': alpha, 'beta': beta}
    else:
        KGEgroup = {'KGE': np.nan, 'r': np.nan,
                    'alpha': np.nan, 'beta': np.nan}
    return KGEgroup  # or output np.vstack((KGE, r, alpha, beta))


def cal_KGE_K2012(obs, est, preprocess=True):
    # Modified Kling-Gupta Efficiency (Kling et al. 2012 -
    # https://doi.org/10.1016/j.jhydrol.2012.01.011)
    # calculate error in timing and dynamics r (Pearson's correlation
    # coefficient)
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        est_mean = np.mean(est, axis=0, dtype=np.float64)
        obs_mean = np.mean(obs, dtype=np.float64)
        r = np.sum((est - est_mean) * (obs - obs_mean), axis=0, dtype=np.float64) / \
            np.sqrt(np.sum((est - est_mean) ** 2, axis=0, dtype=np.float64) *
                    np.sum((obs - obs_mean) ** 2, dtype=np.float64))
        # calculate error in sestad of flow gamma (avoiding cross correlation with
        # bias by dividing by the mean)
        gamma = (np.std(est, axis=0, dtype=np.float64) / est_mean) / \
            (np.std(obs, dtype=np.float64) / obs_mean)
        # calculate error in volume beta (bias of mean discharge)
        beta = np.mean(est, axis=0, dtype=np.float64) / \
            np.mean(obs, axis=0, dtype=np.float64)
        # calculate the modified Kling-Gupta Efficiency KGE'
        KGE = 1 - np.sqrt((r - 1) ** 2 + (gamma - 1) ** 2 + (beta - 1) ** 2)
        KGEgroup = {'KGE': KGE, 'r': r, 'gamma': gamma, 'beta': beta}
    else:
        KGEgroup = {'KGE': np.nan, 'r': np.nan,
                    'gamma': np.nan, 'beta': np.nan}
    return KGEgroup  # or output np.vstack((KGE, r, gamma, beta))


def cal_KGE_P2018(obs, est, preprocess=True):
    # Non-Parametric Kling-Gupta Efficiency (Pool et al. 2018 -
    # https://doi.org/10.1080/02626667.2018.1552002)
    # calculate error in timing and dynamics r (Spearman's correlation
    # coefficient)
    if preprocess:
        obs, est = data_preprocess(obs, est)
    if len(obs) > 1:
        sim_rank = np.argsort(np.argsort(est, axis=0), axis=0)
        obs_rank = np.argsort(np.argsort(obs, axis=0), axis=0)
        r = np.sum((obs_rank - np.mean(obs_rank, axis=0, dtype=np.float64)) *
                   (sim_rank - np.mean(sim_rank, axis=0, dtype=np.float64)), axis=0) / \
            np.sqrt(np.sum((obs_rank - np.mean(obs_rank, axis=0, dtype=np.float64)) ** 2, axis=0) *
                    (np.sum((sim_rank - np.mean(sim_rank, axis=0, dtype=np.float64)) ** 2, axis=0)))
        # calculate error in timing and dynamics alpha (flow duration curve)
        sim_fdc = np.sort(est / (est.shape[0] * np.mean(
            est, axis=0, dtype=np.float64)), axis=0)
        obs_fdc = np.sort(
            obs / (obs.shape[0] * np.mean(obs, axis=0, dtype=np.float64)), axis=0)
        alpha = 1 - 0.5 * np.sum(np.abs(sim_fdc - obs_fdc), axis=0)
        # calculate error in volume beta (bias of mean discharge)
        beta = np.mean(est, axis=0) / \
            np.mean(obs, axis=0, dtype=np.float64)
        # calculate the non-parametric Kling-Gupta Efficiency KGEnp
        KGE = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2)
        KGEgroup = {'KGE': KGE, 'r': r, 'alpha': alpha, 'beta': beta}
    else:
        KGEgroup = {'KGE': np.nan, 'r': np.nan,
                    'alpha': np.nan, 'beta': np.nan}
    return KGEgroup  # or output np.vstack((KGE, r, alpha, beta))


def cal_KGE_T2021(obs, est, preprocess=True):
    # KGE" proposed by Tang et al 2021
    # https://doi.org/10.1175/JCLI-D-21-0067.1
    if preprocess:
        obs, est = data_preprocess(obs, est)
    ind_nan = np.isnan(obs) | np.isnan(est)
    obs = obs[~ind_nan]
    est = est[~ind_nan]
    if len(obs)>2:
        est_mean = np.mean(est, axis=0, dtype=np.float64)
        obs_mean = np.mean(obs, axis=0, dtype=np.float64)
        r = np.sum((est - est_mean) * (obs - obs_mean), axis=0, dtype=np.float64) / \
            np.sqrt(np.sum((est - est_mean) ** 2, axis=0, dtype=np.float64) *
                    np.sum((obs - obs_mean) ** 2, dtype=np.float64))
        # calculate error in sestad of flow alpha
        alpha = np.std(est, axis=0) / np.std(obs, dtype=np.float64)
        # calculate error in volume beta (bias of mean discharge)
        beta = (np.mean(est) - np.mean(obs)) / np.std(obs)
        # calculate the Kling-Gupta Efficiency KGE
        KGE = 1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta) ** 2)
        KGEgroup = {'KGE': KGE, 'r': r, 'alpha': alpha, 'beta': beta}
    else:
        KGEgroup = {'KGE': np.nan, 'r': np.nan,
                    'alpha': np.nan, 'beta': np.nan}
    return KGEgroup  # or output np.vstack((KGE, r, alpha, beta))


def cal_contingency(obs, est, Tre=0.1):
    # Tre: rain/no rain threshold
    # POD(Probability of Detection),FOH(frequency of hit)
    # FAR(False Alarm Ratio), CSI(Critical Success Index)
    # HSS(Heidke skillscore),Ebert et al. [2007]
    if np.isnan(Tre):
        POD = np.nan
        FAR = np.nan
        FOH = np.nan
        CSI = np.nan
        HSS = np.nan
    else:
        if len(obs) > 1:
            n11 = np.sum((obs > Tre) & (est > Tre))
            n10 = np.sum((obs <= Tre) & (est > Tre))
            n01 = np.sum((obs > Tre) & (est <= Tre))
            n00 = np.sum((obs <= Tre) & (est <= Tre))
        try:
            POD = n11 / (n11 + n01)
        except:
            POD = np.nan
        try:
            FOH = n11 / (n11 + n10)
            FAR = n10 / (n11 + n10)
        except:
            FOH = np.nan
            FAR = np.nan
        try:
            CSI = n11 / (n11 + n01 + n10)
        except:
            CSI = np.nan
        try:
            HSS = 2 * (n11 * n00 - n10 * n01) / ((n11 + n01) *
                                                 (n01 + n00) + (n11 + n10) * (n10 + n00))
        except:
            HSS = np.nan
    contingency_group = {'POD': POD, 'FOH': FOH, 'FAR': FAR,
                         'CSI': CSI, 'HSS': HSS}
    return contingency_group


if __name__ == '__main__':
    a = np.random.rand(100)
    b = np.random.rand(100)
    metric_values, metric_names = evaluate(a, b, 0.1)
    print(metric_names, '\n', metric_values)
