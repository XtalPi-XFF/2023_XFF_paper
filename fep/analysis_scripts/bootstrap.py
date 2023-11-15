import numpy as np
import scipy
import sklearn.metrics
from typing import Union

# The code for bootstrap analysis was adapted from https://github.com/OpenFreeEnergy/cinnabar/blob/main/cinnabar/stats.py

def bootstrap_statistic(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    dy_true: Union[np.ndarray, None] = None,
    dy_pred: Union[np.ndarray, None] = None,
    ci: float = 0.95,
    statistic: str = "RMSE",
    nbootstrap: int = 1000,
) -> dict:

    """Compute mean and confidence intervals of specified statistic.

    Parameters
    ----------
    y_true : ndarray with shape (N,)
        True values
    y_pred : ndarray with shape (N,)
        Predicted values
    dy_true : ndarray with shape (N,) or None
        Errors of true values. If None, the values are assumed to have no errors
    dy_pred : ndarray with shape (N,) or None
        Errors of predicted values. If None, the values are assumed to have no errors
    ci : float, optional, default=0.95
        Interval for confidence interval (CI)
    statistic : str
        Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','KTAU','RAE']
    nbootstrap : int, optional, default=1000
        Number of bootstrap samples
    plot_type : str, optional, default='dG'
        'dG' or 'ddG'

    Returns
    -------
    rmse_stats : dict of float
        'mean' : mean RMSE
        'stderr' : standard error
        'low' : low end of CI
        'high' : high end of CI
    """

    def compute_statistic(y_true_sample: np.ndarray, y_pred_sample: np.ndarray, statistic: str):
        """Compute requested statistic.

        Parameters
        ----------
        y_true : ndarray with shape (N,)
            True values
        y_pred : ndarray with shape (N,)
            Predicted values
        statistic : str
            Statistic, one of ['RMSE', 'MUE', 'R2', 'rho','RAE','KTAU']

        """

        def calc_RAE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            MAE = sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
            mean = np.mean(y_true_sample)
            MAD = np.sum([np.abs(mean - i) for i in y_true_sample]) / float(len(y_true_sample))
            return MAE / MAD

        def calc_RRMSE(y_true_sample: np.ndarray, y_pred_sample: np.ndarray):
            rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample))
            mean_exp = np.mean(y_true_sample)
            mds = np.sum([(mean_exp - i) ** 2 for i in y_true_sample]) / float(len(y_true_sample))
            rrmse = np.sqrt(rmse**2 / mds)
            return rrmse

        if statistic == "RMSE":
            return np.sqrt(sklearn.metrics.mean_squared_error(y_true_sample, y_pred_sample))
        elif statistic == "MUE":
            return sklearn.metrics.mean_absolute_error(y_true_sample, y_pred_sample)
        elif statistic == "R2":
            slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(
                y_true_sample, y_pred_sample
            )
            return r_value**2
        elif statistic == "rho":
            return scipy.stats.pearsonr(y_true_sample, y_pred_sample)[0]
        elif statistic == "RAE":
            return calc_RAE(y_true_sample, y_pred_sample)
        elif statistic == "KTAU":
            return scipy.stats.kendalltau(y_true_sample, y_pred_sample)[0]
        else:
            raise Exception("unknown statistic '{}'".format(statistic))

    # not used?
    def unique_differences(x):
        """Compute all unique differences"""
        N = len(x)
        return np.array([(x[i] - x[j]) for i in range(N) for j in range(N) if (i != j)])

    if statistic == "RMSE":
        y_true = np.concatenate(y_true)
        y_pred = np.concatenate(y_pred)
        dy_pred = np.concatenate(dy_pred)

        if dy_true is None:
            dy_true = np.zeros_like(y_true)
        if dy_pred is None:
            dy_pred = np.zeros_like(y_pred)
        assert len(y_true) == len(y_pred)
        assert len(y_true) == len(dy_true)
        assert len(y_true) == len(dy_pred)
        sample_size = len(y_true)
        s_n = np.zeros(
            [nbootstrap], np.float64
        )  # s_n[n] is the statistic computed for bootstrap sample n
        for replicate in range(nbootstrap):
            y_true_sample = np.zeros_like(y_true)
            y_pred_sample = np.zeros_like(y_pred)
            for i, j in enumerate(
                np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)
            ):
                y_true_sample[i] = np.random.normal(loc=y_true[j], scale=np.fabs(dy_true[j]), size=1)
                y_pred_sample[i] = np.random.normal(loc=y_pred[j], scale=np.fabs(dy_pred[j]), size=1)
            s_n[replicate] = compute_statistic(y_true_sample, y_pred_sample, statistic)

        rmse_stats = dict()
        rmse_stats["mle"] = compute_statistic(y_true, y_pred, statistic)
        rmse_stats["stderr"] = np.std(s_n)
        rmse_stats["mean"] = np.mean(s_n)
        # TODO: Is there a canned method to do this?
        s_n = np.sort(s_n)
        low_frac = (1.0 - ci) / 2.0
        high_frac = 1.0 - low_frac
        rmse_stats["low"] = s_n[int(np.floor(nbootstrap * low_frac))]
        rmse_stats["high"] = s_n[int(np.ceil(nbootstrap * high_frac))]

    else:
        mle = 0.0
        ct_sum = 0
        s_n = np.zeros(
            [nbootstrap], np.float64
        )  # s_n[n] is the statistic computed for bootstrap sample n
        n_file = len(y_true)

        for yy_true, yy_pred, dyy_pred in zip(y_true, y_pred, dy_pred):
            sample_size = len(yy_true)
            ct_sum += sample_size
            dyy_true = np.zeros_like(yy_true)
            for replicate in range(nbootstrap):
                y_true_sample = np.zeros_like(yy_true)
                y_pred_sample = np.zeros_like(yy_pred)
                for i, j in enumerate(
                    np.random.choice(np.arange(sample_size), size=[sample_size], replace=True)
                ):
                    y_true_sample[i] = np.random.normal(loc=yy_true[j], scale=np.fabs(dyy_true[j]), size=1)
                    y_pred_sample[i] = np.random.normal(loc=yy_pred[j], scale=np.fabs(dyy_pred[j]), size=1)
                s_n[replicate] += sample_size * compute_statistic(y_true_sample, y_pred_sample, statistic)
            mle += sample_size * compute_statistic(yy_true, yy_pred, statistic)
        s_n = [s/ct_sum for s in s_n]
        mle /= ct_sum
        mean = np.mean(s_n)
        rmse_stats = dict()
        rmse_stats["mle"] = mle
        rmse_stats["mean"] = mean
        # TODO: Is there a canned method to do this?
        s_n = np.sort(s_n)
        low_frac = (1.0 - ci) / 2.0
        high_frac = 1.0 - low_frac
        rmse_stats["low"] = s_n[int(np.floor(nbootstrap * low_frac))]
        rmse_stats["high"] = s_n[int(np.ceil(nbootstrap * high_frac))]

    return rmse_stats

if __name__ == "__main__":
    import pandas as pd
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-files', help='files containing the fep results', nargs='+')
    parser.add_argument('-type', help='analyze type, dg or ddg')
    args = parser.parse_args()
    files = args.files
    tp = args.type

    if tp == "ddg":
        val = "CCC"
        val_err = "CCC_err"
    else:
        val = "Pred dG"
        val_err = "Pred dG err"

    sims = []
    sims_err = []
    exps = []
    for f in files:
        data = pd.read_csv(f)
        try:
            sims.append(data[val].values)
            sims_err.append(data[val_err].values)
        except:
            sims.append(data["FEP"].values)
            sims_err.append(data["FEP_err"].values)
        exps.append(data["EXP"].values)

    statistics = ['RMSE','R2', "KTAU"]
    #statistics = ['RMSE']
    string = []
    for statistic in statistics:
        s = bootstrap_statistic(exps, sims, None, sims_err, statistic=statistic)
        string.append(
            f"{statistic + ':':5s}{s['mle']:5.2f} [95%: {s['low']:5.2f}, {s['high']:5.2f}]")
    print(string)
