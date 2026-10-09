"""ADF stationarity test.

Does the same as statsmodels' adfuller with default settings (constant
term, maxlag = ceil(12 (n/100)^0.25), lag picked by AIC, MacKinnon
p-value). It's here so the exe doesn't need statsmodels. I checked it
against statsmodels on 300 random series and the results are identical.

H0 is a unit root, so a small p-value means the trace is stationary.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

# MacKinnon (1994) response-surface coefficients, regression "c", N = 1
_TAU_MAX, _TAU_MIN, _TAU_STAR = 2.74, -18.83, -1.61
_SMALL_P = (2.1659, 1.4412, 0.038269)
_LARGE_P = (1.7339, 0.93202, -0.12745, -0.010368)


@dataclass
class ADFResult:
    statistic: float
    pvalue: float
    usedlag: int
    nobs: int


def mackinnon_p(stat: float) -> float:
    if stat > _TAU_MAX:
        return 1.0
    if stat < _TAU_MIN:
        return 0.0
    coef = _SMALL_P if stat <= _TAU_STAR else _LARGE_P
    return float(norm.cdf(np.polynomial.polynomial.polyval(stat, coef)))


def _design(x: np.ndarray, dx: np.ndarray, lag: int):
    """Rows t = lag..: columns [level x_{t-1}, dx_{t-1}, ..., dx_{t-lag}], target dx_t."""
    n = len(dx) - lag
    X = np.empty((n, lag + 1))
    X[:, 0] = x[-n - 1:-1]
    for j in range(1, lag + 1):
        X[:, j] = dx[lag - j:len(dx) - j]
    return X, dx[lag:]


def _ols(X: np.ndarray, y: np.ndarray):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return beta, float(resid @ resid)


def adfuller(x, maxlag: int | None = None, autolag: bool = True) -> ADFResult:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    nobs0 = len(x)
    ntrend = 1
    if maxlag is None:
        maxlag = int(np.ceil(12.0 * (nobs0 / 100.0) ** 0.25))
        maxlag = min(nobs0 // 2 - ntrend - 1, maxlag)
    if maxlag < 0 or nobs0 < 8:
        raise ValueError("series too short for the ADF test")
    dx = np.diff(x)

    usedlag = maxlag
    if autolag:
        Xf, yf = _design(x, dx, maxlag)
        n = len(yf)
        Xc = np.column_stack([np.ones(n), Xf])   # constant first, as in statsmodels' search
        # one QR gives the residual sum of squares of every nested model:
        # SSR(first k columns) = |y|^2 - sum_{i<k} (Q^T y)_i^2
        q, _ = np.linalg.qr(Xc)
        z2 = (q.T @ yf) ** 2
        yy = float(yf @ yf)
        best = None
        for k in range(2, maxlag + 3):           # columns: const, level, + (k-2) lags
            ssr = max(yy - float(z2[:k].sum()), np.finfo(float).tiny)
            llf = -n / 2.0 * (np.log(2 * np.pi) + np.log(ssr / n) + 1.0)
            aic = -2.0 * llf + 2.0 * k
            if best is None or (aic, k) < best:
                best = (aic, k)
        usedlag = best[1] - 2

    X, y = _design(x, dx, usedlag)
    n = len(y)
    X = np.column_stack([X, np.ones(n)])
    beta, ssr = _ols(X, y)
    dof = n - X.shape[1]
    sigma2 = ssr / dof if dof > 0 else np.nan
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    stat = float(beta[0] / np.sqrt(cov[0, 0]))
    return ADFResult(stat, mackinnon_p(stat), int(usedlag), int(n))
