"""PCA + clustering of the traces (scikit-learn).

Every trace becomes a small normalised 2D histogram (distance x log G), so
traces with different lengths can be compared. Same idea as Lemmer 2016
and Cabosart 2019.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, field

import numpy as np

from .traces import Trace


@dataclass
class ClusterSettings:
    feature: str = "2D histogram"      # "2D histogram" | "1D histogram" | "resampled trace"
    g_min: float = -6.0
    g_max: float = 0.5
    g_bins: int = 28
    z_min: float = -0.2
    z_max: float = 1.5
    z_bins: int = 28
    resample_points: int = 128
    pca_components: int = 10           # 0 = cluster on raw features
    method: str = "k-means"            # "k-means" | "gaussian mixture" | "agglomerative" | "spectral"
    n_clusters: int = 3
    random_state: int = 0
    sort_by: str = "conductance"       # "conductance" | "size"
    k_scan: tuple = (2, 8)

    def to_dict(self) -> dict:
        return asdict(self)


def features(traces: list[Trace], cs: ClusterSettings) -> np.ndarray:
    ge = np.linspace(cs.g_min, cs.g_max, cs.g_bins + 1)
    if cs.feature == "1D histogram":
        F = np.array([np.histogram(t.logG, ge)[0] for t in traces], float)
    elif cs.feature == "resampled trace":
        zz = np.linspace(cs.z_min, cs.z_max, cs.resample_points)
        F = np.array([np.interp(zz, t.z, np.clip(t.logG, cs.g_min, cs.g_max),
                                left=cs.g_max, right=cs.g_min) for t in traces])
        return F
    else:
        ze = np.linspace(cs.z_min, cs.z_max, cs.z_bins + 1)
        F = np.array([np.histogram2d(t.z, np.clip(t.logG, cs.g_min, cs.g_max - 1e-9), [ze, ge])[0].ravel()
                      for t in traces], float)
    s = F.sum(axis=1, keepdims=True)
    s[s == 0] = 1
    return F / s


@dataclass
class ClusterResult:
    labels: np.ndarray
    scores: np.ndarray                 # PCA scores (n_traces, n_components)
    explained: np.ndarray              # explained variance ratio per component
    n_clusters: int
    silhouette: float = np.nan
    k_scan: dict = field(default_factory=dict)   # k -> silhouette
    feature_shape: tuple = ()
    centroids: np.ndarray | None = None          # mean feature per cluster


def _cluster(X: np.ndarray, cs: ClusterSettings, k: int) -> np.ndarray:
    from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
    from sklearn.mixture import GaussianMixture

    if cs.method == "gaussian mixture":
        return GaussianMixture(k, covariance_type="full", random_state=cs.random_state,
                               n_init=3).fit_predict(X)
    if cs.method == "agglomerative":
        return AgglomerativeClustering(k, linkage="ward").fit_predict(X)
    if cs.method == "spectral":
        return SpectralClustering(k, random_state=cs.random_state, affinity="nearest_neighbors",
                                  n_neighbors=min(15, len(X) - 1)).fit_predict(X)
    return KMeans(k, n_init=10, random_state=cs.random_state).fit_predict(X)


def run(traces: list[Trace], cs: ClusterSettings, scan_k: bool = False) -> ClusterResult:
    from sklearn.decomposition import PCA
    from sklearn.metrics import silhouette_score

    if len(traces) < max(3, cs.n_clusters):
        raise ValueError(f"At least {max(3, cs.n_clusters)} traces are needed; {len(traces)} are available.")
    F = features(traces, cs)
    Fc = F - F.mean(axis=0)
    n_comp = min(max(2, cs.pca_components or 2), Fc.shape[0], Fc.shape[1])
    pca = PCA(n_components=n_comp, random_state=cs.random_state)
    scores = pca.fit_transform(Fc)
    X = scores if cs.pca_components else Fc

    labels = _cluster(X, cs, cs.n_clusters)
    sil = silhouette_score(X, labels) if len(set(labels)) > 1 else np.nan

    scan = {}
    if scan_k:
        lo, hi = cs.k_scan
        for k in range(max(2, lo), min(hi, len(traces) - 1) + 1):
            lab = _cluster(X, cs, k)
            if len(set(lab)) > 1:
                scan[k] = float(silhouette_score(X, lab))

    # relabel clusters 0..k-1 sorted by mean conductance (high -> low) or size
    uniq = np.unique(labels)
    if cs.sort_by == "size":
        key = {u: -np.count_nonzero(labels == u) for u in uniq}
    else:
        key = {u: -np.mean([np.median(traces[i].logG) for i in np.flatnonzero(labels == u)]) for u in uniq}
    order = sorted(uniq, key=lambda u: key[u])
    remap = {u: i for i, u in enumerate(order)}
    labels = np.array([remap[u] for u in labels])
    cent = np.array([F[labels == i].mean(axis=0) for i in range(len(order))])

    shape = (cs.z_bins, cs.g_bins) if cs.feature == "2D histogram" else (F.shape[1],)
    return ClusterResult(labels, scores, pca.explained_variance_ratio_, len(order), float(sil),
                         scan, shape, cent)
