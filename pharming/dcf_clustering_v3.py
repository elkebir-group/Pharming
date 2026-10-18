import argparse
from collections import defaultdict
import pickle
import multiprocessing
import numpy as np
import pandas as pd
from scipy.special import gammaln, xlogy, xlog1py

from . import clonelib
from .utils import timeit_decorator
from .genotype_tree import GenotypeTree

# ---------------------------------------------------------------------------
# Convergence / numerical constants
# ---------------------------------------------------------------------------

NORMALIZED_TOLERANCE = 1e-4     # convergence threshold per SNV
INVALID_LOG_LIKELIHOOD = -1e11  # penalty for an out-of-range (tree, DCF) pair





# ---------------------------------------------------------------------------
# DCF-to-VAF parameter extraction
# ---------------------------------------------------------------------------

def tree_dcf_parameters(tree, cn_prob):
    """Return (alpha, beta, vmin, vplus) for the linear map v = alpha*dcf + beta.

    Matches GenotypeTree.dcf_to_v() exactly:
        alpha = m_star / F
        beta  = (1/F) * sum((mut_copies(g) - m_star) * cn_prob[(g[0],g[1])]
                            for g in desc_genotypes
                            if (g[0],g[1]) != split_geno_cn)

    Note: beta != vmin in general.  vmin sums over gamma (all nodes) with
    weight mut_copies(g), whereas beta sums over desc_genotypes only with
    weight (mut_copies(g) - m_star).

    Returns (0, 0, nan, nan) for trees that have no SNV placement node, which
    marks every DCF value as invalid so those trees are never selected.
    """
    if not hasattr(tree, 'm_star') or tree.m_star is None:
        return 0.0, 0.0, np.nan, np.nan

    F = tree.compute_F(cn_prob)
    if F == 0:
        return 0.0, 0.0, np.nan, np.nan

    split_cn = (tree.split_geno[0], tree.split_geno[1])
    alpha = tree.m_star / F
    beta = sum(
        (tree.mut_copies(g) - tree.m_star) * cn_prob[(g[0], g[1])]
        for g in tree.desc_genotypes
        if (g[0], g[1]) != split_cn
    ) / F

    vmin = tree.v_minus(cn_prob)
    vplus = tree.v_plus(cn_prob, vmin)

    return float(alpha), float(beta), float(vmin), float(vplus)


# ---------------------------------------------------------------------------
# Sufficient-statistics group for one (tree, SNV-list) assignment
# ---------------------------------------------------------------------------

def make_dcf_group(tree, snvs, cn_prob, data):
    """Reduce one (tree, SNV-list) pair to the scalars needed by the optimizer.

    The group stores:
      alpha, beta   — linear map  v = alpha * dcf + beta
      valid_low     — smallest DCF that keeps v in [vmin, vplus]  = (vmin-beta)/alpha
      valid_high    — largest  DCF that keeps v in [vmin, vplus]  = (vplus-beta)/alpha
      alt, ref      — summed read counts across SNVs and cells
      log_choose    — sum of log C(d, a) terms (constant w.r.t. DCF)
      size          — number of SNVs (used for the invalid-tree penalty)
    """
    alt = np.asarray(data.var_marg(snvs), dtype=float)
    total = np.asarray(data.total_marg(snvs), dtype=float)
    ref = total - alt

    alpha, beta, vmin, vplus = tree_dcf_parameters(tree, cn_prob)

    # Solve v = alpha*dcf + beta for the DCF bounds that keep v in [vmin, vplus].
    # beta != vmin in general (see tree_dcf_parameters docstring).
    if np.isnan(vmin) or np.isnan(vplus) or alpha == 0:
        valid_low = np.nan
        valid_high = np.nan
    else:
        valid_low = float((vmin - beta) / alpha)
        valid_high = float((vplus - beta) / alpha)

    return {
        "alpha": float(alpha),
        "beta": float(beta),
        "valid_low": valid_low,
        "valid_high": valid_high,
        "alt": float(alt.sum()),
        "ref": float(ref.sum()),
        "log_choose": float(np.sum(
            gammaln(total + 1) - gammaln(alt + 1) - gammaln(ref + 1)
        )),
        "size": int(alt.size),
    }


# ---------------------------------------------------------------------------
# Group-level likelihood and derivative
# ---------------------------------------------------------------------------

def _dcf_group_is_valid(group, dcf, tol=1e-12):
    lo = group["valid_low"]
    hi = group["valid_high"]
    if np.isnan(lo) or np.isnan(hi):
        return False
    return lo - tol <= dcf <= hi + tol


def _dcf_group_log_likelihood(group, dcf):
    if not _dcf_group_is_valid(group, dcf):
        return group["size"] * INVALID_LOG_LIKELIHOOD

    v = group["alpha"] * dcf + group["beta"]
    # Mirror the [0, 1] check in vectorized_tree_cluster_likelihood so that
    # malformed numerical input cannot produce nan in xlogy / xlog1py and
    # corrupt candidate comparisons.
    if v < 0.0 or v > 1.0:
        return group["size"] * INVALID_LOG_LIKELIHOOD
    # xlogy(0, 0) = 0 and xlog1py(0, -1) = 0, so boundary cases where
    # alt=0 or ref=0 are handled correctly without clipping.
    return float(
        group["log_choose"]
        + xlogy(group["alt"], v)
        + xlog1py(group["ref"], -v)
    )


def _dcf_group_derivative(group, dcf):
    """d/d(dcf) of the binomial log-likelihood for one group.

    Clipping v here is intentional: the bisection solver evaluates the
    derivative strictly inside a validity interval (with a margin), so v
    should never actually reach 0 or 1.  The clip is purely defensive and
    does not affect the likelihood values returned by _dcf_group_log_likelihood.
    """
    v = group["alpha"] * dcf + group["beta"]
    v = np.clip(v, 1e-15, 1 - 1e-15)
    return float(group["alpha"] * (
        group["alt"] / v - group["ref"] / (1.0 - v)
    ))


def _exact_dcf_log_likelihood(groups, dcf):
    return sum(_dcf_group_log_likelihood(g, dcf) for g in groups)


# ---------------------------------------------------------------------------
# Bisection optimizer over one concave validity interval
# ---------------------------------------------------------------------------

def _solve_dcf_interval(groups, left, right, steps=64):
    """Find the maximum of the concave log-likelihood on (left, right)."""
    width = right - left
    margin = min(1e-12, width / 4)
    lo = left + margin
    hi = right - margin

    if lo >= hi:
        return (left + right) / 2.0

    d_lo = sum(_dcf_group_derivative(g, lo) for g in groups)
    d_hi = sum(_dcf_group_derivative(g, hi) for g in groups)

    # Derivative non-positive at left boundary → likelihood already decreasing;
    # optimum is at the left boundary.
    if np.isnan(d_lo) or d_lo <= 0:
        return lo
    # Derivative non-negative at right boundary → likelihood still increasing;
    # optimum is at the right boundary.
    if np.isnan(d_hi) or d_hi >= 0:
        return hi

    for _ in range(steps):
        mid = (lo + hi) / 2.0
        d_mid = sum(_dcf_group_derivative(g, mid) for g in groups)
        if np.isnan(d_mid):
            # v hit a clamp boundary; current interval is degenerate
            break
        if d_mid > 0:
            lo = mid
        else:
            hi = mid

    return (lo + hi) / 2.0


# ---------------------------------------------------------------------------
# Global DCF optimizer across all validity intervals
# ---------------------------------------------------------------------------

def maximize_dcf_groups(groups, initial, lower=0.025, upper=1.0):
    """Find the DCF in [lower, upper] that maximises the total log-likelihood.

    The log-likelihood is piecewise concave: each (tree, SNV) group contributes
    a concave piece on its validity interval and a flat penalty outside it.
    We enumerate all interval endpoints, bisect within each active sub-interval,
    and return the global best.
    """
    if len(groups) == 0:
        return float(np.clip(initial, lower, upper))

    # Collect all validity-interval endpoints that fall inside [lower, upper]
    endpoints = {float(lower), float(upper)}
    for g in groups:
        if not np.isnan(g["valid_low"]) and lower < g["valid_low"] < upper:
            endpoints.add(g["valid_low"])
        if not np.isnan(g["valid_high"]) and lower < g["valid_high"] < upper:
            endpoints.add(g["valid_high"])
    endpoints = sorted(endpoints)

    candidates = list(endpoints)
    for left, right in zip(endpoints[:-1], endpoints[1:]):
        mid = (left + right) / 2.0
        active = [g for g in groups if _dcf_group_is_valid(g, mid)]
        if active:
            candidates.append(_solve_dcf_interval(active, left, right))
        else:
            candidates.append(mid)

    return max(
        (c for c in candidates if lower <= c <= upper),
        key=lambda dcf: _exact_dcf_log_likelihood(groups, dcf),
        default=float(np.clip(initial, lower, upper)),
    )


# ---------------------------------------------------------------------------
# Build sufficient-statistics groups for one cluster
# ---------------------------------------------------------------------------

def dcf_groups_for_cluster(clust_group_map, q, data,
                            cn_prob_cache=None,
                            alt_marginal=None,
                            total_marginal=None):
    """Collect make_dcf_group dicts for every (segment, tree, SNV-list) triple
    assigned to cluster q."""
    groups = []
    for ell in clust_group_map:
        if q not in clust_group_map[ell]:
            continue

        cn_prob = cn_prob_cache[ell] if cn_prob_cache is not None \
            else dict(data.cn_proportions(ell))
        if (1, 1) not in cn_prob:
            cn_prob[(1, 1)] = 0

        for tree, snvs in clust_group_map[ell][q]:
            if len(snvs) == 0:
                continue

            if alt_marginal is not None and total_marginal is not None:
                # Use pre-marginalised counts stored on self
                alt = np.asarray(alt_marginal[snvs], dtype=float)
                total = np.asarray(total_marginal[snvs], dtype=float)
                ref = total - alt
                alpha, beta, vmin, vplus = tree_dcf_parameters(tree, cn_prob)
                if np.isnan(vmin) or np.isnan(vplus) or alpha == 0:
                    valid_low, valid_high = np.nan, np.nan
                else:
                    valid_low = float((vmin - beta) / alpha)
                    valid_high = float((vplus - beta) / alpha)
                groups.append({
                    "alpha": float(alpha),
                    "beta": float(beta),
                    "valid_low": valid_low,
                    "valid_high": valid_high,
                    "alt": float(alt.sum()),
                    "ref": float(ref.sum()),
                    "log_choose": float(np.sum(
                        gammaln(total + 1)
                        - gammaln(alt + 1)
                        - gammaln(ref + 1)
                    )),
                    "size": int(alt.size),
                })
            else:
                groups.append(make_dcf_group(tree, snvs, cn_prob, data))

    return groups


# ---------------------------------------------------------------------------
# Vectorised assignment scoring
# ---------------------------------------------------------------------------

def vectorized_tree_cluster_likelihood(dcfs, tree_parameters, alt, total):
    """Score every (cluster, tree, SNV) triple in one NumPy operation.

    Parameters
    ----------
    dcfs : (K,) array of cluster DCF values
    tree_parameters : list of (alpha, beta, vmin, vplus) tuples, length T
    alt : (M,) array of variant read counts
    total : (M,) array of total read counts

    Returns
    -------
    scores : (K, T, M) array of log-likelihoods
    """
    dcfs = np.asarray(dcfs, dtype=float)           # (K,)
    alt = np.asarray(alt, dtype=float)             # (M,)
    total = np.asarray(total, dtype=float)         # (M,)
    ref = total - alt

    params = np.asarray(tree_parameters, dtype=float)  # (T, 4)
    alpha = params[:, 0]   # (T,)
    beta = params[:, 1]    # (T,)
    vmin = params[:, 2]    # (T,)
    vplus = params[:, 3]   # (T,)

    # v[k, t] = dcfs[k] * alpha[t] + beta[t]
    v = dcfs[:, None] * alpha[None, :] + beta[None, :]   # (K, T)

    # Validity mask: v must be in [vmin, vplus] and in [0, 1].
    # Trees with nan parameters are always invalid.
    valid = (
        (v >= vmin[None, :])
        & (v <= vplus[None, :])
        & (v >= 0.0)
        & (v <= 1.0)
        & ~np.isnan(vmin[None, :])
        & ~np.isnan(vplus[None, :])
    )   # (K, T)

    # Binomial log-likelihood components.
    # xlogy and xlog1py handle zero-boundary cases correctly (xlogy(0,0)=0),
    # so we do not clip v — clipping would silently make impossible observations
    # (e.g. v=0 with alt>0) merely unlikely rather than -inf.
    log_choose = gammaln(total + 1) - gammaln(alt + 1) - gammaln(ref + 1)  # (M,)

    # Broadcast to (K, T, M)
    likelihood = (
        log_choose[None, None, :]
        + xlogy(alt[None, None, :], v[:, :, None])
        + xlog1py(ref[None, None, :], -v[:, :, None])
    )

    return np.where(valid[:, :, None], likelihood, INVALID_LOG_LIKELIHOOD)


# ---------------------------------------------------------------------------
# Main clustering class
# ---------------------------------------------------------------------------

class DCF_Clustering:
    def __init__(self, nrestarts=25, seed=1026, verbose=False,
                 cna_restriction=True, rng=None, iterations=100,
                 nsegs=15, thresh_prop=0.05, max_cn_states=4,
                 tolerance=NORMALIZED_TOLERANCE,
                 assignment_chunk_size=10_000):

        self.nrestarts = nrestarts
        self.max_cn_states = max_cn_states
        self.rng = rng if rng is not None else np.random.default_rng(seed)
        self.iterations = iterations
        self.verbose = verbose
        self.cna_restriction = cna_restriction
        self.nsegs = nsegs
        self.thresh_prop = thresh_prop
        self.tolerance = tolerance
        self.assignment_chunk_size = assignment_chunk_size

        print(f"Verbose on: {self.verbose}")
        print(f"# Sampled Segs: {self.nsegs} # Restarts: {self.nrestarts}")

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def init_cluster_centers(self, k, include_truncal=True):
        """Randomly initialise k cluster DCF values in [0, 1]."""
        init = self.rng.uniform(low=0.0, high=1.0, size=k)
        if include_truncal:
            init[np.argmax(init)] = 0.98
        return init

    # ------------------------------------------------------------------
    # E-step: assign SNVs to (cluster, tree) pairs
    # ------------------------------------------------------------------

    def optimize_cluster_and_tree_assignments(self, snv_trees, dcfs, ell):
        """Assign each SNV in segment ell to the (cluster, tree) pair with
        the highest binomial log-likelihood.

        Returns
        -------
        best_cluster : dict  snv → cluster index
        best_tree    : dict  snv → GenotypeTree
        likelihood   : dict  snv → log-likelihood of best assignment
        clust_group  : dict  cluster → [(tree, [snv, ...])]
        """
        snvs = np.asarray(self.data.seg_to_snvs[ell], dtype=int)
        cn_prob = self.cn_prob_cache[ell]

        # Trees with a loss edge sort first (has_loss()=True → reverse=True puts
        # True before False).  Copy-neutral trees therefore come last and are
        # selected only when no loss tree scores higher.
        snv_trees = sorted(snv_trees, reverse=True, key=lambda T: T.has_loss())

        a_vec = self.alt_marginal    # (M_total,) pre-marginalised
        d_vec = self.total_marginal  # (M_total,)

        tree_parameters = [
            self.tree_parameter_cache[id(tree)] for tree in snv_trees
        ]

        opt_clust = np.empty(len(snvs), dtype=int)
        opt_tree = np.empty(len(snvs), dtype=int)
        opt_like = np.empty(len(snvs), dtype=float)

        # Process SNVs in chunks to cap peak memory (K × T × chunk).
        for start in range(0, len(snvs), self.assignment_chunk_size):
            stop = min(start + self.assignment_chunk_size, len(snvs))
            chunk = snvs[start:stop]

            scores = vectorized_tree_cluster_likelihood(
                dcfs, tree_parameters, a_vec[chunk], d_vec[chunk]
            )  # (K, T, chunk_size)

            flat = scores.reshape(len(dcfs) * len(snv_trees), len(chunk))
            best_flat = np.argmax(flat, axis=0)

            opt_clust[start:stop] = best_flat // len(snv_trees)
            opt_tree[start:stop] = best_flat % len(snv_trees)
            opt_like[start:stop] = np.max(flat, axis=0)

        # Build output dictionaries
        best_cluster = {}
        best_tree = {}
        likelihood = {}
        clust_to_snv_indices = defaultdict(list)

        for i, m in enumerate(snvs):
            q = int(opt_clust[i])
            best_cluster[m] = q
            best_tree[m] = snv_trees[int(opt_tree[i])]
            likelihood[m] = float(opt_like[i])
            clust_to_snv_indices[q].append(i)

        clust_group = defaultdict(list)
        for q, indices in clust_to_snv_indices.items():
            idx_arr = np.asarray(indices, dtype=int)
            for t in np.unique(opt_tree[idx_arr]):
                mask = opt_tree[idx_arr] == t
                clust_group[q].append((
                    snv_trees[int(t)],
                    snvs[idx_arr[mask]].tolist(),
                ))

        # Sanity check
        all_snvs = [s for pairs in clust_group.values() for _, sl in pairs for s in sl]
        assert len(all_snvs) == len(snvs), (
            f"SNV count mismatch after assignment: {len(all_snvs)} != {len(snvs)}"
        )

        return best_cluster, best_tree, likelihood, clust_group

    # ------------------------------------------------------------------
    # M-step: optimise cluster DCF centers
    # ------------------------------------------------------------------

    def _build_dcf_group_cache(self, dcfs, clust_group_map):
        """Pre-compute sufficient-statistics groups for every cluster.

        Stored in self.dcf_group_cache so that compute_cost() can use them
        without re-slicing the data.  Must be called before compute_cost().
        """
        self.dcf_group_cache = {
            q: dcf_groups_for_cluster(
                clust_group_map, q, self.data,
                self.cn_prob_cache, self.alt_marginal, self.total_marginal,
            )
            for q in range(len(dcfs))
        }

    def optimize_cluster_centers(self, dcfs, clust_group_map):
        """Optimise each cluster's DCF center given the current assignments."""
        self._build_dcf_group_cache(dcfs, clust_group_map)
        new_dcfs = np.array([
            maximize_dcf_groups(self.dcf_group_cache[q], dcfs[q])
            for q in range(len(dcfs))
        ])
        return new_dcfs

    # ------------------------------------------------------------------
    # Likelihood / cost evaluation
    # ------------------------------------------------------------------

    def compute_cost(self, dcfs):
        """Return the negative total log-likelihood (lower = better fit).

        Requires self.dcf_group_cache to be populated by a prior call to
        _build_dcf_group_cache (which optimize_cluster_centers always does).
        """
        return sum(
            -_exact_dcf_log_likelihood(self.dcf_group_cache.get(q, []), dcf)
            for q, dcf in enumerate(dcfs)
        )



    # ------------------------------------------------------------------
    # Main EM loop
    # ------------------------------------------------------------------

    def decifer(self, data, dcfs, segments):
        """Run one full EM restart.

        Parameters
        ----------
        data     : Data  (see data.py)
        dcfs     : (k,) array of initial cluster DCF values
        segments : list of segment indices to use

        Returns
        -------
        norm_likelihood, dcfs, CNA_tree, CLUSTER_ASSIGNMENTS,
        TREE_ASSIGNMENTS, raw_likelihood, n_snvs
        """
        init_dcfs = dcfs.copy()
        prev_cost = np.inf
        self.k = len(dcfs)
        self.data = data

        # Marginalise read counts over cells once per restart.
        self.alt_marginal = np.asarray(data.var.sum(axis=0), dtype=float)
        self.total_marginal = np.asarray(data.total.sum(axis=0), dtype=float)

        # Per-segment caches (built once; do not depend on DCF values).
        self.cn_prob_cache = {}
        self.tree_parameter_cache = {}

        m = sum(data.num_snvs(ell) for ell in segments)

        # ---- pre-build candidate trees and parameter cache ----
        cn_props = {}
        candidate_trees = {}   # ell → [(cna_tree_s, [GenotypeTree, ...])]

        for ell in segments:
            cn_props[ell] = data.thresholded_cn_prop(
                ell, thresh=self.thresh_prop, start_state=(1, 1)
            )
            self.cn_prob_cache[ell] = dict(data.cn_proportions(ell))
            if (1, 1) not in self.cn_prob_cache[ell]:
                self.cn_prob_cache[ell][(1, 1)] = 0

            states = set(cn_props[ell].keys())
            cna_trees = clonelib.get_cna_trees(states, 1, 1)

            candidate_trees[ell] = []
            for s in cna_trees:
                if len(s) == 0:
                    seg_states = list(cn_props[ell].keys())
                    assert len(seg_states) == 1
                    cn = seg_states[0]
                    scriptT = [[((*cn, 0, 0), (*cn, 1, 0))]]
                else:
                    scriptT = clonelib.get_genotype_trees(s)

                T_SNVs = [
                    GenotypeTree(edges=edge_list, id=i)
                    for i, edge_list in enumerate(scriptT)
                ]
                for tree in T_SNVs:
                    self.tree_parameter_cache[id(tree)] = tree_dcf_parameters(
                        tree, self.cn_prob_cache[ell]
                    )
                candidate_trees[ell].append((s, T_SNVs))

        # ---- EM iterations ----
        CLUSTER_ASSIGNMENTS = {}
        TREE_ASSIGNMENTS = {}
        CNA_tree = {}
        CLUST_GROUP_MAPPING = {}

        for j in range(self.iterations):
            CLUSTER_ASSIGNMENTS = {}
            TREE_ASSIGNMENTS = {}
            CNA_tree = {}
            CLUST_GROUP_MAPPING = {}

            for ell in segments:
                seg_like = -np.inf

                if self.cna_restriction:
                    # Pick the CNA tree that yields the best SNV assignment.
                    for s, T_SNVs in candidate_trees[ell]:
                        cluster, tree, likelihood, clust_group = \
                            self.optimize_cluster_and_tree_assignments(
                                T_SNVs, dcfs, ell,
                            )
                        new_seg_likelihood = sum(likelihood.values())

                        if new_seg_likelihood > seg_like:
                            CNA_tree[ell] = s
                            CLUSTER_ASSIGNMENTS[ell] = cluster
                            TREE_ASSIGNMENTS[ell] = tree
                            seg_like = new_seg_likelihood
                            CLUST_GROUP_MAPPING[ell] = clust_group

                else:
                    # Pool all genotype trees across all CNA trees.
                    all_T_SNVs = []
                    for s, T_SNVs in candidate_trees[ell]:
                        all_T_SNVs.extend(T_SNVs)
                    # Trees are already registered in tree_parameter_cache.

                    cluster, tree, likelihood, clust_group = \
                        self.optimize_cluster_and_tree_assignments(
                            all_T_SNVs, dcfs, ell,
                        )
                    new_seg_likelihood = sum(likelihood.values())

                    if new_seg_likelihood > seg_like:
                        CNA_tree[ell] = None   # no single CNA tree
                        CLUSTER_ASSIGNMENTS[ell] = cluster
                        TREE_ASSIGNMENTS[ell] = tree
                        CLUST_GROUP_MAPPING[ell] = clust_group
                        seg_like = new_seg_likelihood

            # M-step
            dcfs = self.optimize_cluster_centers(dcfs, CLUST_GROUP_MAPPING)
            dcfs[dcfs > 0.999] = 1.0
            dcfs[dcfs < 1e-3] = 0.0

            # Convergence check (normalised by number of SNVs)
            new_cost = self.compute_cost(dcfs)
            diff = abs(new_cost - prev_cost) / max(m, 1)

            if self.verbose:
                print(
                    f"Iter {j}: cost={new_cost:.4f}  "
                    f"norm_diff={diff:.2e}  dcfs={np.round(dcfs, 3)}"
                )

            if diff < self.tolerance:
                break
            prev_cost = new_cost

        # ---- post-processing ----
        clust_assign = defaultdict(list)
        for ell in CLUSTER_ASSIGNMENTS:
            for snv, q in CLUSTER_ASSIGNMENTS[ell].items():
                clust_assign[q].append(snv)

        if self.verbose:
            print(f"Restart complete. norm_cost={new_cost/max(m,1):.6f}")
            print(f"Initial DCFs : {init_dcfs}")
            print(f"Inferred DCFs: {dcfs}")

        # Preserve original API: result[0] is normalised negative log-likelihood
        # (cost / m), lower = better fit.  External callers rely on this sign.
        norm_cost = new_cost / max(m, 1)

        return (norm_cost, dcfs, CNA_tree,
                CLUSTER_ASSIGNMENTS, TREE_ASSIGNMENTS,
                new_cost, m)

    # ------------------------------------------------------------------
    # Model selection
    # ------------------------------------------------------------------

    @staticmethod
    def elbow_criteria_model_selection(scalar_objs, mink, maxk, ubleft=0.06):
        """Select k using an elbow criterion on normalised likelihoods.

        Parameters
        ----------
        scalar_objs : dict  k → float  (norm likelihood; higher = better)
        """
        if mink == maxk:
            return mink, scalar_objs, {}

        # Enforce monotonicity: if likelihood drops going from k-1 to k,
        # cap it at the k-1 value.
        for k in range(mink + 1, maxk + 1):
            if scalar_objs[k] < scalar_objs[k - 1]:
                scalar_objs[k] = scalar_objs[k - 1]

        chk = lambda v: v if v != 0.0 else 0.01
        left = (lambda k: min(
            (scalar_objs[k] - scalar_objs[k - 1]) / abs(chk(scalar_objs[k - 1])),
            ubleft,
        ) if k > mink else ubleft)
        right = (lambda k: (scalar_objs[k + 1] - scalar_objs[k]) / abs(chk(scalar_objs[k])))

        elbow = {k: left(k) - right(k) for k in range(mink, maxk)}
        selected = max(range(mink, maxk), key=lambda k: elbow[k])
        return selected, scalar_objs, elbow

    # ------------------------------------------------------------------
    # Outer loop: restarts + model selection
    # ------------------------------------------------------------------

    @timeit_decorator
    def run(self, data, k_vals, cores=1, nfull=1, segs=None):
        print(f"Number initial segs: {self.nsegs}, number restarts: {self.nrestarts}")
        print(segs)

        results_by_k = {}
        scalar_likelihoods = {}   # k → float  used for elbow

        for k in k_vals:
            all_results = []

            if cores <= 1:
                for i in range(self.nrestarts):
                    init_dcfs = self.init_cluster_centers(k)
                    if segs is None:
                        init_segs = data.sample_segments(
                            self.nsegs, self.rng,
                            max_cn_states=self.max_cn_states,
                            min_snvs=1, thresh=self.thresh_prop,
                        )
                    else:
                        n = min(self.nsegs, len(segs))
                        init_segs = self.rng.choice(segs, n, replace=False)
                    res = self.decifer(data, init_dcfs, init_segs)
                    all_results.append(res)
                    if self.verbose:
                        print(f"Restart {i}: norm_likelihood={res[0]:.6f}  DCFs={res[1]}")
            else:
                if segs is None:
                    sample_segs = [
                        data.sample_segments(
                            self.nsegs, self.rng,
                            max_cn_states=self.max_cn_states,
                            min_snvs=1, thresh=self.thresh_prop,
                        )
                        for _ in range(self.nrestarts)
                    ]
                else:
                    n = min(self.nsegs, len(segs))
                    sample_segs = [
                        self.rng.choice(segs, n, replace=False)
                        for _ in range(self.nrestarts)
                    ]
                args = [
                    (data, self.init_cluster_centers(k), sample_segs[i])
                    for i in range(self.nrestarts)
                ]
                with multiprocessing.Pool(cores) as pool:
                    all_results = pool.starmap(self.decifer, args)

            # Sort ascending: result[0] is normalised negative log-likelihood,
            # lower = better.
            sorted_results = sorted(all_results, key=lambda x: x[0])
            top_results = sorted_results[:nfull]

            for res in top_results:
                print(f"  Sampled-seg norm cost: {res[0]:.6f}  DCFs: {res[1]}")

            # Re-run top solutions on the full segment set
            if segs is None:
                all_segments = [
                    ell for ell in data.seg_to_snvs
                    if data.num_snvs(ell) >= 1
                    and data.num_cn_states(ell, self.thresh_prop) <= self.max_cn_states
                ]
            else:
                all_segments = list(segs)

            full_results = [
                self.decifer(data, res[1], all_segments) for res in top_results
            ]
            full_results = sorted(full_results, key=lambda x: x[0])

            for res in full_results:
                print(f"  Full-seg norm cost:    {res[0]:.6f}  DCFs: {res[1]}")

            best = full_results[0]
            results_by_k[k] = best
            # Negate for elbow: elbow expects higher = better likelihood
            scalar_likelihoods[k] = -best[0]

        print("Starting model selection using elbow criteria...")
        selected, updated_obj, elbow = self.elbow_criteria_model_selection(
            scalar_likelihoods, min(k_vals), max(k_vals)
        )

        return results_by_k[selected]


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main(args):
    data = pd.read_pickle(args.pickle_path)
    print(data)

    dec = DCF_Clustering(
        nrestarts=args.num_restarts,
        seed=args.seed,
        cna_restriction=args.restrict_CNA_trees,
        nsegs=args.nsegs,
        verbose=args.verbose,
        thresh_prop=args.thresh_prop,
        max_cn_states=args.max_cn_states,
    )

    if args.segments is not None:
        # Caller supplied an explicit list — use it as-is (they own the filter).
        with open(args.segments) as f:
            segs = [int(line.strip()) for line in f]
    else:
        # Mirror the quality filter that run() applies when segs=None, so the
        # CLI and the programmatic path see the same segment universe.
        segs = [
            ell for ell in data.seg_to_snvs
            if data.num_snvs(ell) >= 1
            and data.num_cn_states(ell, dec.thresh_prop) <= dec.max_cn_states
        ]

    if args.excl_segments is not None:
        excl = set(args.excl_segments)
        segs = [s for s in segs if s not in excl]

    if args.ground_truth is not None:
        ground_truth = read_ground_truth_text_file(args.ground_truth)
        res = dec.decifer(data, np.array(ground_truth), data.segments)
        print(f"GT norm likelihood: {res[0]}")
        print(f"DCFs: {res[1]}")
        if args.output_path is not None:
            with open(args.output_path, 'wb') as f:
                pickle.dump(res, f)
        if args.dcfs is not None:
            with open(args.dcfs, 'w') as f:
                for d in res[1]:
                    f.write(f"{d}\n")
        return
    elif args.clusters is not None:
        k_vals = [args.clusters]
    else:
        mink = min(args.mink, args.maxk)
        maxk = max(args.mink, args.maxk)
        # One extra k so maxk can itself be an elbow point
        k_vals = list(range(mink, maxk + 2))

    best = dec.run(data, k_vals=k_vals, cores=args.cores,
                   nfull=args.nfull, segs=segs)

    print("Best DCFs:")
    print(best[1])
    print(f"Norm likelihood: {best[0]}")

    dcfs = best[1]

    if args.output_path is not None:
        with open(args.output_path, 'wb') as f:
            pickle.dump(best, f)

    if args.dcfs is not None:
        with open(args.dcfs, 'w') as f:
            for d in dcfs:
                f.write(f"{d}\n")


def read_ground_truth_text_file(file_path):
    with open(file_path) as f:
        return [float(line.strip()) for line in f]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform DCF clustering on data from a pickled object"
    )
    parser.add_argument("-d", "--pickle_path", type=str)
    parser.add_argument("-g", "--ground_truth", type=str)
    parser.add_argument("-o", "--output_path", type=str)
    parser.add_argument("-r", "--num-restarts", type=int, default=25)
    parser.add_argument("--nsegs", type=int, default=15)
    parser.add_argument("--max-cn-states", type=int, default=4)
    parser.add_argument("--nfull", type=int, default=1)
    parser.add_argument(
        "--no-cna-restriction",
        dest="restrict_CNA_trees",
        action="store_false",
        help="allow genotype trees from different CNA histories in the same segment "
             "(biologically incoherent; off by default)",
    )
    parser.set_defaults(restrict_CNA_trees=True)
    parser.add_argument("-k", "--clusters", type=int)
    parser.add_argument("--thresh-prop", type=float, default=0.05)
    parser.add_argument("--mink", type=int, default=2)
    parser.add_argument("--maxk", type=int, default=5)
    parser.add_argument("-D", "--dcfs", type=str)
    parser.add_argument("-s", "--seed", type=int, default=21)
    parser.add_argument("-j", "--cores", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--segments", type=str)
    parser.add_argument("--excl-segments", type=int, nargs='+')

    args = parser.parse_args()
    main(args)