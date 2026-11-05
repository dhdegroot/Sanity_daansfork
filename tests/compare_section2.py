"""Check Sanity's output against an independent NumPy implementation of Section 2 of the error-bar note.

"note" is E. van Nimwegen, "Getting proper error-bars in Sanity" (update_final_proper_error_bars.pdf); equation numbers
below are its Section 2. "SI" is the Supplementary Information of Breda, Zavolan & van Nimwegen, Nat. Biotechnol. 2021.

For every gene of a small count matrix this script solves q_g at every grid point by bisection (with its own Lambert W
iteration), computes L(v) from Eqs. 36-38 with the actual mean <x*>, the estimates at the selected grid point (Eqs. 32,
33, 39, 41, 43, 44, 46) and the MARG averages (Eqs. 47-50), and compares them with the -e 1 -npy output of the Sanity
binary for all four variance methods and both values of --measured-expansion (posterior-mode: Eqs. 44 and 46 for every
cell; likelihood-maximum: log(n_c/N_c) +- 1/sqrt(n_c) for nonzero counts, Eqs. 44 and 46 for zero counts), checking
that the two expansions differ only in the nonzero entries of the measured files. It also checks the identities the
note relies on:
sum_c f*_c = 1, <x*> = log(n) - q, Eq. 41 against the diagonal of the inverted precision matrix of Eq. 35, Eq. 40 and
Eq. 43 against the full inverse, and that Eqs. 36-38 equal the SI's L(v) up to a constant over the grid.

Usage: python3 compare_section2.py [--sanity ../bin/Sanity] [--counts matrix.tsv] [--keep]
Without --counts a synthetic matrix of 24 genes x 60 cells is generated (deterministically). The script needs NumPy.
"""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np


# ----------------------------------------------------------------------------------------------------------------------
# Section 2 in NumPy
# ----------------------------------------------------------------------------------------------------------------------

def lambert_w_of_exp(s):
    """W(exp(s)) for the principal branch, by Newton's method on w + log(w) = s."""
    s = np.asarray(s, dtype=float)
    w = np.where(s > 1.0, s - np.log(np.maximum(s, 1.0)), np.exp(np.minimum(s, 1.0)))
    for _ in range(100):
        residual = w + np.log(w) - s
        step = residual * w / (w + 1.0)
        w = w - step
        if np.all(np.abs(step) <= 1e-15 * np.maximum(np.abs(w), 1e-300)):
            break
    return w


def fractions_at(q, counts, n, v, log_N):
    """f*_c of Eq. 31 for a given q, with -q = <log f> - <log N> enforced by the caller through sum_c f*_c = 1."""
    return lambert_w_of_exp(np.log(n * v) + counts * v + log_N - q) / (n * v)


def solve_q(counts, n, v, log_N, sum_N):
    """q at one v by bisection on R(q) = 1 - sum_c f*_c(q), which increases with q."""
    residual = lambda q: 1.0 - fractions_at(q, counts, n, v, log_N).sum()
    lo = hi = np.log(sum_N) + 0.5 * v
    step = max(0.5 * v, 1.0)
    while residual(lo) > 0:
        lo -= step
        step *= 2
    step = max(0.5 * v, 1.0)
    while residual(hi) < 0:
        hi += step
        step *= 2
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) < 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-15 * max(1.0, abs(mid)):
            break
    return 0.5 * (lo + hi)


def log_likelihood(counts, n, v, log_N, f):
    """L(v) of Eqs. 36-38 (note), with x*_c of Eq. 32 and the actual mean <x*>."""
    C = counts.size
    x = np.log(f) + np.log(n) - log_N                               # Eq. 32
    deviation = x - x.mean()
    log_posterior = (-0.5 * (C - 1) * np.log(v) + (counts * x).sum() - (n * f).sum()
                     - (deviation ** 2).sum() / (2 * v))             # Eq. 36
    nfv = n * f * v
    log_det = -C * np.log(v) + np.log1p(nfv).sum() + np.log((nfv / (nfv + 1)).mean())  # Eq. 38
    return log_posterior - 0.5 * log_det                            # Eq. 37


def log_likelihood_si(counts, n, v, log_N, f, q):
    """The previous L(v) (SI Eqs. 20, 32, 33 without the C log(n) constant), for the equivalence check."""
    C = counts.size
    delta = np.log(f) - log_N + q
    L = -0.5 * C * np.log(v) + (counts * delta - delta ** 2 / (2 * v)).sum() - n * q
    a = f + 1 / (n * v)
    log_det = np.log(a).sum() + np.log(1 - (f ** 2 / a).sum())
    return L - 0.5 * log_det


def estimates_at(counts, n, v, log_N, f):
    """Everything Section 2 reports at one v: x*, epsilon^2, x^m, its error, mu and sigma_mu^2."""
    C = counts.size
    x = np.log(f) + np.log(n) - log_N                 # Eq. 32
    mu = x.mean()                                     # Eq. 33
    v_c = v / (n * f * v + 1)                         # Eq. 39
    mean_v = v_c.mean()
    eps_sq = v_c + v_c ** 2 / (C * (v - mean_v))      # Eq. 41
    sigma_mu_sq = (mean_v / C) / (1 - mean_v / v)     # Eq. 43
    x_m = x + (x - mu) / (n * f * v)                  # Eq. 44
    x_m_err_sq = 1 / (n * f)                          # Eq. 46
    # --measured-expansion likelihood-maximum: for a nonzero count the Poisson log-likelihood n_c x - N_c e^x is
    # expanded around its own maximum log(n_c/N_c), where the curvature is -n_c; zero counts keep Eqs. 44 and 46.
    nonzero = counts > 0
    x_m_lm = np.where(nonzero, np.log(np.where(nonzero, counts, 1.0)) - log_N, x_m)
    x_m_lm_err_sq = np.where(nonzero, 1 / np.where(nonzero, counts, 1.0), x_m_err_sq)
    return dict(x=x, eps_sq=eps_sq, x_m=x_m, x_m_err_sq=x_m_err_sq, x_m_lm=x_m_lm, x_m_lm_err_sq=x_m_lm_err_sq,
                mu=mu, sigma_mu_sq=sigma_mu_sq, v_c=v_c)


def precision_matrix(n, v, f):
    """M of Eq. 35."""
    C = f.size
    return np.diag(n * f + 1 / v) - 1 / (C * v)


def fit_gene(counts, log_N, sum_N, grid, method):
    """The fit of one gene as Sanity does it, with the reference formulas; returns a dict of the reported quantities."""
    C = counts.size
    n = counts.sum()
    K = grid.size
    q_k = np.empty(K)
    L_k = np.empty(K)
    for k, v in enumerate(grid):
        q_k[k] = solve_q(counts, n, v, log_N, sum_N)
        f = fractions_at(q_k[k], counts, n, v, log_N)
        L_k[k] = log_likelihood(counts, n, v, log_N, f)
    P_k = np.exp(L_k - L_k.max())
    P_k /= P_k.sum()
    if method == "MARG":
        # Eqs. 47-50, with Eq. 41 and Eq. 43 at every grid point.
        x = np.zeros(C); eps_sq = np.zeros(C); x_m = np.zeros(C); x_m_err_sq = np.zeros(C)
        mu_k = np.empty(K); sig_k = np.empty(K)
        per_bin = []
        for k, v in enumerate(grid):
            f = fractions_at(q_k[k], counts, n, v, log_N)
            est = estimates_at(counts, n, v, log_N, f)
            per_bin.append(est)
            mu_k[k], sig_k[k] = est["mu"], est["sigma_mu_sq"]
        x = sum(P * e["x"] for P, e in zip(P_k, per_bin))
        eps_sq = sum(P * (e["eps_sq"] + (e["x"] - x) ** 2) for P, e in zip(P_k, per_bin))
        x_m = sum(P * e["x_m"] for P, e in zip(P_k, per_bin))
        x_m_err_sq = sum(P * (e["x_m_err_sq"] + (e["x_m"] - x_m) ** 2) for P, e in zip(P_k, per_bin))
        x_m_lm = sum(P * e["x_m_lm"] for P, e in zip(P_k, per_bin))
        x_m_lm_err_sq = sum(P * (e["x_m_lm_err_sq"] + (e["x_m_lm"] - x_m_lm) ** 2) for P, e in zip(P_k, per_bin))
        mu = (P_k * mu_k).sum()
        sigma_mu_sq = (P_k * (sig_k + (mu_k - mu) ** 2)).sum()
        v_report = (P_k * grid).sum()
        return dict(x=x, eps_sq=eps_sq, x_m=x_m, x_m_err_sq=x_m_err_sq, x_m_lm=x_m_lm, x_m_lm_err_sq=x_m_lm_err_sq,
                    mu=mu, sigma_mu_sq=sigma_mu_sq, v=v_report, P=P_k, q_k=q_k, L_k=L_k)
    if method == "MLE":
        k = int(np.argmax(L_k))
    elif method == "MAP":
        k = int(np.argmax(P_k / grid))
    elif method == "EAP":
        k = int(np.argmin(np.abs(grid - (P_k * grid).sum())))
    else:
        raise ValueError(method)
    v = grid[k]
    f = fractions_at(q_k[k], counts, n, v, log_N)
    est = estimates_at(counts, n, v, log_N, f)
    est.update(v=v, P=P_k, q_k=q_k, L_k=L_k, k=k, f=f)
    return est


# ----------------------------------------------------------------------------------------------------------------------
# Identity checks at one grid point
# ----------------------------------------------------------------------------------------------------------------------

def check_identities(counts, log_N, sum_N, grid, report):
    """The identities of Section 2 at the MAP grid point of one gene, and the L(v) equivalence over the grid."""
    n = counts.sum()
    C = counts.size
    est = fit_gene(counts, log_N, sum_N, grid, "MAP")
    v, f, q = est["v"], est["f"], est["q_k"][est["k"]]
    report("sum_c f*_c - 1", abs(f.sum() - 1), 1e-12)
    report("<x*> - (log n - q)", abs(est["mu"] - (np.log(n) - q)), 1e-12)
    # Eq. 30 at the optimum: x*_c - <x*> = v (n_c - n f*_c).
    report("Eq. 30 residual", np.abs((est["x"] - est["mu"]) - v * (counts - n * f)).max(), 1e-9)
    # Eqs. 40, 41, 43 against the inverse of M (Eq. 35); the tolerance allows for the conditioning of a C x C inverse.
    cov = np.linalg.inv(precision_matrix(n, v, f))
    v_c = est["v_c"]
    cov_40 = np.diag(v_c) + np.outer(v_c, v_c) / (C * (v - v_c.mean()))
    report("Eq. 40 vs inv(M), max abs", np.abs(cov - cov_40).max(), 1e-8 * np.abs(cov).max())
    report("Eq. 41 vs diag inv(M), max rel", (np.abs(np.diag(cov) - est["eps_sq"]) / np.diag(cov)).max(), 1e-8)
    report("Eq. 43 vs mean of inv(M), rel", abs(cov.sum() / C ** 2 - est["sigma_mu_sq"]) / est["sigma_mu_sq"], 1e-8)
    # Eqs. 36-38 minus the SI's L(v) is constant over the grid (docs/likelihood_equivalence.md).
    diff = []
    for k, vk in enumerate(grid):
        fk = fractions_at(est["q_k"][k], counts, n, vk, log_N)
        diff.append(est["L_k"][k] - log_likelihood_si(counts, n, vk, log_N, fk, est["q_k"][k]))
    diff = np.array(diff)
    predicted = n * np.log(n) - n - 0.5 * (C + 1) * np.log(n) + 0.5 * np.log(C)
    report("L_new - L_old spread over grid", np.ptp(diff), 1e-8 * max(1.0, abs(predicted)))
    report("L_new - L_old vs predicted constant", abs(diff.mean() - predicted), 1e-8 * max(1.0, abs(predicted)))


# ----------------------------------------------------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------------------------------------------------

def synthetic_counts():
    """24 genes x 60 cells: expression levels from very low to high, a few special rows, fractional counts allowed."""
    rng = np.random.default_rng(2024)
    C = 60
    library = rng.lognormal(np.log(2000), 0.5, C)
    means = np.geomspace(0.02, 200, 20)
    counts = np.empty((24, C))
    for g, mean in enumerate(means):
        ltq = rng.normal(np.log(mean / library.mean()), rng.uniform(0.1, 1.5))
        counts[g] = rng.poisson(library * np.exp(ltq + rng.normal(0, 0.6, C)))
    counts[20] = 0; counts[20, 7] = 1                      # a single count
    counts[21] = 1                                         # one count everywhere
    counts[22] = rng.poisson(50, C) * 0.25                 # fractional counts
    counts[23] = np.round(rng.poisson(5000, C) * rng.uniform(0.5, 1.5, C))  # high counts
    return counts


def read_tsv(path):
    rows = [line.rstrip("\n").split("\t") for line in open(path)]
    header = rows[0]
    data = rows[1:]
    names = [r[0] for r in data]
    counts = np.array([[float(x) for x in r[1:]] for r in data])
    cells = header[1:] if len(header) == counts.shape[1] + 1 else header
    return names, cells, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sanity", type=Path, default=Path(__file__).resolve().parents[1] / "bin/Sanity")
    parser.add_argument("--counts", type=Path, help="A TSV count matrix instead of the synthetic one")
    parser.add_argument("--keep", action="store_true", help="Keep the temporary directory")
    parser.add_argument("--nbin", type=int, default=40)
    args = parser.parse_args()

    if args.counts:
        names, cells, counts = read_tsv(args.counts)
    else:
        counts = synthetic_counts()
        names = [f"gene_{g}" for g in range(counts.shape[0])]
        cells = [f"cell_{c}" for c in range(counts.shape[1])]
    counts = counts[counts.sum(1) > 0]
    N = counts.sum(0)
    log_N = np.log(N)
    sum_N = N.sum()
    vmin, vmax, nbin = 0.001, 50.0, args.nbin
    grid = vmin * np.exp(np.log(vmax / vmin) / (nbin - 1) * np.arange(nbin))

    failures = []

    def report(label, value, tolerance, gene=None):
        ok = value <= tolerance
        print(f"  {'ok  ' if ok else 'FAIL'} {label:42s} {value:.3e} (tolerance {tolerance:.1e})")
        if not ok:
            failures.append((gene, label, value))

    temporary = tempfile.mkdtemp(prefix="sanity-section2-")
    directory = Path(temporary)
    matrix = directory / "counts.tsv"
    with matrix.open("w") as stream:
        stream.write("GeneID\t" + "\t".join(cells) + "\n")
        for name, row in zip(names, counts):
            stream.write(name + "\t" + "\t".join(f"{x:.17g}" for x in row) + "\n")

    print("Identity checks at the MAP grid point (a few genes)")
    for g in range(min(counts.shape[0], 6)):
        print(f" gene {names[g]} (n = {counts[g].sum():g}, zeros = {(counts[g] == 0).sum()} of {counts.shape[1]})")
        check_identities(counts[g], log_N, sum_N, grid, lambda l, v, t: report(l, v, t, names[g]))

    quantities = {  # Sanity file -> (reference key, kind)
        "log_transcription_quotients": ("x", "value"),
        "ltq_error_bars": ("eps_sq", "sqrt"),
        "measured_ltq": ("x_m", "value"),
        "measured_ltq_error_bars": ("x_m_err_sq", "sqrt"),
        "mu": ("mu", "value"),
        "d_mu": ("sigma_mu_sq", "sqrt"),
        "variance": ("v", "value"),
        "likelihood": ("P", "value"),
        "delta": ("delta", "value"),
        "d_delta": ("eps_sq", "sqrt"),
    }
    # --measured-expansion likelihood-maximum only changes the two measured files.
    quantities_lm = dict(quantities, measured_ltq=("x_m_lm", "value"), measured_ltq_error_bars=("x_m_lm_err_sq", "sqrt"))
    keys = sorted({key for key, _ in quantities.values()} | {key for key, _ in quantities_lm.values()})
    for method in ["MAP", "MLE", "EAP", "MARG"]:
        reference = {key: [] for key in keys}
        for g in range(counts.shape[0]):
            est = fit_gene(counts[g], log_N, sum_N, grid, method)
            est["delta"] = est["x"] - est["mu"]
            for key in keys:
                reference[key].append(est[key])
        outputs = {}
        for expansion, files in [("posterior-mode", quantities), ("likelihood-maximum", quantities_lm)]:
            out = directory / f"out_{method}_{expansion}"
            outputs[expansion] = out
            command = [str(args.sanity), "-f", str(matrix), "-d", str(out), "-v_m", method, "-e", "1", "-npy",
                       "-n", "2", "-nbin", str(nbin), "-vmin", str(vmin), "-vmax", str(vmax),
                       "--measured-expansion", expansion]
            process = subprocess.run(command, capture_output=True, text=True)
            if process.returncode:
                raise RuntimeError(process.stdout + process.stderr)
            print(f"Method {method}, --measured-expansion {expansion}: Sanity versus the reference ({counts.shape[0]} genes)")
            for file_name, (key, kind) in files.items():
                path = out / f"{file_name}.npy"
                if not path.exists():
                    print(f"  (no {file_name}.npy)")
                    continue
                sanity = np.load(path)
                expected = np.array(reference[key])
                if kind == "sqrt":
                    expected = np.sqrt(expected)
                expected = expected.reshape(sanity.shape)
                # x^m divides by n f* v, which amplifies the difference between the two solvers of q; measure it relative.
                scale = 1.0 + np.abs(expected)
                report(file_name + " max |diff|/(1+|ref|)", (np.abs(sanity - expected) / scale).max(), 1e-7, method)
        # The expansion must leave everything else untouched: every other file bit-identical, and the zero-count
        # entries of the measured files too.
        print(f"Method {method}: the two expansions")
        zero = counts == 0
        for file_name in quantities:
            a = np.load(outputs["posterior-mode"] / f"{file_name}.npy")
            b = np.load(outputs["likelihood-maximum"] / f"{file_name}.npy")
            if file_name.startswith("measured"):
                differ = float(np.any(a[zero] != b[zero]))
                report(file_name + " zero-count entries differ", differ, 0.0, method)
            else:
                report(file_name + " bit-identical (0 = yes)", float(not np.array_equal(a, b)), 0.0, method)

    if failures:
        print(f"FAILED: {len(failures)} checks:")
        for gene, label, value in failures:
            print(f"  {gene}: {label} = {value:.3e}")
    else:
        print("PASSED: Sanity matches the Section 2 reference for MAP, MLE, EAP and MARG with both measured expansions, "
              "and the identities hold.")
    if args.keep:
        print("Outputs kept in", directory)
    else:
        shutil.rmtree(directory)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
