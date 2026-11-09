import numpy as np
from numba import njit
from scipy.optimize import minimize

from hmm_functions import initB, MIN_COALESCENCE_TIME

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------
# EM for the demographic parametrization of the HMM (see DemoStructure in hmm_functions.py).
# Emissions are Poisson with rate mut_rate * window_size * (A[state][k] . theta), so the expected complete-data
# log-likelihood only depends on G[s] (expected number of windows in state s) and S[s][k] (expected count of
# observation class k in state s). Transitions and starting probabilities are not estimated: they stay as the
# demography JSON defines them.
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------

def flatten_observations(obs_list):
    """obs dicts from Load_observations (chrom -> (ploidy, windows, states, outgroups)) -> list of sequences."""
    seqs = []
    for obs in obs_list:
        for chrom in obs:
            for z in range(obs[chrom].shape[0]):
                seqs.append(obs[chrom][z])
    return seqs


def obs_indices(seq, max_obs):
    """Index of each (window, state) observation in the emission table built by initB (same encoding as obs_to_ind)."""
    powers = max_obs ** np.arange(seq.shape[2])
    return (seq * powers).sum(axis=2)


@njit
def forward_backward(E, start, trans):
    """Scaled forward-backward. E[t, i] = emission probability. Returns loglik and gamma (T x N)."""
    T, N = E.shape
    alpha = np.zeros((T, N))
    beta = np.ones((T, N))
    c = np.zeros(T)
    alpha[0] = start * E[0]
    c[0] = alpha[0].sum()
    if c[0] == 0:
        c[0] = 1e-300
    alpha[0] /= c[0]
    for t in range(1, T):
        for j in range(N):
            s = 0.0
            for i in range(N):
                s += alpha[t - 1, i] * trans[i, j]
            alpha[t, j] = s * E[t, j]
        c[t] = alpha[t].sum()
        if c[t] == 0:
            c[t] = 1e-300
        alpha[t] /= c[t]
    for t in range(T - 2, -1, -1):
        for i in range(N):
            s = 0.0
            for j in range(N):
                s += trans[i, j] * E[t + 1, j] * beta[t + 1, j]
            beta[t, i] = s / c[t + 1]
    gamma = alpha * beta
    for t in range(T):
        g = gamma[t].sum()
        if g > 0:
            gamma[t] /= g
    return np.log(c).sum(), gamma


def e_step(seqs, hmm, max_obs, mut_rate):
    B = initB(hmm.emissions, mut_rate, max_obs)
    n_states, n_out = seqs[0].shape[1], seqs[0].shape[2]
    loglik = 0.0
    G = np.zeros(n_states)
    S = np.zeros((n_states, n_out))
    states = np.arange(n_states)
    for seq in seqs:
        E = B[states[None, :], obs_indices(seq, max_obs)]
        ll, gamma = forward_backward(E, hmm.starting_probabilities, hmm.transitions)
        loglik += ll
        G += gamma.sum(axis=0)
        S += np.einsum('ts,tsk->sk', gamma, seq)
    return loglik, G, S


def m_step_emissions(structure, G, S, theta0, mut_rate, window_size=1000, min_time=MIN_COALESCENCE_TIME):
    """Maximise sum_{s,k} S[s,k] log(lambda_sk) - G[s] lambda_sk over theta (concave), lambda = mut_rate*L*A.theta."""
    c = mut_rate * window_size
    rows = [(s, k) for s in range(structure.A.shape[0]) for k in range(structure.A.shape[1]) if structure.A[s, k].any()]
    R = np.array([structure.A[s, k] for s, k in rows])
    Sr = np.array([S[s, k] for s, k in rows])
    Gr = np.array([G[s] for s, k in rows])
    scale = 1000.0  # optimise theta in thousands of generations

    def f(u):
        t = np.maximum(R @ (u * scale), 1e-12)
        return -(Sr * np.log(c * t) - Gr * c * t).sum()

    def grad(u):
        t = np.maximum(R @ (u * scale), 1e-12)
        return -((Sr / t - Gr * c) @ R) * scale

    cons = [{'type': 'ineq', 'fun': lambda u: R @ (u * scale) - min_time, 'jac': lambda u: R * scale}]
    res = minimize(f, theta0 / scale, jac=grad, method='SLSQP', constraints=cons,
                   bounds=[(min_time / scale, None)] * len(theta0), options={'maxiter': 500, 'ftol': 1e-12})
    return res.x * scale


def fit_em(structure, seqs, max_obs, max_iter=50, tol=0.1,
           mut_rate=1.25e-8, rec_rate=1.2e-8, window_size=1000, log=print):
    theta, rec_times, props = structure.theta.copy(), structure.rec_times.copy(), structure.props.copy()
    log('iteration\tloglikelihood\t' + '\t'.join(structure.param_keys))
    prev = -np.inf
    first = None
    decreased = False
    for it in range(max_iter + 1):
        hmm = structure.hmm(theta, rec_times, props, rec_rate, window_size)
        loglik, G, S = e_step(seqs, hmm, max_obs, mut_rate)
        if first is None:
            first = loglik
        log(f'{it}\t{loglik:.4f}\t' + '\t'.join(f'{t:.1f}' for t in theta))
        if loglik < prev - 1e-6 * abs(prev):
            log(f'# warning: log-likelihood decreased by {prev - loglik:.4f}')
            decreased = True
        if it == max_iter or loglik - prev < tol:
            break
        prev = loglik
        theta = m_step_emissions(structure, G, S, theta, mut_rate, window_size)
    info = dict(iterations=it, converged=bool(it < max_iter), loglikelihood=round(float(loglik), 2),
                loglikelihood_initial=round(float(first), 2), monotone=not decreased,
                sequences=len(seqs), windows=int(sum(len(s) for s in seqs)))
    return structure.hmm(theta, rec_times, props, rec_rate, window_size), info


def em_report(structure, hmm, info, mut_rate=1.25e-8, window_size=1000):
    """A readable summary of what -em estimated, written to <out>/em.json.

    Says per parameter what it is, which lineage pairs share it, its value in the demography JSON and its fitted
    value, and keeps the quantities that were NOT estimated in a separate block so they are not mistaken for
    estimates. The full HMM is kept under "hmm" for reuse.
    """
    names, outs, A = structure.state_names, structure.outgroup_names, structure.A
    initial = dict(zip(structure.param_keys, structure.theta))

    params = []
    for p, key in enumerate(structure.param_keys):
        pop, _, entry = key.partition('@')
        pairs = [f'{names[s]} vs {outs[structure.order[k]]}'
                 for s in range(A.shape[0]) for k in range(A.shape[1]) if A[s, k, p] == 1]
        fitted = float(hmm.theta[p])
        params.append(dict(parameter=key, shared_population=pop, entry_time=float(entry),
                           lineage_pairs=pairs,
                           initial=round(float(initial[key]), 1), fitted=round(fitted, 1),
                           ratio=round(fitted / initial[key], 3) if initial[key] else None))

    rates = []
    for s, state in enumerate(names):
        for k in range(len(outs)):
            terms = ' '.join(f'{"+" if c > 0 else "-"} {key}'
                             for c, key in zip(A[s, k], structure.param_keys) if c).lstrip('+ ')
            rates.append(dict(state=state, outgroup=outs[structure.order[k]], expression=terms,
                              coalescence_time=round(float(hmm.emissions[s][k]) / window_size, 1),
                              expected_variants_per_window=round(float(hmm.emissions[s][k]) * mut_rate, 4)))

    report = dict(
        fitted='emission coalescence times',
        units='generations',
        run={k: (round(info[k], 2) if k.startswith('loglikelihood') else info[k])
             for k in ('individuals', 'sequences', 'windows', 'iterations', 'converged', 'monotone',
                       'loglikelihood_initial', 'loglikelihood') if k in info},
        emissions='conditional' if structure.conditional else 'independent',
        coalescence_times=params,
        emission_rates=rates,
        not_fitted=dict(
            note='-em estimates the emission coalescence times only. The values below come from the demography '
                 'JSON and were not estimated.',
            admixture_times=dict(zip(names, [round(float(t), 1) for t in hmm.rec_times])),
            ancestry_proportions=dict(zip(names, [round(float(p), 6) for p in hmm.starting_probabilities]))),
        hmm={k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in vars(hmm).items()})
    return report
