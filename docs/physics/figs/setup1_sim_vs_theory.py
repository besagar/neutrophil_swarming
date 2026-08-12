"""Setup 1: exact Boltzmann fraction vs Laplace closed form vs Euler-Maruyama simulation.

SDE (nondim, d=2):  dP_a = [(L-1)P_a + lam(|P|^2 P_a - |P|^4 P_a)] dt~ + sqrt(2 vt L) dW_a
Integrator: Euler-Maruyama, dt~ = 0.005 (explicit; see CLAUDE.md rule 5).

Each (L, vartheta, init) cell gets N independent walkers; all cells are packed into
one flat array so the whole sweep is a single vectorized time loop.

Two initial conditions per cell -- all walkers at P=0, and all at |P|=r_+ -- so that
a gap between the two markers diagnoses NON-equilibration (Kramers time > run time)
rather than a failure of the theory.
"""
import numpy as np
import matplotlib.pyplot as plt

OUT = __file__.rsplit('/', 1)[0]
LAM = 2.0
DT = 0.005
NSTEPS = 80000          # t~ = 400
BURN = 0.4              # discard this fraction before accumulating
SAMPLE_EVERY = 25
N = 400                 # walkers per (L, vartheta, init) cell
THETAS = [0.005, 0.02, 0.05]
SEED = 20260721


def mu_of(L, lam=LAM):
    disc = 1.0 + 4.0 * (L - 1.0) / lam
    return np.sqrt(disc) if disc > 0 else None


def frac_exact(L, theta, d=2, lam=LAM, n=200000):
    mu = mu_of(L, lam)
    if mu is None:
        return 0.0
    if L >= 1.0:
        return 1.0
    r_b = np.sqrt(max(0.5 * (1.0 - mu), 0.0))
    r_plus = np.sqrt(0.5 * (1.0 + mu))
    T = theta * L
    rmax = 1.6 * max(1.0, r_plus) + (240.0 * T / lam) ** (1.0 / 6.0)
    r = np.linspace(0.0, rmax, n)
    F = -0.5 * (L - 1.0) * r**2 - 0.25 * lam * r**4 + lam * r**6 / 6.0
    w = r ** (d - 1) * np.exp(-F / T - (-F / T).max())
    return np.trapezoid(np.where(r >= r_b, w, 0.0), r) / np.trapezoid(w, r)


def frac_laplace(L, theta, d=2, lam=LAM):
    mu = mu_of(L, lam)
    if mu is None:
        return 0.0
    if L >= 1.0:
        return 1.0
    T = theta * L
    a = 1.0 - L
    Fpp = lam * mu * (1.0 + mu)
    dF = lam * (1.0 + mu) ** 2 * (1.0 - 2.0 * mu) / 48.0
    r_plus = np.sqrt(0.5 * (1.0 + mu))
    if d == 1:
        ratio = 0.5 * np.sqrt(Fpp / a) * np.exp(dF / T)
    else:
        ratio = np.sqrt(T * Fpp / (2 * np.pi)) / (a * r_plus) * np.exp(dF / T)
    return 1.0 / (1.0 + ratio)


# ---------------- build the packed parameter arrays ----------------
Lgrid = np.linspace(1 - LAM / 4 + 0.01, 0.99, 16)
cells = []                                   # (L, theta, init_tag)
for th in THETAS:
    for L in Lgrid:
        for tag in (0, 1):
            cells.append((L, th, tag))

L_c = np.repeat([c[0] for c in cells], N)
th_c = np.repeat([c[1] for c in cells], N)
tag_c = np.repeat([c[2] for c in cells], N)
a_c = L_c - 1.0
sig_c = np.sqrt(2.0 * th_c * L_c * DT)
mu_c = np.sqrt(np.maximum(1.0 + 4.0 * a_c / LAM, 0.0))
rb_c = np.sqrt(np.maximum(0.5 * (1.0 - mu_c), 0.0))      # barrier radius
rp_c = np.sqrt(np.maximum(0.5 * (1.0 + mu_c), 0.0))      # minimum radius

M = L_c.size
rng = np.random.default_rng(SEED)
px = np.where(tag_c == 0, 0.0, rp_c).astype(np.float64)
py = np.zeros(M)
px += rng.normal(0.0, 1e-3, M)               # break the exact-zero fixed point
py += rng.normal(0.0, 1e-3, M)

acc = np.zeros(M)
nacc = 0
noise = np.empty(M)
burn_steps = int(BURN * NSTEPS)

print(f'{M} walkers, {NSTEPS} steps, dt={DT} (t~={NSTEPS*DT:g})')
for step in range(NSTEPS):
    s = px * px + py * py
    g = a_c + LAM * (s - s * s)              # dP/dt = g * P
    rng.standard_normal(M, out=noise)
    px += g * px * DT + sig_c * noise
    rng.standard_normal(M, out=noise)
    py += g * py * DT + sig_c * noise
    if step >= burn_steps and step % SAMPLE_EVERY == 0:
        acc += (px * px + py * py) >= rb_c * rb_c
        nacc += 1
    if step % 20000 == 0:
        print(f'  step {step}')

f_sim_flat = (acc / nacc).reshape(len(cells), N).mean(axis=1)
f_sim = {(c[0], c[1], c[2]): v for c, v in zip(cells, f_sim_flat)}
# binomial-ish error bar from walker-to-walker spread
f_err_flat = (acc / nacc).reshape(len(cells), N).std(axis=1) / np.sqrt(N)
f_err = {(c[0], c[1], c[2]): v for c, v in zip(cells, f_err_flat)}

# ---------------- plot ----------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), constrained_layout=True, sharey=True)
Ldense = np.linspace(1 - LAM / 4 + 1e-4, 1.0 - 1e-4, 500)
xd = -np.log(Ldense)

for ax, th in zip(axes, THETAS):
    ax.plot(xd, [frac_exact(L, th) for L in Ldense], color='#1b4965', lw=2.4,
            label='exact (quadrature)', zorder=2)
    ax.plot(xd, [frac_laplace(L, th) for L in Ldense], color='#c96a3c', lw=1.4,
            ls='--', label='Laplace (closed form)', zorder=3)
    xs = -np.log(Lgrid)
    for tag, mk, lab in [(0, 'o', r'sim, init $P=0$'), (1, '^', r'sim, init $|P|=r_+$')]:
        ys = [f_sim[(L, th, tag)] for L in Lgrid]
        es = [f_err[(L, th, tag)] for L in Lgrid]
        ax.errorbar(xs, ys, yerr=es, fmt=mk, ms=5.5, mfc='none', lw=0,
                    elinewidth=1, color='#2b2b2b' if tag == 0 else '#8c2f39',
                    label=lab, zorder=4)
    ax.axvline(-np.log(1 - 3 * LAM / 16), color='0.8', lw=0.9, ls=':')
    ax.set_title(rf'$\vartheta={th}$', fontsize=11)
    ax.set_xlabel(r'$-\log\mathcal{L}$')
    ax.spines[['top', 'right']].set_visible(False)

axes[0].set_ylabel(r'$f$  (fraction with $|P|>r_-$)')
axes[0].set_ylim(-0.05, 1.08)
axes[0].legend(frameon=False, fontsize=8.5, loc='lower left')
fig.suptitle(rf'Setup 1, $d=2$, $\lambda={LAM}$ — Euler–Maruyama ($N={N}$/point, '
             rf'$\Delta\tilde t={DT}$, $\tilde t={NSTEPS*DT:g}$) vs theory', fontsize=11)
fig.savefig(f'{OUT}/setup1_sim_vs_theory.png', dpi=160)
print('saved', f'{OUT}/setup1_sim_vs_theory.png')

print(f'\n{"vt":>7} {"L":>7} {"exact":>8} {"laplace":>8} {"sim(0)":>8} {"sim(r+)":>8} {"gap":>7}')
for th in THETAS:
    for L in Lgrid[::3]:
        s0, s1 = f_sim[(L, th, 0)], f_sim[(L, th, 1)]
        print(f'{th:7.3f} {L:7.3f} {frac_exact(L, th):8.4f} {frac_laplace(L, th):8.4f} '
              f'{s0:8.4f} {s1:8.4f} {abs(s0-s1):7.4f}')
