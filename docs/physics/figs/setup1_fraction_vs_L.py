import numpy as np
import matplotlib.pyplot as plt

# Setup 1 free energy (nondim):  F(r) = -1/2 (L-1) r^2 - 1/4 lam r^4 + 1/6 lam r^6
# Stationary density  rho(r) ~ r^(d-1) exp(-F/T),  T = theta * L
# f = fraction beyond the barrier radius r_- = sqrt((1-mu)/2), mu = sqrt(1+4(L-1)/lam)

def frac(L, lam, theta, d, n=200000, rmax_pad=1.6):
    a = L - 1.0
    disc = 1.0 + 4.0 * a / lam
    if disc <= 0.0:
        return 0.0                      # no non-trivial extrema
    mu = np.sqrt(disc)
    s_minus = 0.5 * (1.0 - mu)
    if a >= 0.0:
        return 1.0                      # P=0 unstable: single polarized state
    r_b = np.sqrt(max(s_minus, 0.0))
    T = theta * L

    r_plus = np.sqrt(0.5 * (1.0 + mu))
    rmax = rmax_pad * max(1.0, r_plus) + (6.0 * 40.0 * T / lam) ** (1.0 / 6.0)
    r = np.linspace(0.0, rmax, n)
    F = -0.5 * a * r**2 - 0.25 * lam * r**4 + lam * r**6 / 6.0
    e = -F / T
    w = r ** (d - 1) * np.exp(e - e.max())

    tot = np.trapezoid(w, r)
    top = np.trapezoid(np.where(r >= r_b, w, 0.0), r)
    return top / tot


fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)

# ---- panel A: vary lambda at fixed noise ----
ax = axes[0]
theta = 5e-3
ax.set_xlim(0, 1.5)
for lam, c in zip([0.5, 1.0, 2.0, 3.0], ['#1b4965', '#3d8ea8', '#c96a3c', '#8c2f39']):
    Ls = np.linspace(max(1e-6, 1 - lam / 4 - 0.02 * lam), 1.0 - 1e-9, 600)
    f = np.array([frac(L, lam, theta, 2) for L in Ls])
    ax.plot(-np.log(Ls), f, color=c, lw=2, label=rf'$\lambda={lam}$')
    Lmax = 1 - 3 * lam / 16          # equal-depth (Maxwell) point
    ax.plot(-np.log(Lmax), frac(Lmax, lam, theta, 2), 'o', color=c, ms=5,
            mec='white', mew=1.0, zorder=5)
ax.set_title(rf'$d=2$,  $\vartheta={theta}$   (dots: $\mathcal{{L}}=1-3\lambda/16$)', fontsize=10)
ax.legend(frameon=False, fontsize=9)

# ---- panel B: vary noise at fixed lambda ----
ax = axes[1]
lam = 2.0
for th, c in zip([2e-4, 1e-3, 5e-3, 2e-2], ['#1b4965', '#3d8ea8', '#c96a3c', '#8c2f39']):
    Ls = np.linspace(max(1e-6, 1 - lam / 4 - 0.02 * lam), 1.0 - 1e-9, 600)
    f2 = np.array([frac(L, lam, th, 2) for L in Ls])
    f1 = np.array([frac(L, lam, th, 1) for L in Ls])
    ax.plot(-np.log(Ls), f2, color=c, lw=2, label=rf'$\vartheta={th}$')
    ax.plot(-np.log(Ls), f1, color=c, lw=1.1, ls='--')
ax.axvline(-np.log(1 - 3 * lam / 16), color='0.6', lw=0.9, ls=':')
ax.set_title(rf'$\lambda={lam}$;  solid $d=2$, dashed $d=1$', fontsize=10)
ax.legend(frameon=False, fontsize=9)

for ax in axes:
    ax.set_xlabel(r'$-\log\mathcal{L}$')
    ax.set_ylabel(r'$f$  (fraction in polarized minimum)')
    ax.set_ylim(-0.02, 1.02)
    ax.axhline(0.5, color='0.85', lw=0.8, zorder=0)
    ax.spines[['top', 'right']].set_visible(False)

out = '/private/tmp/claude-502/-Users-romang-Desktop-codes-GL-motility/50a5fd07-7dbb-4bb1-9dd5-95f2d7ca20e6/scratchpad/frac_vs_L.png'
fig.savefig(out, dpi=160)
print('saved', out)

# a few numbers
print('\nlam=2, theta=1e-3, d=2:')
for L in [0.51, 0.55, 0.625, 0.7, 0.8, 0.9, 0.99]:
    print(f'  L={L:.3f}  -logL={-np.log(L):.4f}  f={frac(L,2.0,1e-3,2):.4f}')
