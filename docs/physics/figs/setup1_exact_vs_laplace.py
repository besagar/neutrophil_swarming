"""Setup 1: exact Boltzmann fraction (quadrature) vs closed-form Laplace estimate.

Exact:    f = int_{r_-}^inf r^(d-1) e^{-F/T} dr  /  int_0^inf ...
Laplace:  f = [1 + Z_0/Z_+]^{-1} with Gaussian wells at r=0 and r=r_+.
Both use the same exact F, s_pm, dF = lam (1+mu)^2 (1-2mu)/48, F''_+ = lam mu (1+mu).
"""
import numpy as np
import matplotlib.pyplot as plt

OUT = __file__.rsplit('/', 1)[0]


def _mu(L, lam):
    disc = 1.0 + 4.0 * (L - 1.0) / lam
    return np.sqrt(disc) if disc > 0 else None


def frac_exact(L, lam, theta, d, n=200000):
    mu = _mu(L, lam)
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


def frac_laplace(L, lam, theta, d):
    mu = _mu(L, lam)
    if mu is None:
        return 0.0
    if L >= 1.0:
        return 1.0
    T = theta * L
    a = 1.0 - L                       # |F''(0)|
    Fpp = lam * mu * (1.0 + mu)       # F''(r_+)
    dF = lam * (1.0 + mu) ** 2 * (1.0 - 2.0 * mu) / 48.0
    r_plus = np.sqrt(0.5 * (1.0 + mu))
    if d == 1:
        ratio = 0.5 * np.sqrt(Fpp / a) * np.exp(dF / T)
    else:
        ratio = np.sqrt(T * Fpp / (2 * np.pi)) / (a * r_plus) * np.exp(dF / T)
    return 1.0 / (1.0 + ratio)


lam = 2.0
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), constrained_layout=True)

for ax, d in zip(axes, [1, 2]):
    for th, c in zip([2e-4, 1e-3, 5e-3], ['#1b4965', '#c96a3c', '#8c2f39']):
        Ls = np.linspace(1 - lam / 4 + 1e-4, 1.0 - 1e-4, 500)
        x = -np.log(Ls)
        ax.plot(x, [frac_exact(L, lam, th, d) for L in Ls],
                color=c, lw=2.2, label=rf'$\vartheta={th}$')
        ax.plot(x, [frac_laplace(L, lam, th, d) for L in Ls],
                color=c, lw=1.2, ls='--')
    # spinodal (mu=0) and Maxwell (mu=1/2) points
    ax.axvline(-np.log(1 - lam / 4), color='0.75', lw=0.9, ls=':')
    ax.axvline(-np.log(1 - 3 * lam / 16), color='0.75', lw=0.9, ls=':')
    ax.set_title(rf'$d={d}$,  $\lambda={lam}$;  solid = exact, dashed = Laplace', fontsize=10)
    ax.set_xlabel(r'$-\log\mathcal{L}$')
    ax.set_ylabel(r'$f$')
    ax.set_ylim(-0.03, 1.05)
    ax.legend(frameon=False, fontsize=9, loc='lower left')
    ax.spines[['top', 'right']].set_visible(False)

fig.savefig(f'{OUT}/setup1_exact_vs_laplace.png', dpi=160)
print('saved', f'{OUT}/setup1_exact_vs_laplace.png')

print(f'\nlam={lam}, d=2, theta=1e-3   (mu -> 0 is the spinodal)')
print(f"{'L':>8} {'mu':>7} {'exact':>9} {'laplace':>9}")
for L in [0.505, 0.53, 0.56, 0.60, 0.625, 0.66, 0.75, 0.90, 0.98, 0.999]:
    print(f'{L:8.3f} {_mu(L, lam):7.3f} {frac_exact(L, lam, 1e-3, 2):9.4f} '
          f'{frac_laplace(L, lam, 1e-3, 2):9.4f}')
