---
title: Spherical harmonics and Wigner D matrices
---

Working notes for `mostlyharmless`. Notation is pinned to the code: coefficients are
indexed by $m+l$ with $m$ **ascending**, and $D^l_{m'm}(R) = \langle lm'|U(R)|lm\rangle$.

## Spherical harmonics

Orthonormal, Condon–Shortley phase carried inside $P_l^m$, polar $\theta$, azimuthal $\varphi$:

$$Y_l^m(\theta,\varphi) = \sqrt{\frac{2l+1}{4\pi}\,\frac{(l-m)!}{(l+m)!}}\; P_l^m(\cos\theta)\; e^{im\varphi}.$$

Negative order, conjugate symmetry, and the reality condition for a real field:

$$P_l^{-m} = (-1)^m \frac{(l-m)!}{(l+m)!} P_l^{\,m},
\qquad
Y_l^{-m} = (-1)^m \overline{Y_l^{\,m}},
\qquad
a_{l,-m} = (-1)^m \overline{a_{lm}}.$$

Orthonormality and the addition theorem:

$$\int Y_l^m \,\overline{Y_{l'}^{m'}}\; d\Omega = \delta_{ll'}\delta_{mm'},
\qquad
\sum_{m=-l}^{l} Y_l^m(\hat n_1)\,\overline{Y_l^m(\hat n_2)}
= \frac{2l+1}{4\pi}\,P_l(\hat n_1 \cdot \hat n_2).$$

Low orders, for hand checks:

$$Y_0^0 = \frac{1}{\sqrt{4\pi}},
\qquad
Y_1^0 = \sqrt{\frac{3}{4\pi}}\cos\theta,
\qquad
Y_1^{\pm 1} = \mp\sqrt{\frac{3}{8\pi}}\,\sin\theta\,e^{\pm i\varphi},
\qquad
Y_2^0 = \sqrt{\frac{5}{16\pi}}\,(3\cos^2\theta - 1).$$

The real basis replaces $e^{im\varphi}$ by $\cos m\varphi$ for $m>0$ and by $\sin|m|\varphi$
for $m<0$, with a factor $\sqrt{2}$ on the $m \neq 0$ terms.

## Wigner $D$ matrices

A rotation mixes spherical harmonics **within** each degree and never across degrees:

$$Y_l^m(R^{-1}\hat n) = \sum_{m'=-l}^{l} D^l_{m'm}(R)\; Y_l^{m'}(\hat n),
\qquad
D^l_{m'm}(R) = \langle lm'|U(R)|lm\rangle .$$

For $R = R_z(\alpha)R_y(\beta)R_z(\gamma)$ the matrix factorises into two phase ramps
around one real matrix:

$$D^l_{m'm}(R) = e^{-im'\alpha}\; d^l_{m'm}(\beta)\; e^{-im\gamma},
\qquad
d^l_{m'm}(\beta) = \langle lm'|e^{-i\beta J_y}|lm\rangle \in \mathbb{R}.$$

The tilt is real because $J_y$ is purely imaginary; the $z$-rotations are diagonal because
$Y_l^m$ is an eigenfunction of $J_z$ with eigenvalue $m$. Properties:

$$D(R_1)D(R_2) = D(R_1R_2),
\qquad
D(R)^{\dagger} = D(R^{-1}),
\qquad
d^l_{m'm}(0) = \delta_{mm'},
\qquad
d^l_{00}(\beta) = P_l(\cos\beta),$$

$$d^l_{m'm} = (-1)^{m-m'} d^l_{mm'},
\qquad
d^l_{m'm}(-\beta) = d^l_{mm'}(\beta).$$

The $l=1$ case, rows and columns ordered $(1,0,-1)$:

$$d^1(\beta) =
\begin{pmatrix}
\frac{1+\cos\beta}{2} & -\frac{\sin\beta}{\sqrt{2}} & \frac{1-\cos\beta}{2} \\
\frac{\sin\beta}{\sqrt{2}} & \cos\beta & -\frac{\sin\beta}{\sqrt{2}} \\
\frac{1-\cos\beta}{2} & \frac{\sin\beta}{\sqrt{2}} & \frac{1+\cos\beta}{2}
\end{pmatrix}.$$

The explicit sum used in the code:

$$\begin{aligned}
d^l_{m'm}(\beta) &= \sqrt{(l+m)!\,(l-m)!\,(l+m')!\,(l-m')!} \\
&\quad \times \sum_s \frac{(-1)^{m'-m+s}}{s!\,(l+m-s)!\,(l-m'-s)!\,(m'-m+s)!} \\
&\quad \times \left(\cos\tfrac{\beta}{2}\right)^{2l+m-m'-2s}
\left(\sin\tfrac{\beta}{2}\right)^{m'-m+2s},
\end{aligned}$$

with $s$ ranging over all values that keep every factorial argument non-negative.

## Coefficients under a frame change

If a field has body-frame coefficients $a_{lm}$, its observer-frame coefficients are

$$b_{lm'} = \sum_{m=-l}^{l} D^l_{m'm}(R)\; a_{lm},$$

contracting the **second** index of $D$; $l$ is preserved. In the conventions of
`spot_test.py` and `spot_pattern.py` — observer axes $(u,v,w) = (\text{East}, \text{North},
\text{line of sight})$, body $z$ the spin axis, and

$$R(\alpha,\beta,\mathrm{PA}) = R_z\!\left(\tfrac{\pi}{2} - \mathrm{PA}\right) R_y(\beta) R_z(\alpha)$$

mapping body to observer — this becomes

$$b_{lm'} = e^{-im'(\pi/2-\mathrm{PA})} \sum_{m=-l}^{l} d^l_{m'm}(\beta)\; e^{-im\alpha}\; a_{lm},$$

that is: a diagonal phase ramp in $m$, a dense real orthogonal $(2l+1)\times(2l+1)$ block,
and a diagonal phase ramp in $m'$. Because $D$ is a representation, two successive
rotations compose as $D(R)D(R_0) = D(RR_0)$.

## Traps

- **Ordering.** $d^l$ is not symmetric ($d_{m'm} = (-1)^{m-m'}d_{mm'}$), so an unstated
  row/column order silently transposes it. Ascending $m$ throughout here.
- **Which index is contracted.** $b_{lm'} = \sum_m D^l_{m'm}a_{lm}$, not
  $\sum_m D^l_{mm'}a_{lm}$.
- **Condon–Shortley.** Drop the $(-1)^m$ and every odd-$m$ coefficient flips sign, and the
  reality condition changes with it.
- **Real versus complex.** They differ by $\sqrt{2}$ per $m \neq 0$ term, and by a sign.

## How these were checked

The harmonic normalisation and phase agree with SymPy's `Ynm` to machine precision, and the
exact quadrature reproduces orthonormality to $\sim10^{-15}$. The rotation law above was
checked against a direct geometric rotation of the function on a Gauss–Legendre grid
(agreement $1.1\times10^{-15}$), and the explicit $d^l$ sum against $e^{-i\beta J_y}$
($9.99\times10^{-15}$ for $l \le 10$). The addition theorem and the symmetry relations are
standard and were not re-derived here; prose is not machine-checked.

## References

- NIST *Digital Library of Mathematical Functions*, chapters 14 (Legendre and spherical
  harmonics) and 34 ($3j$, $6j$, $9j$ symbols and $D$ matrices).
- D. A. Varshalovich, A. N. Moskalev and V. K. Khersonskii, *Quantum Theory of Angular
  Momentum* (World Scientific, 1988) — the standard reference for $D$-matrix conventions.
