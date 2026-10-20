# Roadmap — `mostlyharmless`

Analytic interferometric visibilities of spherical-harmonic surface maps.

**Status: living document.** Expected to be modified as the work proceeds. Sections 1–6
are the plan; section 0 records decisions already taken and section 7 collects what is
still open.

---

## 0. Decisions recorded

| # | decision |
|---|---|
| D1 | **Banner name:** `mostlyharmless.visylm` — short for "visibility of Ylm" — holding spherical harmonics, Wigner matrices, Zernike polynomials and the visibility chain under one roof. Chosen in place of `harmonics`, which was rejected as too bland. Names checked and found already in use elsewhere: `spherical` (which has a `spherical.modes` submodule), `harmonix` (the Dholakia & Pope code), `modvis` / `modevis`. A search turned up no use of `visylm`; the check was not exhaustive. |
| D2 | **No notebooks.** They work against lemma compression. The certificate is pytest; the publication is LaTeX. |
| D3 | **Emit data plus LaTeX fragments, never generated code.** The general implementation is hand-written; the notebook/script only produces tables. |
| D4 | **uv-tracks and baseline geometry are out of scope** — done by other software. Visibility functions receive `(u, v)` as given; `q` and `psi` are derived internally and do not appear in a caller's code. |
| D5 | **Coefficients are exact rationals**, stored inspectably as integer numerator/denominator pairs. |
| D6 | **`make all` checks before it documents.** The default goal runs the certificates and stops if they fail, so the PDF is never built from a red suite: its verification register asserts what the certificates establish. `make docs` on its own stays the plain LaTeX path, needing no Python unless a generated file is missing, which is the escape hatch for a typo or a caption. |
| D7 | **Measure:** planar over the projected disk, `dx_E dx_N`. The exact solid-angle alternative (`dOmega = dx_E dx_N / n`) is stated as a remark, not implemented. |
| D8 | **Conventions:** Condon–Shortley phase inside `P_l^m`; complex harmonics `Y_l^m` with `m >= 0` only, used in the derivations; maps use the real basis `R_l^0 = Y_l^0`, `R_l^{m,c} = sqrt(2) Re Y_l^m`, `R_l^{m,s} = sqrt(2) Im Y_l^m`, expanded with real coefficients `c_lm, s_lm`; **negative orders are never formed**, so `\|m\|` is never needed and there is no signed-order convention to get wrong; `u, v` in wavelengths; the `2*pi` inside the exponent; observer frame `(East, North, line of sight)`. |
| D9 | **Notation:** `rho = sin(theta)`; `x_E, x_N` direction cosines; `q = sqrt(u^2+v^2)`; `psi = arctan(v/u)`; `a = 2*pi*q`; `k = (n-m)/2` for the Jacobi degree. `\|m\|` is never used as a symbol: orders are always written signed, and the parity condition is stated as `l+m` odd or even. |
| D10 | **Verification rungs** are labelled on every claim: `exact` (free parameters), `range-limited` (exactly zero for all integer cases in a stated range — a proof of the finite statement), `numerical` (agreement to a stated tolerance). Prefer exact expansion and coefficient comparison to `simplify`, so a pass is a proof of the instance rather than trust in a heuristic. |
| D11 | **Non-spherical sources: perturb in a deformation field.** The surface is `r(n) = R[1 + sum eps_lm Y_lm(n)]`, with an independent brightness perturbation, rather than a single flattening scalar; the `eps_lm` are the fitting parameters. The first-order visibility correction has **two** pieces — the interior integrand change and a **limb boundary term** from the moving silhouette — and both are required. Shape and gravity darkening enter at the same order and must be included together. Validate by the `eps^2` scaling of the residual, not merely by agreement. The reference implementation is a direct sky-plane ray-trace, chosen because it shares no algebra with the perturbation. |
| D12 | **Triaxiality is required** — spinning stars in a binary — so a single flattening parameter will not do. Rotational and tidal quadrupoles about different axes are the natural input; the deformation-field parameterisation of D11 covers this. |
| D13 | **Coefficients are half-complex and negative-order-free.** A real map is stored as `a_l0 = c_l0` and `a_lm = (c_lm - i s_lm)/sqrt(2)` for `m > 0`, one complex array per degree of length `l+1`, equivalently `I = sum_l [a_l0 Y_l0 + sum_{m>0} 2 Re(a_lm Y_lm)]`. The map is real by construction — there is no second set of stored coefficients to disagree with the first — and no unphysical direction for an optimiser to populate. Fitting goes through a flat real vector of exactly `(l+1)^2` numbers. The order is non-negative wherever it occurs, so the magnitude appears nowhere. |
| D14 | **Negative orders are indices, never coefficients.** `Y_l^-m` is *defined* as `(-1)^m conj(Y_l^m)`, which coincides with the standard convention, and is used only in the rotation machinery, where the order indexes the representation. Coefficients are formed only for `m >= 0`; where the rotation would call for a negative-order coefficient its value is substituted, so the transform reads `b_lm' = D^l_{m'0} a_l0 + sum_{m=1..l} [ D^l_{m'm} a_lm + (-1)^m D^l_{m',-m} conj(a_lm) ]` with `m' >= 0`. The `m = 0` term appears once, not twice, because `a_l0` is real. The classical Wigner `d` and `D` therefore stay in the manual, rather than being relegated to the certificates. |
| D15 | **Flux, mean intensity, and normalisation.** With the convention of D8 (no prefactor), `V(0,0)` is the total flux — a linear functional of *all* the coefficients, not of `a_00` alone, because the disk integral weights each harmonic by the foreshortening factor `cos(theta)` rather than being orthogonal as the sphere integral is. Normalising the visibility to unity at the origin therefore means dividing by the instantaneous flux, which for a rotating spotted star is time-dependent and discards the light-curve information; the unnormalised convention is kept. Separately, `a_00` is `sqrt(4*pi)` times the mean intensity over the **sphere** — not the mean over the visible disk, and not the flux. Fixing `a_00` is a choice of overall brightness scale and does not by itself make `V(0,0) = 1` unless the map is constant. Where a flux constraint is wanted, its row is the mode visibilities at `q = 0`, which Stage 2 produces anyway. |
| D16 | **No convenience features the user can add trivially, and `a_00 = 1`.** Two assumptions do the work that parameters, helper methods or wrappers otherwise would: the sphere has **unit radius**, and maps are normalised so that **`a_00 = 1`**, which makes every other coefficient a mode amplitude relative to uniform brightness. Both are documented assumptions rather than code constraints or conveniences such as `normalised()`, because a wrapper reduces the readability of the code for no real gain. The convention is preserved by rotation, the `l = 0` block being untouched by it, so the free coefficients of a map are those with `l >= 1`. The container itself stays general — it still holds arbitrary maps, including the zero map — and it is the modelling convention, not the container, that fixes `a_00`. |
| D17 | **No `special/` package.** Bessel functions come straight from scipy and sympy, with no convention of ours to protect; the classical anchors — the uniform disk, and later the limb-darkened disk — are visibility functions and live in `visylm/visibility.py`. A module with one customer and no convention of its own is folded in rather than kept. |
| D18 | **Negative orders are extended, not computed.** The visibility is extended to negative orders by `V_l^-m = conj(V_l^m)`, the analogue of the harmonic extension `Y_l^-m = (-1)^m conj(Y_l^m)`; the sign is absent because the `(-1)^m` of the harmonic extension is cancelled by the `(-1)^m` that reversing the baseline contributes. A stored `a_lm` therefore contributes `a_lm V_l^m + (-1)^m conj(a_lm) conj(V_l^m) = 4 pi (-i)^m N_lm R_lm(q) Re(a_lm exp(i m psi))`, which is equation (35) of the manual. There is consequently no need for two real branches, in the derivation or in the code: one mode function suffices, and the real branches become a remark for comparing with the literature. |
| D19 | **Coefficients are computed just in time, and the emitter is separate.** The even-family coefficients are exact rationals computed on demand and memoised, not stored in a table: there is then no maximum `l` to choose, and no data file for the package to depend on. A short program, `symbolic/emit.py`, calls the library and writes the LaTeX fragment the manual inputs, so the only code that knows about LaTeX is the emitter; `symbolic/figures.py` does the same for the figures, importing matplotlib lazily so that it stays an optional extra. `visylm/tables.py` is dropped from the planned layout, nothing having to load a table. The manual's printed table goes to `l <= 3`. |
| D20 | **The application is intensity (HBT) interferometry.** The measured quantity is `|V|^2`, so the phase of `V` is not observable — but it is retained throughout, both for visualisation and because discarding it early would make the intermediate results harder to check, and squared only at the end. Two design consequences: the deliverable for deformed sources is a *diagnostic* evaluable on a sparse set of `(u,v)` points rather than an image, since a fitting loop over orientations makes a numerical transform per model too expensive; and the (u,v) sampling being sparse means the high-frequency behaviour of the transform is not the design driver. |
| D21 | **Brightness is proportional to `g`, for now.** The surface brightness of the deformed star is set by the bolometric von Zeipel law, `T_eff^4 ∝ g` and hence `I ∝ g`, which is the simplest law that is also the standard one for radiative envelopes. The literature favours a weaker dependence — `T_eff ∝ g^beta` with `beta ≈ 0.25` radiative but nearer 0.08 for convective envelopes, so `I ∝ g^(4 beta)` and about `g^0.3` in the convective case, a much flatter law. This is recorded rather than acted on: the brightness law enters the renderer as a single exponent, so changing it later is a parameter and not a rewrite, and the figures will say which exponent they use. |
| D22 | **The observable is normalised at the origin.** Intensity interferometry measures `|V|^2` in units of the zero-baseline correlation, so what is modelled is `V/V(0,0)`, which is 1 at the origin by definition. Every difference is therefore taken between *normalised* visibilities, so `Delta V(0,0) = 0` by construction, and the absolute flux — which the normalisation removes — is not modelled at all. This also settles the normalisation rung 1 had left open: the affine model is normalised by its own origin, and the unobservable flux difference between the star and its reference is divided out rather than carried. The *squared* difference `|V_tri|^2 - |V_ref|^2` is deliberately not used: the phase is retained, and the difference is taken on `V` itself. Worth recording why the normalisation is needed even against a volume-equivalent reference: the flux is the *projected area*, not the volume, so the two agree only for a sphere, and the deformed star's projected area depends on the viewing direction. |
| D23 | **The azimuth reference is North, through East.** Position angles are measured from the observer's North toward East, the astronomical convention, which fixes both the reference axis and the sense. With `u` along North and `v` along East (D24), `psi = arctan2(v, u)` *is* the position angle, so `PA = psi` and no offset appears anywhere. Had the older convention been kept — `u` along East — the consequence would have been `PA = pi/2 - psi`, the other sign corresponding to North-through-West, and the two differ by a *reflection* rather than a rotation, flipping every odd-`m` sign and changing the meaning of the mirror symmetry of section 8. The body's longitude zero remains the body `x` axis, which is already the convention for the harmonic azimuth. The real-basis branch labelling (section 7 item 1) is a separate convention and is untouched by this. |
| D24 | **The sky's orientation is shared by the ground and the image, departing from the radio convention.** `u` is the **North** component of the baseline and `v` the **East** component, so `V(u,v) = ∫ I(x_E, x_N) exp(-2 pi i (u x_N + v x_E))`. The third axis points from the sky *toward* the observer — into the ground for an observer looking up — so (North, East, toward the observer) is right-handed *given* the sky convention that North is up and East is to the left. The point of the swap is that `psi` becomes the position angle directly (D23), and that the coordinates on the sky and on the ground have the same orientation. The cost is deliberate: radio interferometry takes `u` as East-West and `v` as North-South, so this is a departure, justified only because our `(u,v)` are used throughout as *sky-plane* coordinates rather than as baseline hardware coordinates. Every image is therefore drawn with North vertical and East horizontal increasing leftward, images and `(u,v)` planes alike. *Status: decided, not yet implemented.* The obstacle is that the mode reduction's `psi` is tied to the harmonics' `phi`, which is measured from the body `x` axis; relabelling `(u,v)` alone rotates `psi` by `pi/2` while `phi` stays put, and V15 detects exactly that as a 0.14 residual. Doing it properly means re-referencing the body `x` axis to North as well — so the display maps the body `x` to the vertical — together with the transform and the V15 pairing. |
| D25 | **The interactive GUI is a separate, optional front end.** Matplotlib widgets, since there is no Jupyter here, and nothing added to the library. Two radio buttons --- image or `(u,v)`, analytic or numeric --- between them select the four combinations the Stage 6 figures lay out. The spin's position angle is *pinned to zero*, because for a single observation the sky's orientation is a gauge; the model keeps it for real data. The light is **fixed at the von Zeipel law**, `I ∝ g` (D21), with the exponent to be exposed as a slider later; and the analytic row keeps **only `l <= 2`**, discarding `l = 4` and above. Note what that is and is not: it is *not* "first order in `eps`". The `l = 4` content first appears at `eps^2` from `P_2^2`, so dropping it *is* the `eps^2` truncation — but the projection keeps the `l = 2` coefficient to all orders, including its own `eps^2` correction, so the retained set is "`l <= 2`, resummed in `eps`". The difference from a true first-order truncation is itself `O(eps^2)`, below the model's error, so the simpler and more accurate choice is `l <= 2`: `SPHERE_L_MAX = 2`, with no further decision to make. The free parameters are `eps_tide`, `eps_cent`, `i`, `lambda_tide` and `phi_tide`. The spin phase `alpha` is deliberately *absent*: this front end is for a spin and tide alone, and with no starspots the phase is degenerate with the tide's longitude — only their difference is observable — so a slider for it would move a gauge. The **latitude slider is constrained** rather than the amplitudes: changing `lambda` can only tighten the closure bound, so the slider stops at the largest latitude for which the current amplitudes are still feasible (closed form, a quadratic in `cos^2 theta`), instead of letting `eps` move underneath the user. Updates are **coalesced, not cancelled** --- numpy cannot be interrupted --- so a worker thread computes, a latest-values slot holds the newest parameters, and results computed from stale ones are discarded; the GUI thread only draws. Preview settings (8000 samples, 120 pixels, 20 bisection steps) keep a redraw in tens of milliseconds, and the level radii and brightness pattern are cached across orientation changes because they depend on the deformation alone. |

---

## 1. Formula inventory

Status: **derived** = worked out and checked during design discussion; **new** = still to do.

| # | formula | status | section |
|---|---|---|---|
| F1 | `Y_l^m` for `m >= 0`: definition, CS phase; the real basis `R_l^{m,c}`, `R_l^{m,s}` and the map expansion | derived | 02 |
| F2 | Orthonormality; addition theorem | derived | 02 |
| F3 | Real basis (cos/sin) and its normalisation | derived | 02 |
| F4 | Rotation law `Y_lm(R^-1 n) = sum D^l_{m'm}(R) Y_lm'(n)` | derived | 03 |
| F5 | zyz factorisation `D = exp(-i m' alpha) d(beta) exp(-i m gamma)`; explicit sum for `d` | derived | 03 |
| F6 | Properties: unitarity, `d(0)=I`, `d_00 = P_l(cos beta)`, symmetries | derived | 03 |
| F7 | Coefficient transform `b_lm' = sum_m D^l_{m'm}(R) a_lm`, body → observer | derived | 03 |
| F8 | Parity factorisation of `i^m` into a family phase times a real sign | derived | 03 |
| F9 | Jacobi–Anger reduction: `V = 2 pi (-i)^m e^{i m psi} int f(rho) J_m(2 pi q rho) rho drho` | derived | 04 |
| F10 | Planar measure convention, and the `1/n` alternative | derived | 04 |
| F11 | Parity dichotomy: `l+m` odd ⇒ factor `sqrt(1-rho^2)`; even ⇒ polynomial in `rho^2` | derived | 04 |
| F12 | Odd family: parity doubling, Legendre recurrence, the `K` identity, the collapse | derived | 05 |
| F13 | Odd result: single spherical Bessel `j_l(2 pi q)/q` with double-factorial coefficient | derived | 05 |
| F14 | Zernike `R_n^m` from Jacobi; orthogonality on the disk | derived | 06 |
| F15 | Zernike projection coefficients for the even family | structure derived; table to build | 06 |
| F16 | Zernike–Hankel identity: `int R_n^m J_m rho drho = (-1)^((n-m)/2) J_{n+1}(a)/a` | derived | 06 |
| F17 | Even result: finite sum of odd-order ordinary Bessels | derived | 06 |
| F18 | Full-map visibility by linearity over modes | trivial | 07 |
| F19 | Limb-darkening profile and its `m = 0` harmonic coefficients | **new** | 08 |
| F20 | Product identity / Gaunt coupling; limb darkening raised in order | **new** | 08 |
| F21 | Anchors: uniform disk `2 J_1(x)/x`; linearly limb-darkened disk | known classically | 09 |

---

## 2. Package layout and public API

```
src/mostlyharmless/
├── __init__.py
├── visylm/
│   ├── __init__.py
│   ├── spherical_harmonics.py    # F1-F3
│   ├── coefficients.py           # the Coeffs container (D13)
│   ├── wigner.py                 # F4-F8
│   ├── zernike.py                # F14
│   ├── visibility.py             # F9-F13, F15-F18, plus the classical anchors
│   └── deformation.py            # the potential, its level surface, the gravity there
└── symbolic/
    ├── __init__.py
    ├── derivations.py            # one function per certificate, returning a residual
    ├── render.py                 # the reference ray renderer of Stage 6
    ├── affine.py                 # rung 1: the projected stretch and its visibility
    ├── emit.py                   # LaTeX fragments for the manual; run as a module
    └── figures.py                # the figures of docs/generated/; run as a module
```

### `mostlyharmless.visylm`

```python
# --- F1-F3 ---------------------------------------------------------------
sph_harm(l, m, theta, phi)                  -> complex   # CS phase, theta polar
real_sph_harm(l, m, theta, phi)             -> float     # cos/sin basis
sph_harm_addition(l, th1, ph1, th2, ph2)    -> complex   # theorem LHS, for checks

# --- F4-F8 ---------------------------------------------------------------
wigner_d(l, beta)                           -> ndarray   # d^l, index [m'+l, m+l]
wigner_D(l, alpha, beta, gamma)             -> ndarray   # full D
rotate_coefficients(coeffs, alpha, beta, pa) -> Coeffs   # body -> observer, our (alpha, beta, PA)

# --- F14-F15 -------------------------------------------------------------
zernike_radial(n, m, rho)                   -> ndarray   # R_n^m, numeric, via Jacobi
zernike_polynomial(n, m)                    -> sympy expr # R_n^m, exactly
mode_polynomial(l, m)                       -> sympy expr # P_l^m(sqrt(1-rho^2)) = rho^m A(rho^2)

# --- F9-F18 --------------------------------------------------------------
visibility_mode(l, m, u, v)                 -> ndarray   # complex, Y_l^m, m >= 0
visibility(coeffs, u, v)                    -> ndarray   # full map, by linearity
radial_integral(l, m, q)                    -> ndarray   # the quadrature reference
zernike_coefficients(l, m)                  -> tuple     # even-family c_k, exact rationals

# --- the deformation field (Stage 6) -------------------------------------
Deformation(tide, tide_direction, centrifugal, spin_direction)
    .potential(x) / .gradient(x)             # Phi and grad Phi, analytic
    .level_radius(direction)                 # the inner crossing of Phi = -1
    .surface_gravity(direction)              # |grad Phi| there
    .axis_lengths() / .volume() / .equivalent_radius()
    .first_order_radii(direction)            # for comparison with the solver

# --- classical anchors ---------------------------------------------------
visibility_uniform_disk(u, v, theta_d)      -> ndarray   # 2 J_1(x)/x
visibility_limb_darkened(u, v, theta_d, ld_coeffs) -> ndarray   # F21
```

Public entry points take `(u, v)`; `q` and `psi` are derived internally (D9) and never appear in a
caller's code. The intermediate `radial_integral(l, m, q)` is the workhorse and is deliberately
public, because it dispatches on the parity of `l + m` and so keeps the odd and even paths
separately visible and separately tested. `visibility_mode` applies the closed-form azimuthal and
phase factors; `visibility` sums over modes.

A `Coeffs` container holds the ragged `(l, m)` coefficient set so that the ordering and the
parity split are not the caller's problem. Kept inspectable (D5).

**Naming collision to settle.** The limb-darkening literature uses `u`, `u_1`, `u_2` for the
coefficients of the law, which clashes with `u` as a baseline coordinate. The signatures above use
`ld_coeffs`, and the manual will need the same distinction stated explicitly.

The classical anchors live in `visylm/visibility.py` rather than in a special-functions module.
The Bessel functions come straight from scipy and sympy, and the anchors are visibility functions,
so a separate module would have a single customer and no convention of its own (D17).

### `mostlyharmless.symbolic`

```python
check_odd_collapse(l_max)                   -> residual
check_zernike_jacobi(n_max)                 -> residual
check_zernike_hankel(n_max)                 -> residual
check_rotation_operator(l_max)              -> residual
odd_case_expression(l, m)                   -> sympy expression
even_case_coefficients(l, m)                -> list[Fraction]
generate_table(l_max, path)                 -> None
emit_latex_table(l_max)                     -> str
```

Plus a CLI so the Makefile can drive emission:
`python -m mostlyharmless.symbolic emit --lmax 10`.

---

## 3. Data artifacts

| artifact | contents | notes |
|---|---|---|
| `docs/generated/coefficients.tex` | the even-family coefficients for the manual | emitted by `symbolic/emit.py` from the just-in-time coefficients (D19); printed to `l <= 3` |
| `docs/generated/verification.tex` | register of certificates: tag, statement, function, rung, scope | emitted from the test list |
| `docs/generated/*.pdf` | the starspot on both hemispheres, Re and Im V(u,v) normalised, the l_max convergence, and for Stage 6 the triaxial star with its gravity darkening and the visibility difference from its volume-equivalent sphere | emitted by `symbolic/figures.py`; rendered by the test suite at small size |

There is no `data/` directory and no coefficient data file. Coefficients are computed just in
time in exact rational arithmetic (D19), so nothing has to read a table, and the emitters are the
only code that knows about LaTeX or matplotlib.

All generated files are committed. The LaTeX fragments have a staleness test that regenerates them
in memory and fails on mismatch. The figures do not: byte-comparing matplotlib output across
versions would fail for reasons that mean nothing, so they are covered by a rendering test at small
size instead.

---

## 4. Certificates

One tag, one function, one row in the register, one citation in the text.

| tag | claim | rung |
|---|---|---|
| V1 | `Y_l^m` orthonormal on the exact quadrature grid | exact, `l <= L` |
| V2 | harmonic extension `Y_l^-m = (-1)^m conj(Y_l^m)`, and that the implementation reproduces it | exact |
| V3 | addition theorem | numerical: 40-digit SymPy, plus float at higher `l` (the symbolic route tangles in the sign of `sin theta`) |
| V4 | explicit `d` sum equals `exp(-i beta J_y)` | numerical, `l <= 6` |
| V5 | `d^T d = I`; `d(0) = I`; `d_00 = P_l(cos beta)`; symmetry | numerical, `l <= 6` |
| V6 | rotation equals a direct geometric rotation; and the non-negative-order form equals the full-range transform | numerical, `l <= 4` |
| V7 | composition of two rotations | numerical, `l <= 3` |
| V8 | Jacobi–Zernike relation against the explicit series | exact, `n <= 8` |
| V9 | Zernike orthogonality | exact rational, `n <= 8` |
| V10 | Zernike–Hankel identity, by power series | exact, `n <= 6`, powers to `a^18` |
| V11 | the collapse coefficients coincide | exact, `l <= 8` |
| V12 | projection coefficients equal coefficient matching | exact rational, `l <= 6` |
| V13 | odd closed form against the quadrature, `q = 0` included | numerical, `l <= 5` |
| V14 | even Bessel sum against the quadrature, `q = 0` included | numerical, `l <= 6` |
| V15 | the whole map against a 2-D FFT of the rendered image, low-frequency band | numerical |
| V16 | uniform disk reduces to `2 J_1(x)/x`, and the flux is `sqrt(pi)/2` | numerical, `l = 0` exactly |
| V17 | Jacobi–Anger azimuthal factors | numerical, roundoff level |
| V18 | the `K` identity, in our Condon–Shortley convention | numerical, `l <= 6` |
| V19 | the deformation field against its first-order relations | numerical, residual/`eps^2` below 50 |
| V20 | the level surface is level, and the analytic gradient matches a difference | numerical, `1e-8` |

Two habits, from D10: prefer exact expansion and coefficient comparison to `simplify`; and
state the rung and the scope next to every claim, in the register and in the text.

---

## 5. Documentation structure

```
docs/
  mostlyharmless.tex
  sections/
    01-conventions.tex      symbol table, measures, phase conventions   (D8, D9)
    02-harmonics.tex        F1-F3
    03-wigner.tex           F4-F8
    04-visibility.tex       F9-F11, the master reduction
    05-odd.tex              F12-F13, lemmas 1-4
    06-even.tex             F14-F17, the Zernike route
    07-full-map.tex         F18
    08-triaxial.tex         the deformation field, the affine model, rung 2
    09-limb-darkening.tex   F19-F20
    10-numerics.tex         accuracy, range of validity
  generated/                emitted fragments
  refs.bib
```

Every derivation section opens with a conventions preamble pulled from section 01, and each
checked step carries its tag. Section 05 is the template for the compressed lemma style:
structural steps in full, mechanical algebra delegated to a named certificate.

---

## 6. Build order, with gates

*Status: Stages 0 and 1 (a and b) closed — 52 tests green, manual 8 pages, 13 certificates.*

**Stage 0 — skeleton.** Repo layout, `pyproject.toml`, pytest harness, section 01 with the
symbol table. **Done.** Gate passed: 9 smoke tests green under Thonny's Python 3.10.11
(numpy 2.2.6, scipy 1.15.3, sympy 1.13.3, pytest 9.1.1), and the manual compiles to 3 pages with
cross-references resolved.

**Stage 1 — harmonics and Wigner.** F1–F8, V1–V7. Largely transcription of material already
verified during design, plus the certificate workflow itself. Split in two, because only the
second needs the `Coeffs` container and the real/complex interface settled (section 7).

**Stage 1a — spherical harmonics.** F1–F3, V1–V3; `visylm/spherical_harmonics.py`,
`symbolic/derivations.py`, section 02. **Done.** Gate passed: 24 tests green, manual 5 pages,
cross-references resolved. It also caught a real error in section 01 on the way — the harmonic
normalisation had been written with `|m|` instead of the signed `m`, which leaves every negative
order unnormalised by `[(l-|m|)!/(l+|m|)!]^2`, invisible for `m >= 0`. Found by writing the
reality relation as certificate V2.

**Stage 1b — Wigner.** F4–F8, V4–V7; `visylm/wigner.py`, section 03. **Done.** Gate passed: 52
tests green, manual 8 pages, sections 02–03 compile. Both blocking decisions were taken first:
half-complex `Coeffs` (D13) and negative orders as indices only (D14).
*Gate for both:* all certificates pass; sections 02–03 compile.

**Stage 2 — master reduction.** F9–F11, plus V16 and V17. `visylm/visibility.py`,
`docs/sections/04-visibility.tex`. **Written**; gate pending the test run.
*Gate:* a uniform disk comes out as `2 J_1(x)/x`, and the flux as `sqrt(pi)/2`. Note this is
*exact*, not a convergence statement: the map of a uniformly bright star is the `l = 0` harmonic
alone, so the series terminates after one term, and the truncation ringing that afflicts
sharp-edged features does not arise.

**Stage 3 — odd family.** F12–F13, V11, V13, V18. `docs/sections/05-odd.tex`, the closed form in
`visylm/visibility.py`. **Done.** 70 tests pass (62 + 8). The closed form is private and
reached through `visibility_mode`, which dispatches on the parity of `l+m`; `radial_integral` stays
pure quadrature, since a reference that dispatches to the formula under test would make the
comparison vacuous.
*Gate:* the closed form matches direct quadrature across a grid of `(l, m, q)`, including `q = 0`,
where the limit is `1/3` for `l = 1` and zero for `l >= 2`.

**Stage 4 — even family.** F14–F17, V8–V10, V12, V14; the library coefficients and the LaTeX emitter.
**Done.** 82 tests pass (70 + 12). Coefficients computed just in time in exact rational arithmetic and memoised (D19), with no
data file; `symbolic/emit.py` writes `docs/generated/coefficients.tex` for the manual, printed to
`l <= 3`; the even closed form is dispatched in `visibility_mode`, so both families are then closed.
*Gate:* the coefficients are exact and self-consistent (V9 orthogonality, V12 coefficient matching);
the Bessel sums match quadrature (V14) with `q = 0` included; the generated fragment is committed,
passes its staleness test, and the manual compiles with it.

**Stage 5 — full map.** F18, V15. `docs/sections/07-full-map.tex`, `check_full_map_fft`.
**Done.** 90 tests pass (87 + 3). The comparison is over a low-frequency band, since the
disk's rim is a discontinuity and its Gibbs ringing contaminates high frequencies at any
resolution; Hermitian symmetry is checked separately and tightly.
*Gate:* FFT agreement on a rendered map. This is also where truncation ringing shows up: a map
that is discontinuous *on the sphere* — a sharp-edged spot or a polar cap — has a slowly
converging expansion, and the reference image and the harmonic sum will differ near the edge at
any finite `l_max`. Worth a stated convergence check in `l_max`, since it is the one place the
truncation error is visible.

**Stage 6 — triaxial deformation.** Section 8's method, in this order:

1. **The potential and its level surface.** **Done.** `visylm/deformation.py` carries `phi_1` = a
   tidal quadrupole plus a centrifugal one, both `l = 2`, the latter with an `l = 0` piece that
   renormalises the effective monopole. `level_radius` finds the *inner* crossing of `Phi = -1` and
   `surface_gravity` takes `|grad Phi|` there. Certificates V19 (the first-order relations) and V20
   (the surface really is level, and the analytic gradient matches a central difference). The oblate
   case reproduces the equatorial bulge *and* the reduced equatorial gravity from the one field.
   *Range of validity:* with `E = eps_tide + 1.5 eps_cent`, the surface closes only for `E < 0.148`,
   because the quadrupole grows as `r^2` while the monopole falls as `1/r`; beyond that it opens, as
   a Roche lobe does. About 0.15 for the tide alone, 0.10 for the spin alone, 0.059 for two equal
   amplitudes. So the twenty per cent first asked for is *not* attainable in this model, which is
   itself the useful result.
2. **The physical renderer.** **Done.** `symbolic/render.py` samples the level surface by
   direction and forward-maps it onto the sky. For modest deformation the body is convex, so that
   map is one-to-one; and since the surface is an equipotential the outward normal *is* the
   gradient, so visibility is the sign of `N . e_obs`, with no z-buffer and no nested root find.
   Each sample carries the *sky-plane* area it projects to, `(N.e_obs) r^2/(n.N) dOmega`, which is
   what our planar visibility integrates. Level radii come from a vectorised bisection on a common
   bracket. Products: `surface_samples`, `sky_image`, and `direct_visibility`, a quadrature, so no
   Gibbs ringing at the limb. Tests include the decisive null case — a zero deformation rendering
   as the uniform disk, flux `pi` and `2 J_1(x)/x`. Brightness is `g^power` with the exponent a
   parameter (D21); limb darkening is absent by construction until Stage 7, and the captions must
   say so.
3. **The two figures.** Both are **done**: `figure_triaxial` and `figure_delta_visibility` in
   `symbolic/figures.py`. The render shows the star down the spin axis and inclined by 60
   degrees, each in *both hemispheres*, through the same sky axes per orientation — so the far
   panel is the near view **rotated by half a turn**, the far observer's offset `(a, b)` being the
   antipodal surface point rather than the mirrored one. A convex body's silhouette is the same
   from either side, so the two panels of a view share an outline and differ only in the shading,
   which is where the gravity darkening lives and so where the triaxiality shows.
   The default is `eps_tide = 0.09`, `eps_cent = 0.03` — a 3:1 pair, which makes the perturbation
   vanish at the tide's perpendicular, putting that semi-axis at exactly the reference radius — at
   orbital `phase = 45` degrees and inclination 60, with principal semi-axes near
   `(1.28, 1.00, 0.962)`. This is essentially the last usable step. *Ceiling:* closure needs
   `eps_tide + 1.5 eps_cent < 0.148`, which at 3:1 means `eps_tide < 0.099`, so 0.09/0.03 is 91
   per cent of the way and the surface is close to pinching along the tide axis. Two consequences
   worth expecting: the ellipsoid approximation is poor enough that the *outlines* should
   disagree, and the projection weight `r^2/(n.N)` grows sharply toward the bulge tip, so the
   renderer's own sampling is coarsest exactly where the model is worst — artefact and model
   error have to be told apart there. Beyond the ceiling the star opens toward the companion,
   which is Roche-lobe overflow rather than a numerical failure.
   `figure_stretched_sphere` shows the actual rendering against the affine prediction — the same
   potential's brightness evaluated on the unit sphere and stretched in the sky plane, samples
   moved to `Lambda x`, their area grown by `|det Lambda|` and their brightness divided by it so
   each carries its flux. Both rows are shown with their **mean removed**, which is D22 in image
   space: the affine model conserves total flux while the render does not (it sets `I ∝ g` with no
   renormalisation), so the two would otherwise differ by a uniform factor of order the area ratio
   that swamps the pattern entirely. What remains is the departure from uniform brightness, which
   is what the model actually claims to reproduce.
   *The phase is not cosmetic.* With the tidal axis *in* the plane of the spin and the line of
   sight, reflection through that plane fixes both quadrupoles, so the image keeps the mirror
   symmetry of an axisymmetric star however triaxial the body is, and the null test has nothing to
   see. Only a phase away from that plane breaks it — which is what "two orbital phases" meant.
   Setting the tide to zero gives the pure-spin case, whose far panel is then an exact mirror and
   whose cut is a pure cosine series: the base case to compare against. Brightness is `g^power` in
   the potential's units, where the undeformed `g` is 1, so all four panels share one colour scale.
   The sampler splats each sample bilinearly over four pixels, a nearest-pixel deposit having
   speckled at the earlier sampling.
4. **Rung 1, the affine model.** **Done.** `symbolic/affine.py` has `body_axes`,
   `projected_stretch` and `affine_visibility`. `Lambda = sqrt(A A^T)` with `A = P diag(a, b, c)`
   and `P` the two rows (east, north): the projection is exact for an ellipsoid, so the only
   approximation is that the deformed body *is* one. `V_affine(u) = V_sphere(Lambda^T u)` carries
   **no flux factor**, the `1/|det Lambda|` that the change of variables puts on the image
   cancelling the Jacobian's `|det Lambda|` — which settles the normalisation this roadmap had
   left open. `V_sphere` is a sphere carrying the *same* brightness pattern, so the gravity
   darkening is carried through and only the shape is modelled. Three tests, the decisive one
   being that the sampled surface's projected area equals `pi |det Lambda|`.
5. **Rung 2**, documented in section 8, implemented only if rung 1's error demands it.

*Gate:* the affine model against the ray-trace over a sweep of `eps`, inclinations and orientations,
with the error *measured* and reported, and the oblate certificate of item 1. The reference for the
certificate is the **direct Fourier sum over the sampled surface**, not an FFT: an FFT imposes the
periodicity of the picture and so Gibbs ringing at the limb, which no refinement removes. The FFT
stays for filling the (u,v) plane in the figures, where low `(u,v)` — where HBT operates — is
trustworthy.

**Stage 7 — limb darkening.** F19–F20. Genuinely new; this is where Clebsch–Gordan / Gaunt
enters. Not to be started before stages 1–5 are closed.

**Stage 8 — the interactive GUI.** D25. A separate module, matplotlib widgets, no new
dependency, and nothing added to the library. Two radio buttons select the four combinations the
Stage 6 figures lay out; sliders cover the five parameters, with the latitude slider constrained
by the closure bound and the position angle pinned as a gauge. Brightness is fixed at von Zeipel
and the analytic row keeps only `l <= 2`, which is the `eps^2` truncation of the higher
multipoles with the `l = 2` coefficient resummed — no further truncation to choose.
Caching separates the cheap orientation sweeps from the expensive deformation changes, and stale
updates are discarded rather than queued.
*Gate:* the compute path is factored out of the widgets and tested headlessly — the four
combinations render at preview settings within a stated time budget, and the closed-form closure
bound agrees with the numerical scan in `level_radii`, which is a certificate worth having
regardless of the GUI.

---

## 7. Open items

1. **`real_sph_harm` labelling.** The manual labels the real basis `(l, m, c)` and `(l, m, s)`
   with `m > 0`, while the code reaches the sine branch as `real_sph_harm(l, -m, ...)`. The two
   agree under D14. D18 has much reduced the stakes: the real basis is no longer in the visibility
   path, so the call sites are the real-basis certificates, the tests and the figures. Worth doing
   for the call sites' sake, not urgent.
2. **Triaxial ellipsoid sources.** Now Stage 6, planned in sections 6 and 8 rather than open as a
   question: the dipole flux weight, the coupling of the families and the two rungs are recorded
   there. What remains undecided about it is the target case (item 5) and the diagnostic (item 4).
3. **The reference direction of the baseline azimuth.** `psi = arctan2(v, u)` is measured from
   `+u` = East towards `+v` = North, whereas a sky position angle is measured from North towards
   East. Two variants have been raised and neither is settled. The first is simply to shift by
   `pi/2`: the packed phase of a complex mode, `(-i)^m exp(i m psi) = exp(i m (psi - pi/2))`, then
   has no stray factor. The second is to adopt the astronomical reference throughout, `psi` as a
   position angle, which achieves the same tidying. Both cost something: the body azimuth `phi` is
   a longitude, and the mode formulas pair `cos(m phi)` with `cos(m psi)`, so re-expressing `psi`
   as a position angle makes the cosine branch contribute `cos(m pi/2 - m PA)`, which is
   `+/- sin(m PA)` for odd `m` and `+/- cos(m PA)` for even `m` --- the branch names would
   alternate meaning with `m`. Deliberately left open.
4. **The form of the diagnostic for triaxiality.** The mirror-symmetry null of section 8 is the
   candidate: the `psi`-antisymmetric part of `|V|^2` about the projected spin axis, zero for an
   oblate star and `O(eps)` for a triaxial one. Alternatives — the phase of the azimuthal structure,
   or a fit of the low-order azimuthal harmonics — have not been compared against it, and the
   comparison should account for what the sparse sampling can actually constrain.
5. **The target case for Stage 6.** A synchronous circular binary is the clean case, the tidal and
   spin quadrupoles being perpendicular. A general oblique, eccentric or non-synchronous system
   gives an arbitrary angle between them and, in the spin frame, a deformation that rotates with the
   phase. To be settled before the figures, since it fixes what they show.
6. **Spectral scope.** Gravity darkening is temperature-dependent, and intensity interferometry can
   observe in more than one band, so a *ratio* of `|V|^2` between bands would isolate the
   temperature variation from the shape. In scope now or later.

### Closed since this list was written

* **`make docs` strictness.** Settled by D6 as amended: `make all` checks before it documents and
  stops if the certificates fail, while `make docs` stays the plain LaTeX path.
* **A bibliography.** `docs/refs.bib` with six entries, natbib in author-year style and `bibtex`
  folded into `make docs`: no undefined citations, no bibtex warnings.
* **`Coeffs` opacity.** Settled as a frozen dataclass with a public, read-only `.blocks` and the
  ordering documented: inspectable at the REPL, and still safe to share.
* **Real versus complex interface.** Settled by D13 and D18 — complex storage, a real parameter
  vector for fitting, and no real branches in the visibility.
* **`special/bessel.py`.** Does not survive; folded into `visylm/visibility.py` (D17).

---

## 8. Non-spherical sources: triaxiality and the deformation field

Approach agreed (D11, D12); not yet scheduled. The analysis below is retained in full because it is
what the approach rests on.

**Parameterisation.** A triaxial ellipsoid is an affine image of the sphere: `P = L n` with
`|n| = 1`, where `L` is the stretch (symmetric positive definite in the body axes). The surface
map is defined on the unit sphere and inherited, which keeps the spherical-harmonic expansion
meaningful.

**What survives.** Three structural features carry over unchanged:

- the visible region in `n` is still a **hemisphere**, with axis the stretched line of sight
  `L^{-1} e_los`, so the domain does not degenerate into an awkward region;
- the kernel is still a **plane wave** in `n`, with a linearly stretched frequency
  `w = L^T u_2D`. This is the part the "stretch u,v" intuition gets right;
- the flux element is still **linear in `n`**:
  `dx_E dx_N = |det L| (n . q) dOmega` with `q = L^{-1} e_los`. For `L = I` this reduces to the
  familiar `cos(theta)`, so the formulation is a strict generalisation of what we already have.

Note the consistency check that comes free: the weight vanishes exactly where the visible region
ends, which is where the limb must be.

**What breaks.**

- The entry-point map is the larger root of a quadratic, hence **nonlinear**, so the projected
  image is a nonlinear warp of the spherical image. There is no exact relation of the form
  `V_ellipsoid(u,v) = V_sphere(M u, M v)`.
- The flux weight is a **dipole** (a pure `l = 1` harmonic) rather than `n_z`. Multiplying by it
  couples `l` to `l +- 1` through the spherical-harmonic product rule (Gaunt / Clebsch-Gordan),
  which **flips the parity** and therefore couples the odd and even families. Expect finite sums
  of spherical Bessels, not a single term.

**Cases where stretching is exact.** Uniform brightness — the silhouette is an ellipse, so the
visibility is the elliptical-disk transform, which is literally a stretched uniform disk. And any
model defined directly on the sky plane rather than on the surface, where the Fourier scaling
theorem applies. The second is a modelling choice worth weighing.

**Related physics.** The surface normal is no longer radial, so limb darkening — a function of the
angle between the normal and the line of sight — changes meaning and must be handled consistently
or the law will be wrong. Gravity darkening is a separate multiplicative factor; for a rotating
star it is axisymmetric about the spin axis, but a tidal term perturbs the surface gravity too, so
the brightness field needs the same generality as the shape field (D11).

**Prior art.** Interferometry of rapidly rotating, oblate, gravity-darkened stars is a mature
area — Domiciano de Souza and collaborators, notably the Achernar work. To be surveyed before
anything is built, since the value added here is closed form and speed, not feasibility.

**Measured, not guessed.** A first-order perturbation leaves an `O(eps^2)` residual: about 1% at
`eps = 0.1` and 4% at `eps = 0.2`, against interferometric precisions often better than 1%. The
usable range is therefore to be *measured* by sweeping `eps` at fixed `(u,v)` and finding where the
residual crosses the accuracy target, rather than assumed.

**Agreed approach** (D11, D12), cheapest first.

1. **Perturbation in a deformation field** about the sphere, carrying both a shape field and a
   brightness field, each expanded in spherical harmonics. Covers oblateness, tides and general
   small asphericity. This is Stage 6, ahead of limb darkening: it is cheap to evaluate pointwise,
   which is what a fitting loop over orientations needs, whereas limb darkening is the harder
   consistency problem.
2. **Reference implementation: direct sky-plane ray-trace.** For each node in the silhouette, solve
   the line–ellipsoid quadratic, take the entry point, evaluate the map there, and accumulate the
   Fourier sum. Chosen over the hemisphere formulation precisely because it shares no algebra with
   the derivation it validates. Converge its grid before drawing conclusions from any comparison.
   Sum the series directly at the chosen `(u,v)` rather than transforming a rendered image: an FFT
   imposes periodicity on a picture whose limb is a discontinuity, and the Gibbs ringing that
   follows does not shrink with the grid, whereas direct summation converges with the sampling.
3. **Exact treatment via the dipole / Gaunt coupling** — the honest generalisation, held in reserve
   in case the measured range from item 1 proves too small to be useful.

**The deformation field, and where it comes from.** The shape and the brightness are not
independent. For a body whose surface is an equipotential, both follow from one potential
perturbation `phi_1`, because on a level surface the normal is along `grad Phi`:

    g_eff = |grad Phi| ,   delta r = -phi_1 / g_0 ,   delta g = d phi_1 / d r ,

to first order in `phi_1`, with `phi_1` defined so that `-grad phi_1` is the extra acceleration.
For uniform rotation `phi_1 = -(1/2) Omega^2 varpi^2`, which gives the equatorial bulge and the
reduced equatorial gravity, both as they should be. Gravity darkening then enters multiplicatively,
`T_eff ∝ g^beta` with `beta = 0.25` for radiative envelopes and near `0.08` for convective ones, so
the brightness perturbation has the same angular dependence as `phi_1`: a separable
`phi_1 = f(r) Y_l^m` gives a brightness perturbation that is also a pure `Y_l^m`. Rotation
contributes `l = 2, m = 0` about the spin axis. The tidal quadrupole's axis is the companion's
direction, generally at an angle to the spin, so *in the spin frame* it carries `m = 0, +-1, +-2`;
only in the body's principal-axis frame is the triaxial figure diagonal, hence `m = 0, +-2` alone.
Which frame matters, since our coefficients live in the spin frame, where `m = +-1` is what a tilt
looks like. Three unequal axes need the two quadrupoles non-parallel; in a synchronous circular
binary they are perpendicular, synchronous rotation putting the companion's direction in the orbital
plane at right angles to the spin, and that is the clean triaxial case. Note also that the
centrifugal term is not a pure quadrupole: its `l = 0` piece renormalises the effective monopole.
That is why truncating the deformation at `l = 2` is physical rather than crude. A smooth,
gravity-darkened, deformed star is then the monopole plus a handful of modes, and its visibility is a
short sum of the mode visibilities of section 7. Caveat: the surface is assumed an equipotential,
which holds for uniform rotation and no meridional circulation; under differential rotation the shape
and the gravity are independent fields.

**The two rungs.** Rung 1 is the **affine model**: the deformed brightness pattern is the linear
image of the undeformed one, so `V(u) = V_sphere(Lambda^T u)` with `Lambda` the stretch. It is exact
for the projected shape of a uniform ellipsoid, one line at the call site, pointwise, and cheap
enough for a million orientations. It is a *model*, not a controlled approximation of a
surface-mapped ellipsoid: its difference from the ray-trace is first order in the deformation, and
that difference is to be *measured*, not assumed. Whether the stretch preserves flux or surface
brightness is left unspecified, the applications being insensitive to it. Rung 2 is the first-order
surface perturbation in the stretched parameterisation: the frequency stretch is a directional
derivative of the `V` we already have, and the only genuinely new object is the transverse part of
the dipole weight, which needs the solid-angle sibling of our radial integrals. Documented now,
implemented only if rung 1's measured error demands it.

**What the application wants.** The measured quantity is `|V|^2`, so the phase is unavailable (D20)
and the target is a *diagnostic* for a sparse set of baselines, not an image. The candidate is a
null test on mirror symmetry: an axisymmetric star's visible-hemisphere image is symmetric under
reflection in the plane containing the line of sight and the spin axis, so `|V|^2` has only cosine
harmonics in the azimuth measured from the projected spin axis. A triaxial figure is not invariant
under that reflection, so the antisymmetric part appears at `O(eps)`. That needs only pairs of
baselines mirrored about an axis, not a filled plane, and a null is easier to measure than an
amplitude. Rotation or orbital phase modulates it, which is the rotation-synthesis idea. For the
difference against a spherical reference, use the *volume-equivalent* sphere, `R = (abc)^(1/3)`:
same volume, so the same flux for uniform brightness, hence `Delta V(0,0) = 0` and only asphericity
shows.

**Cost, and its mitigation.** A numerical silhouette integral is milliseconds per `(u,v)`, which is
hopeless inside a fitting loop (`10^4` data points x `10^6` samples). But the image does not depend
on `(u,v)`: build it once per model evaluation and FFT it to obtain the whole uv plane, then
interpolate to the sampled points. Not as fast as the closed form, which costs `O(number of modes)`
Bessel evaluations per point, but it keeps the numerical route viable as a fallback.

**Binary-specific consequences.** Orbit and spin introduce structure the spherical case lacks.

- **Two independent angles.** The surface map rotates with the spin phase; the tidal bulge is
  oriented by the orbital phase. For non-synchronous rotation these differ, and conflating them
  would be wrong.
- **Orientation is a frame-change problem.** The deformation is naturally defined in the co-rotating
  frame, so presenting it to the observer is exactly the rotation already handled by the Wigner D
  machinery of section 03. That machinery is not only for surface maps.
- **Orbital phase becomes a fitting axis.** The ellipsoidal term varies with orbital phase, and the
  classical photometric ellipsoidal modulation runs at twice the orbital frequency; the same
  harmonic content appears in the visibility.
- **Reflex motion.** The star's wobble about the system barycentre is a sky-plane offset, hence a
  linear phase ramp in `(u,v)`. Trivial to include, but it must not be mistaken for structure in the
  visibility amplitude.

**Prior art for the binary case.** Tidally and rotationally distorted stars in close binaries are
the classical ellipsoidal variables, where this same expansion in spherical harmonics is the
standard framework in the photometric domain. Worth surveying before deriving: the harmonic
conventions and the standard coefficients may already exist in a form we can adopt, and it provides
a well-defined comparison for the photometric limit (`q -> 0`).
