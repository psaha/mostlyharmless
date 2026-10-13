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
| D6 | **`make docs` is not strict for now:** certificates run, failures are reported, the PDF still builds. Tighten once the stages are closed. |
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
│   ├── visibility.py             # F9-F13, F15-F18
│   └── tables.py                 # load / generate / verify the rational table
├── special/
│   ├── __init__.py
│   └── bessel.py                 # thin documented wrappers; the classical anchors
└── symbolic/
    ├── __init__.py
    ├── derivations.py            # one function per certificate, returning a residual
    └── emit.py                   # LaTeX fragments and the data file
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

# --- F14 -----------------------------------------------------------------
zernike_radial(n, m, rho)                   -> ndarray   # R_n^m, built from Jacobi

# --- F9-F18 --------------------------------------------------------------
visibility_mode(l, m, u, v)                 -> ndarray   # complex, one mode
visibility(coeffs, u, v)                    -> ndarray   # full map, by linearity
radial_integral(l, m, q)                    -> ndarray   # the intermediate radial part
zernike_coefficients(l, m)                  -> list[Fraction]
```

Public entry points take `(u, v)`; `q` and `psi` are derived internally (D9) and never appear in a
caller's code. The intermediate `radial_integral(l, m, q)` is the workhorse and is deliberately
public, because it dispatches on the parity of `l + m` and so keeps the odd and even paths
separately visible and separately tested. `visibility_mode` applies the closed-form azimuthal and
phase factors; `visibility` sums over modes.

A `Coeffs` container holds the ragged `(l, m)` coefficient set so that the ordering and the
parity split are not the caller's problem. Kept inspectable (D5).

### `mostlyharmless.special`

```python
j_l(l, x)                                   -> ndarray   # spherical Bessel j_l
J_m(m, x)                                   -> ndarray   # cylindrical Bessel J_m
visibility_uniform_disk(u, v, theta_d)      -> ndarray   # anchor: 2 J_1(x)/x
visibility_limb_darkened(u, v, theta_d, ld_coeffs) -> ndarray   # anchor: F21
```

**Naming collision to settle.** The limb-darkening literature uses `u`, `u_1`, `u_2` for the
coefficients of the law, which clashes with `u` as a baseline coordinate. The signatures above use
`ld_coeffs`, and the manual will need the same distinction stated explicitly.
```

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
| `data/zernike_coefficients.csv` | columns `l, m, k, numerator, denominator` | exact rationals (D5); a few hundred rows for `l <= 10`; committed |
| `docs/generated/coefficients.tex` | the same table for the manual | emitted via `sympy.latex` |
| `docs/generated/verification.tex` | register of certificates: tag, statement, function, rung, scope | emitted from the test list |

All generated files are committed, with a staleness test that regenerates in memory and
fails on mismatch.

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
| V8 | Jacobi–Zernike relation against the explicit series | exact, `n <= L` |
| V9 | Zernike orthogonality | exact rational |
| V10 | Zernike–Hankel identity | exact, `n <= L` |
| V11 | odd collapse residual zero | exact, `l <= L`, both parities of `m` |
| V12 | projection coefficients equal coefficient matching | exact rational |
| V13 | odd closed form against direct quadrature in `q` | numerical, stated tolerance |
| V14 | even Bessel sum against direct quadrature in `q` | numerical, stated tolerance |
| V15 | full-map visibility against a 2-D FFT of a rendered map | numerical |
| V16 | uniform disk reduces to `2 J_1(x)/x` | numerical, tight |

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
    08-limb-darkening.tex   F19-F20
    09-numerics.tex         the table, accuracy, range of validity
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

**Stage 2 — master reduction.** F9–F11, plus V16 as the end-to-end anchor that the
conventions are right.
*Gate:* a uniform disk comes out as `2 J_1(x)/x`.

**Stage 3 — odd family.** F12–F13, V11, V13.
*Gate:* the closed form matches direct quadrature across a grid of `(l, m, q)`.

**Stage 4 — even family.** F14–F17, V8–V10, V12, V14; the emitted table and LaTeX fragment.
*Gate:* the table is exact and self-consistent; the Bessel sums match quadrature.

**Stage 5 — full map.** F18, V15.
*Gate:* FFT agreement on a rendered map.

**Stage 6 — limb darkening.** F19–F20. Genuinely new; this is where Clebsch–Gordan / Gaunt
enters. Not to be started before stages 1–5 are closed.

---

## 7. Still open

1. **`Coeffs` container**: opaque type with accessors, or a plain list of arrays with the
   ordering documented? Opaque is safer; a plain structure is easier to inspect at the REPL.
   Given D5, the plain route is defensible.
2. **Strictness of `make docs`** (D6): revisit once stages 1–5 are closed, and decide between
   a warning banner and a hard failure.
3. **Real versus complex interface.** D8 says the real cos/sin basis for the interface, but the
   coefficient container's natural home is the complex basis. Worth deciding when `Coeffs` is
   designed, since it fixes the public signatures.
4. **Whether `special/bessel.py` survives.** If the Bessel wrappers end up used only inside
   `visylm/visibility.py`, fold them in rather than keeping a one-customer module.
5. **Triaxial ellipsoid sources.** Raised after the first version of this roadmap. Not a small
   extension: the linear-stretch intuition gets the plane-wave argument right
   (`w = L^T u_2D`), but the flux weight becomes a dipole instead of `cos(theta)`, which
   couples the odd and even families and destroys the single-spherical-Bessel result. To be
   scoped as its own stage once the approach is agreed — see section 8.
6. **`real_sph_harm` labelling.** The manual labels the real basis as `(l, m, c)` and
   `(l, m, s)` with `m > 0`, while the code reaches the sine branch as `real_sph_harm(l, -m, ...)`.
   The two agree under D14, but an explicit branch argument would read better at the call site.
   Note that `sph_harm`'s acceptance of negative orders is now *correct* — it implements the
   extension (6) — so V2 has gone back to being a certificate rather than a code-level check.

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
   small asphericity. Scheduled to follow Stage 6, because the dipole / Gaunt coupling it needs is
   introduced there for limb darkening.
2. **Reference implementation: direct sky-plane ray-trace.** For each node in the silhouette, solve
   the line–ellipsoid quadratic, take the entry point, evaluate the map there, and accumulate the
   Fourier sum. Chosen over the hemisphere formulation precisely because it shares no algebra with
   the derivation it validates. Converge its grid before drawing conclusions from any comparison.
3. **Exact treatment via the dipole / Gaunt coupling** — the honest generalisation, held in reserve
   in case the measured range from item 1 proves too small to be useful.

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
