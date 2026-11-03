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
| D8 | **Conventions:** Condon–Shortley phase inside `P_{lm}`; complex harmonics `Y_{lm}` with `m >= 0` only, used in the derivations; maps use the real basis `R_{l,0} = Y_{l0}`, `R_{l,m,c} = sqrt(2) Re Y_{lm}`, `R_{l,m,s} = sqrt(2) Im Y_{lm}`, expanded with real coefficients `c_lm, s_lm`; **negative orders are never formed**, so `\|m\|` is never needed and there is no signed-order convention to get wrong; `u, v` in wavelengths; the `2*pi` inside the exponent; observer frame `(East, North, line of sight)`. |
| D9 | **Notation:** `rho = sin(theta)`; `x_E, x_N` direction cosines; `q = sqrt(u^2+v^2)`; `psi = arctan(v/u)`; `a = 2*pi*q`; `k = (n-m)/2` for the Jacobi degree. `\|m\|` is never used as a symbol: orders are always written signed, and the parity condition is stated as `l+m` odd or even. |
| D10 | **Verification rungs** are labelled on every claim: `exact` (free parameters), `range-limited` (exactly zero for all integer cases in a stated range — a proof of the finite statement), `numerical` (agreement to a stated tolerance). Prefer exact expansion and coefficient comparison to `simplify`, so a pass is a proof of the instance rather than trust in a heuristic. |
| D11 | **Non-spherical sources: perturb in a deformation field.** The surface is `r(n) = R[1 + sum eps_lm Y_lm(n)]`, with an independent brightness perturbation, rather than a single flattening scalar; the `eps_lm` are the fitting parameters. The first-order visibility correction has **two** pieces — the interior integrand change and a **limb boundary term** from the moving silhouette — and both are required. Shape and gravity darkening enter at the same order and must be included together. Validate by the `eps^2` scaling of the residual, not merely by agreement. The reference implementation is a direct sky-plane ray-trace, chosen because it shares no algebra with the perturbation. |
| D12 | **Triaxiality is required** — spinning stars in a binary — so a single flattening parameter will not do. Rotational and tidal quadrupoles about different axes are the natural input; the deformation-field parameterisation of D11 covers this. |
| D13 | **Coefficients are half-complex and negative-order-free.** A real map is stored as `a_l0 = c_l0` and `a_lm = (c_lm - i s_lm)/sqrt(2)` for `m > 0`, one complex array per degree of length `l+1`, equivalently `I = sum_l [a_l0 Y_l0 + sum_{m>0} 2 Re(a_lm Y_lm)]`. The map is real by construction — there is no second set of stored coefficients to disagree with the first — and no unphysical direction for an optimiser to populate. Fitting goes through a flat real vector of exactly `(l+1)^2` numbers. The order is non-negative wherever it occurs, so the magnitude appears nowhere. |
| D14 | **Negative orders are indices, never coefficients.** `Y_{l,-m}` is *defined* as `(-1)^m conj(Y_{lm})`, which coincides with the standard convention, and is used only in the rotation machinery, where the order indexes the representation. Coefficients are formed only for `m >= 0`; where the rotation would call for a negative-order coefficient its value is substituted, so the transform reads `a_lm' = D^l_{m'0} A_l0 + sum_{m=1..l} [ D^l_{m'm} A_lm + (-1)^m D^l_{m',-m} conj(A_lm) ]` with `m' >= 0`. The `m = 0` term appears once, not twice, because `A_l0` is real. (The letters were renamed to match the manual: capital `A` is body-frame, lower case `a` observer-frame, and the visibility always takes the latter. This entry originally wrote `a` for body-frame and `b` for observer-frame, which collided with the rest of the document.) The classical Wigner `d` and `D` therefore stay in the manual, rather than being relegated to the certificates. |
| D15 | **Flux, mean intensity, and normalisation.** With the convention of D8 (no prefactor), `V(0,0)` is the total flux — a linear functional of *all* the coefficients, not of `a_00` alone, because the disk integral weights each harmonic by the foreshortening factor `cos(theta)` rather than being orthogonal as the sphere integral is. Normalising the visibility to unity at the origin therefore means dividing by the instantaneous flux, which for a rotating spotted star is time-dependent and discards the light-curve information; the unnormalised convention is kept. Separately, `a_00` is `sqrt(4*pi)` times the mean intensity over the **sphere** — not the mean over the visible disk, and not the flux. Fixing `a_00` is a choice of overall brightness scale and does not by itself make `V(0,0) = 1` unless the map is constant. Where a flux constraint is wanted, its row is the mode visibilities at `q = 0`, which Stage 2 produces anyway. |
| D16 | **No convenience features the user can add trivially, and `a_00 = 1`.** Two assumptions do the work that parameters, helper methods or wrappers otherwise would: the sphere has **unit radius**, and maps are normalised so that **`a_00 = 1`**, which makes every other coefficient a mode amplitude relative to uniform brightness. Both are documented assumptions rather than code constraints or conveniences such as `normalised()`, because a wrapper reduces the readability of the code for no real gain. The convention is preserved by rotation, the `l = 0` block being untouched by it, so the free coefficients of a map are those with `l >= 1`. The container itself stays general — it still holds arbitrary maps, including the zero map — and it is the modelling convention, not the container, that fixes `a_00`. |
| D17 | **No `special/` package.** Bessel functions come straight from scipy and sympy, with no convention of ours to protect; the classical anchors — the uniform disk, and later the limb-darkened disk — are visibility functions and live in `visylm/visibility.py`. A module with one customer and no convention of its own is folded in rather than kept. |
| D18 | **Negative orders are extended, not computed.** The visibility is extended to negative orders by `V_l^-m = conj(V_{lm})`, the analogue of the harmonic extension `Y_{l,-m} = (-1)^m conj(Y_{lm})`; the sign is absent because the `(-1)^m` of the harmonic extension is cancelled by the `(-1)^m` that reversing the baseline contributes. A stored `a_lm` therefore contributes `a_lm V_{lm} + (-1)^m conj(a_lm) conj(V_{lm}) = 4 pi (-i)^m N_lm R_lm(q) Re(a_lm exp(i m psi))`, which is equation (35) of the manual. There is consequently no need for two real branches, in the derivation or in the code: one mode function suffices, and the real branches become a remark for comparing with the literature. |
| D19 | **Coefficients are computed just in time, and the emitter is separate.** The even-family coefficients are exact rationals computed on demand and memoised, not stored in a table: there is then no maximum `l` to choose, and no data file for the package to depend on. A short program, `symbolic/emit.py`, calls the library and writes the LaTeX fragment the manual inputs, so the only code that knows about LaTeX is the emitter; `symbolic/figures.py` does the same for the figures, importing matplotlib lazily so that it stays an optional extra. `visylm/tables.py` is dropped from the planned layout, nothing having to load a table. The manual's printed table goes to `l <= 3`. |
| D20 | **The application is intensity (HBT) interferometry.** The measured quantity is `|V|^2`, so the phase of `V` is not observable — but it is retained throughout, both for visualisation and because discarding it early would make the intermediate results harder to check, and squared only at the end. Two design consequences: the deliverable for deformed sources is a *diagnostic* evaluable on a sparse set of `(u,v)` points rather than an image, since a fitting loop over orientations makes a numerical transform per model too expensive; and the (u,v) sampling being sparse means the high-frequency behaviour of the transform is not the design driver. |
| D21 | **Brightness is proportional to `g`, for now.** The surface brightness of the deformed star is set by the bolometric von Zeipel law, `T_eff^4 ∝ g` and hence `I ∝ g`, which is the simplest law that is also the standard one for radiative envelopes. The literature favours a weaker dependence — `T_eff ∝ g^beta` with `beta ≈ 0.25` radiative but nearer 0.08 for convective envelopes, so `I ∝ g^(4 beta)` and about `g^0.3` in the convective case, a much flatter law. **Closed as a decision, not a task.** The exponent stays unexposed. In the numeric and fitted rows it would be trivial --- `surface_samples` already takes `power` --- but the first-order row is too approximate to carry another parameter: its brightness is a linearisation of the pattern, so an exponent other than 1 would have to be linearised again and the row would gain a knob whose effect is smaller than its own error. Exposing it for two rows out of three would also invite reading the third as a like-for-like comparison, which it is not. |
| D22 | **The observable is normalised at the origin.** Intensity interferometry measures `|V|^2` in units of the zero-baseline correlation, so what is modelled is `V/V(0,0)`, which is 1 at the origin by definition. Every difference is therefore taken between *normalised* visibilities, so `Delta V(0,0) = 0` by construction, and the absolute flux — which the normalisation removes — is not modelled at all. This also settles the normalisation rung 1 had left open: the affine model is normalised by its own origin, and the unobservable flux difference between the star and its reference is divided out rather than carried. The *squared* difference `|V_tri|^2 - |V_ref|^2` is deliberately not used: the phase is retained, and the difference is taken on `V` itself. Worth recording why the normalisation is needed even against a volume-equivalent reference: the flux is the *projected area*, not the volume, so the two agree only for a sphere, and the deformed star's projected area depends on the viewing direction. |
| D23 | **The azimuth reference is North, through East.** Position angles are measured from the observer's North toward East, the astronomical convention, which fixes both the reference axis and the sense. With `u` along North and `v` along East (D24), `psi = arctan2(v, u)` *is* the baseline's position angle, so no offset appears anywhere between the baseline and the spin's `gamma` (D29). Had the older convention been kept — `u` along East — the consequence would have been `PA = pi/2 - psi`, the other sign corresponding to North-through-West, and the two differ by a *reflection* rather than a rotation, flipping every odd-`m` sign and changing the meaning of the mirror symmetry of section 8. The body's longitude zero remains the body `x` axis, which is already the convention for the harmonic azimuth. The real-basis branch labelling (T2) is a separate convention and is untouched by this. |
| D24 | **The sky's orientation is shared by the ground and the image, departing from the radio convention.** `u` is the **North** component of the baseline and `v` the **East** component, so `V(u,v) = ∫ I(x_E, x_N) exp(-2 pi i (u x_N + v x_E))`. The third axis points from the sky *toward* the observer — into the ground for an observer looking up — so (North, East, toward the observer) is right-handed *given* the sky convention that North is up and East is to the left. The point of the swap is that `psi` becomes the position angle directly (D23), and that the coordinates on the sky and on the ground have the same orientation. The cost is deliberate: radio interferometry takes `u` as East-West and `v` as North-South, so this is a departure, justified only because our `(u,v)` are used throughout as *sky-plane* coordinates rather than as baseline hardware coordinates. Every image is therefore drawn with North vertical and East horizontal increasing leftward, images and `(u,v)` planes alike. *Status: decided, not yet implemented.* The obstacle is that the mode reduction's `psi` is tied to the harmonics' `phi`, which is measured from the body `x` axis; relabelling `(u,v)` alone rotates `psi` by `pi/2` while `phi` stays put, and V15 detects exactly that as a 0.14 residual. Doing it properly means re-referencing the body `x` axis to North as well — so the display maps the body `x` to the vertical — together with the transform and the V15 pairing. |
| D25 | **The interactive GUI is a separate, optional front end.** Matplotlib widgets, since there is no Jupyter here, and nothing added to the library. Two radio buttons: the view is image or `(u,v)`, and the model has **four** settings, chosen so that the shape model, the brightness approximation and the truncation can be varied one at a time --- `numeric` (exact level surface), `first-order` (ellipsoid shape, brightness linearised as in the manual's equation 89), and two fitted orders, `l <= 4` and `l <= 8`, whose coefficients are chosen by least squares in the image plane rather than taken from the brightness pattern. The readings are these: `numeric` against `first-order` is the whole analytic stack at `O(eps)`; `first-order` against the exact pattern at the same order isolates the *brightness* linearisation at `O(eps^2)`; and `l <= 4` against `l <= 8` isolates the *truncation* from the shape model, which the earlier rows conflated. The top row was at `l <= 6` until a visual pass showed that it still carried a visible hard edge at the ceiling in the tilt check, where the brightness falls to zero; at 8 it does not. The `affine` row --- ellipsoid shape with the pattern's coefficients projected rather than fitted --- was in the interface until the fitted row existed and left it under T4, being strictly worse in the same norm; its code is still there pending that item. The spin's position angle is *pinned to zero*, because for a single observation the sky's orientation is a gauge; the model keeps it for real data. The light is **fixed at the von Zeipel law**, `I ∝ g` (D21), with the exponent to be exposed as a slider later; and the analytic row keeps **only `l <= 2`**, discarding `l = 4` and above. Note what that is and is not: it is *not* "first order in `eps`". The `l = 4` content first appears at `eps^2` from `P_{2,2}`, so dropping it *is* the `eps^2` truncation — but the projection keeps the `l = 2` coefficient to all orders, including its own `eps^2` correction, so the retained set is "`l <= 2`, resummed in `eps`". The difference from a true first-order truncation is itself `O(eps^2)`, below the model's error, so the simpler and more accurate choice is `l <= 2`: `SPHERE_L_MAX = 2`, with no further decision to make. The free parameters are `eps_tide`, `eps_spin`, `i`, `lambda_tide` and `phi_tide`. The spin phase `alpha` is deliberately *absent*: this front end is for a spin and tide alone, and with no starspots the phase is degenerate with the tide's longitude — only their difference is observable — so a slider for it would move a gauge. The **latitude slider is constrained** rather than the amplitudes, and in the direction that is easy to get backwards: raising `lambda` *loosens* the bound, since at the equator the centrifugal quadrupole deepens the potential by `1.5 eps_spin` while along the spin it contributes nothing, so the infeasible region is near `lambda = 0` and the limit on the slider is a **lower** one. The eps sliders get upper limits from the same function. Clamping is applied after every change rather than by restricting ranges, because the feasible set is a polygon in the amplitudes rather than a box. Sliders show **degrees**; the conversion lives at the widget boundary, so the model and everything outside the GUI stay in radians. *Implemented so far:* single-threaded redraws, with a **numerical resolution control** of three labelled radio buttons in a column --- coarse, medium, fine. Discrete rather than continuous because every change costs a full redraw. Each level sets the surface samples, the image pixels and the `(u,v)` grid together (10 000/140/17, 30 000/220/23, 150 000/360/31), keeping the sample density per pixel roughly constant; only the samples are physics, so the control trades speed against noise, which is what makes it the way to tell a sampling artefact from a model error. The finest level matches the Stage 6 figures --- 150 000 samples for a rendering, a 31-point grid for a visibility map --- so it is publication quality at the cost of a visible pause, since the numeric `(u,v)` row costs `grid^2 x samples`, about fifty times the coarsest level. Widget references are held on the figure so that they are not collected; and clamping is written **back into the widgets**, so the thumb tracks the state rather than sitting where it was dragged. That last one matters because the clamp often moves a parameter other than the one being dragged --- raising `eps_tide` pushes `lambda` up --- so every slider is refreshed, with a guard flag to stop `Slider.set_val` re-firing its own callback. *Planned, not yet done:* coalescing updates --- numpy cannot be interrupted, so the intent is a worker thread, a latest-values slot holding the newest parameters, and stale results discarded, with the GUI thread only drawing --- and caching the level radii and brightness pattern across orientation changes, since they depend on the deformation alone. |
| D29 | **The third Euler angle is `gamma`, not `PA`.** The rotation is `R_z(gamma) R_y(beta) R_z(alpha)` in the zyz convention, and section 3 already wrote it that way in one sentence while writing `PA` in the next --- three names for one angle (`gamma`, `PA`, `pa`). The manual now uses `gamma` throughout; the code keeps `pa` as the argument name, which names the observable, and the GUI slider says `gamma`. The baseline's position angle stays `psi` and is a different object. |
| D26 | **The linearised level radii are not offered.** Taking the ellipsoid's radii from `delta r = -phi_1` instead of the level condition would remove the last root-find from the analytic rows and make them fully closed form. It is exact only to first order, and the next term is not small where it matters: `delta r = E + 3E^2 + ...` diverges as the closure bound is approached, giving a bulge of 0.144 where the cubic gives 0.366 at the ceiling. The saving is three bisections per parameter set --- some fifty scalar evaluations against the 10^5 the replaced quadrature cost --- so a shortcut that is safe only in a corner and buys nothing measurable is a hazard rather than an option. The brightness is linearised instead, which is exact at first order with no such divergence, and `principal_axes` always solves the level condition. |
| D28 | **Limb darkening is a viewing effect, so it lives in the sampler.** The GUI applies linear limb darkening, `I = I (1 - u_limb (1 - mu))` with `mu = sqrt(1 - rho^2)` the cosine to the line of sight, and every row carries the same law: `render.surface_samples` and `render.sphere_samples` multiply the sampled brightness, `affine.first_order_coefficients` adds it in closed form as a constant plus a dipole about the observer's axis --- `mu = sqrt(4 pi / 3) Y_10`, one coefficient at `m = 0` --- and the fitted rows inherit it from the samples they fit. Both samplers keep only visible points, so `mu > 0` throughout and the law is monotone across the disc. A quadratic law in `rho^2` was tried first and rejected as too drastic at mid-disc: it gives `1 - 0.75 u_limb` there against the linear law's `1 - 0.5 u_limb`. The numeric row's law is a placeholder for what a stellar-atmosphere model would give, so its exact shape is provisional; using the same `rho^2` law everywhere means the rows differ in the shape and the brightness pattern, which is what they are for, and not additionally in the limb. The consequence to keep in mind is that a fit in the *body* frame absorbs the darkening as though it were painted onto the surface, so its coefficients are meaningful only at the orientation they were fitted at --- a limitation of the surrogate's frame, and one of the things the disk-mode route of section 9 would remove. This is the GUI-level half of T1; the analytic visibility still has no limb term. |
| D27 | **The surrogate is fitted in `(x,y)`, not in `(u,v)`.** By Plancherel an image-space fit is an `L2` fit over *all* `(u,v)`, so the answer no longer depends on which baselines were sampled --- which is what left the sparse-`(u,v)` least-squares underdetermined, since the coefficient count at `l <= 4` is comparable to a sparse baseline count. The procedure: unstretch the numerical image, `I'(y) = |det Lambda| I_num(Lambda y)` (the inverse of the manual's eq. 81), then project `I'` onto the *projected* harmonics on the unit disk --- the section 6 objects, with `exp(i m phi)` angular dependence and a Jacobi radial part. Their Gram matrix depends only on the mode count: not on `Lambda` and not on the node, so it is computed and inverted once, and a node costs one quadrature over the samples the numeric model already produced. Three consequences: `A_lm` becomes a linear functional of the image samples and therefore smooth in the parameters, which is what the interpolant needs; the sampling-driven part of the conditioning disappears, leaving only the intrinsic mode near-dependence, now a single `Lambda`-independent condition number checked once; and because the unstretched image should be supported inside the unit disk, the overshoot is a free diagnostic on `Lambda`. *Caveat recorded:* this is not a noise-weighted fit to real baselines. If the surrogate is used inside a fit to actual interferometric data, the `(u,v)` metric with those weights matches the use better; for interpolating the model, the image metric is the more robust choice because it is uniform in baseline rather than in sampling. |
| D30 | **Style: the manual is mathematics, not prose about mathematics.**  (i) Punctuation in a display only where it resolves an ambiguity.  (ii) Concise notation in preference to explanation in words.  (iii) No editorial text --- motivation, judgements about which approximation is good enough, advice --- which belongs in `expos/`, the human-written document. |

---

## 1. Formula inventory

Status: **derived** = worked out and checked during design discussion; **new** = still to do.

| # | formula | status | section |
|---|---|---|---|
| F1 | `Y_{lm}` for `m >= 0`: definition, CS phase; the real basis `R_{l,m,c}`, `R_{l,m,s}` and the map expansion | derived | 02 |
| F2 | Orthonormality | derived | 02 |
| F2b | the `l = 2` zonal corollary, from the rotation law | derived | 03 |
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
| F14 | Zernike `R_{n,m}` from Jacobi; orthogonality on the disk | derived | 06 |
| F15 | Zernike projection coefficients for the even family | derived, table built (V12, V14) | 06 |
| F16 | Zernike–Hankel identity: `int R_{n,m} J_m rho drho = (-1)^((n-m)/2) J_{n+1}(a)/a` | derived | 06 |
| F17 | Even result: finite sum of odd-order ordinary Bessels | derived | 06 |
| F18 | Full-map visibility by linearity over modes | trivial | 07 |
| F19 | Limb-darkening profile and its `m = 0` harmonic coefficients | **new** | 11 |
| F20 | Product identity / Gaunt coupling; limb darkening raised in order | **new** | 11 |
| F21 | Anchors: the uniform disk `2 J_1(x)/x`, and the linearly limb-darkened disk | uniform disk derived (V16); limb-darkened **new** | 04, 11 |
| F22 | Deformation field: a tidal plus centrifugal quadrupole with its `l = 0` piece, the level condition, and the closure bound | derived | 08 |
| F23 | First-order relations `delta r = -phi_1` and `delta g = 4 phi_1` | derived | 08 |
| F24 | First-order brightness and its coefficients in closed form, body and observer frames | derived | 08 |
| F25 | Rung 1: `Lambda = sqrt(A A^T)`, `V = V_sphere(Lambda^T u)`, and the image-space counterpart | derived | 08 |
| F26 | The half-complex pairing of the map and of the visibility, and the closed-form per-mode contribution | derived | 10 |
| F27 | The disk-mode surrogate: `V ~ sum_M c_M Zhat_M(Lambda^T u)`, with the projection for `c_M` | proposed, not implemented | 09 |
| F28 | The fitted rung: the disk harmonics `Phi_lm = Y_lm(n)(n.e)` with `n` the lifted direction, their Gram matrix, the unstretched-image fit and its two-way residual | derived and implemented | 08 |

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
│   ├── visibility.py             # F9-F18, plus the classical anchor
│   └── deformation.py            # the potential, its level surface, the gravity there
└── symbolic/
    ├── __init__.py
    ├── derivations.py            # one function per certificate, returning a residual
    ├── render.py                 # the reference ray renderer of Stage 6
    ├── affine.py                 # rung 1, and the closed-form first-order row
    ├── fitted.py                 # rung 1 with fitted coefficients: disk harmonics, the fit
    ├── emit.py                   # the LaTeX fragment; run as a module
    ├── maps.py                   # the renderer: sky_map, sky_map_batch
    └── gui.py                    # the interactive model of D25; imports only the library
```

### `mostlyharmless.visylm`

```python
# --- F1-F3 ---------------------------------------------------------------
sph_harm(l, m, theta, phi)                  -> complex   # CS phase, theta polar
real_sph_harm(l, m, theta, phi)             -> float     # cos/sin basis
sph_harm_addition(l, th1, ph1, th2, ph2)    -> complex   # theorem LHS, for the checks
normalisation(l, m) / associated_legendre(l, m, x)       # the pieces, public for the checks
m_values(l) / factorial_ratio(l, m)

# --- F4-F8 ---------------------------------------------------------------
wigner_d(l, beta)                           -> ndarray   # d^l, index [m'+l, m+l]
wigner_D(l, alpha, beta, gamma)             -> ndarray   # full D
rotation_matrix(alpha, beta, pa)            -> ndarray   # R = R_z(PA) R_y(beta) R_z(alpha)
rotate_coefficients(coeffs, alpha, beta, pa) -> Coeffs   # body -> observer

# --- F14-F15 -------------------------------------------------------------
zernike_radial(n, m, rho)                   -> ndarray   # R_{n,m}, numeric, via Jacobi
zernike_polynomial(n, m)                    -> sympy expr # R_{n,m}, exactly
mode_polynomial(l, m)                       -> sympy expr # P_{lm}(sqrt(1-rho^2)) = rho^m A(rho^2)
zernike_coefficients(l, m)                  -> tuple     # even-family c_k, exact rationals

# --- F9-F18 --------------------------------------------------------------
visibility_mode(l, m, u, v)                 -> ndarray   # complex, Y_{lm}, m >= 0, dispatches on parity
visibility(coeffs, u, v)                    -> ndarray   # full map, by linearity
radial_integral(l, m, q, n_nodes=None)      -> ndarray   # the quadrature reference
visibility_uniform_disk(u, v, theta_d)      -> ndarray   # the anchor, 2 J_1(x)/x

# --- the deformation field (Stage 6) -------------------------------------
Deformation(tide, tide_direction, centrifugal, spin_direction)
    .axes()                                  # (tide axis, spin axis)
    .potential(x) / .gradient(x)             # Phi and grad Phi, analytic
    .level_radius(direction)                 # the inner crossing of Phi = -1
    .level_radii(directions)                 # vectorised, on a common bracket
    .surface_gravity(direction)              # |grad Phi| there
    .axis_lengths() / .volume() / .equivalent_radius()
    .first_order_radii(direction)            # for comparison with the solver
```

Public entry points take `(u, v)`; `q` and `psi` are derived internally (D9) and never appear in a
caller's code.

The parity dispatch lives in `visibility_mode`, which is why the odd and even closed forms stay
separately visible and separately tested. `radial_integral` is the **quadrature reference** those
closed forms are checked against, not a production path: nothing in the library calls it, and the
certificates do.

A `Coeffs` container holds the ragged `(l, m)` coefficient set so that the ordering and the parity
split are not the caller's problem. Kept inspectable (D5).

**Naming collision to settle.** The limb-darkening literature uses `u`, `u_1`, `u_2` for the
coefficients of the law, which clashes with `u` as a baseline coordinate. The manual will need that
distinction stated explicitly when limb darkening (T1) lands, and the anchor function will have to
carry it.

The classical anchor lives in `visylm/visibility.py` rather than in a special-functions module.
The Bessel functions come straight from scipy and sympy, and the anchor is a visibility function, so
a separate module would have a single customer and no convention of its own (D17). The
limb-darkened anchor is F21 and is not written yet.

### `mostlyharmless.symbolic`

```python
# --- render.py, the reference ray renderer of Stage 6 --------------------
sphere_directions(n_samples) / sky_basis(line_of_sight)
surface_samples(field, eye, n_samples=..., power=...)   -> points, brightness, weights
sphere_samples(field, eye, ..., linear=False)           -> the same, on the reference sphere
first_order_brightness(field, directions, power) / sphere_brightness(field, power, ...)
sky_image(points, brightness, weights, east, north, half_width, n)  -> array
direct_visibility(points, brightness, weights, east, north, u, v)   -> complex

# --- affine.py, rung 1 and the closed-form row ---------------------------
principal_axes(field)                       -> (frame, radii)
projected_stretch(field, eye)               -> ndarray   # 2x2 Lambda, sky order (east, north)
sphere_coefficients(field, l_max=..., power=..., n_samples=..., eye=..., linear=False) -> Coeffs
first_order_coefficients(field, eye)        -> Coeffs    # the closed form of section 8.7
affine_visibility(sphere_visibility, stretch, u, v) -> ndarray

# --- fitted.py, the fitted rung (manual section 8.8) --------------------
disk_harmonic(l, m, east, north)            -> ndarray   # Phi_lm, the projection onto the disk
design(l_max=4) / gram(l_max=4) / block_condition(l_max=4)
fit_image(samples, stretch, east, north, l_max=4) -> Fit  # coefficients, residual, conditioning
residual_of(coefficients, samples, stretch, east, north, l_max=4) -> (residual, norm)
model_image(coefficients, stretch, half_width=..., n=...)  -> ndarray  # [north, east]

# --- derivations.py -----------------------------------------------------
check_*                                     -> residual  # thirty functions; tags in section 4

# --- emit.py, run as a module -------------------------------------------
coefficients_table(l_max=...) / write(directory=None)
sky_map / sky_map_batch       # the renderer the panels use

# --- gui.py, D25 --------------------------------------------------------
Parameters [field]                          # the five sliders
closure_depth / largest_spin / smallest_latitude / clamped
compute_image(parameters, model=..., ...) / compute_visibility(parameters, model=..., ...)
build(parameters) / main()
```

`emit` and `figures` are run as modules by the `data` target --- `python -m
mostlyharmless.symbolic.emit` and `... .figures` --- and neither takes arguments: the defaults *are*
the committed artifacts, which is what keeps the staleness test meaningful.

---

## 3. Data artifacts

| artifact | contents | notes |
|---|---|---|
| `docs/generated/coefficients.tex` | the even-family coefficients for the manual | emitted by `symbolic/emit.py` from the just-in-time coefficients (D19); printed to `l <= 3` |

There is no `data/` directory and no coefficient data file. Coefficients are computed just in
time in exact rational arithmetic (D19), so nothing has to read a table, and the emitters are the
only code that knows about LaTeX or matplotlib.

All generated files are committed. The LaTeX fragments have a staleness test that regenerates them
in memory and fails on mismatch. The figures do not: byte-comparing matplotlib output across
versions would fail for reasons that mean nothing, so they are covered by a rendering test at small
size instead.

---

## 4. Certificates

One tag, one claim, one citation in the text, and one register row per check. A tag may cover
a small group when the checks share one claim --- V1 is three rows --- but no check should be
left outside every tag. The table below lists tags, not functions.

| tag | claim | rung |
|---|---|---|
| V1 | `Y_{lm}` orthonormal on the exact quadrature grid | exact, `l <= L` |
| V2 | harmonic extension `Y_{l,-m} = (-1)^m conj(Y_{lm})`, and that the implementation reproduces it | exact |
| V3 | the `l = 2` zonal corollary, against the general `sph_harm_addition` | numerical: 40-digit SymPy, plus float at higher `l` (the symbolic route tangles in the sign of `sin theta`) |
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
| V21 | the closed-form first-order harmonics against the quadrature that projects the same pattern | numerical, residual `~1e-5` |
| V22 | the Zernike transform identity, against the library's `R_{n,m}`, by Gauss-Legendre | numerical, `n <= 8`, 400 nodes |
| V23 | the mode at opposite baselines, `V_{lm}(-u,-v) = (-1)^m V_{lm}(u,v)` | algebraic, `1e-12` |
| V24 | the first-order coefficients in the body and observer frames, through the rotation law | algebraic, `1e-12` |
| V25 | the real basis is orthonormal, and matches the complex one up to the phase | numerical, `l <= 4` |
| V26 | the fit returns a synthetic image's own coefficients | algebraic, `1e-12` |
| V27 | the Gram matrix: analytic entry `1/4`, block structure, midpoint rule | numerical, `5e-3` |
| V28 | the row at `eps = 0` is the classical limb-darkened disk, `V(q)/V(0)` | numerical |

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
    08-triaxial.tex         the deformation field, the affine model, rung 1, the closed form
    09-disk-modes.tex       (retired) a disk-mode surrogate, and a possible route to limb darkening
    10-models.tex           the two interactive models, and the geometry checks
    11-limb-darkening.tex   F19-F20
    12-numerics.tex         accuracy, range of validity
  generated/                emitted fragments
  refs.bib
```

Sections 01--10 exist; 11 and 12 are planned. Every derivation section opens with a conventions
preamble pulled from section 01, and each checked step carries its tag. Section 05 is the template for the compressed lemma style:
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
that is discontinuous *on the sphere* — a coastline, the edge of a dark patch — has a slowly
converging expansion, and the reference image and the harmonic sum will differ near the edge at
any finite `l_max`. The answer taken here is the taper rather than a convergence study, since the
truncation error is visible.

**Stage 6 — triaxial deformation.** Section 8's method, in this order:

1. **The potential and its level surface.** **Done.** `visylm/deformation.py` carries `phi_1` = a
   tidal quadrupole plus a centrifugal one, both `l = 2`, the latter with an `l = 0` piece that
   renormalises the effective monopole. `level_radius` finds the *inner* crossing of `Phi = -1` and
   `surface_gravity` takes `|grad Phi|` there. Certificates V19 (the first-order relations) and V20
   (the surface really is level, and the analytic gradient matches a central difference). The oblate
   case reproduces the equatorial bulge *and* the reduced equatorial gravity from the one field.
   *Range of validity:* with `E = eps_tide + 1.5 eps_spin`, the surface closes only for `E < 0.148`,
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
   parameter (D21); limb darkening is absent by construction until T1, and the captions must
   say so.
3. **The two figures.** Both are **done**: `figure_triaxial` and `figure_delta_visibility` in
   `symbolic/maps.py`. The render shows the star down the spin axis and inclined by 60
   degrees, each in *both hemispheres*, through the same sky axes per orientation — so the far
   panel is the near view **rotated by half a turn**, the far observer's offset `(a, b)` being the
   antipodal surface point rather than the mirrored one. A convex body's silhouette is the same
   from either side, so the two panels of a view share an outline and differ only in the shading,
   which is where the gravity darkening lives and so where the triaxiality shows.
   The default is `eps_tide = 0.096`, `eps_spin = 0.032` — a 3:1 pair, which makes the perturbation
   vanish at the tide's perpendicular, putting that semi-axis at exactly the reference radius — at
   orbital `phase = 45` degrees and inclination 60. This is essentially the last usable step.
   *Ceiling:* closure needs `eps_tide + 1.5 eps_spin < 0.148`, which at 3:1 means
   `eps_tide < 0.099`, so 0.096/0.032 is 97 per cent of the way and the surface is close to
   pinching along the tide axis. Two consequences
   worth expecting: the ellipsoid approximation is poor enough that the *outlines* should
   disagree, and the projection weight `r^2/(n.N)` grows sharply toward the bulge tip, so the
   renderer's own sampling is coarsest exactly where the model is worst — artefact and model
   error have to be told apart there. Beyond the ceiling the star opens toward the companion,
   which is Roche-lobe overflow rather than a numerical failure.
   Two further figures, `stretched-sphere.pdf` and `visibility-maps.pdf`, were written and then
   **dropped**: both compared the rendering against the affine prediction, which the GUI now shows
   interactively in its image pane, and neither was ever included in the manual. The comparison
   they made survives as a *measurement* --- the affine residual test prints the worst
   `|V_affine - V_star|` at the working amplitudes and at the ceiling --- which is where I2's
   two numbers come from. One point they made is worth keeping in view: the two must be compared
   with their **mean removed**, which is D22 in image space, because the affine model conserves
   total flux while the render does not (it sets `I ∝ g` with no renormalisation), so otherwise
   they differ by a uniform factor of order the area ratio that swamps the pattern entirely.
   *The phase is not cosmetic.* With the tidal axis *in* the plane of the spin and the line of
   sight, reflection through that plane fixes both quadrupoles, so the image keeps the mirror
   symmetry of an axisymmetric star however triaxial the body is, and the null test has nothing to
   see. Only a phase away from that plane breaks it — which is what "two orbital phases" meant.
   Setting the tide to zero gives the pure-spin case, whose far panel is then an exact mirror and
   whose cut is a pure cosine series: the base case to compare against. Brightness is `g^power` in
   the potential's units, where the undeformed `g` is 1, so all four panels share one colour scale.
   The sampler splats each sample bilinearly over four pixels, a nearest-pixel deposit having
   speckled at the earlier sampling.
4. **Rung 1, the affine model.** **Done.** `symbolic/affine.py` has `principal_axes`,
   `projected_stretch`, `sphere_coefficients`, `first_order_coefficients` and `affine_visibility`. `Lambda = sqrt(A A^T)` with `A = P diag(a, b, c)`
   and `P` the two rows (east, north): the projection is exact for an ellipsoid, so the only
   approximation is that the deformed body *is* one. `V_affine(u) = V_sphere(Lambda^T u)` carries
   **no flux factor**, the `1/|det Lambda|` that the change of variables puts on the image
   cancelling the Jacobian's `|det Lambda|` — which settles the normalisation this roadmap had
   left open. `V_sphere` is a sphere carrying the *same* brightness pattern, so the gravity
   darkening is carried through and only the shape is modelled. Three tests, the decisive one
   being that the sampled surface's projected area equals `pi |det Lambda|`.
5. **Rung 2**, documented in section 8, implemented only if rung 1's error demands it.

*Gate:* the affine model against the ray-trace over a sweep of `eps`, inclinations and orientations,
with the error *measured* and reported, and the oblate certificate of its own item 1. The reference
for the
certificate is the **direct Fourier sum over the sampled surface**, not an FFT: an FFT imposes the
periodicity of the picture and so Gibbs ringing at the limb, which no refinement removes. The FFT
stays for filling the (u,v) plane in the figures, where low `(u,v)` — where HBT operates — is
trustworthy.

**Stage 7 — the interactive GUI.** D25. **Done.** A separate module, matplotlib widgets only, no
new dependency, and nothing added to the library: `compute_image` and `compute_visibility` are the
compute path, `build` is the interface, and `Parameters` holds the five sliders, in degrees. Three
models are offered --- numerical, affine and first-order --- over either the image or a sparse
`(u,v)` set, at three resolutions. The latitude slider is clamped by the closure bound and the
position angle is a gauge, so neither can be set to something the model cannot render; brightness
is von Zeipel, and the analytic rows keep only `l <= 2`.

The interface is deliberately not described further here. D25 records what was built and the
layout is obvious from running it (`make gui`); prose that restates widget state is prose that goes
stale, which is how this section came to describe a version with two radio buttons.

*Gate:* the compute path is tested headlessly by eleven tests, including that the closed-form
closure bound agrees with the numerical scan in `level_radii` --- a certificate worth having
regardless of the GUI --- and that the linearised-brightness row differs from the affine one only
at `O(eps^2)`.

---

## 7. Work outstanding

Three lists, by kind: what will be done and whose design is settled; what is deliberately undecided;
and what may or may not be pursued. Entries are labelled rather than numbered --- `T` to be
implemented, `O` for an open decision, `I` for an idea --- so that inserting one cannot break a
cross-reference elsewhere.

### Yet to be implemented

**T9 — Tikhonov regularisation for the fitted rows.** The fitted rows ring at the rim, where
the brightness falls to zero, and near the cusp at the equator; with the grayscale rendering the
ripples are the most conspicuous thing in the panel. The truncation is real --- the textbook
power-law decay of section 7, and `earth.taper` is the same remedy already applied to real data ---
so this is about presentation, not about the model being wrong.

The classical remedy is tapering the coefficients with Lanczos σ-factors,
`sin(pi l/L)/(pi l/L)`, which Arfken and Weber introduce for exactly this in Fourier series,
with Fejér (`1 - l/L`) and raised-cosine as the other members of the family: all are positive
kernels, so they convolve the image with something smooth and overshoot cannot occur. But that
theory assumes an orthonormal series, and our basis is orthonormal on the *sphere*, not on the
*disc* --- the fit is a least-squares projection onto a non-orthogonal span. Tapering would
still damp the ripples and blur the edge, but without the guarantee. Hence Tikhonov: it is
data-aware and degrades gracefully instead of blurring everything equally.

The code is small, because the fit already solves normal equations: replace
`gram(l_max)` by `gram(l_max) + lam * penalty`, with a diagonal `P` rising with `l` (say
`(l/l_max)^(2p)`, or an indicator above a threshold), and `lam` scaled to the trace of the Gram
so that a slider from 0 to 1 means a fraction of the data term rather than a number with units.
Then a `lam` argument through `fit_image`, a slider, and the certificate untouched: at
`lam = 0` both of its residuals reduce to what they are now.

One property to state rather than discover: the residual the row prints is the *data* residual,
`norm - rhs . values`, which does not include the penalty. Regularising therefore *raises* it,
correctly --- the fit stops being optimal in its own norm --- so the row becomes "cleaner
picture, worse number", and the number is the honest one.

**T17 — keep the manual and the human-written exposition in step.** The two documents have
different jobs and are now diverging by design in one respect and by accident in another.

*By design:* the exposition, `expos/`, is a paper --- "AI sans Slop" --- whose order is motivation,
background, machinery, results, and whose notation is abbreviated: $\epsilon_\text{s}$,
$\epsilon_\text{t}$, $\hat{\mathbf e}_\text{s}$, $\hat{\mathbf e}_\text{t}$ where the manual writes
$\epsilon_\text{spin}$, $\epsilon_\text{tide}$, $\hat{\mathbf e}_\text{spin}$,
$\hat{\mathbf e}_\text{tide}$. That abbreviation is right for a journal and wrong for the manual,
whose job is to be unambiguous rather than compact, so the difference should stay and be written
down rather than discovered.

*By accident:* the exposition's last three sections are placeholders of one line each --- "Introduce
$\Lambda$", "SH coefficients", "Two versions of $\mu$" --- and those three are exactly the material
the manual now carries in finished form (the three-step stretch, the hemispheric basis, the two
cosines). The risk is not that the exposition is unwritten; it is that when it is written it will be
written from an older reading of the manual, which has changed substantially.

Three things worth doing. A short notation map --- manual symbol, exposition symbol, and the one or
two cases where the object itself differs rather than the label --- kept where both documents can
see it, not in either one's prose. A rule that the exposition quotes the manual's results by section
and equation rather than restating them, so a change in the manual shows up as a dangling reference
rather than as a quiet contradiction. And a pass over the exposition's sections 1--5 for statements
that the manual has since overtaken: the gauge of section 8.2, the limb's two cosines, and the
retirement of the affine row all postdate most of what is written there.

### Closed since this list was written

* **The five-part document structure (T11), superseded rather than done.** The parts were: spherical
  harmonics and Wigner D; the $V_{lm}$ structure; simplifying $X_{lm}$; spin, tide and limb
  darkening; and the models. They are a description of the material, not a good order for a document
  whose job is to be checkable in place --- the manual's ordering, where each certificate sits under
  the object it certifies, is the one that serves that job, and reordering the sections to match the
  parts would have put the certificates at a distance from what they certify.

  The narrative order the item was reaching for belongs to the exposition, `expos/`, where it now
  exists: motivation, background, machinery, physics, results. What survives of the item is the
  coordination between the two documents, which is T17.

* **Rendering and projection in the stored half-complex basis (T13).** The last two real-basis
  call sites are gone. `maps.sky_map_batch` -- and `sky_map`, which calls it -- now sums the stored
  orders with the $2\,\mathrm{Re}$ pairing, $2N_{lm}P_{lm}[\mathrm{Re}(a)\cos m\varphi -
  \mathrm{Im}(a)\sin m\varphi]$, and `earth.channel` projects onto $\overline{Y_{lm}}$ with the
  $1/4\pi$ that makes the projection an average rather than an integral.

  Two things turned out simpler than the entry expected. The emitted `earth_levels.py` needed no
  regeneration: the old path divided by $4\pi$ and by $\sqrt2$ against $\sqrt2\,\mathrm{Re}\,Y$
  entries, so the new $a_{lm} = \frac{1}{4\pi}\sum w f\overline{Y}$ gives the same numbers, and the
  emitter already wrote complex blocks. And `test_maps` is *stronger* than the entry feared, not
  weaker: its per-mode reference is still written with `real_sph_harm`, so the comparison is now one
  complex sum against one real sum rather than two complex sums.

  `real_sph_harm`, `to_real` and `from_real` all survive, with narrower jobs: the first is V25's
  subject and the independent reference in `test_maps`, and the last two are how `to_vector` and
  `from_vector` work, which is the fitting interface in starry's order. Section 2.2 is now purely
  the translation to other conventions, which is what its opening sentence always claimed it was
  for.

* **The spot model retired (T10).** The code went first: the RGB map GUI kept at tag
  `rgb-mapgui`, the renderer and cap tests moved to `test_maps.py`, the view sweep to the
  surviving GUI's tests, and the three-channel data deleted because `make earth` regenerates it.
  The manual was rewritten around the survivor in the same pass --- section 10.2's subject is the
  Earth map, its sliders the three angles a fixed map leaves --- so the "still describes the
  retired GUI" sentence in this item was stale by the time the item closed.

  What survives of the spot is one retrospective paragraph, "The retired starspot", saying why the
  Earth replaced it: real data with the same sharp edge, structure at every degree rather than one
  analytic profile, and a taper that answers the ringing a discontinuity provokes. That is history
  rather than mathematics, so it is a candidate for trimming under D30 rather than a requirement
  of the model.

* **The affine model retired (T4).** Its radio button had gone long before; what remained was the
  code path, and it is now gone. `MODEL_TITLES` no longer carries an `affine` entry, the fallback
  branch of `compute_visibility` raises on an unknown model name instead of quietly building one,
  and `SPHERE_L_MAX`, which only that branch used, went with it.

  Three things survive deliberately. `affine.sphere_coefficients` stays as the quadrature
  reference for V21, so the certificate still compares the closed form against an independent
  projection rather than against itself. `affine.affine_visibility` stays, since the first-order
  and fitted rows both reach their visibility through it. And the first-order row stays in the
  interface, being closed form with no fit in it at all.

  The tests that pinned the retired row now use survivors: the finite-panel sweep covers
  `numeric` and `first-order`, the analytic-against-numeric comparison uses `first-order`, the
  order test takes the fitted row as its exact comparison, and the imaginary-part sign test drops
  its affine column. `test_the_fitted_row_beats_the_projected_one_in_its_own_norm` keeps
  `sphere_coefficients`, because there it is a projection to compare against rather than a model
  claim --- the property, that the fit is optimal in the image norm, is why the row is a diagnostic
  rather than another approximation.

* **`real_sph_harm` labelling (T2), closed as documented rather than unified.** The manual labels the
  real basis `(l, m, c)` and `(l, m, s)` with $m>0$; the code reaches the sine branch as
  `real_sph_harm(l, -m, ...)`; and D14 records the correspondence along with the rule that negative
  orders are indices and never coefficients. Once the real basis left the main line for section
  2.2, the call sites the item was written for numbered two, both documented.

* **The map GUI (T5).** `symbolic/map_gui.py`, `make mapgui`, with the scaffolding extracted into
  `symbolic/gui_core.py` and used by both GUIs. It replaced two of the manual's figures, carrying
  both panels rather than the figures' four: the spin-frame against observer-frame distinction is
  the inclination slider. Having a colatitude is what makes the rotation phase mean anything ---
  the cap is zonal, so a rotation about its own axis does nothing, which the first version's slider
  discovered by doing nothing and the tests now pin as a property. `gui_core` holds the coalescing
  timer, the colour-scale rule, the East-leftward flip, the colourbar layout, the slider factory
  and the widget-retention line, each written twice before, with the pure parts tested. Both GUIs
  carry a preset box of named parameter sets, labelled by what should be seen.

* **The manual's figures retired (T6).** With T5 in place the three figures went, with the
  `figure_*` functions, `write()`, the `FIGURES` list, the figure half of the `data` target, three
  rendering tests and the `\includegraphics` calls. The convergence figure's number became a
  certificate, since retired because its claim is textbook. The manual is figure-free, and section
  1's "North up, East to the left" is a statement about the display rather than about a picture.

* **The in-place certificate format (T7), settled on 6.2 and 6.3.** Those two are the pattern:
  each certificate its own item, a $\checkmark$ in place of the word "Certificate", the Python
  function's name where the tag used to be, the claim and scope in the sentence, and the tag and
  rung on a second line at the item's margin. The rest of the manual is still in three styles ---
  certification folded into the prose, scope in the same sentence, or register tables at the ends
  of sections 8, 9 and 10 --- so this closed as the format decision with the per-section conversion
  left to do rather than as a finished conversion.

* **All-subscript indices (T8).** `R_{n,m}`, `V_{lm}`, `Y_{l,-m}`, `R_{l,m,c}`, `R_{l,m,s}`
  throughout, so orders sit on the same level as degrees. Powers and squares untouched, since those
  superscripts are exponents; the Wigner matrices keep $D^l_{m'm}$ and $d^l_{m'm}$ because the
  degree is a superscript in every convention that uses them; Jacobi parameters stay a parenthesis
  pair. No certificate covers a notation change, so the check was that the PDF builds and that the
  sections reconciling our conventions with published ones say which form we use.

* **The limb-darkened disk, as an anchor for the closed-form row (T1).** The route was settled long
  before the item closed, and the manual now carries it: the limb enters as a dipole about the line
  of sight, applied in the observer frame, so no coefficient set is rotated after the fact; the
  azimuthal variation of $\mu$ on a deformed star enters only through the dropped $O(u\epsilon)$
  terms; and exactness for a deformed star is the numerical row's job, per sample.

  The classical case is proved rather than quoted. Section 8.8 derives
  $\frac{V(q)}{V(0)} = 2[(1-u)J_1(x)/x + u\,j_1(x)/x]/[(1-u)+\tfrac23u]$ from the standard
  $(1-\rho^2)^{1/2}J_0$ integral, notes that its two terms are the $l = 0$ and $l = 1$ rows at
  $\epsilon = 0$, and that the $u = 0$ limit is the $2J_1(x)/x$ anchor of section 4. Certificate
  V28 evaluates the row at $\epsilon = 0$ against that expression, over a range of $q$ and against
  scipy's Bessel functions rather than ours, with the value at the origin checked first so the
  ratio has a fixed normalisation.

* **The mean-radius gauge for the deformation (T14), implemented in both text and code.** The level
  value is $a = 1 - \epsilon_\text{spin}$, the level surface sitting at $-1/a$, so the first-order
  displacement is the pure quadrupole $E - \langle E\rangle$ and the $l = 0$ part of the perturbation
  distorts nothing. Three consequences, all in the text: the first-order surface is an ellipsoid
  whose semi-axes sum to the mean radius; $\det\Lambda = 1 + O(\epsilon^{2})$ follows from the
  displacement being traceless, so the flux factor is dropped rather than cancelled; and the
  brightness keeps its $l = 0$ part, $g = 1 - 2\epsilon_\text{spin} - 4q$, which is physical and
  cannot be gauged away.

  The code follows: `Deformation.level_value` returns $-1/(1-\epsilon_\text{spin})$,
  `first_order_radii` returns the quadrupole displacement and the gauged $\delta g$, the first-order
  $a_{00}$ is $1 - 2\epsilon_\text{spin}$, and the two certificates that pin them were restated.
  Two guard tests were added: the angular mean of the level radii is 1 to $O(\epsilon^{2})$, and
  $\det\Lambda = 1 + O(\epsilon^{2})$. One bug was found and fixed on the way --- the closure test
  in `level_radii` took its bracket from the gauged level rather than from the depth, which widened
  the bracket past the outer root and reported an open surface.

* **One account of limb darkening, and one meaning of $\mu$ per row (T16).** The brightness model
  moved to the numerical row, where it belongs: section 8.4 declares $I(x,y)$ general --- a
  stellar-atmosphere model with no simple form is what the fit actually expects --- and gives the
  instance every figure uses, (86), gravity times a linear limb factor written with $z/r$ rather
  than $\mu$. Equation (85) is then visibly its $u_\text{limb} = 0$ case.

  Section 8.8 owns the limb and nothing else does. Its opening no longer says the limb is absent by
  construction, which contradicted the machinery below it; it says the numerical row has applied it
  sample by sample since (86), and that what remains is the projection onto the disk basis. Both
  cosines are defined there and only there: the surface's $z/r$ in (103) and the model rows'
  $\sqrt{1-|\Lambda^{-1}(x,y)|^{2}}$ in (104), with the agreement to first order and the counting of
  the difference. The section then says which rows carry the factor exactly --- the closed-form one,
  in its coefficients --- and which absorb it, the fitted row, whose basis $\Yup_{lm}$ is limb-free.

  The fitted row now refers to section 8.4 for its image instead of restating it, and (115) is
  written with $\tilde I$ and $z/r$. The appendix's duplicate of the $\mu$ definition is the one
  piece of the item left standing.

* **The translation to Dholakia & Pope's conventions is written, and needs no table.**
  Section 4.1 carries the basis: $\Yup_{lm}$ is their $\tilde y_{\mathrm{HSH}}$ for $l+m$ even and
  their $\tilde y_{\mathrm{CHSH}}$ for $l+m$ odd, and (38) gives both cases at once, with the
  normalisation $N_{lm}$ against their unnormalised functions and their index spelled out as
  $\tfrac12(l(l+2)\mp m)$. Section 4.2 carries the transform: the same sign and normalisation,
  with $u_{\mathrm{DP}} = 2\pi u$, and the letter collisions $(\rho,\phi) = (q,\psi)$,
  $\theta = \varphi$. The $l,m \to j$ direction is what the comparison needs, and (38) states it
  in full; the inverse is not needed and is not written.

  One dependence remains stated in the text rather than resolved: their paper does not say where
  $\theta$ is measured from, so (38) holds for a shared azimuth origin, and a different one would
  rotate every $m$ by $e^{im\delta}$ --- a convention, not a physical difference.

* **The negative order is plain conjugation, and the convention now lives in one function.**
  The extension is `Y_{l,-m} = conj(Y_lm)`, with no phase, and the `(-1)^m` has moved to the
  visibility, `V_{l,-m} = (-1)^m conj(V_lm)`. The two arguments for it are the ones that should
  have decided it the first time: `m -> -m` is `phi -> -phi`, so conjugation is what one expects
  and nothing else; and the parity asymmetry belongs on the visibility, whose pairing equation
  already splits odd from even `m`, rather than on the map, whose expansion does not.

  (9) and (13) are unchanged, the real-branch relations (30) are unchanged, and the rotation law
  (24) is unchanged in form --- the sign is conserved, not destroyed, and it lands where the
  pairing needs it. What changed is the definitions (7), (8), the Legendre negative-order
  relation (now `P_{l,-m} = (l-m)!/(l+m)! P_lm`, with no phase), section 4's list of where the
  angular parity appears, and the code.

  The interesting part is not the convention but the **five copies**. `unfolded`,
  `from_unfolded`, `Coeffs.__getitem__`, the visibility sum and V6's check each computed the
  extension independently, and two changes of convention in two days produced failures from the
  copies rather than from the mathematics --- one silent (`__getitem__` was missed entirely, and
  a map came out complex), one caught by a docstring that had gone stale. All five now call
  `coefficients.negative_order`, so the next change is one edit. The audit that confirms it is
  `grep -rn 'conjugate(block\|conjugate(row\|conjugate(self.blocks' src/`, which should find only
  the helper itself; a bare `grep np.conjugate` does not work, because conjugating a *harmonic* is
  a genuine operation in the addition theorem and the zonal identity and is not the extension.

  Also worth recording: `make data` and `make earth-levels` are independent, and both write
  *coefficients*. Either one forgotten after a convention change leaves a stale cache that no
  test notices --- `earth_levels.py` did exactly that for twenty minutes, rendering the Earth at
  the wrong peak while the suite was green.

* **The Condon--Shortley phase is dropped, and T12 was found to have a mistake in it.** The
  phase put a `(-1)^m` into the real basis, into the Zernike polynomial and into both closed
  forms; dropping it removes three of those from the page and costs one factor in the code,
  where `rotate_coefficients` now sandwiches the standard (phase-carrying) `d` matrices with
  `diag((-1)^m)`. The trade is the same one T12 made, in the same direction: a simpler
  equation, a more complicated implementation.

  What did *not* change, contrary to the proposal: the extension relation keeps its `(-1)^m`,
  because both `Y_lm` and `Y_{l,-m}` pick up the same factor. And the real basis is now the
  negative of what it was at odd `m`, since it is `sqrt(2) Re Y_lm`.

  Found on the way, and more important than the change itself: **T12 had broken the mode
  formula**. Section 4 was given a prefactor of `sqrt(pi)` in place of `2 pi`, on the reasoning
  that the `sqrt(4 pi)` had to go somewhere. It did not: `2 pi` comes from the angular integral
  alone and does not move with the normalisation. The code was never changed, so every test
  still passed; the only disagreement was between section 4 and the appendix, which had kept
  an untouched copy of the same display. Two copies of one formula caught what 157 tests could
  not, which is the argument for the duplication the summary appendix used to provide --- and
  a warning about the merge that removed it.

  Two follow-ups the green suite did not catch on its own. `earth_levels.py` had been emitted
  under the phase-carrying basis, so the Earth rendered at a different peak (1.263 instead of
  1.351) while every test passed; `make earth-levels` fixed it, and the lesson is that any
  committed artefact holding *coefficients* is a cache of the convention and has to be
  regenerated with it. And section 3's rotation law is the one equation the change made longer,
  which is where the trade is visible: three `(-1)^m` off the basis and the closed forms, one
  equation and one code factor on.

* **The two appendices are one, and it is a brief rather than a summary.** Appendix A (the
  results) and Appendix B (the same displays with the certificates beside them) had 16 displays
  each, verified identical in order and content --- so B was a strict superset and A duplicated
  it. They are now a single appendix, keeping B's checkmark items and A's opening, its
  conventions sentence and its title. Four pages lighter, and the drift hazard of editing two
  copies is gone.

  The merge follows a change of destination rather than a change of mind about the content.
  Appendix A was meant to become section 0, but a good section 0 has to say something about the
  certificates, and that is not knowable until the certificates have been looked at; and the
  human-written document in `expos/` is an *exposition of the manual*, not a second source. So
  the appendix is now the **brief**: the results in the order a paper would present them, each
  beside the certificate that backs it, terse and complete. `expos/` may reorder and reword
  freely; where the two disagree, the appendix is right.

  What the merge added, beyond deleting a copy, is the paragraph the appendix never had: what a
  certificate *is* --- a differential test, two independent computations, a stated tolerance,
  a residual returned --- the three classes, and the two limits, that independence is only as
  good as the routes and that the trusted base is numpy rather than a kernel.

* **The polar cap is retired, and V28 with it.** The cap survived T10 only because the
  certificates used it, and the argument for keeping it was that V28 needed a discontinuous map.
  That argument fails on inspection: the power-law decay of a truncated discontinuous map is
  textbook, nothing in this work turns on the exponent, and V28 measured a property of a
  synthetic map rather than a claim of the method. The two-level Earth had already replaced the
  cap in the GUI, being real data with the same kind of sharp edge plus a taper that answers the
  ringing. So `cap_coefficients`, its constants, `_zonal_starspot` and `check_convergence_rate`
  are gone, and section 7 now cites the known result and describes the taper as the answer
  instead of certifying a rate.

  V15 stays, because its claim is about *our* code --- one rotation and both parity families
  against a direct transform --- and it is the only check that would catch an error common to
  both families. It now runs on a truncation of the Earth's own coefficients, which keeps it on
  an artefact the package renders rather than a synthetic profile. The appendix loses its
  truncation-rate subsection, and V24 is left as the only *independent* certificate.

* **T12 closed: the harmonic normalisation is now `Y_00 = 1`.** The orthonormal convention,
  `int |Y|^2 dOmega = 1`, was replaced by the un-normalised one, so the *average* of `|Y_lm|^2`
  over the sphere is unity and the integral is `4 pi`. The motive was readability: the `sqrt(4 pi)`
  was written nine times over and is now written once, as `sqrt(pi)` in the mode formula, plus
  `1/4 pi` in the projection. Two things came out better than expected. In `sphere_coefficients`
  the `1/4 pi` of the projection and the `4 pi/N` of the golden-angle samples **cancel**, leaving a
  plain average over samples; and `a_00 = 1` is now the statement that a map has mean one, so
  `MEAN_FACTOR` is 1 and the display chain is trivial. `a_00 = 1` for an undeformed map follows
  from the convention rather than being imposed on it.

  What did not change: `X_lm` is built from `P_lm` and Bessel functions, so the closed forms, the
  radial integrals and every radial certificate were untouched --- which is most of the risk, and
  it did not materialise.

  What did change: `N_lm`, the mode formula's prefactor (`2 pi` to `sqrt(pi)`), the addition
  theorem (`4 pi/(2l+1)` to `1/(2l+1)`, `4 pi/5` to `1/5`), the five numbers, the analytic
  first-order coefficients, the two projections, the display chain, three test expectations that
  were arithmetic, and four certificates whose *statements* changed rather than their tolerances:
  V1 (`int |Y|^2 = 4 pi`), V3 (the theorem without its `4 pi`), V16 (flux `pi`), V25 (mean square
  one). V1's comparison against SymPy's `Ynm` now writes the `sqrt(4 pi)` out, which makes it the
  one place where the convention difference is asserted numerically rather than described.

  Two lessons, both about the sequence rather than the mathematics. First, the change was done in
  four stages with `make check` between them, and the failure list was the checklist; doing it in
  one pass would have produced twenty failures with no way to tell a stale expectation from a real
  slip. Second, three of the six defects found were **no-op edits**: `str.replace` calls whose
  pattern no longer matched, which silently did nothing and were reported as failures by the user
  rather than by me. Every edit that carried an `assert` either worked or stopped.

* **Five loose threads closed in one pass.** (i) Section 1 said negative orders are "wanted in
  exactly one place"; it now says they are *formed* in one place (the rotation law) and elsewhere
  appear only through the definition, which is true of the signed-order relabelling and of the real
  branches. (ii) I2's measurements were all pre-dating the `sky_basis` fix and are replaced: the
  affine error is **second order** in `eps` (0.0016 at `E = 0.033`, ratio 0.23 on halving), which is
  what linearising a stretch should give --- the "first order" reading was the reflected frame.
  (iii) Section 7 said limb darkening would buy accuracy "when it comes"; it has come, in section 10.
  (iv) `derivations._zonal_starspot` was made to delegate to `maps.cap_coefficients` instead of
  being a second implementation, the original reason for the duplication having expired with the
  plotting code; the two were verified identical to 0.0 before the merge. Both have since been
  retired with the cap itself. (v) D21 is closed as a decision:
  the brightness exponent stays unexposed, trivial in two rows and too approximate in the third.

* **The GUI scaffolding is shared, and the map GUI exists.** `gui_core.py` holds the coalescing
  timer (`Coalescer`), the colour-scale cache, the East-leftward flip, the panel layout, the slider
  factory and the widget-retention line; `gui.py` and `map_gui.py` keep only their models and their
  slider specifications. The move exposed a double flip: `map_gui.compute_map` and
  `compute_visibility` flipped their own columns *and* `build` now flips again, which would have
  mirrored the image back to where it started. No test could have caught it --- they compare shapes
  and statistics, not orientation --- which is the argument for keeping the flip in one place and
  doing it once. Six new tests cover the core's pure parts. Was T5.

* **All-subscript indices.** `Y_{lm}`, `V_{lm}`, `P_{lm}`, `R_{n,m}`, `K_{lm}`, `Q_{lm}`, with
  `Y_{l,-m}` for a negative order and `R_{l,m,c}` / `R_{l,m,s}` / `R_{l,m,\sigma}` for the real
  basis. About 240 substitutions across the manual, the roadmap's own formulas and the docstrings of
  `spherical_harmonics.py`, `zernike.py` and `visibility.py` --- three times the estimate, because the
  braced and primed variants (`Y_l^{m'}`, `R_{l'}^{m',\sigma'}`) were not in the count. Since then the Zernike radial has gone back to `R_n^m` --- the superscript is the universal external form, and the subscript looked wrong beside every source that uses it --- and our radial transform, which had been `R_{lm}(q)`, is now `X_{lm}(q)`: one letter could not mean both, and `X` was free.  What stayed
  superscripted is exactly what should be: powers, the Wigner degrees `D^l` and `d^l`, the Jacobi
  parameters `P_k^{(m,0)}`, and transposes. Section 5's convention note needed no change, since it
  accounts for phase and normalisation rather than for the placement of the index.

* **The manual is figure-free, and the figures' machinery is gone.** `figures.py` became `maps.py`,
  keeping only `sky_map` and its batch form after the cap and V28 were retired;
  `figure_map`, `figure_uv`, `figure_convergence` and `write()` are deleted, with
  `tests/test_figures.py`, the `FIGURES` list, the figure half of the `data` target, and the three
  committed PDFs. The two figure environments in sections 1 and 4 became sentences: the display
  convention they carried (North up, East to the left, D23) and the reading of the uv panels (real
  even, imaginary odd) now stand in the text, pointing at `make mapgui` for the live version. The
  manual is 32 pages, down from 34. The convergence rate they plotted became certificate V28, in
  section 7, and has since been retired along with the cap.

* **The fitted rung was implemented (was T3).** `symbolic/fitted.py` carries the disk harmonics
  `Phi_lm = Y_lm(n)`, the design plan and its Gram matrix, `fit_image` with both residual routes,
  `residual_of` for comparing against a projection, and `model_image` for display.  It is the model's
  fourth row, labelled *fitted, `l = 0,2,4`*, with the residual in the panel titles; the Gram
  conditioning is printed once at startup, being a property of the basis and not of the parameters.
  Certificates V26 and V27.
* **O1 was not an open decision.** "The reference direction of the baseline azimuth" presented `u`
  along East as the status quo and the task as a `pi/2` shift. D23 had already fixed the opposite ---
  `u` along North, `v` along East --- so `psi = arctan2(v, u)` *is* the position angle with no offset
  anywhere; and D23 records the North-through-West alternative as *rejected*, which is what the item
  had mistaken for the current convention. Nothing else needed changing: the manual is consistent
  throughout (section 1 on the axes and the azimuth reference, section 3 on PA, section 4 on the
  baseline direction, section 8 on the reordering), and so is the code --- every axis label in the
  figures and in the GUI reads `u` (North) and `v` (East).
* **Two open items retired as stale.** "Triaxial ellipsoid sources" was superseded by Stage 6 and
  said so in its own text, and "the target case for Stage 6" was to be settled "before the figures"
  --- figures that have since been deleted. Both were describing future work that had happened.
* **Limb darkening left the build order.** It was Stage 7, with a question rather than a design and
  a precondition that is now satisfied; it is T1, where the obligation and the two candidate routes
  are set out. The GUI moves from Stage 8 to Stage 7, so the build order now ends with a completed
  stage. Section 7 is split into three labelled lists --- T to be implemented, O open decisions,
  I ideas --- so that a cross-reference cannot be broken by inserting an entry.
* **Section 8 was a stale design annex.** It opened "not yet scheduled" though the work was done,
  quoted an `O(eps^2)` residual of 1 per cent where the measurement is 13, and left the rung-1 flux
  question open that D22 settled. The duplicated parts are gone; what remains is the cost argument
  (build the image once, FFT once, interpolate), the binary-specific consequences --- two independent
  angles, orientation as a frame change, orbital phase as a fitting axis, reflex motion --- and an
  honest statement of which prior art is and is not surveyed. The subsection that held the cleanup
  list is retired with it, the list being empty. Was cleanup C10.
* **The Jupiter caricature was a relic.** `figure_jupiter`, `jupiter_coefficients`, the
  `JUPITER_*` constants, the `_banded_moments` helper and the rendering test are gone; it was in
  `__all__` and in the Makefile's figure list but written by nothing and included nowhere, and it
  never earned its place. `eval_legendre` stayed for `cap_coefficients`, which has itself now
  gone; it survives because `maps` and `derivations` still import it. Was cleanup C12.
* **Two more figures were tested but written nowhere.** `figure_triaxial` and
  `figure_delta_visibility` are gone, with `triaxial_field`, the `TRIAXIAL_*` constants and the
  `DELTA_*` defaults that only they used, and with their rendering tests. Superseded by the GUI in
  the same way as the pair C1 removed: they showed the deformed star against the affine prediction,
  which the GUI now does interactively. `SPHERE_L_MAX` stays, since the GUI imports it. Was cleanup
  C11.
* **The `data` target listed figures nothing wrote.** `FIGURES` in the Makefile named six PDFs;
  `figures.write()` emits three. `triaxial.pdf` and `delta-visibility.pdf` went with C11, and
  `jupiter.pdf` was never written at all, so the list is now the three that exist. `figure_jupiter`
  remains in the module, tested and unused --- logged as C12.
* **`docs/generated/verification.tex` was promised and never built.** Dropped rather than
  implemented: the registers now live at the end of the sections that make the claims, so a
  generated global register would be a second copy --- and a copy that can drift, which is what
  cleanup C9 was about. Two sentences elsewhere assumed such a thing existed, in section 1.5 and in
  the preamble comment, and now say what is true: each tag resolves to a named function, and the
  register of the section that cites it lists it. Was cleanup C6.
* **`V2` was cited nowhere, and part of section 2 was untagged.** `V2` now has its citation where the
  extension is defined, in section 1, and its row in section 2's register; the real basis gets a tag
  of its own, `V25`, covering both of its checks; and `V3`'s row now says that it covers the complex
  and real forms of the addition theorem, since two functions sit behind it. All thirty checks are
  now behind a tag, and every tag is both registered and cited. Was cleanup C9.
* **The formula inventory was stale and incomplete.** F15 now records that its table is built and
  certified by V12 and V14; F19--F21 point at the sections that exist, with the anchors split
  between the uniform disk (derived, V16, section 4) and the limb-darkened one (new, section 11);
  and six entries were added for material the inventory claimed to cover but did not --- the
  deformation field and its closure bound, the first-order relations, the closed-form coefficients,
  rung 1 with its image-space counterpart, the pairing, and the disk-mode surrogate. Was cleanup C8.
* **Stage 8 (now Stage 7) described a superseded GUI.** It spoke of two radio buttons and four combinations; the
  interface has three models, a radio for image against `(u,v)` and one for resolution, with the
  sliders in degrees. Rather than chase the current layout, the paragraph now records only what is
  stable --- the compute path, the three models, the clamps, the mode count --- and points at D25
  and at `make gui` for the rest. Was cleanup C5.
* **Stage 6 named a renamed function and quoted superseded defaults.** It said `body_axes`,
  which is `principal_axes`, and gave `eps_tide = 0.09`, `eps_spin = 0.03` with everything
  downstream of that arithmetic. The defaults are now 0.096 and 0.032, the fraction of the
  closure bound recomputed at 97 per cent, and the three principal semi-axes dropped: they were
  tied to the old amplitudes and cannot be re-derived without running `axis_lengths()`. Was
  cleanup C4.
* **The documentation plan numbered the sections wrongly.** It listed `10-limb-darkening.tex`
  and `11-numerics.tex`, sketching a manual without the pairing section. Section 5 now reads
  10-limb-darkening, 11-numerics; the pairing section it sketched has since been dissolved into
  sections 1 and 4, so the manual now has ten sections and two still to come. Was cleanup C2.
* **The API listing described a package that did not exist.** Section 2 now names what the tree
  contains: `gui.py` is in the layout, the five functions that were never written are gone
  (`visibility_limb_darkened`, `odd_case_expression`, `even_case_coefficients`,
  `generate_table`, `emit_latex_table`), and the `render`, `affine`, `emit`, `figures` and
  `gui` entry points are listed. Two prose claims went with them: that `radial_integral`
  dispatches on parity (it is the quadrature reference; `visibility_mode` dispatches), and
  that there is a `--lmax` CLI (both modules are run bare, and their defaults are the
  committed artifacts). Was cleanup C3.
* **The certificates table stopped at V20.** V21--V24 are now in section 4, matching the manual's
  registers word for word, and the rule above the table was corrected to describe what is
  actually done: a tag may cover a group of checks, as V1 does. Was cleanup C7.
* **The two superseded Stage 6 figures.** `stretched-sphere.pdf` and `visibility-maps.pdf` are
  gone, with their renderers and their rendering tests: the GUI shows the same comparison
  interactively, neither figure was ever included in the manual, and the number they printed is
  now asserted by the affine residual test. Was cleanup C1.
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

## 8. Non-spherical sources: what is left of the design notes

The approach (D11, D12) is now **implemented**, and the account of it is section 8 of the manual:
the deformation field, the level surface, the first-order relations, the closed-form coefficients and
rung 1. The measured error is I2 above. What follows is only what the manual does not carry, the
rest having been removed as duplicated.

**Cost, and its mitigation.** A numerical silhouette integral is milliseconds per `(u,v)`, which is
hopeless inside a fitting loop (`10^4` data points x `10^6` samples). But the image does not depend on
`(u,v)`: build it once per model evaluation and FFT it to obtain the whole uv plane, then interpolate
to the sampled points. Not as fast as the closed form, which costs `O(number of modes)` Bessel
evaluations per point, but it keeps the numerical route viable as a fallback --- and it is the thing
any surrogate of I3 would be checked against.

**Binary-specific consequences.** Orbit and spin introduce structure the spherical case lacks.

- **Two independent angles.** The surface map rotates with the spin phase; the tidal bulge is
  oriented by the orbital phase. For non-synchronous rotation these differ, and conflating them would
  be wrong.
- **Orientation is a frame-change problem.** The deformation is naturally defined in the co-rotating
  frame, so presenting it to the observer is exactly the rotation the Wigner D machinery of manual
  section 3 already handles. That machinery is not only for surface maps.
- **Orbital phase becomes a fitting axis.** The ellipsoidal term varies with orbital phase, and the
  classical photometric ellipsoidal modulation runs at twice the orbital frequency; the same harmonic
  content appears in the visibility.
- **Reflex motion.** The star's wobble about the system barycentre is a sky-plane offset, hence a
  linear phase ramp in `(u,v)`. Trivial to include, but it must not be mistaken for structure in the
  visibility amplitude.

**Prior art, one half done.** The single-star literature is covered to the extent the manual's
`dholakia2025` citation reaches it --- closed-form visibilities for a rotating, spotted, spherical
star. Two things are not surveyed: rapidly rotating oblate gravity-darkened stars (Domiciano de
Souza and collaborators, notably the Achernar work), and the ellipsoidal variables of close
binaries, where this same expansion in spherical harmonics is the standard photometric framework and
the conventions may already exist in a form we could adopt. In both cases the value added here is
closed form and speed, not feasibility.
