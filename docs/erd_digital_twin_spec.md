# ERD Digital Twin Technical Specification

## Purpose

This document is the concrete technical specification for the toy electrodeless-ring-discharge digital twin implemented in this repository. It has two jobs:

1. define the exact PDE plant, reduction pipeline, and controller now used by the codebase;
2. serve as a clean source document for figures, schematics, and report text.

The design target is a control-oriented surrogate problem with the following qualitative structure:

- an annular ring of density-like material is maintained by a source and wall-loss balance,
- the ring is only marginally stable in open loop,
- small stochastic perturbations can grow into visible symmetry breaking, swirl, and transport,
- a reduced-order model captures the dominant low-dimensional dynamics,
- a forecast-based controller attempts to preserve ring cohesion and symmetry.

Throughout, the saved high-dimensional state remains
\[
  x(t) = \begin{bmatrix}\operatorname{vec}(n(\cdot,t)) \\ \operatorname{vec}(\omega(\cdot,t))\end{bmatrix},
\]
so that the plant, dataset, POD, SINDy, and MPC contracts remain aligned.

## 1. Plant PDE System

### 1.1 Domain and coordinates

The plant evolves on a rectangularized annulus
\[
(r,\phi) \in [R_{\min}, R_{\max}] \times [0, 2\pi),
\]
with:

- uniform radial grid in \(r\),
- uniform periodic grid in \(\phi\),
- code convention `x := phi`, `y := r` inside the mesh implementation.

The rectangularized geometry is used only for the numerical operators. Diagnostic visuals may remap snapshots back onto an annular polar plot for presentation.

### 1.2 State variables

The evolving plant state contains exactly two fields:

- \(n(r,\phi,t)\): density-like ring field,
- \(\omega(r,\phi,t)\): vorticity-like field.

A diagnostic potential \(\varphi(r,\phi,t)\) is stored in the code as `psi` and is recomputed from \(\omega\) at each RK stage.

### 1.3 Potential-vorticity closure

The diagnostic Poisson problem is
\[
\nabla_\perp^2 \varphi = -\omega,
\]
with periodic azimuthal boundary conditions and radial Dirichlet conditions on \(\varphi\).

The velocity field is computed from the potential as
\[
 u_r = -\partial_\phi \varphi,
 \qquad
 u_\phi = \partial_r \varphi.
\]

This sign convention is fixed by the existing Poisson solver implementation and is preserved so downstream model and controller assumptions do not break.

### 1.4 Zonal / nonzonal decomposition

The zonal mean and fluctuation operators are
\[
\langle f \rangle(r,t) = \frac{1}{2\pi}\int_0^{2\pi} f(r,\phi,t)\,d\phi,
\qquad
\widetilde f = f - \langle f \rangle.
\]

In the code these appear as `_zonal_mean()` and `_tilde()`.

### 1.5 Nonlinear transport operator

The nonlinear advection is written in Poisson-bracket form
\[
\{A, B\} = \partial_r A\,\partial_\phi B - \partial_\phi A\,\partial_r B.
\]

The implementation uses an Arakawa-type discrete Jacobian instead of direct pointwise velocity advection. This is important because it is less numerically diffusive and better preserves the transport character needed for visible swirl and reduced-order structure.

### 1.6 Plant equations currently implemented

The code now supports two compatible closures under the same saved state and
artifact contracts:

- `mhw`: the thresholded modified-Hasegawa-Wakatani-style branch,
- `wave_landau`: the main experimental branch with explicit packet
  propagation, burst-speed amplification, and thresholded growth/saturation.

The current annular modified Hasegawa-Wakatani-like plant is
\[
\boxed{
\begin{aligned}
\partial_t n &= C\bigl(\widetilde\varphi - \widetilde n\bigr)
              - \{\varphi, n\}
              - \kappa_{\mathrm{eff}}(r,t)\,\partial_\phi \varphi
              - V_{\mathrm{ph}}(r,t)\,\partial_\phi \widetilde n
              - T_{\mathrm{sc}}(r,t)\,\{\varphi,n\}
              + \nu_n \, \nabla_\perp^{2p} n
              + S_{\mathrm{bal}}(n)
              - \kappa_w(r)\,n
              - D_E(r,t)\,\widetilde n,
\\[0.5em]
              \partial_t \omega &= C\bigl(\widetilde\varphi - \widetilde n\bigr)
              - T_{\mathrm{sc}}(r,t)\,\{\varphi, \omega\}
              - V_{\mathrm{ph}}(r,t)\,\partial_\phi \widetilde\omega
              + \nu_\omega \, \nabla_\perp^{2p} \omega
              - \gamma_\omega\,\omega
              + \Gamma_{\mathrm{curv}}(n)
              + \Gamma_{\mathrm{baro}}(n)
              + \Gamma_{\mathrm{flux}}(n,\varphi)
              + \Gamma_{\mathrm{landau}}(\varphi)
              + \Gamma_{\mathrm{fb},n\leftrightarrow\varphi}
              + f_\omega(r,\phi,t;u)
              + d_\omega(r,\phi,t)
              - D_E^\omega(r,t)\,\widetilde\omega.
\end{aligned}}
\]

Each term has a concrete purpose:

- \(C(\widetilde\varphi - \widetilde n)\): adiabaticity-like coupling between density and potential fluctuations,
- \(-\{\varphi, \cdot\}\): nonlinear transport,
- \(-\kappa_{\mathrm{eff}} \partial_\phi \varphi\): drift-wave drive in the density equation,
- \(\nu_n \nabla_\perp^{2p} n\), \(\nu_\omega \nabla_\perp^{2p}\omega\): Laplacian or biharmonic dissipation,
- \(-\gamma_\omega \omega\): vorticity drag,
- \(S_{\mathrm{bal}}(n)\): mass-balanced fueling that maintains a marginal operating point,
- \(-\kappa_w n\): wall loss,
- \(-D_E \widetilde n\), \(-D_E^\omega \widetilde \omega\): zonal-shear decorrelation of nonzonal activity,
- \(-V_{\mathrm{ph}}\partial_\phi \widetilde n\), \(-V_{\mathrm{ph}}\partial_\phi \widetilde \omega\): supercritical diamagnetic phase advection of unstable structures,
- \(-T_{\mathrm{sc}}\{\varphi,n\}\), \(-T_{\mathrm{sc}}\{\varphi,\omega\}\): supercritical transport amplification of nonlinear advection,
- \(\Gamma_{\mathrm{curv}}(n)\): linear density-to-vorticity curvature coupling,
- \(\Gamma_{\mathrm{baro}}(n)\): nonlinear baroclinic / interchange-style density-to-vorticity coupling,
- \(\Gamma_{\mathrm{flux}}(n,\varphi)\): flux-balanced / bHW-style vorticity source from radial transport,
- \(\Gamma_{\mathrm{landau}}(\varphi)\): large-scale potential feedback that keeps the ring lively near threshold,
- \(\Gamma_{\mathrm{fb},n\leftrightarrow\varphi}\): local positive fluctuation feedback that only turns on above threshold,
- \(f_\omega\): control forcing,
- \(d_\omega\): stochastic + weak coherent disturbance forcing.

### 1.7 Thresholded drift drive and zonal-shear regulation

The baseline drift-wave drive profile is
\[
\kappa(r) = \kappa_0 \exp\!\left(-\frac{(r-r_\star)^2}{2\sigma_\kappa^2}\right).
\]

The leading branch now uses a thresholded, state-dependent drive. Define the
dimensionless zonal density-gradient ratio
\[
 g(r,t) = \frac{-\partial_r \langle n\rangle(r,t)}{n_{\mathrm{amp}}/\sigma_\star},
\]
the normalized zonal ExB shear
\[
 s(r,t) = \frac{|\partial_r \langle u_\phi\rangle(r,t)|}{s_{\mathrm{ref}}},
\]
and the instability margin
\[
 m(r,t) = g(r,t) - g_{\mathrm{crit}} - \chi_s s(r,t).
\]

The code uses a smooth positive part
\[
 \mathcal A(m) = \frac12\left(m + \sqrt{m^2 + \delta_m^2}\right),
\]
where \(g_{\mathrm{crit}}\), \(\chi_s\), \(s_{\mathrm{ref}}\), and \(\delta_m\)
are controlled by

- `critical_gradient_ratio`,
- `shear_suppression_gain`,
- `shear_ref`,
- `threshold_width`.

The threshold activation uses the gradient magnitude, but the actual drift
drive keeps the signed flank direction. Writing
\[
 g_{\mathrm{sgn}}(r,t) = \frac{-\partial_r \langle n\rangle(r,t)}{n_{\mathrm{amp}}/\sigma_\star},
\]
the effective drift drive on the active threshold branch is
\[
\kappa_{\mathrm{eff}}(r,t) = \kappa(r) + \kappa_0 G_\kappa \mathcal A(m(r,t))\,g_{\mathrm{sgn}}(r,t),
\]
where `gradient_drive_gain = G_kappa`.

The same state is also used to set a supercritical phase-advection speed
\[
V_{\mathrm{ph}}(r,t) = V_{\mathrm{ph},0}\,\mathcal A(m(r,t))\,g_{\mathrm{sgn}}(r,t),
\]
with `phase_advection_gain = V_{ph,0}`. This term is included only on the
nonsymmetric components of \(n\) and \(\omega\), so coherent unstable patches
move faster around the ring as they become more supercritical.

The nonlinear advection itself is also amplified above threshold:
\[
T_{\mathrm{sc}}(r,t) = 1 + G_T \mathcal A(m(r,t)),
\]
with `supercritical_transport_gain = G_T`. This is a surrogate for the jump
from coherent drift-wave propagation to faster convective transport and
filamentation once the ring becomes locally supercritical.

Interpretation:

- the ring can remain coherent when the zonal profile is below threshold or when
  zonal shear is strong enough to quench growth,
- stochastic perturbations matter because they can push the profile/shear balance
  across this margin,
- once the margin turns positive, drift activity and velocity increase together,
- unstable structures gain their own azimuthal phase speed rather than waiting
  to be moved only by externally imposed forcing,
- inner and outer flanks of the ring can propagate in opposite azimuthal
  directions because the signed gradient changes across the ring.

The legacy additive zonal-gradient drive remains available when the threshold
parameters are zero, but it is no longer the preferred tuning branch.

### 1.7.1 Explicit `wave_landau` branch

When
\[
\texttt{pde.dynamics_model} = \texttt{"wave\_landau"},
\]
the code switches to a more explicit packet model:
\[
\boxed{
\begin{aligned}
\partial_t n
&=
C(\widetilde\varphi-\widetilde n)
-T_{\mathrm{sc}}\{\varphi,n\}
-V_{\mathrm{ph}}\partial_\phi \widetilde n
+G_n\left(c_\phi^{\mathrm{w}}\partial_\phi \widetilde\omega + c_r^{\mathrm{w}}\partial_r \widetilde\omega\right)
-\kappa_{\mathrm{eff}}\partial_\phi \varphi
+\Gamma_{n,\mathrm{Landau}}
+\Gamma_{n,\mathrm{vort}}
+\Gamma_{n,\mathrm{pack}}
+\Gamma_{n,\varphi}
+\nu_n\nabla_\perp^{2p}n
+S_{\mathrm{bal}}(n)
-\kappa_w n
-D_E\widetilde n,
\\[0.5em]
\partial_t \omega
&=
C(\widetilde\varphi-\widetilde n)
-T_{\mathrm{sc}}\{\varphi,\omega\}
-V_{\mathrm{ph}}\partial_\phi \widetilde\omega
+G_\omega\left(c_\phi^{\mathrm{w}}\partial_\phi \widetilde n + c_r^{\mathrm{w}}\partial_r \widetilde n\right)
+\Gamma_{\omega,\mathrm{Landau}}
+\Gamma_{\omega,n}
+\Gamma_{\omega,\varphi}
+\nu_\omega\nabla_\perp^{2p}\omega
-\gamma_\omega\omega
+f_\omega(r,\phi,t;u)
+d_\omega(r,\phi,t).
\end{aligned}}
\]

This is the branch intended to realize the current qualitative target:

- coherent and visibly dynamic ring below breakup,
- packets that propagate azimuthally and radially from the state itself,
- velocity that increases with instability activation,
- seed-dependent breakup once the same packets grow past threshold.

The self-generated packet speeds are
\[
c_\phi^{\mathrm{w}} = G_{\phi}^{\mathrm{w}}\,B_{\mathrm{w}}
\tanh\!\left(
0.55\,\frac{\widetilde\omega}{\omega_{\mathrm{ref}}}
+0.30\,\frac{\partial_\phi \widetilde n}{n_{\mathrm{ref}}}
+0.15\,\frac{\partial_\phi \widetilde\omega}{\omega_{\mathrm{ref}}}
\right),
\]
\[
c_r^{\mathrm{w}} = G_{r}^{\mathrm{w}}\,B_{\mathrm{w}}
\tanh\!\left(
0.55\,\frac{\partial_r \widetilde n}{g_{\mathrm{ref}}}
-0.25\,\frac{\partial_\phi \widetilde\omega}{\omega_{\mathrm{ref}}}
+0.20\,\frac{\partial_r \widetilde\omega}{g_{\omega,\mathrm{ref}}}
\right),
\]
with burst-speed amplification
\[
B_{\mathrm{w}}
=
1 + G_{\mathrm{burst}}^{\mathrm{w}}\mathcal A(m)
+0.35\tanh\!\left(\frac{|\widetilde\omega|}{\omega_{\mathrm{ref}}}\right).
\]

The density-side growth and corrugation terms are
\[
\Gamma_{n,\mathrm{Landau}}
=
\alpha_n \mathcal A(m)
\left(
\widetilde n - c_n\frac{\widetilde n^3}{n_{\mathrm{ref}}^2}
\right),
\]
\[
\Gamma_{n,\mathrm{vort}}
=
\beta_{n\omega}\,\bigl(0.35+\mathcal A(m)\bigr)\,
n_{\mathrm{ref}}\tanh\!\left(\frac{\widetilde\omega}{\omega_{\mathrm{ref}}}\right),
\]
and
\[
\Gamma_{n,\mathrm{pack}}
=
\beta_{n,\mathrm{pack}}\,
S(r)\,B_{\mathrm{w}}\,
\bigl(0.12+\mathcal A(m)\bigr)\,
\tanh\!\left(
\frac{\widetilde\omega + 0.35\,\widetilde\varphi}{\omega_{\mathrm{ref}}}
\right),
\]
with its area-weighted mean removed so it redistributes density rather than
changing the total mass,
while the vorticity-side growth is
\[
\Gamma_{\omega,\mathrm{Landau}}
=
\alpha_\omega \mathcal A(m)
\left(
\widetilde\omega - c_\omega\frac{\widetilde\omega^3}{\omega_{\mathrm{ref}}^2}
\right).
\]

The density equation also includes a weak state-aware inward confinement flux
\[
\Gamma_{n,\mathrm{conf}}
=
-\partial_r\!\bigl(v_{\mathrm{conf}}\,n\bigr),
\]
where the confinement speed is
\[
v_{\mathrm{conf}}
=
G_{\mathrm{conf}}\,
\frac{1}{1 + G_{\mathrm{rel}}\mathcal A(m)}
\left[
-\frac{r-r_\star}{\sigma_\star}
\exp\!\left(
-\frac12\left(\frac{(r-r_\star)/\sigma_\star}{1.6}\right)^2
\right)
\right].
\]

Operationally:

- below threshold it holds the ring radially coherent instead of allowing an
  immediate wall strike,
- above threshold it weakens, so the same packets that propagate azimuthally
  can also peel radially away and degrade the ring.

In practice this term is now deliberately weak on the active branch. It is no
longer the main mechanism that holds the ring together, because using a strong
full-field inward pinch made the plant look artificially squeezed and pushed it
toward inward collapse instead of a genuine near-threshold operating point.

The primary operating-point support is now an axisymmetric zonal profile term
\[
\Gamma_{n,\mathrm{zonal}}
=
G_{\mathrm{zp}}\,
\mathcal R_{\mathrm{zp}}(r,t)\,
\bigl(n_{\mathrm{eq}}(r)-\langle n\rangle(r,t)\bigr),
\]
with release factor
\[
\mathcal R_{\mathrm{zp}}(r,t)
=
\frac{1}{
1 + G_{\mathrm{zp,rel}}\mathcal A(m)
+ G_{\mathrm{zp,dam}}
\frac{|\langle n\rangle-n_{\mathrm{eq}}|}{n_{\mathrm{amp}}}
}.
\]

This is the key structural change in the current branch:

- it acts only on the zonal mean, not directly on nonzonal wave packets,
- it keeps the mean ring near the marginally stable operating profile,
- it releases once stochastic perturbations genuinely grow and damage the ring,
- it is therefore much more compatible with the intended control story than a
  hard inward pinch on the full two-dimensional density field.

Two numerical safeguards now matter to the qualitative behavior of this
branch:

- Density RK stages are rate-limited so one substep cannot remove more than a
  fixed fraction of the local density. This prevents the old failure mode
  where low-density pockets hit the positivity floor and the mass controller
  rebuilt a diffuse axisymmetric ring.
- Vorticity RK stages are also rate-limited, but more weakly. This keeps the
  active branch lively while preventing late-time \(\omega\) blow-up that
  would otherwise destroy ROM-friendliness.

In code, the main tuning knobs for this branch are:

- `wave_density_coupling_gain`,
- `wave_vorticity_coupling_gain`,
- `wave_packet_phi_gain`,
- `wave_packet_r_gain`,
- `wave_burst_speed_gain`,
- `omega_landau_gain`,
- `density_landau_gain`,
- `density_vortex_gain`,
- `density_packet_drive_gain`,
- `radial_confinement_gain`,
- `radial_release_gain`,
- `zonal_profile_restore_gain`,
- `zonal_profile_release_gain`,
- `zonal_profile_damage_release_gain`.

### 1.8 Density-to-vorticity couplings

The deeper physics branch adds explicit density-to-vorticity coupling:
\[
\Gamma_{\mathrm{curv}}(n) = \beta_{\mathrm{curv}} \, \kappa_{\mathrm{eff}}(r,t)\,\partial_\phi n.
\]

In code this is controlled by `pde.curvature_omega_gain`.

The deeper nonlinear branch adds a baroclinic-style source:
\[
\Gamma_{\mathrm{baro}}(n)
= \beta_{\mathrm{baro}} \, \kappa_{\mathrm{eff}}(r,t)
\left(\frac{\partial_r \widetilde n}{G_r}\right)
\left(\frac{\partial_\phi n}{G_\phi}\right),
\]
where \(G_r\) and \(G_\phi\) are reference gradient scales derived from the equilibrium-ring amplitude and width. In code this is controlled by `pde.baroclinic_omega_gain`.

The flux-balanced source is
\[
\Gamma_{\mathrm{flux}}(n,\varphi)
= \beta_{\mathrm{flux}}
\left(
\partial_r \langle n\rangle \,\widetilde u_r
- \partial_r(\widetilde u_r \widetilde n)
\right),
\qquad \widetilde u_r = -\partial_\phi \widetilde\varphi,
\]
with `pde.flux_balance_omega_gain = beta_flux`. On the current threshold branch
this source is multiplied by the same activation \(\mathcal A(m)\), so it
stays weak below threshold and grows only after the ring becomes locally
supercritical.

The large-scale potential feedback is
\[
\Gamma_{\mathrm{landau}}(\varphi) = \lambda_\varphi \widetilde\varphi,
\]
with `pde.landau_phi_gain = lambda_phi`.

The newest threshold-only fluctuation-growth loop is
\[
\Gamma_{\mathrm{fb},n\leftrightarrow\varphi}
=
\beta_{\mathrm{fb}} \mathcal A(m)
\begin{cases}
0.35\left(\widetilde n - c_n \widetilde n^3 / n_{\mathrm{ref}}^2\right), & \text{in the density equation},\\
\left(\widetilde\omega - c_\omega \widetilde\omega^3 / \omega_{\mathrm{ref}}^2\right), & \text{in the vorticity equation},
\end{cases}
\]
with `pde.supercritical_feedback_gain = beta_fb`.

This is a deliberate surrogate for the fact that, in the truly unstable
regime, nonzonal fluctuations do not merely drift faster; they also experience
local negative damping and then saturate nonlinearly. It is the minimum extra
ingredient needed to realize the intended qualitative picture:

- below threshold: coherent, fast, but non-breaking ring turbulence;
- above threshold: stochastic perturbations trigger rapid local growth and
  breakup.

Finally, the shear-decorrelation surrogate is
\[
D_E(r,t) = \mu_E\,s(r,t),
\qquad
D_E^\omega(r,t) = 0.75\,\mu_E\,s(r,t),
\]
with `pde.shear_damping_gain = mu_E`.

The intent is:

- density asymmetries should be able to generate vorticity locally,
- stochastic perturbations can grow through the plant dynamics instead of only through externally scripted forcing,
- the nonlinear branch can generate local positive and negative vortices rather than only amplifying one separable carrier pattern,
- the open-loop operating point can be made marginally stable rather than immediately breaking.

These terms are the current serious-physics extensions being evaluated against the earlier disturbance-dominated branches.

### 1.9 Dissipation operator

The code supports:
\[
\nabla_\perp^{2p} f =
\begin{cases}
\nabla_\perp^2 f, & p=1,\\
-\nabla_\perp^4 f, & p=2.
\end{cases}
\]

Operationally:

- `hyper_p = 1` gives Laplacian diffusion,
- `hyper_p = 2` gives dissipative biharmonic regularization.

The biharmonic branch is used in most visually active runs because it suppresses small-scale blow-up while allowing larger-scale transport to remain lively.

### 1.10 Source, equilibrium ring, wall loss, and mass correction

The reference ring profile is
\[
 n_{\mathrm{eq}}(r) = n_{\mathrm{bg}} + n_{\mathrm{amp}} \exp\!\left(-\frac{(r-r_\star)^2}{2\sigma_\star^2}\right).
\]

The static source envelope is
\[
 S(r) = S_0 \exp\!\left(-\frac{(r-r_\star)^2}{2\sigma_S^2}\right).
\]

The wall-loss layer is
\[
 \kappa_w(r) = \kappa_{w0} \exp\!\left(-\frac{(R_{\max}-r)^2}{2\delta_w^2}\right).
\]

The implemented fueling term is not simply \(S(r)\). It contains a fixed floor
plus dynamic feedback so the ring profile can actually build toward threshold:
\[
 S_{\mathrm{bal}}(n) = \alpha_{\mathrm{src}}(n)\,S(r),
\]
with
\[
\alpha_{\mathrm{src}}(n)
= \alpha_{\mathrm{floor}}
+ \alpha_{\mathrm{wall}}\frac{L_w(n)}{\int S(r)\,dA}
+ \alpha_{\mathrm{mass}}\,\varepsilon_M(n),
\]
where \(L_w(n)\) is the instantaneous wall-loss rate and \(\varepsilon_M\) is
the normalized total-mass error.

After the RK update, the code also applies a bounded additive mass relaxation
along an axisymmetric source/ring profile. This is intentionally *not* an exact
projection to the target mass every time step. The rationale is:

- exact projection was numerically clean but too strongly re-axisymmetrized the
  density field,
- the present branch needs the ring mass to remain realistic without erasing
  late-time nonzonal structure,
- a relaxed correction keeps the inventory in range while preserving the
  predictive disturbance-growth behavior needed for POD/SINDy/MPC.

### 1.11 Boundary conditions

The effective discrete boundary conditions are:

- periodic in \(\phi\) for all fields,
- radial Neumann for \(n\) and \(\omega\),
- radial Dirichlet for the diagnostic potential \(\varphi\) in the Poisson solve.

### 1.12 Control forcing

The control input is
\[
 u(t) = [u_0(t), u_1(t), u_2(t), u_3(t), u_4(t)]^\top.
\]

It enters only through the vorticity equation via
\[
 f_\omega(r,\phi,t;u) = u_0 b_0(r) + u_1 b_1(r,\phi) + u_2 b_2(r,\phi) + u_3 b_3(r,\phi) + u_4 b_4(r,\phi),
\]
with:
\[
 b_1 \sim g_1(r)\cos\phi,\quad
 b_2 \sim g_1(r)\sin\phi,\quad
 b_3 \sim g_2(r)\cos 2\phi,\quad
 b_4 \sim g_2(r)\sin 2\phi.
\]

This preserves the required indirect mapping
\[
 u \to f_\omega \to \omega \to \varphi \to \mathbf u \to n \to \text{metrics},
\]
which is essential for the ROM/MPC story.

### 1.13 Disturbance model

The current disturbance is intentionally hybrid:

1. weak coherent low-mode carriers,
2. seeded stochastic band-limited forcing in azimuthal Fourier space,
3. drifting localized dipole eddies,
4. stage-wise amplitude scheduling over the run.

The disturbance is reproducible by seed, but meant to look less scripted than the earlier deterministic-only carriers.

The code precomputes:

- correlated modal coefficients for the stochastic band,
- correlated phase-drift processes for low modes,
- correlated drifting eddy trajectories, amplitudes, and optionally orientations.

For localized eddies, the code now supports two geometries:

- `fixed_pair`: the stronger historical branch based on stochastic positive/negative Gaussian pairs,
- `rotating_pair`: a stochastic dipole whose orientation evolves in time, used mainly to probe less scripted operating-point behavior.

The forcing is projected to zero net mean each step so that it does not inject spurious net mass or a constant directional bias.

### 1.14 Time integration

The plant uses RK4. Each RK stage:

1. clips fields to safe finite bounds,
2. enforces Neumann projection on \(n\) and \(\omega\),
3. solves Poisson for \(\varphi\),
4. evaluates the full right-hand side,
5. advances to the next RK stage.

After the RK4 update:

1. density is mass-projected,
2. vorticity is clipped and projected onto radial Neumann boundaries,
3. the diagnostic potential and velocity are recomputed for output and metrics.

## 2. Plant Metrics and Diagnostics

### 2.1 Core scalar metrics used in modeling and control

The active scalar metrics are:

- profile deviation \(J_{\mathrm{prof}}\),
- low-mode nonaxisymmetric energy \(E_{\mathrm{low}}\),
- wall-loss rate \(L_w\),
- ring thickness \(\sigma_r\).

#### Profile deviation

Let
\[
 \bar n(r_i,t) = \frac{1}{N_\phi}\sum_{j=0}^{N_\phi-1} n(r_i,\phi_j,t).
\]
Then
\[
 J_{\mathrm{prof}}(t) = \frac12 \sum_i r_i\,\bigl(\bar n_i - n_{\mathrm{eq}}(r_i)\bigr)^2\,\Delta r.
\]

This measures degradation of the axisymmetric ring profile.

#### Low-mode energy

Using the azimuthal FFT
\[
 \hat n_{i,m} = \frac{1}{N_\phi}\sum_j n_{i,j} e^{-im\phi_j},
\]
with cutoff \(m_c = 4\),
\[
 E_{\mathrm{low}}(t) = \frac12 \sum_i
 \left(r_i e^{-\frac{(r_i-r_\star)^2}{2\sigma_w^2}}\right)
 \left(\sum_{m=1}^{4} |\hat n_{i,m}|^2\right)
 \Delta r.
\]

This measures low-order symmetry breaking around the ring.

#### Wall-loss rate

Using annular cell weights \(\Delta A_{ij} = r_i\Delta r\Delta\phi\),
\[
 L_w(t) = \sum_{i,j} \kappa_w(r_i)\,n_{i,j}(t)\,\Delta A_{ij}.
\]

#### Ring thickness

Define
\[
 M(t) = \sum_i \bar n_i r_i\Delta r,
\qquad
 r_{\mathrm{mean}}(t) = \frac{1}{M}\sum_i r_i\bar n_i r_i\Delta r,
\]
\[
 \sigma_r^2(t) = \frac{1}{M}\sum_i (r_i-r_{\mathrm{mean}})^2\bar n_i r_i\Delta r,
\qquad
 \sigma_r = \sqrt{\sigma_r^2}.
\]

### 2.2 Derived degradation metrics

The code also computes:
\[
L_{w,\mathrm{cum}}(t_k) = \sum_{\ell \le k} L_w(t_\ell)\,\Delta t,
\]
\[
J_{\mathrm{prof},\mathrm{excess}} = \max(J_{\mathrm{prof}} - J_{\mathrm{prof}}(0), 0),
\]
\[
E_{\mathrm{low},\mathrm{excess}} = \max(E_{\mathrm{low}} - E_{\mathrm{low}}(0), 0),
\]
\[
\sigma_{r,\mathrm{excess}} = \max(\sigma_r - \sigma_r(0), 0).
\]

These feed a composite badness score:
\[
\mathrm{badness}(t) = w_J J_{\mathrm{prof},\mathrm{excess}} + w_E E_{\mathrm{low},\mathrm{excess}} + w_L L_{w,\mathrm{cum}} + w_S \sigma_{r,\mathrm{excess}}.
\]

### 2.3 Visual-tuning diagnostics

To prevent regression toward boring or obviously scripted runs, the plant also logs diagnostics that are not used directly by MPC.

#### Variation diagnostics

These include:

- ratio of nonaxisymmetric energy to axisymmetric energy during the first 20 percent of the run,
- band-power ratios between low, mid, and high azimuthal mode ranges,
- relative frame-to-frame density change.

These are used to gate whether a run is visibly active at all.

#### Transport diagnostics

These include:

- median and maximum relative changes in \(n\) and \(\omega\),
- ring-band velocity magnitude statistics,
- separate statistics for \(|u_r|\) and \(|u_\phi|\),
- ratio \(|u_\phi|/|u_r|\) to diagnose radial vs azimuthal transport dominance,
- sign changes and sign-mixing metrics for the mean azimuthal velocity,
- low-mode phase drift in the \(m=1\) density coefficient,
- ring-band acceleration statistics for \(|u|\),
- total mass drift and mass span,
- radial row-correlation metrics used to detect stripe-like structures.
- spectral entropy of the ring-band azimuthal spectrum to detect collapse back to a single smooth mode.
- ring-mean gradient ratio, zonal shear, instability margin, and activation,
- supercritical fraction and threshold-burst counts,
- radial particle-flux burst metrics based on the ring-band \(q90\) of \(|n u_r|\),
- confinement and zonal-support magnitudes and release factors,
- inner-edge, outer-edge, and ring-core mass fractions,
- peak-time fractions for profile-deviation and low-mode-energy metrics, used to
  detect the undesirable pattern “instant spike then settle”.

These diagnostics are the main tuning tools for deciding whether a run is:

- too quiet,
- too radial,
- too one-directional,
- too coherent across radius,
- below threshold and therefore too quiet,
- supercritical but not actually transporting,
- numerically unphysical.

## 3. POD and Controlled SINDy Model

### 3.1 Data contract

Each plant snapshot is saved as the stacked field
\[
 x_k = \begin{bmatrix}\operatorname{vec}(n_k) \\ \operatorname{vec}(\omega_k)\end{bmatrix}
 \in \mathbb{R}^{2N},
 \qquad N = N_r N_\phi.
\]

The dataset loader preserves this contract across train/validation/test runs.

### 3.2 POD reduction

Given centered snapshots
\[
 X = [x_0-\bar x,\; x_1-\bar x,\; \dots,\; x_{K-1}-\bar x],
\]
compute the SVD
\[
 X = U\Sigma V^\top.
\]

A rank \(r\in\{4,5,6,7,8\}\) is selected subject to configuration and energy retention. The reduced coordinates are
\[
 a_k = U_r^\top (x_k - \bar x).
\]

### 3.3 Reduced derivative estimation

For each contiguous training trajectory, the reduced derivative is estimated by second-order central difference:
\[
 \dot a_k \approx \frac{a_{k+1}-a_{k-1}}{2\Delta t}.
\]

Endpoint samples are excluded per trajectory so time alignment with the control sequence remains exact.

### 3.4 Controlled SINDy structure

The reduced model fits
\[
 \dot a = f(a,u).
\]

The feature structure is split into:

- autonomous reduced-state features: polynomial library of degree 2 in \(a\), no bias,
- control residual features: linear terms in \(u\) and bilinear terms \(a_i u_j\).

Pure quadratic control terms \(u_i u_j\) are intentionally excluded.

### 3.5 TrappingSR3 hybrid fit

The implemented controlled model is hybrid:
\[
 \dot a = f_{\mathrm{trap}}(a) + g_{\mathrm{ctrl}}(a,u).
\]

Here:

- \(f_{\mathrm{trap}}\) is fitted with PySINDy TrappingSR3 to improve stability of the autonomous reduced dynamics,
- \(g_{\mathrm{ctrl}}\) is a sparse residual model in \([u, a_i u_j]\).

This preserves the required control-ready structure while still allowing the autonomous part to be regularized toward stable local behavior.

### 3.6 Model validation artifacts

The model run writes:

- POD basis and mean state,
- serialized controlled SINDy model,
- coefficient matrices and feature names,
- fit diagnostics including optimizer and fallback behavior,
- held-out rollout metrics,
- true-vs-predicted movies and contact sheets.

The held-out rollout selection is now activity-aware so the representative artifact is not a visually dead case.

## 4. MPC Controller

### 4.1 Controller structure

The controller is random-shooting MPC over the learned reduced model.

At time step \(k\):

1. observe the current full state \(x_k\),
2. project to reduced coordinates \(a_k\),
3. sample candidate control sequences,
4. propagate the reduced model forward over the horizon,
5. lift predicted reduced states back to full fields,
6. evaluate plant-style metrics from the lifted density field,
7. apply the first control from the minimum-cost candidate.

### 4.2 Horizon parameterization

The horizon length is \(H\). Candidate sequences are piecewise constant:

- `shoot_segments = S`,
- `shoot_seg_len = L`,
- `H = S L`.

Each candidate samples \(S\) control vectors and repeats each for \(L\) steps. Candidate 0 is always the hold-previous-control baseline.

### 4.3 Lifted-field cost

For each predicted step \(h\), the stage cost is
\[
\ell_h =
 w_J J_{\mathrm{prof}}(h)
 + w_{Jg}\,\max\bigl(0, J_{\mathrm{prof}}(h)-J_{\mathrm{prof,ref}}\bigr)^2
 + w_E E_{\mathrm{low}}(h)^2
 + w_{Eg}\,\max\bigl(0, E_{\mathrm{low}}(h)-E_{\mathrm{low,ref}}\bigr)^2
 + w_L L_w(h)
 + w_\sigma \sigma_r(h)^2
 + w_u \|u_h\|_2^2
 + w_{\Delta u}\|u_h-u_{h-1}\|_2^2.
\]

The horizon objective is
\[
 J = \sum_{h=1}^{H} \ell_h.
\]

The reference values \(J_{\mathrm{prof,ref}}\) and \(E_{\mathrm{low,ref}}\) are taken from the current state at the start of the MPC solve.

### 4.4 Metric evaluation inside MPC

The controller does not trust reduced coordinates directly as proxies for report metrics. Instead, each predicted reduced state is lifted back to full-order fields, and the cost is computed from the lifted density field using the same discrete formulas as the plant metrics.

This is important because the demonstration story depends on regulating physically meaningful quantities rather than an arbitrary reduced norm.

## 5. End-to-End Workflow

### 5.1 Data generation

`erd_fipy/scripts/generate_training_runs.py` produces reproducible train/validation/test trajectory folders and a manifest. The training trajectories use persistently exciting piecewise-constant controls so the controlled SINDy model sees enough variation in \(u\).

### 5.2 Model fit

`model/scripts/run_pipeline.py` loads the manifest, stacks snapshots, computes POD, fits the controlled SINDy model, validates on held-out trajectories, and writes a self-contained model run folder.

### 5.3 Closed-loop comparison

`control/scripts/run_closed_loop.py` loads the fitted model and a plant configuration, then produces:

- an open-loop baseline run,
- a closed-loop MPC run,
- comparison plots,
- relative-delta summaries,
- acceptance checks.

## 6. Run-Folder Artifact Contract

Every stage writes a timestamped run folder:
\[
\texttt{outputs/runs/<timestamp>_<tag>/}.
\]

Plant runs contain:

- `config.yaml`
- `fields/snapshots.h5`
- `controls/control_timeseries.h5`
- `metrics/curves.json`
- `metrics/variation_diagnostics.json`
- `metrics/transport_diagnostics.json`
- `plots/*.png`
- `movies/*.gif`
- `logs/run.log`

Model runs additionally contain:

- `model/basis_U.npy`
- `model/mean_state.npy`
- `model/sindy_model.pkl`
- `model/Xi.npy`
- `model/feature_names.json`
- `model/lineage.json`

Control runs additionally contain:

- `stages/open_loop/*`
- `stages/closed_loop/*`
- `metrics/relative_deltas.json`
- `metrics/acceptance.json`
- copied model artifacts under `model/`

## 7. Current Tuning Goals

The current tuning target is:

- marginally stable open-loop operating point,
- dynamic propagation around the ring even without breakup,
- breakup driven by stochastic perturbations that start small and grow,
- seed-dependent qualitative behavior,
- less visibly scripted forcing,
- stronger azimuthal transport and less radial striping,
- preserved ROM compatibility.

Two branches are especially important:

1. `wave_landau_marginal`
   - dynamic but still coherent,
   - intended as the operating-point debug mode,
   - should show packet propagation without full ring breakup.
2. `wave_landau_burst`
   - active near-threshold branch,
   - intended to let stochastic packets accelerate and grow into ring degradation.

The older threshold-mHW branches remain available for comparison, but they are
no longer the main development path.

## 8. Literature Inspirations

The current plant branch is not a literal implementation of any one plasma
edge model. It is a two-field control-oriented surrogate assembled from a small
set of ideas that are consistent with the drift-wave / zonal-flow literature:

1. **Modified / balanced Hasegawa-Wakatani dynamics**
   - The base two-field structure comes from modified Hasegawa-Wakatani models,
     where nonzonal density and potential fluctuations couple through an
     adiabaticity term and are regulated by self-generated zonal flows.
   - Key use here: keep the digital twin ROM-friendly while still allowing
     turbulence-zonal-flow competition.

2. **Threshold behavior near the Dimits-like transition**
   - The current branch explicitly models a critical-gradient margin and zonal
     shear suppression because the target operating point is not fully stable
     and not fully broken; it sits near a turbulence onset threshold.
   - Key use here: open-loop should be coherent-but-risky, with stochastic
     disturbances sometimes crossing a margin and growing.

3. **Flux-driven profile build and bursty relaxation**
   - The source floor plus wall-loss balance and radial-flux feedback are meant
     to emulate a flux-driven system whose profile slowly steepens and is then
     relaxed by transport events, rather than a purely initial-value problem.
   - Key use here: breakup should emerge from the evolving state, not from a
     fully scripted disturbance movie.

4. **Interchange / curvature-style vorticity generation**
   - The density-to-vorticity curvature and baroclinic terms are included so
     density perturbations can seed local vortices and streamer-like transport
     without adding new state variables.
   - Key use here: unstable regions should accelerate and generate motion
     through the PDE itself.

Representative references used for direction:

- Numata, Ball, and Dewar, “Bifurcation in electrostatic resistive drift wave turbulence” (Physics of Plasmas, 2007).
- St-Onge, Zhu, and Diamond, “A Flux-balanced Hasegawa-Wakatani Model” (arXiv:1706.02266, 2017).
- Zhu et al., “Phase Transition in the Modified Hasegawa-Wakatani Model” (arXiv:2406.13312, 2024).
- Camargo, Biskamp, and Scott, “Intermittent transport in drift-wave turbulence” / related drift-wave burst literature.

## 9. Figure / Schematic Prompt Notes

A paper-style workflow figure should include these boxes and arrows:

1. **Annular PDE Plant**
   - inputs: `u(t)`, disturbance seed/config,
   - states: `n(r,phi,t)`, `omega(r,phi,t)`, `psi(r,phi,t)`
2. **Diagnostics and Metrics**
   - `J_prof`, `E_low`, `L_w`, `sigma_r`
   - variation / transport diagnostics
3. **Snapshot Stacking**
   - `x = [vec(n), vec(omega)]`
4. **POD Reduction**
   - `x -> a`
5. **Controlled SINDy / TrappingSR3**
   - `dot(a) = f(a,u)`
6. **ROM Forecast**
   - horizon rollout in reduced space
7. **Random-Shooting MPC**
   - candidate control sequences
   - lifted-field metric cost
8. **Closed-Loop Plant**
   - apply first action
   - repeat

The visual story should emphasize the indirect control path:
\[
 u \rightarrow \omega\text{-forcing} \rightarrow \omega \rightarrow \varphi \rightarrow \mathbf u \rightarrow n \rightarrow \text{metrics}.
\]

## 10. Maintenance note

This document should be updated whenever any of the following change:

- the plant equations,
- forcing or disturbance parameterization,
- the set of primary metrics,
- the reduced-model library or optimizer strategy,
- the MPC cost or sampling policy,
- the run-folder contract.
