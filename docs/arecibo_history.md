# Driven Arecibo electron fluid

The two-way Maxwell coupling is in `../coupled_sim`. Its driver calls this
library's momentum-continuity adapter and the existing FDTD field steps at each
timestep. See `../coupled_sim/README.md`. The prescribed-field calculations
below remain as earlier controls; they are not the two-way result.

## Directional plasma-line asymmetry

Native paired spectra at 17:03:50 UTC, 117.75 km, give a fitted enhanced-line
SNR of 89.98 and an opposite-line excess of -0.00624. The empirical off-range
standard deviation of the opposite-line template amplitude is 0.0186.
At 17:05:40 UTC, 117.60 km, the corresponding values are 13.83, 0.0102 and
0.0126. Both opposite-branch fits are consistent with no enhancement.
These fits use the same 44 microsecond sinc-squared template, mirrored about
zero frequency at 5.34003 MHz, with a separately fitted constant baseline.
The baseline is not clipped at zero. Three times the empirical scatter lies
32.1 and 25.6 dB below the respective enhanced amplitudes. Those are sensitivity
scales, not formal confidence bounds: adjacent range gates are correlated and
the noise background is estimated from the same frame. This check concerns
the 5.34 MHz feature at its matched height/frequency, not every frequency.

The recovered plotting convention reverses the archived frequency axis.
The original manuscript identifies the enhanced branch as RF upshifted; the
power arrays alone do not independently verify the receiver's RF/IQ sign.
For a vertical upward-looking radar, a downward positive-frequency Langmuir
wave produces the upshifted return. An upward Langmuir wave produces the
downshifted return. Upgoing plasma waves and upshifted plasma lines must not
be used interchangeably.

In the unmagnetized warm-fluid illustration,

\[
k_z^2(z)=\frac{\omega^2-\omega_{pe}^2(z)}{3v_{Te}^2}-k_x^2.
\]

At fixed frequency, this becomes negative on the denser side of the turning
level. There is an evanescent disturbance there, rather than a freely
propagating upward Langmuir channel. Below that level both signs of real kz
are locally possible; the source and radiation conditions select the
outgoing downward wave. An upward wave launched below the turning level
would approach it and reflect, so an increasing density gradient alone does
not prove that no upward waves can ever exist. The absence of an incident
Langmuir wave from below and of appreciable backscatter are model assumptions.
In the present magnetized vertical fluid closure, replace omega squared in
the numerator by D = omega^2 (omega^2-Omega^2)/(omega^2-Omega_z^2).
Its local kz squared changes sign at resonance; the 430 MHz match occurs
4.05 m below it. This is a local bottomside argument, not a prohibition on
topside excitation through another electromagnetic channel.

Applying the existing coherent radar operator to the density phasor and its
conjugate gives an opposite return 76--87 dB below the enhanced return as the
integration window changes. The small residual is window sensitive and is
not a converged prediction of the suppression factor. It supports strong
directionality, not an exactly zero return. The calculation retains the
shared approximation kB=2kr for both sidebands and the current collisional
fluid model; it does not add kinetic damping. Propagation-related plasma-line
asymmetry was discussed by [Muldrew (1978)](https://doi.org/10.1029/JA083iA11p05207);
the asymmetry alone does not uniquely identify the excitation mechanism.

```sh
conda run -n base python check_plasma_line_asymmetry.py
cargo run --release --manifest-path plotter/Cargo.toml --bin plasma_asymmetry
```

## Vertical electric-field gradients

For the measured-profile case normalized to the observed enhancement using the
coherent-beam ion-line calibration, within 10 m of resonance:

| Quantity | Model value |
|---|---:|
| Peak cold Maxwell vertical field | 0.0397 V/m |
| Peak magnitude of complex cold-field derivative | 0.267 V/m^2 |
| Peak absolute derivative of cold-field magnitude | 0.143 V/m^2 |
| Peak derivative of the imposed vertical source | 2.16e-5 V/m^2 |
| Peak derivative of the effective vertical drive, including magnetic coupling | 4.54e-5 V/m^2 |
| Warm total vertical-field derivative at the 430 MHz Bragg height | 0.0119 V/m^2 |

The complex derivative includes spatial phase changes. It is different from
the derivative of field magnitude. The cold total field includes its induced
polarization response and therefore is not the independent source in the
fluid calculation. Its large resonant gradient must not be added again while
independently generating a Poisson field. The Langmuir wave has its own field
gradient: in 1D, d_z E_self = -e n1/epsilon0. The fitted 1.90 ppm density wave
at the Bragg point consequently requires approximately 0.012 V/m^2 in this
wave's self-field. These numbers share the previous normalization assumptions;
they are not direct electric-field measurements.

For this fixed spatial shape, linear scaling gives a conditional required
cold-field gradient

\[
g_{\rm HF}=0.267\sqrt{\mathrm{SNR}_{\star}/89.98}\ \mathrm{V\,m^{-2}}.
\]

SNR here is the on-frequency matched-filter power ratio, not dB. Choosing a
target SNR gives a gradient for this particular field profile. It does not
define an instability threshold or establish the experiment's detection limit.

Continuity involves the flux gradient d_z(n0 u_z). In the simplest
unmagnetized forcing expression, the electric term is proportional to
d_z(n0 E_source,z) = n0 d_z E_source,z + E_source,z d_z n0. A spatially uniform
imposed field can compress a nonuniform electron background. For the actual
magnetized calculation the relevant forced velocity is

\[
\widetilde u_{z,\rm forced}=-\frac{e}{m_e}R_{zz}\widetilde E_{\rm eff},\quad
R_{zz}=\frac{a^2+\Omega_z^2}{a(a^2+\Omega_x^2+\Omega_z^2)},\quad a=\nu-i\omega,
\]
\[
\widetilde E_{\rm eff}=\widetilde E_{z,\rm src}
+\frac{\Omega_x\Omega_z}{a^2+\Omega_z^2}\widetilde E_{x,\rm src}
+\frac{a\Omega_x}{a^2+\Omega_z^2}\widetilde E_{y,\rm src}.
\]

The compressional source thus involves d_z(n0 Rzz E_eff), including phase and
magnetic polarization, rather than the maximum of d_z Ez alone. Eliminating
the horizontal source algebraically reproduces the computed density response.

A prescribed-source control smooths all complex components over 0.02--20 m,
renormalizing each to the same peak effective drive within 10 m of resonance.
It holds n0, temperature, magnetic field and collisions fixed. Smoother fields
still excite waves, and the radar response can increase: removing local phase
cancellation can outweigh the reduction of peak gradient. A constant-vector
source gives zero vertical source gradient in the interaction region yet a
nonzero density response. Moving its boundaries from [-240,160] m to
[-400,260] m leaves returned power unchanged to numerical precision. With
constant density and collision rate as well, its interior density response
falls from 7.43e-5 to 5.74e-19. Thus the effect in this control is not caused by
the distant source taper.

These smoothing and constant-field controls are not new self-consistent HF
Maxwell solutions. They refute a gradient-only criterion for the prescribed
fluid drive; they do not prove that arbitrary flat pump fields can be generated
by a transmitter in the measured layer. The density scale length at resonance
in this case is 1.013 km. A physical occurrence boundary must constrain the
density profile, incident HF field, polarization, damping and radar response
together. Linear mode-conversion treatments retain this coupled structure;
see [Gondarenko et al. (2003)](https://doi.org/10.1029/2003JA009985).

```sh
conda run -n base python scan_hf_field_gradient.py
cargo run --release --manifest-path plotter/Cargo.toml --bin hf_gradient
```

## HF power from the measured enhancement

The 17:03:50 UTC frame gives an excess-power SNR of 73.68 (18.67 dB) in
the brightest native bin at 5.33447 MHz and 117.75 km. A narrow line convolved
with the 44 us rectangular-envelope response fits a center of 5.34003 MHz and
an on-frequency SNR of 89.98 (19.54 dB). Letting the duration vary gives
43.13 us, consistent with the archived 44 us processing window. The fit uses
the signed raw excess relative to off-range noise, not the deconvolved peak
or display contrast. It excludes the separate 4--5 MHz feature.

The experimenter supplied 1 MW radar peak transmit power and 150 K system
temperature. The matched-filter noise bandwidth is 1/(44 us) = 22.727 kHz;
FFT bin spacing (12.207 kHz) is not this bandwidth. Thus the narrow-line
received power estimate is SNR times k_B T_sys / tau = 4.235e-15 W. Averaging
pulse powers reduces estimator uncertainty; it does not multiply this SNR
by the number of pulses.

### Ion-line normalization and cross-check

The ion line is integrated over +/-2 MHz, which captures its code-spread
power while excluding the plasma lines. Integrating over Doppler removes code
phase dependence by Parseval. For a constant-amplitude 44 us pulse, the
remaining range weight is a triangle with half-width c tau/2 = 6.596 km.
This is an integrated-power check, not a recreation of the unsaved waveform.
Widening to +/-2.5 MHz changes the integral by +0.85%; narrowing to +/-1.5 MHz
changes it by -1.13% over the calibration gates.

The prediction uses the native-cycle fitted ne, Te and Ti, with the density
scale anchored by broad plasma-line frequencies, not enhanced-line power.
The integrated ion-source weight is approximately
F_i = alpha_e^2 / [(1+alpha_e)(1+alpha_e+alpha_i)], where
alpha_s = n_s e^2/(epsilon_0 k_B T_s k_Bragg^2). Direct numerical integration
of the user's ISR spectrum implementation, including the measured magnetic
geometry, agrees with this approximation within 0.9% at three check heights.

With a nominal 61.1 dBi Gaussian radar beam, the observed ion-line integral is
0.3189 times the nominal prediction at 125--135 km. Therefore the nominal
absolute radar normalization does NOT pass unchanged: the discrepancy is
4.96 dB. Apply this single empirically determined multiplier to the modeled
received power while holding beam geometry fixed. On withheld 110--125 km
gates the measured/scaled-predicted ratio is 1.0057 in the median, with 0.89%
relative RMS. This checks profile consistency after calibration, not the
unscaled absolute radar model. The 0.9% scatter is not total absolute accuracy;
beam geometry, pulse-envelope assumptions and cross-frequency calibration
remain systematic uncertainties.

### Coherent-density radar response

For a normalized Gaussian radar mode u, integral |u|^2 dA = 1, the one-sideband
Born prediction is

\[
 P_{\rm PL}=P_{\rm radar}r_e^2\lambda_r^2
 \left|\int dz\,\frac{\widetilde n_1(z)}{2}
 e^{2ik_rz}\int dA\,u^2(x,y,z)e^{ik_xx}\right|^2.
\]

The factor 1/2 is essential because the saved complex density is a peak
phasor: n1 = Re[n_tilde exp(-i omega t)]. The radar beam is filled. The
horizontal pump phase, finite-distance beam spreading, and round-trip Fresnel
phase are included. Omitting these phases and coherently summing a beam volume
would overestimate the returned power. The vertical integral includes the
whole complex density solution, rather than assigning an arbitrary slab
thickness at the Bragg point. Changing its integration window changes power
by 0.016%; numerical transverse quadrature reproduces the analytic overlap.

After the ion-line normalization, the coherent baseline requires an incident
complex-vector peak field of 0.213 mV/m at the model reference height. The
Bragg density amplitude is 1.90 ppm. The corresponding incoming vertical HF
flux is 4.35e-11 W/m^2. This is the current 5.33 MHz, 300 K fluid response;
the field calculation has not been rerun at the fitted 5.340 MHz line center.

For a 125 km lossless free-space-equivalent path and the model's 16.75 degree
vacuum-equivalent incidence, this is 8.9 W EIRP. The equivalent vacuum field
is 0.131 mV/m RMS. It is a conditional result for the coherent one-way fluid
model, not an identified transmitter power or a unique consequence of beam
filling. The HF antenna gain toward the layer and any additional path loss
remain explicit:

\[
P_{\rm HF}=\frac{\mathrm{EIRP}_{125}}{G_{\rm HF}}
 \left(\frac{r}{125\,\mathrm{km}}\right)^2 10^{L_{\rm extra}/10}.
\]

### Transverse phase-coherence sensitivity

Beam filling specifies occupied area, not phase correlation. Keep unit filling
but multiply the transverse covariance by exp[-d_perp^2/(2 ell_perp^2)]. With
a = 2/w^2 and b = k_r/R_curvature, the horizontal response factor is

\[
H(\ell_\perp)=\frac{a^2}{a^2+b^2+a/\ell_\perp^2}
\exp\left[-\frac{k_x^2a}{2(a^2+b^2+a/\ell_\perp^2)}\right].
\]

Relative to the coherent case, required power scales as H(infinity)/H(ell).
The calculation is a radar-response sensitivity for the same local vertical
wave amplitude, not a self-consistent simulation of lateral phase disorder.
It gives approximately 49 W EIRP at ell=10 m and 4.54 kW at ell=1 m, compared
with 8.9 W for the coherent limit. The slight nonmonotonicity at large ell is
real within this assumed covariance: broadening the angular spectrum can
improve overlap with the vertical radar. This means the coherent limit need
not be a strict lower bound. Very short ell comparable to the Bragg wavelength
would also require revisiting the 1D fluid approximation.

Sources for the radar baseline and ion-line normalization approach:
[Mathews (2013)](https://doi.org/10.5194/hgss-4-19-2013) and
[Vierinen et al. (2017)](https://doi.org/10.1002/2017GL073141).

```sh
conda run -n base python estimate_hf_power.py
cargo run --release --manifest-path plotter/Cargo.toml --bin hf_power
```

## Bragg match and downward continuation

The extended stationary calculation uses the same imposed Maxwell field and
smooth measured density, with a 1 mm mesh from 160 m below to 80 m above the
vertical fluid resonance. The displayed interval ends at 100 m below resonance,
before the absorbing buffer. The earlier Rust time-domain result is preserved;
this extension solves its frequency-domain equations, not a longer transient.

For the local collisionless limit of this particular 1D3V fluid closure,

\[
c_e^2=3k_BT_e/m_e,\qquad
D(\omega)=\frac{\omega^2(\omega^2-\Omega_e^2)}{\omega^2-\Omega_{ez}^2},
\qquad k_z=-\sqrt{\frac{D-\omega_{pe}^2(z)}{c_e^2}}.
\]

Here altitude increases upward, so negative real k_z is the downward branch.
The frequency is fixed at 5.33 MHz. Decreasing density increases |k_z| while
shortening the wavelength. At 300 K:

| Height relative to resonance | Wavelength | Downward wave-number magnitude |
|---|---|---|
| -4.0526 m: 430 MHz Bragg match | 0.34860 m | 18.024 rad/m |
| -10 m | 0.22100 m | 28.431 rad/m |
| -20 m | 0.15523 m | 40.476 rad/m |
| -40 m | 0.10843 m | 57.946 rad/m |
| -100 m | 0.06666 m | 94.265 rad/m |

The resonance is at 117.9416 km above the radar. The formal fluid bound is
sqrt(D)/c_e = 284.10 rad/m (minimum wavelength 2.2116 cm) as density tends to
zero. Over the supplied density profile down to 100.2 km, the local expression
reaches only 232.54 rad/m. These are dispersion diagnostics, not evidence that
the wave survives to those heights. At large k times Debye length the fluid
closure cannot supply the missing kinetic/Landau-damping constraint. There is
no sharply defined damping cutoff without specifying an amplitude criterion.

The full driven solution eventually contains a dominant smooth forced response.
Its phase gradient must not be interpreted as the outgoing Langmuir wavelength.
A separate homogeneous continuation below -10 m is normalized to the complex
density at its uppermost cell. It removes further local HF forcing and follows
the outgoing solution into the same density and collision profiles. This is a
diagnostic continuation, not an exact decomposition of the driven solution.

For a local collision-only attenuation estimate, set a = nu - i omega and
R_zz = (a^2 + Omega_ez^2) / [a (a^2 + Omega_e^2)]. The downward spatial root is
k_z = -sqrt[(i omega/R_zz - omega_pe^2)/c_e^2], with negative imaginary part.
The amplitude attenuation is exp[-integral Im(k_z) dz] along the downward path.
Starting at the Bragg height, its attenuation reaches -20 dB at -11.91 m,
-40 dB at -24.04 m, and -107.48 dB at -100 m. These decibels express the damping
factor (20 log amplitude, equivalently 10 log squared amplitude), excluding
inhomogeneous amplitude transport and renewed forcing. They are not an ISR
detection threshold. Landau damping is absent from both extended calculations.

Halving the mesh from 2 to 1 mm changes the Bragg amplitude by 0.0072%, the
continued phase wave number at -100 m by 0.112%, and its very small amplitude by
1.46%. Enlarging the domain leaves those results unchanged to better than
1e-8 relative. The 1 mm grid Nyquist wave number is 3142 rad/m, well above the
formal fluid bound; the displayed limit is not a grid cutoff.

```sh
conda run -n base python extend_fluid_kz.py
cargo run --release --manifest-path plotter/Cargo.toml --bin fluid_kz
```

The HDF5 output retains the full complex solution, source, background, outgoing
continuation, convergence checks and local-dispersion arrays. Figures show the
Bragg match, 100 m continuation, and the wider-profile fluid bound separately.

## Recorded HF field: current default

The default now uses the saved complex Maxwell fields and current. The source
amplitude is the **incident complex-vector field norm**, initially 0.5 mV/m.
The vertical drive must be matched to the actual one-dimensional Poisson solver:

\[
\widetilde{\mathbf E}_{\rm src}=
(\widetilde E_{x,c},\widetilde E_{y,c},
 \widetilde E_{z,c}+i\widetilde J_{z,c}/(\omega\epsilon_0)).
\]

Here c denotes the recorded cold Maxwell solution and the time convention is
exp(-i omega t). Vertical continuity and Poisson give
`d_t E_self = -J_z/epsilon0`, hence the cold vertical polarization field is
`-i J_z,c/(omega epsilon0)`. Subtracting it avoids adding the cold restoring
response again as an external force. Horizontal fields are held at their
recorded values; the full magnetic force couples all three velocity components.
The pressure-free model recovers the recorded cold vertical field to 0.35%,
including interpolation between the Maxwell and fluid staggered meshes.

The Fourier-x Helmholtz transverse field is NOT the correct input to a strictly
vertical continuity equation: dropping horizontal compression while retaining
that split produced a spurious response and roughly 598 times the incident
power. That failed control is retained separately in
`output/fluid/fdtd_driven_vertical_hodge_control.h5`. The corrected stationary
source work is 0.874 times the incident vertical flux for this case. This passes
a necessary power check but is not a self-consistent conversion efficiency;
the prescribed pump still cannot respond to the warm-plasma loading.

This is one-way coupling with vertical derivatives only. Horizontal fluid
compression and feedback on the imposed HF wave are omitted. The recorded-field
mode fixes the measured profile, frequency and magnetic geometry to that solved
HF case. Temperature, fluid collisions, amplitude and pulse duration can be
varied as response controls. Airy and Gaussian modes remain available separately.

The native reference run uses a 40 m window from 24 m below to 16 m above
resonance, 2048 cells, 300 K, the recorded collision profile, a 50 µs ramp, and a
pulse ending at 1.5 ms. It starts from rest and runs to 4 ms. The browser window
is centered on resonance; doubling the native domain changes the stationary
Bragg amplitude by less than 0.00008%. Halving its cell spacing changes that
amplitude by 0.69%. The native and WASM solvers each agree with an independent
harmonic solution to about 0.002% over 10 µs. The existing fluid benchmarks and
new complex-source tests pass. These checks validate this reduced calculation,
not a fully coupled kinetic conversion model.

Results: the local fluid 430 MHz wavelength match is 4.053 m below resonance.
The harmonic density-wave phase propagates downward and its wavelength shortens.
At the chosen incident normalization, the stationary fractional amplitude at
that height is 0.000442%. The wavelength-averaged pressure flux is downward.
The corrected pulse loses 70.3 dB of energy between the first
post-switch-off sample and 4 ms. Physical collisions remove 89.62% of source
work; numerical buffers remove 10.38%. The energy residual is
0.0038%. Landau damping
is absent.

Run from the workspace root (the first command recreates the result file):

```sh
cargo build --release --manifest-path fluid_browser/Cargo.toml
conda run -n base python drive_fluid_from_fdtd.py
conda run -n base python drive_fluid_from_fdtd.py --check
conda run -n base python drive_fluid_from_fdtd.py --run
conda run -n base python prepare_fdtd_fluid_figures.py
cargo run --release --manifest-path plotter/Cargo.toml --bin fdtd_fluid -- --animation
```

`output/fluid/fdtd_driven.h5` retains the full source, smooth equilibrium,
collisions, harmonic solution, convergence cases, native time-height evolution,
energy ledger and final state. The movie animates one stationary HF period
slowed to four seconds; the time-height plot is the separate native turn-on and
switch-off calculation. Neither is calibrated radar power.

Current scope: validate the simple fluid model first. Landau damping and a
velocity-space kinetic solver are deferred. The measured-background browser
model remains a fluid–Poisson calculation with a prescribed electric-field source.

Rust advances the continuity and momentum equations in the browser. The background
is the existing smooth fit to calibrated Arecibo electron density. A prescribed
electric field is supplied by the user; Poisson's equation supplies the restoring
field. The source angle rotates E in the vertical–B plane, not an HF ray.

```sh
conda run -n base python prepare_fluid_browser.py
cd fluid_browser
rustup run stable cargo test --release
RUSTFLAGS='-C target-feature=+simd128' rustup run stable cargo build --release --target wasm32-unknown-unknown
cp target/wasm32-unknown-unknown/release/arecibo_fluid.wasm ../web/fluid/
```

Serve the workspace root with the existing local server, then open `/web/fluid/`.
The numerical solver runs inside a Web Worker; no solve request goes to Python.
HDF5 export runs locally using vendored h5wasm 0.10.3 (license included).

The optional Airy source is Ai((z-zc)/ellA), with cutoff and Airy scale evaluated
from the smooth measured density. Its signed phase reversals are retained,
while a green magnitude envelope remains visible before/after the pulse.
The amplitude slider specifies the Airy peak field, not the incident field.
The default E angle is 90 degrees (transverse to vertical); other angles are
controlled forcing experiments, not self-consistent oblique HF propagation.
This scalar, locally linear, unmagnetized reference is imposed on the magnetized
fluid. It must not be described as a solved Maxwell mode-conversion problem.
Reference: Fitzpatrick, Plasma Physics, Cutoffs, equations 6.19–6.21,
https://farside.ph.utexas.edu/teaching/plasma/Plasma/node76.html.

The Airy-mode window is 640 m and uses 32768 cells, preserving 1.953 cm spacing.
The unbuffered interval is 448 m, with more than four Airy lengths on each
side of the scalar cutoff for the default observation. Airy mode rejects
windows with fewer than three Airy lengths on either side before the buffers.
The source is smoothly tapered only in those buffers. This is a local reference
window, not a claim that the Airy tail ends there or that the locally linear
density approximation remains exact across the domain. The Gaussian control
remains available and is required for the uniform-density counterfactual.
Source shape and all settings are included in the HDF5 export.

`node check_airy_wasm.mjs` tests the delivered 32768-cell WASM, the applied
signed source, finite evolution, and its energy ledger. This startup check
is not a long-time domain-convergence or kinetic validation.

The linear model has one spatial coordinate and three velocity components.
Its pressure acceleration is `-3 kB Te/me * d_z(n1/n0)`, an effective barotropic
closure with the background supported by an external equilibrium force.
It gives the uniform Bohm–Gross dispersion but is not a complete stratified
kinetic pressure tensor. Ions are fixed. Number is conserved by staggered fluxes;
all time derivatives use RK4. Poisson integration sets the lower exterior
induced field to zero; charge conservation keeps the upper exterior field zero.
A polarization-induced potential difference is allowed. Face velocities vanish
at the ends. A cubic velocity-drag buffer
occupies each outer 15%; it is distinct from physical collisions.

Only nu_ei + nu_en enters physical momentum drag. Electron–electron collisions
conserve total electron momentum, and their transport effects are not included.
Landau damping is collisionless and is not included in this first fluid model.
See [collision conservation](https://farside.ph.utexas.edu/teaching/plasma/lectures1/node34.html)
and [Landau damping](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node109.html).

The default now evaluates physical rates on every velocity face. Electron–ion
drag uses the NRL (2023) Maxwellian momentum frequency with a density- and
temperature-dependent Coulomb logarithm, singly charged ions, and Te < 10 eV.
Electron–neutral drag uses Schunk & Nagy (2009), Table 4.6 (N2, O2 and O), with
the saved MSIS 2.1 neutral prior interpolated in time and log density in height.
The neutral prior preserves the original range-as-altitude convention, without
a station-elevation correction. It is a climatological prior, not a measurement.
The selected electron temperature is still an assumption. Gradient controls
stretch electron density only, leaving neutral densities at their real heights.
Manual constant rates remain available for controlled comparisons.

At the default 5.33 MHz resonance, Te=300 K, the rates are 2698/s (ei) and
3066/s (en). The weakly damped magnetized fluid branch gives a 341 microsecond
amplitude e-fold time and 170 microsecond energy e-fold time. This differs
slightly from nu/2 because of the magnetized velocity polarization. A local
ray estimate gives 28.2 dB collisional energy loss over 1.106 ms to the 430 MHz
matching height, 4.05 m below resonance. It estimates post-launch propagation;
it is not a conversion efficiency, radar power, or a substitute for resolving
the turning region. Losses for half/quarter gradients are 56.4/112.8 dB at
their respective matching heights under the same assumptions.

The new panels show collision profiles, ray attenuation, and a post-pulse
exponential reference. The physical-loss ledger uses the actual spatial drag
at every RK stage; numerical edge losses remain separate. All three neutral
species and the electron-ion/electron-neutral rates are included in HDF5 exports.
Neither the reference nor the ray diagnostic feeds back into the solver.

`conda run -n base python validate_fluid_damping.py` writes
`data/fluid_damping_prior.h5`. It independently compares the diagnostic's
magnetized decay rate with exact eigenvalues of the homogeneous four-state
fluid system. Tested angles are 0, 46.816 and 90 degrees; maximum relative
error is below 4e-11 at nu=6000/s. This is a fluid-model check, not validation
of the Maxwellian collision operator against a kinetic calculation.

The 1 µs exponentially averaged density RMS is a display diagnostic and does
not feed back into the dynamics. History columns retain the largest amplitude
in each height bin. The instantaneous Bragg projection is a spatial diagnostic,
not a spectral power density or calibrated ISR return. Neither EM propagation
nor conversion efficiency is calculated. A fine wave mesh does not improve
the measurement's density-profile resolution.

Tests check equilibrium, continuity/charge conservation, Poisson's law, homogeneous
dispersion and collisional ringdown, magnetic energy exchange, and the absence
of longitudinal compression from a horizontal field in the unmagnetized 1D limit.
These checks do not establish convergence or physical validation for the event.

For constant density, B = 0, and constant collision drag, eliminating velocity
and the self-field gives

\[
\partial_t^2 n_1+\nu\partial_t n_1+\omega_{pe}^2n_1
-3v_{te}^2\partial_z^2n_1
=\frac{en_0}{m_e}\partial_z E_{\rm src},
\quad v_{te}^2=\frac{k_BT_e}{m_e},\quad \nu=\nu_{ei}+\nu_{en}.
\]

For an unforced Fourier mode, the oscillation frequency is
`sqrt(wpe² + 3 vte² k² − ν²/4)` and its amplitude envelope is `exp(−νt/2)`.
The momentum-drag rate is therefore twice the weakly damped wave-amplitude
decay rate in this reference problem.

The Rust benchmark at the 430 MHz Bragg wavelength uses n0 = 3.5e11 m⁻³,
Te = 300 K, B = 0, and ν = 20000 s⁻¹. At 16, 32 and 64 cells per wavelength,
frequency errors relative to the continuum formula are 25.42, 6.41 and 1.61 ppm,
respectively, demonstrating second-order spatial convergence. At 64 cells per
wavelength the amplitude after approximately 100 µs is 0.36787818 of its initial
value, compared with 0.36787874 analytically at the actual integration time.
This benchmark excludes boundary damping and uses a standing mode compatible
with the closed boundaries. It checks the fluid equations, not Landau damping.

```sh
rustup run stable cargo test --release bragg_scale_frequency_and_decay_benchmark -- --nocapture
```

`node check_wasm.mjs` exercises the delivered WASM on the measured background:
at 100 µs, doubling the mesh from 2048 to 4096 cells in an 80 m window changes
total wave energy by 0.00027%. This is an early-time mesh check, not late-time
convergence of the Bragg-scale wave.

The default 17:03:50 UTC profile places the cold electrostatic resonance at
117.94159 km and the warm-fluid 430 MHz match 4.05 m below it (5.33 MHz,
300 K, 35.065 µT, B/vertical 46.816°). A local fluid-ray estimate gives about
1.106 ms for propagation from resonance to that height. Runs of a few hundred
microseconds can be dominated by turn-on phase mixing and need not yet show
the steady, fixed-frequency wavelength curve. The earlier Airy demonstration used
8 ms, with a source pulse ending at 1.5 ms, and compares actual zero-crossing
spacings against the local dispersion.
That comparison is diagnostic; a superposition of transient waves has no unique
local wavelength. The stop time can be increased without resetting the state.

The scale-length table uses analytic derivatives of the saved Hermite density
profile. For the default vertical fluid branch: Ln=1012.69 m at resonance,
the plasma-frequency scale is 2 Ln, the electrostatic turning scale is
0.23237 m, and the scalar unmagnetized EM Airy reference is 44.4628 m at its
different cutoff height. These two Airy scales must not be interchanged.
At the Bragg height, the group speed is 7339 m/s and the local amplitude/energy
damping lengths are 2.5037/1.2519 m. They are local values, not constant
attenuation lengths from resonance. The Debye length is 2.032 mm at resonance;
the gyroradius is 10.934 mm using vte=sqrt(kB Te/me). The Bragg k lambdaD
is 0.03670. The default grid has 17.85 cells per Bragg wavelength but does
not resolve Debye-scale kinetics. The prescribed 1 m Gaussian sigma is a
source control, not an inferred HF Airy scale; amplitude/intensity FWHM are
2.3548/1.6651 m. Measured density sampling remains 300 m.

Startup reports import/download/worker failures instead of staying at the
loading message. The solver and density profile are fetched without cache;
the GUI and worker imports carry a matching revision to avoid mixed assets.

The altitude window can be 40–1280 m, with up to 65536 cells and 100 ms duration.
At least 12 cells per 430 MHz Bragg wavelength are required. Rust advances in
mesh-dependent batches (512 steps at 4096 cells), yielding after about 70 ms of computation; the display updates
four times per second by default. Offscreen plots are not redrawn. Fixed amplitude
scales reveal fading waves; optional rescaling is explicitly labeled.

Source work, physical collision losses, and numerical-buffer losses are integrated
with the same RK4 stages as the state. The energy residual measures their balance.
The decay chart and HDF5 export retain energy and amplitude histories.

`../validate_fluid_gradient.py` independently solves the stationary harmonic
fluid–Poisson equations for the same prescribed Gaussian source. Results are in
`../data/fluid_gradient_validation.h5`. Stretching the profile about resonance
makes the Bragg distance 4.05, 8.11, and 16.21 m for gradient factors 1, 1/2,
and 1/4. At an assumed momentum-drag rate of 3000/s, density amplitude squared
at the respective Bragg heights falls by 11.6 and 37.5 dB relative to the measured
gradient. This is a local response proxy, not calibrated radar power. Without
collisions, a gentler gradient does not suppress this response. Uniform density
still permits a forced disturbance but no downward refraction to a Bragg match.

For gradient factor 1 and 3000/s drag, doubling the harmonic grid changes Bragg
amplitude by 0.018%; doubling the domain at fixed spacing changes it only at
roundoff. Across the six displayed nonuniform cases (no collisions and 3000/s),
grid doubling changes Bragg amplitudes by at most 0.52%; domain doubling changes
them by less than 3e-9 fractionally. This comparison holds the
local source fixed and does not model incident HF-field intensification.
