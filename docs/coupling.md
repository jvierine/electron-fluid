# Coupling to a field model

`coupling::FieldCoupledFluid` is an object-safe host interface. The provided
`YeeFluidCoupler` implements it without depending on a particular Maxwell code.
Other discretizations can implement the same trait; they must supply their own
compatible interpolation, divergence, units and boundary treatment.

## Inputs and units

The current adapter resolves z and one horizontal Fourier mode exp(i kx x).
Each complex value is `[real, imaginary]`, with three components ordered x,y,z.
All array lengths equal the number of cells. x/y fields and density are at z[i];
z current and electric field are at z[i]+width[i]/2. Indices wrap periodically.
Physical end boundaries require absorbing buffers or a compatible host treatment.

| Quantity | Stored value |
|---|---|
| widths | metres |
| step length | c times physical timestep, metres |
| electric | V/m |
| plasma_frequency | omega_pe/c, inverse metres, on each component grid |
| collision_rate | physical drag frequency/c, inverse metres |
| gyro | e B/(m_e c), signed vector, inverse metres |
| thermal_speed | sqrt(3 k_B T_e/m_e)/c |
| sqrt_density | sqrt(n0 in m^-3), on each component grid |
| density_buffer | numerical density damping rate/c; separate from physical collisions |
| current j | -sqrt(m_e n0/epsilon0) u, in V/m |
| density q | sqrt(m_e c_e^2/epsilon0) n1/sqrt(n0), in V/m |

Thus J_SI = epsilon0*c*plasma_frequency*j, and
n1 = q*sqrt(n0)*sqrt(epsilon0/(m_e*c_e^2)). Do not interpret j as A/m².
The pressure-free limit is supported, but density recovery from q then becomes
singular; use the cold current only, without claiming a warm density response.

## One step

1. Start with E and q at n, magnetic field and current at n-1/2.
2. Advance the host magnetic field from E^n.
3. Call `fluid.advance(&electric)`: momentum advances current to n+1/2,
   then conservative continuity advances q to n+1.
4. Advance the host electric field with `fluid.current()` at n+1/2.

This current **replaces** any existing cold electron current. Do not add a
second electron population or a Poisson field computed from the same charge.
Initialize Gauss's law consistently, and use the same charge divergence as the
fluid adapter. The pressure/divergence pair is energy-adjoint on its Yee grid.
Physical collision rates are nonnegative; electron-electron collisions do not
act as bulk drag on the electron population.

Construction checks array shapes, positivity and the gyro fixed-point
contraction. The host remains responsible for the combined electromagnetic
and pressure-wave CFL condition. Invalid electric inputs fail before state is
modified. Geometry validation is performed once; field validation is performed
at each exchange. The original low-level callback adapter is still available
in `maxwell_coupling` for hosts that supply their own time-centering kernels.

See `examples/field_host.rs` for a complete independent Ampere host and
`tests/coupling.rs` for analytic oscillation, drag and nonuniform-charge checks.
