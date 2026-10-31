//! Reusable field-to-fluid coupling interface, independent of any Maxwell solver.
//!
//! Geometry is a periodic nonuniform Yee line with one Fourier-x mode. Physical
//! endpoints need absorbing layers supplied by the host. All lengths are metres;
//! internal time is c*t, rates are divided by c, and thermal speed is c_e/c.
//! Electric field is in V/m. Current is energy-normalized, not A/m²:
//! j = -sqrt(m_e*n0/epsilon0)*u, so J_SI = epsilon0*c*(omega_pe/c)*j.
//! Density q = sqrt(m_e*c_e²/epsilon0)*n1/sqrt(n0) has field units.
//! E and q are at integer times, current at half times. No Poisson solve occurs.
use crate::maxwell_coupling::{Complex, CouplingGeometry, MaxwellFluid};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CouplingError(pub &'static str);
impl std::fmt::Display for CouplingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { f.write_str(self.0) }
}
impl std::error::Error for CouplingError {}

/// Host-facing contract. Advance the fluid once using E^n, then use the returned
/// current in the host's Ampere step. Do not also advance a cold electron current.
/// The host must preserve the same divergence, staggering and initial Gauss law.
pub trait FieldCoupledFluid {
    fn advance(&mut self, electric: &[Vec<Complex>]) -> Result<(), CouplingError>;
    fn current(&self) -> &[Vec<Complex>];
    fn density(&self) -> &[Complex];
    /// Internal step c*dt, in metres.
    fn step_length(&self) -> f64;
}

/// Validated adapter around the existing momentum-continuity implementation.
/// Construction validates static inputs; advance validates the supplied field
/// before changing state. It does not validate the host's electromagnetic CFL.
pub struct YeeFluidCoupler<'a> {
    state: MaxwellFluid,
    geometry: CouplingGeometry<'a>,
    dt: f64,
}
fn components(a: &[Vec<f64>], n: usize, positive: bool) -> bool {
    a.len() == 3 && a.iter().all(|v| v.len() == n && v.iter().all(|x|
        x.is_finite() && if positive { *x > 0. } else { *x >= 0. }))
}
fn field_valid(a: &[Vec<Complex>], n: usize) -> bool {
    a.len() == 3 && a.iter().all(|v| v.len() == n && v.iter().flatten().all(|x| x.is_finite()))
}
impl<'a> YeeFluidCoupler<'a> {
    pub fn new(geometry: CouplingGeometry<'a>, dt: f64,
               current: Vec<Vec<Complex>>, density: Vec<Complex>) -> Result<Self, CouplingError> {
        let n = geometry.widths.len();
        if n < 3 || !geometry.widths.iter().all(|x| x.is_finite() && *x > 0.) {
            return Err(CouplingError("grid needs at least three finite positive cell widths"));
        }
        if !dt.is_finite() || dt <= 0. || !geometry.kx.is_finite()
            || !geometry.thermal_speed.is_finite() || geometry.thermal_speed < 0.
            || !geometry.gyro.iter().all(|x| x.is_finite()) {
            return Err(CouplingError("invalid timestep, wave number, thermal speed or gyro vector"));
        }
        let contraction = 0.5*dt*geometry.gyro.iter().map(|v| v*v).sum::<f64>().sqrt();
        if contraction >= 1. || !contraction.is_finite() {
            return Err(CouplingError("gyro iteration requires dt*|gyro|/2 < 1"));
        }
        if !components(geometry.plasma_frequency,n,false)
            || !components(geometry.collision_rate,n,false)
            || !components(geometry.sqrt_density,n,true)
            || geometry.density_buffer.len()!=n
            || !geometry.density_buffer.iter().all(|x| x.is_finite() && *x>=0.) {
            return Err(CouplingError("invalid component geometry, density, collision or buffer array"));
        }
        if !field_valid(&current,n) || density.len()!=n
            || !density.iter().flatten().all(|x| x.is_finite()) {
            return Err(CouplingError("invalid initial fluid arrays"));
        }
        Ok(Self {state:MaxwellFluid::new(current,density),geometry,dt})
    }
}
impl FieldCoupledFluid for YeeFluidCoupler<'_> {
    fn advance(&mut self, electric: &[Vec<Complex>]) -> Result<(), CouplingError> {
        if !field_valid(electric,self.geometry.widths.len()) {
            return Err(CouplingError("electric field needs three finite grid-sized components"));
        }
        self.state.advance_momentum(&self.geometry,self.dt,electric,
            |old, drive, gyro, dt, damping| (1.-damping)*old+dt*drive+0.5*dt*gyro,
            |rhs, gyro, dt, damping| (rhs+0.5*dt*gyro)/(1.+damping));
        self.state.advance_continuity(&self.geometry,self.dt);
        Ok(())
    }
    fn current(&self) -> &[Vec<Complex>] { &self.state.current }
    fn density(&self) -> &[Complex] { &self.state.density }
    fn step_length(&self) -> f64 { self.dt }
}
