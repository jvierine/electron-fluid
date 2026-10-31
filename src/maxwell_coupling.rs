//! Adapter for the existing momentum-continuity model on a Maxwell Yee grid.
//! The standalone solver uses RK4/Poisson. This adapter exposes momentum and
//! continuity separately for the Maxwell leapfrog schedule, with Poisson OFF.
//! Current j=-sqrt(m*n0/eps0)*u and q=sqrt(m*cs2/eps0)*dn/sqrt(n0) are merely
//! energy-normalized fluid variables. Internal time c*t matches the FDTD API.
//! Fourier-x terms extend the existing vertical model for oblique incidence.
use crate::{gyro_acceleration, spatial_difference};
pub type Complex = [f64; 2];
#[inline(always)]
fn add(a: Complex, b: Complex) -> Complex {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline(always)]
fn mul(a: Complex, b: f64) -> Complex {
    [a[0] * b, a[1] * b]
}
#[inline(always)]
fn ik(a: Complex, k: f64) -> Complex {
    [-k * a[1], k * a[0]]
}

pub struct MaxwellFluid {
    pub current: Vec<Vec<Complex>>,
    pub density: Vec<Complex>,
    rhs: Vec<Vec<Complex>>,
    next: Vec<Vec<Complex>>,
}
pub struct CouplingGeometry<'a> {
    pub widths: &'a [f64],
    pub kx: f64,
    pub gyro: [f64; 3],
    pub plasma_frequency: &'a [Vec<f64>],
    pub collision_rate: &'a [Vec<f64>],
    pub sqrt_density: &'a [Vec<f64>],
    pub thermal_speed: f64,
    pub density_buffer: &'a [f64],
}
impl MaxwellFluid {
    pub fn new(current: Vec<Vec<Complex>>, density: Vec<Complex>) -> Self {
        Self {
            rhs: current.clone(),
            next: current.clone(),
            current,
            density,
        }
    }
    fn gyro(&self, g: &CouplingGeometry, q: usize, i: usize) -> Complex {
        let n = g.widths.len();
        let im = (i + n - 1) % n;
        let ip = (i + 1) % n;
        let mut out = [0.; 2];
        for r in 0..2 {
            // Adjoint interpolation of the energy-normalized current components.
            let v = if q < 2 {
                [
                    self.current[0][i][r],
                    self.current[1][i][r],
                    (self.current[2][i][r] * g.widths[i] + self.current[2][im][r] * g.widths[im])
                        / (g.widths[i] + g.widths[im]),
                ]
            } else {
                [
                    (self.current[0][i][r] + self.current[0][ip][r]) / 2.,
                    (self.current[1][i][r] + self.current[1][ip][r]) / 2.,
                    self.current[2][i][r],
                ]
            };
            out[r] = gyro_acceleration(g.gyro, v)[q];
        }
        out
    }
    /// E^n is supplied by Maxwell. Advance J^(n-1/2) to J^(n+1/2).
    /// Reuse FDTD's established trapezoidal gyro/drag time-centering kernels.
    pub fn advance_momentum<R, N>(
        &mut self,
        g: &CouplingGeometry,
        dt: f64,
        electric: &[Vec<Complex>],
        rhs_kernel: R,
        next_kernel: N,
    ) where
        R: Fn(f64, f64, f64, f64, f64) -> f64,
        N: Fn(f64, f64, f64, f64) -> f64,
    {
        let n = g.widths.len();
        for i in 0..n {
            let ip = (i + 1) % n;
            let pressure = [
                ik(self.density[i], g.thermal_speed * g.kx),
                [0., 0.],
                std::array::from_fn(|r| {
                    g.thermal_speed
                        * g.sqrt_density[2][i]
                        * spatial_difference(
                            self.density[ip][r] / g.sqrt_density[0][ip],
                            self.density[i][r] / g.sqrt_density[0][i],
                            1. / g.widths[i],
                        )
                }),
            ];
            for q in 0..3 {
                let gyro = self.gyro(g, q, i);
                let damping = 0.5 * dt * g.collision_rate[q][i];
                for r in 0..2 {
                    self.rhs[q][i][r] = rhs_kernel(
                        self.current[q][i][r],
                        g.plasma_frequency[q][i] * electric[q][i][r] + pressure[q][r],
                        gyro[r],
                        dt,
                        damping,
                    );
                }
            }
        }
        let contraction = 0.5 * dt * g.gyro.iter().map(|v| v * v).sum::<f64>().sqrt();
        let iterations = if contraction > 0. {
            (1e-14_f64.ln() / contraction.ln()).ceil().max(2.) as usize
        } else {
            1
        };
        for _ in 0..iterations {
            for q in 0..3 {
                for i in 0..n {
                    let gyro = self.gyro(g, q, i);
                    let damping = 0.5 * dt * g.collision_rate[q][i];
                    for r in 0..2 {
                        self.next[q][i][r] = next_kernel(self.rhs[q][i][r], gyro[r], dt, damping);
                    }
                }
            }
            std::mem::swap(&mut self.current, &mut self.next);
        }
    }
    /// Advance density with exactly the J^(n+1/2) passed to Maxwell's Ampere step.
    pub fn advance_continuity(&mut self, g: &CouplingGeometry, dt: f64) {
        let n = g.widths.len();
        for i in 0..n {
            let im = (i + n - 1) % n;
            let dual = 0.5 * (g.widths[i] + g.widths[im]);
            let div = add(
                ik(self.current[0][i], g.kx),
                std::array::from_fn(|r| {
                    spatial_difference(
                        g.sqrt_density[2][i] * self.current[2][i][r],
                        g.sqrt_density[2][im] * self.current[2][im][r],
                        1. / (dual * g.sqrt_density[0][i]),
                    )
                }),
            );
            let d = 0.5 * dt * g.density_buffer[i];
            self.density[i] = mul(
                add(mul(self.density[i], 1. - d), mul(div, dt * g.thermal_speed)),
                1. / (1. + d),
            );
        }
    }
}
