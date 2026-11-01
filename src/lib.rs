//! Linear 1D, three-velocity electron fluid, with a prescribed localized E source.
//! State: dn/n_ref at cell centers, ux,uy,uz on faces. SI throughout.
//! Continuity: dn_t = -d_z(n0 uz).
//! Momentum: u_t = -e/me (E_src + E_self zhat + u x B)
//!                  -3 kB Te/me d_z(dn/n0) zhat -(nu_ei+nu_en+sponge)u.
//! Poisson: d_z E_self=-e dn/eps0, E_self=0 at the lower boundary.
//! Charge-conserving closed end faces also keep the upper exterior field zero.
//! This is an effective barotropic pressure closure; equilibrium forces hold n0
//! fixed. It is not a kinetic closure in a stratified plasma. No Landau damping,
//! ee transport, Maxwell propagation, nonlinear saturation or ion motion here.
//! Staggered conservative differences; RK4; closed end faces and drag buffers.
use std::f64::consts::PI;
pub mod maxwell_coupling;
pub mod coupling;

/// Electron gyro acceleration, shared by standalone and Maxwell-coupled fluid.
#[inline(always)]
pub fn gyro_acceleration(omega: [f64; 3], u: [f64; 3]) -> [f64; 3] {
    [omega[1]*u[2]-omega[2]*u[1],
     omega[2]*u[0]-omega[0]*u[2],
     omega[0]*u[1]-omega[1]*u[0]]
}

/// Conservative flux divergence and enthalpy gradient use this same difference.
#[inline(always)]
pub fn spatial_difference(right: f64, left: f64, inverse_width: f64) -> f64 {
    (right-left)*inverse_width
}
const QE: f64 = 1.602176634e-19;
const ME: f64 = 9.1093837139e-31;
const EPS: f64 = 8.8541878188e-12;
const KB: f64 = 1.380649e-23;
const NR: f64 = 1e11;

#[derive(Clone)]
pub struct Source {
    pub amplitude: f64,
    pub frequency: f64,
    pub angle: f64,
    pub center: f64,
    pub width: f64,
    pub ramp: f64,
    pub duration: f64,
}
impl Source {
    fn envelope(&self, t: f64) -> f64 {
        let rise = (t / self.ramp.max(1e-12)).clamp(0., 1.);
        let fall = if self.duration > 0. {
            ((self.duration - t) / self.ramp.max(1e-12)).clamp(0., 1.)
        } else {
            1.
        };
        (0.5 - 0.5 * (PI * rise).cos()) * (0.5 - 0.5 * (PI * fall).cos())
    }
    fn temporal(&self, t: f64) -> f64 {
        self.amplitude * self.envelope(t) * (2. * PI * self.frequency * t).sin()
    }
}

pub struct Geometry {
    pub n: usize,
    pub lo: f64,
    pub dz: f64,
    pub n0: Vec<f64>,
    nf: Vec<f64>,
    flux_coefficient: Vec<f64>,
    pressure_coefficient: Vec<f64>,
    cs2: f64,
    bx: f64,
    bz: f64,
    nu: Vec<f64>,
    sponge: Vec<f64>,
    self_field: bool,
}
pub struct Simulation {
    pub g: Geometry,
    pub source: Source,
    pub time: f64,
    pub dt: f64,
    pub y: Vec<f64>,
    k1: Vec<f64>,
    k2: Vec<f64>,
    k3: Vec<f64>,
    k4: Vec<f64>,
    tmp: Vec<f64>,
    field: Vec<f64>,
    source_shape: Vec<f64>,
    source_phasor: Option<(Vec<f64>, Vec<f64>)>,
    rms2: Vec<f64>,
    output: Vec<f64>,
    energy_ledger: [f64; 3],
    initial_energy: Option<f64>,
}

/// Cubic Hermite interpolation of the saved smooth fit and its analytic slope.
/// This preserves sub-cell gradients without interpolating the noisy observations.
pub fn density_at(z: f64, zs: &[f64], ns: &[f64], ds: &[f64]) -> f64 {
    let i = zs
        .partition_point(|v| *v <= z)
        .saturating_sub(1)
        .min(zs.len() - 2);
    let h = zs[i + 1] - zs[i];
    let t = ((z - zs[i]) / h).clamp(0., 1.);
    (2. * t * t * t - 3. * t * t + 1.) * ns[i]
        + (t * t * t - 2. * t * t + t) * h * ds[i]
        + (-2. * t * t * t + 3. * t * t) * ns[i + 1]
        + (t * t * t - t * t) * h * ds[i + 1]
}

impl Simulation {
    pub fn new(
        n: usize,
        lo: f64,
        hi: f64,
        n0: Vec<f64>,
        te: f64,
        b: f64,
        b_angle: f64,
        nu: f64,
    ) -> Self {
        let dz = (hi - lo) / n as f64;
        let mut nf = vec![0.; n + 1];
        nf[0] = n0[0];
        nf[n] = n0[n - 1];
        for j in 1..n {
            nf[j] = 0.5 * (n0[j - 1] + n0[j]);
        }
        let omega = QE * b / ME;
        let max_wp2 = n0.iter().copied().fold(0., f64::max) * QE * QE / (EPS * ME);
        let cs2 = 3. * KB * te / ME;
        let dt = 0.06 / (max_wp2 + 4. * cs2 / (dz * dz) + omega * omega).sqrt();
        let sponge = (0..=n)
            .map(|j| {
                let d = (j.min(n - j) as f64) * dz;
                2e6 * (1. - d / (0.15 * (hi - lo))).clamp(0., 1.).powi(3)
            })
            .collect();
        let len = 4 * n + 3;
        let g = Geometry {
            n,
            lo,
            dz,
            pressure_coefficient: n0.iter().map(|v| cs2 * NR / (v * dz)).collect(),
            flux_coefficient: nf.iter().map(|v| v / (NR * dz)).collect(),
            n0,
            nf,
            cs2,
            bx: omega * b_angle.sin(),
            bz: omega * b_angle.cos(),
            nu: vec![nu; n + 1],
            sponge,
            self_field: true,
        };
        let source = Source {
            amplitude: 0.,
            frequency: 5.33e6,
            angle: 0.,
            center: (lo + hi) / 2.,
            width: 1.,
            ramp: 5e-6,
            duration: 0.,
        };
        let mut s = Self {
            g,
            source,
            time: 0.,
            dt,
            y: vec![0.; len],
            k1: vec![0.; len],
            k2: vec![0.; len],
            k3: vec![0.; len],
            k4: vec![0.; len],
            tmp: vec![0.; len],
            field: vec![0.; n + 1],
            source_shape: vec![0.; n + 1],
            source_phasor: None,
            rms2: vec![0.; n],
            output: vec![0.; 13 * n + 12],
            energy_ledger: [0.; 3],
            initial_energy: None,
        };
        s.update_shape();
        s
    }
    fn update_shape(&mut self) {
        for (j, v) in self.source_shape.iter_mut().enumerate() {
            let z = self.g.lo + j as f64 * self.g.dz;
            *v = (-0.5 * ((z - self.source.center) / self.source.width).powi(2)).exp();
        }
    }
    fn electric(g: &Geometry, y: &[f64], out: &mut [f64]) {
        out[0] = 0.;
        for j in 0..g.n {
            out[j + 1] = out[j] - QE * NR / EPS * y[j] * g.dz;
        }
        for v in out {
            *v = if g.self_field { *v } else { 0. };
        }
    }
    fn rhs(
        g: &Geometry,
        source: &Source,
        shape: &[f64],
        phasor: &Option<(Vec<f64>, Vec<f64>)>,
        t: f64,
        y: &[f64],
        dy: &mut [f64],
        field: &mut [f64],
    ) -> [f64; 3] {
        let n = g.n;
        let ux = n;
        let uy = 2 * n + 1;
        let uz = 3 * n + 2;
        Self::electric(g, y, field);
        let drive = source.temporal(t);
        let sx = source.angle.sin();
        let sz = source.angle.cos();
        let (sin, cos) = (2. * PI * source.frequency * t).sin_cos();
        let envelope = source.envelope(t);
        let mut powers = [0.; 3];
        for j in 0..n {
            dy[j] = -spatial_difference(g.flux_coefficient[j+1]*y[uz+j+1],
                                        g.flux_coefficient[j]*y[uz+j], 1.);
        }
        for j in 0..=n {
            if j == 0 || j == n {
                dy[ux + j] = 0.;
                dy[uy + j] = 0.;
                dy[uz + j] = 0.;
                continue;
            }
            let drag = g.nu[j] + g.sponge[j];
            let e = drive * shape[j];
            let ef = if let Some((re, im)) = phasor {
                std::array::from_fn::<_, 3, _>(|c| envelope * (re[c*(n+1)+j]*cos + im[c*(n+1)+j]*sin))
            } else { [e*sx, 0., e*sz] };
            let gyro=gyro_acceleration([g.bx,0.,g.bz],[y[ux+j],y[uy+j],y[uz+j]]);
            dy[ux + j] = -QE / ME * ef[0] + gyro[0] - drag * y[ux + j];
            dy[uy + j] = -QE / ME * ef[1] + gyro[1] - drag * y[uy + j];
            dy[uz + j] = -QE / ME * (ef[2] + field[j]) + gyro[2]
                - drag * y[uz + j]
                - spatial_difference(g.pressure_coefficient[j]*y[j],g.pressure_coefficient[j-1]*y[j-1],1.);
            let kinetic_factor = ME
                * g.nf[j]
                * (y[ux + j] * y[ux + j] + y[uy + j] * y[uy + j] + y[uz + j] * y[uz + j]);
            powers[0] -= QE * g.nf[j] * (ef[0]*y[ux+j] + ef[1]*y[uy+j] + ef[2]*y[uz+j]);
            powers[1] += kinetic_factor * g.nu[j];
            powers[2] += kinetic_factor * g.sponge[j];
        }
        powers.map(|v| v * g.dz)
    }
    pub fn step(&mut self) {
        let h = self.dt;
        if self.initial_energy.is_none() {
            self.initial_energy = Some(self.energy());
        }
        let power1 = Self::rhs(
            &self.g,
            &self.source,
            &self.source_shape,
            &self.source_phasor,
            self.time,
            &self.y,
            &mut self.k1,
            &mut self.field,
        );
        for i in 0..self.y.len() {
            self.tmp[i] = self.y[i] + h * 0.5 * self.k1[i];
        }
        let power2 = Self::rhs(
            &self.g,
            &self.source,
            &self.source_shape,
            &self.source_phasor,
            self.time + h * 0.5,
            &self.tmp,
            &mut self.k2,
            &mut self.field,
        );
        for i in 0..self.y.len() {
            self.tmp[i] = self.y[i] + h * 0.5 * self.k2[i];
        }
        let power3 = Self::rhs(
            &self.g,
            &self.source,
            &self.source_shape,
            &self.source_phasor,
            self.time + h * 0.5,
            &self.tmp,
            &mut self.k3,
            &mut self.field,
        );
        for i in 0..self.y.len() {
            self.tmp[i] = self.y[i] + h * self.k3[i];
        }
        let power4 = Self::rhs(
            &self.g,
            &self.source,
            &self.source_shape,
            &self.source_phasor,
            self.time + h,
            &self.tmp,
            &mut self.k4,
            &mut self.field,
        );
        for i in 0..self.y.len() {
            self.y[i] += h / 6. * (self.k1[i] + 2. * self.k2[i] + 2. * self.k3[i] + self.k4[i]);
        }
        self.time += h;
        for j in 0..3 {
            self.energy_ledger[j] +=
                h / 6. * (power1[j] + 2. * power2[j] + 2. * power3[j] + power4[j]);
        }
        let a = 1. - (-h / 1e-6).exp();
        for j in 0..self.g.n {
            let eta = self.y[j] * NR / self.g.n0[j];
            self.rms2[j] += a * (eta * eta - self.rms2[j]);
        }
    }
    pub fn energy(&self) -> f64 {
        let n = self.g.n;
        let mut field = vec![0.; n + 1];
        Self::electric(&self.g, &self.y, &mut field);
        let mut sum = 0.;
        for j in 0..n {
            let dn = self.y[j] * NR;
            sum += 0.5 * ME * self.g.cs2 * dn * dn / self.g.n0[j];
        }
        for j in 0..=n {
            let w = if j == 0 || j == n { 0.5 } else { 1. };
            let u2 = self.y[n + j].powi(2)
                + self.y[2 * n + 1 + j].powi(2)
                + self.y[3 * n + 2 + j].powi(2);
            sum += w * (0.5 * ME * self.g.nf[j] * u2 + 0.5 * EPS * field[j] * field[j]);
        }
        sum * self.g.dz
    }
    fn snapshot(&mut self) -> *const f64 {
        let n = self.g.n;
        Self::electric(&self.g, &self.y, &mut self.field);
        let d = self.source.temporal(self.time);
        let sx = self.source.angle.sin();
        let sz = self.source.angle.cos();
        let mut max_eta: f64 = 0.;
        let mut br = 0.;
        let mut bi = 0.;
        let mut den = 0.;
        let kb = 4. * PI * 430e6 / 299792458.;
        for j in 0..n {
            let z = self.g.lo + (j as f64 + 0.5) * self.g.dz;
            let eta = self.y[j] * NR / self.g.n0[j];
            let field = 0.5 * (self.field[j] + self.field[j + 1]);
            let src = 0.5 * d * (self.source_shape[j] + self.source_shape[j + 1]);
            let mut src_x = src*sx;
            let mut src_z = src*sz;
            if let Some((re, im)) = &self.source_phasor {
                let (sn,cs) = (2.*PI*self.source.frequency*self.time).sin_cos();
                let field = |c: usize| 0.5*self.source.envelope(self.time)*
                    ((re[c*(n+1)+j]+re[c*(n+1)+j+1])*cs+(im[c*(n+1)+j]+im[c*(n+1)+j+1])*sn);
                src_x=field(0);src_z=field(2);
            }
            let grad = if j > 0 && j + 1 < n {
                (self.g.n0[j + 1] - self.g.n0[j - 1]) / (2. * self.g.dz * self.g.n0[j])
            } else {
                0.
            };
            let values = [
                z,
                self.g.n0[j],
                eta,
                src_z,
                field,
                src_x,
                0.5 * (self.y[3 * n + 2 + j] + self.y[3 * n + 3 + j]),
                0.5 * (self.y[n + j] + self.y[n + j + 1]),
                0.5 * (self.y[2 * n + 1 + j] + self.y[2 * n + 2 + j]),
                grad,
                0.5 * (self.g.nu[j] + self.g.nu[j + 1]),
                self.g.sponge[j],
                (2. * self.rms2[j]).sqrt(),
            ];
            for (a, v) in values.iter().enumerate() {
                self.output[a * n + j] = *v;
            }
            max_eta = max_eta.max(eta.abs());
            // Spatial Bragg projection of instantaneous density; not an ISR PSD.
            let w = (PI * (j as f64 + 0.5) / n as f64).sin().powi(2);
            br += w * self.y[j] * NR * (kb * (z - self.source.center)).cos();
            bi += w * self.y[j] * NR * (kb * (z - self.source.center)).sin();
            den += w * self.g.n0[j];
        }
        let energy = self.energy();
        let diag = [
            self.time,
            self.dt,
            max_eta,
            energy,
            self.y[..n].iter().sum::<f64>() * NR * self.g.dz,
            (br * br + bi * bi) / (den * den),
            self.g.dz,
            self.source.envelope(self.time),
            self.energy_ledger[0],
            self.energy_ledger[1],
            self.energy_ledger[2],
            self.initial_energy.unwrap_or(energy) + self.energy_ledger[0]
                - energy
                - self.energy_ledger[1]
                - self.energy_ledger[2],
        ];
        self.output[13 * n..].copy_from_slice(&diag);
        self.output.as_ptr()
    }
}

// Small C ABI: identical native/WASM solver, no JavaScript numerical stepping.
#[no_mangle]
pub extern "C" fn alloc_f64(n: usize) -> *mut f64 {
    Box::into_raw(vec![0.; n].into_boxed_slice()) as *mut f64
}
#[no_mangle]
pub unsafe extern "C" fn free_f64(p: *mut f64, n: usize) {
    drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(p, n)));
}
#[no_mangle]
pub unsafe extern "C" fn create(
    n: usize,
    lo: f64,
    hi: f64,
    te: f64,
    b: f64,
    angle: f64,
    nu_ei: f64,
    nu_en: f64,
    p: *const f64,
    count: usize,
) -> *mut Simulation {
    if n < 32
        || n > 65536
        || count < 2
        || !([lo, hi, te, b, angle, nu_ei, nu_en]
            .iter()
            .all(|x| x.is_finite()))
        || hi <= lo
        || te <= 0.
        || b < 0.
        || nu_ei < 0.
        || nu_en < 0.
    {
        return std::ptr::null_mut();
    }
    let a = std::slice::from_raw_parts(p, 3 * count);
    let (z, rest) = a.split_at(count);
    let (ne, d) = rest.split_at(count);
    if !a.iter().all(|x| x.is_finite())
        || z.windows(2).any(|x| x[1] <= x[0])
        || lo < z[0]
        || hi > z[count - 1]
    {
        return std::ptr::null_mut();
    }
    let density = (0..n)
        .map(|j| density_at(lo + (j as f64 + 0.5) * (hi - lo) / n as f64, z, ne, d))
        .collect::<Vec<_>>();
    if density.iter().any(|v| *v <= 0.) {
        return std::ptr::null_mut();
    }
    Box::into_raw(Box::new(Simulation::new(
        n,
        lo,
        hi,
        density,
        te,
        b,
        angle,
        nu_ei + nu_en,
    )))
}
#[no_mangle]
pub unsafe extern "C" fn collision_profile(p: *mut Simulation, rates: *const f64, count: usize) -> u32 {
    let s = &mut *p;
    if count != s.g.n + 1 || s.time != 0. { return 0; }
    let rates = std::slice::from_raw_parts(rates, count);
    // The explicit RK step must resolve drag as well as plasma oscillations.
    if rates.iter().any(|v| !v.is_finite() || *v < 0. || *v * s.dt > 0.1) { return 0; }
    s.g.nu.copy_from_slice(rates);
    1
}
#[no_mangle]
pub unsafe extern "C" fn destroy(p: *mut Simulation) {
    drop(Box::from_raw(p));
}
#[no_mangle]
pub unsafe extern "C" fn source(
    p: *mut Simulation,
    amp: f64,
    f: f64,
    angle: f64,
    z: f64,
    width: f64,
    ramp: f64,
    duration: f64,
) -> u32 {
    if ![amp, f, angle, z, width, ramp, duration]
        .iter()
        .all(|x| x.is_finite())
        || amp < 0.
        || f <= 0.
        || width <= 0.
        || ramp <= 0.
        || duration < 0.
    {
        return 0;
    }
    let s = &mut *p;
    s.source = Source {
        amplitude: amp,
        frequency: f,
        angle,
        center: z,
        width,
        ramp,
        duration,
    };
    s.source_phasor = None;
    s.update_shape();
    1
}
#[no_mangle]
pub unsafe extern "C" fn source_profile(p: *mut Simulation, values: *const f64, count: usize) -> u32 {
    let s = &mut *p;
    if count != s.g.n + 1 { return 0; }
    let values = std::slice::from_raw_parts(values, count);
    if values.iter().any(|v| !v.is_finite()) { return 0; }
    s.source_shape.copy_from_slice(values);
    1
}
/// Complex physical E amplitude on velocity faces; component-major x,y,z.
/// Re[Etilde exp(-i omega t)] times the existing pulse envelope. No Gaussian
/// reshaping, angle rotation, or amplitude renormalization is applied.
#[no_mangle]
pub unsafe extern "C" fn source_phasor(p: *mut Simulation, re: *const f64, im: *const f64, count: usize) -> u32 {
    let s = &mut *p;
    if count != 3*(s.g.n+1) { return 0; }
    let re=std::slice::from_raw_parts(re,count);
    let im=std::slice::from_raw_parts(im,count);
    if re.iter().chain(im).any(|x| !x.is_finite()) { return 0; }
    s.source_phasor=Some((re.to_vec(),im.to_vec()));
    1
}
#[no_mangle]
pub unsafe extern "C" fn advance(p: *mut Simulation, steps: usize) -> u32 {
    let s = &mut *p;
    for _ in 0..steps.min(10000) {
        s.step();
    }
    if s.y.iter().any(|x| !x.is_finite()) {
        0
    } else {
        1
    }
}
#[no_mangle]
pub unsafe extern "C" fn snapshot(p: *mut Simulation) -> *const f64 {
    (&mut *p).snapshot()
}
#[no_mangle]
pub unsafe extern "C" fn state_ptr(p: *mut Simulation) -> *const f64 {
    (&*p).y.as_ptr()
}
#[no_mangle]
pub unsafe extern "C" fn time_s(p: *mut Simulation) -> f64 {
    (&*p).time
}
#[no_mangle]
pub unsafe extern "C" fn max_fraction(p: *mut Simulation) -> f64 {
    let s = &*p;
    s.y[..s.g.n]
        .iter()
        .zip(&s.g.n0)
        .fold(0_f64, |a, (q, n)| a.max((q * NR / n).abs()))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn uniform(n: usize) -> Simulation {
        let mut s = Simulation::new(n, 0., 2., vec![3.5e11; n], 300., 0., 0., 0.);
        s.g.sponge.fill(0.);
        s
    }
    #[test]
    fn complex_vector_source_matches_sine_forcing_and_retains_quadrature() {
        let mut a = uniform(64);
        a.source.amplitude = 1e-4;
        a.source.ramp = 0.1e-6;
        a.source.angle = 0.4;
        let mut b = uniform(64);
        b.source = a.source.clone();
        let count=3*(a.g.n+1);
        let mut re=vec![0.;count];let mut im=vec![0.;count];
        for j in 0..=a.g.n {
            im[j]=a.source.amplitude*a.source_shape[j]*a.source.angle.sin();
            im[2*(a.g.n+1)+j]=a.source.amplitude*a.source_shape[j]*a.source.angle.cos();
        }
        assert_eq!(unsafe {source_phasor(&mut b,re.as_ptr(),im.as_ptr(),count)},1);
        for _ in 0..2000 { a.step();b.step(); }
        for (x,y) in a.y.iter().zip(&b.y) {assert!((x-y).abs()<1e-11);}
        // An out-of-plane cosine force does work but cannot compress B=0 plasma.
        let mut q=uniform(64);q.source.ramp=1e-7;
        for j in 0..=q.g.n {re[q.g.n+1+j]=1e-4*q.source_shape[j];}
        im.fill(0.);
        assert_eq!(unsafe {source_phasor(&mut q,re.as_ptr(),im.as_ptr(),count)},1);
        for _ in 0..2000 {q.step();}
        assert!(q.energy()>0.);
        assert!(q.y[..q.g.n].iter().all(|x|x.abs()<1e-20));
        assert!((q.energy()-q.energy_ledger[0]).abs()/q.energy()<1e-5);
    }
    #[test]
    fn tabulated_source_preserves_sign_and_energy_balance() {
        let mut s = uniform(96);
        s.source.amplitude = 1e-4;
        s.source.angle = PI/2.;
        s.source.ramp = 1e-9;
        let shape: Vec<f64> = (0..=s.g.n).map(|j| (2.*PI*j as f64/s.g.n as f64).sin()).collect();
        let count = shape.len();
        assert_eq!(unsafe {source_profile(&mut s, shape.as_ptr(), count)},1);
        assert!(s.source_shape.iter().any(|v| *v < -0.5));
        assert_eq!(unsafe {source_profile(&mut s, shape.as_ptr(), count-1)},0);
        assert_eq!(s.source_shape,shape);
        for _ in 0..1000{s.step();}
        assert!(s.energy()>0.);
        assert!((s.energy()-s.energy_ledger[0]).abs()/s.energy()<1e-4);
        assert!(s.y[..s.g.n].iter().all(|v|v.abs()<1e-15)); // transverse B=0 limit
    }
    #[test]
    fn spatial_collision_profile_removes_only_accounted_energy() {
        let mut s = uniform(96);
        let count = s.g.n + 1;
        let rates: Vec<f64> = (0..count).map(|j| 1e4 * (1. + j as f64 / count as f64)).collect();
        assert_eq!(unsafe { collision_profile(&mut s, rates.as_ptr(), count) }, 1);
        let invalid = vec![-1.; count];
        assert_eq!(unsafe { collision_profile(&mut s, invalid.as_ptr(), count) }, 0);
        assert_eq!(s.g.nu, rates);
        for j in 0..s.g.n {
            s.y[j] = 1e-5 * (2. * PI * (j as f64 + 0.5) / s.g.n as f64).cos();
        }
        let initial = s.energy();
        for _ in 0..10000 { s.step(); }
        let remaining = s.energy();
        assert!(remaining < initial);
        assert_eq!(s.energy_ledger[2], 0.); // no numerical buffer
        assert!((initial - remaining - s.energy_ledger[1]).abs() / initial < 1e-5);
        assert_eq!(unsafe { collision_profile(&mut s, rates.as_ptr(), count) }, 0);
        s.snapshot();
        for j in 0..s.g.n {
            assert_eq!(s.output[10*s.g.n+j], 0.5*(rates[j]+rates[j+1]));
        }
    }
    #[test]
    fn pulse_energy_ledger_separates_collisions_and_buffer() {
        let mut s = uniform(64);
        s.g.nu.fill(20000.);
        s.g.sponge.fill(30000.);
        s.source.amplitude = 1e-3;
        s.source.ramp = 0.1e-6;
        s.source.duration = 1e-6;
        for _ in 0..10000 {
            s.step();
        }
        assert_eq!(s.source.envelope(s.time), 0.);
        let [work, collision, buffer] = s.energy_ledger;
        assert!(work > 0. && collision > 0. && buffer > 0.);
        assert!((buffer / collision - 1.5).abs() < 1e-10);
        assert!((work - s.energy() - collision - buffer).abs() / work < 1e-5);
        let before = s.energy();
        for _ in 0..10000 {
            s.step();
        }
        assert!(s.energy() < before * 0.7);
    }
    #[test]
    fn quiet_equilibrium_stays_quiet() {
        let mut s = uniform(64);
        for _ in 0..100 {
            s.step();
        }
        assert!(s.y.iter().all(|v| *v == 0.));
    }
    #[test]
    fn gauss_law_and_charge_conservation() {
        let mut s = uniform(128);
        s.source.amplitude = 1e-3;
        s.source.ramp = 1e-9;
        for _ in 0..200 {
            s.step();
        }
        let n = s.g.n;
        s.snapshot();
        assert!(s.y[..n].iter().sum::<f64>().abs() < 1e-16);
        for j in 0..n {
            assert!(
                ((s.field[j + 1] - s.field[j]) / s.g.dz + QE * NR / EPS * s.y[j]).abs() < 1e-13
            );
        }
    }
    #[test]
    fn homogeneous_dispersion_and_collisional_ringdown() {
        for nu in [0., 2e4] {
            let mut s = uniform(128);
            s.g.nu.fill(nu);
            let k = 4. * PI;
            let amp = 1e-5;
            for j in 0..s.g.n {
                s.y[j] = amp * (k * (j as f64 + 0.5) * s.g.dz).cos();
            }
            let kd = 2. * (0.5 * k * s.g.dz).sin() / s.g.dz;
            let omega = (3.5e11 * QE * QE / (EPS * ME) + s.g.cs2 * kd * kd - nu * nu / 4.).sqrt();
            let initial_energy = s.energy();
            let steps = (10e-6 / s.dt) as usize;
            for _ in 0..steps {
                s.step();
            }
            let envelope = (-nu * s.time / 2.).exp();
            let expected = amp
                * envelope
                * ((omega * s.time).cos() + nu / (2. * omega) * (omega * s.time).sin());
            let actual = 2.
                * s.y[..s.g.n]
                    .iter()
                    .enumerate()
                    .map(|(j, v)| v * (k * (j as f64 + 0.5) * s.g.dz).cos())
                    .sum::<f64>()
                / s.g.n as f64;
            assert!(
                (actual - expected).abs() / amp < 2e-4,
                "{actual} {expected}"
            );
            if nu == 0. {
                assert!((s.energy() / initial_energy - 1.).abs() < 1e-5);
            }
        }
    }
    /// Validate the continuum dispersion as dz -> 0, independently of the
    /// discrete-k comparison above. Extract phase and amplitude from the actual
    /// density and continuity-derived density derivative, not from a fitted plot.
    #[test]
    fn bragg_scale_frequency_and_decay_benchmark() {
        let k = 4. * PI * 430e6 / 299792458.;
        let length = 8. * PI / k; // four 430 MHz Bragg wavelengths
        let nu: f64 = 2e4;
        let amplitude = 1e-5;
        let mut previous_error = f64::INFINITY;
        println!("cells points_per_wavelength frequency_Hz continuum_error_ppm amplitude_at_100us exact_amplitude");
        for n in [64, 128, 256] {
            let mut s = Simulation::new(n, 0., length, vec![3.5e11; n], 300., 0., 0., nu);
            s.g.sponge.fill(0.);
            let mode: Vec<f64> = (0..n)
                .map(|j| (k * (j as f64 + 0.5) * s.g.dz).cos())
                .collect();
            for (j, v) in mode.iter().enumerate() {
                s.y[j] = amplitude * v;
            }
            let physical_omega =
                (3.5e11 * QE * QE / (EPS * ME) + s.g.cs2 * k * k - nu * nu / 4.).sqrt();
            let kd = 2. * (0.5 * k * s.g.dz).sin() / s.g.dz;
            let discrete_omega =
                (3.5e11 * QE * QE / (EPS * ME) + s.g.cs2 * kd * kd - nu * nu / 4.).sqrt();
            let quadratures = |s: &Simulation| {
                let q = 2. * s.y[..n].iter().zip(&mode).map(|(a, b)| a * b).sum::<f64>() / n as f64;
                let qdot = -2.
                    * (0..n)
                        .map(|j| {
                            (s.g.nf[j + 1] * s.y[3 * n + 3 + j] - s.g.nf[j] * s.y[3 * n + 2 + j])
                                / (NR * s.g.dz)
                                * mode[j]
                        })
                        .sum::<f64>()
                    / n as f64;
                let p = -(qdot + nu * q / 2.) / discrete_omega;
                (p.atan2(q), q.hypot(p))
            };
            let (mut old_phase, initial_envelope) = quadratures(&s);
            let mut phase = 0.;
            let steps = (100e-6 / s.dt).round() as usize;
            for step in 0..steps {
                s.step();
                if step % 8 == 0 || step + 1 == steps {
                    let (new_phase, _) = quadratures(&s);
                    phase += (new_phase - old_phase + PI).rem_euclid(2. * PI) - PI;
                    old_phase = new_phase;
                }
            }
            let frequency = phase / (2. * PI * s.time);
            let error = (frequency / (physical_omega / (2. * PI)) - 1.).abs();
            let envelope = quadratures(&s).1 / initial_envelope;
            let exact = (-nu * s.time / 2.).exp();
            println!(
                "{n} {} {frequency:.6} {:.6} {envelope:.8} {exact:.8}",
                n / 4,
                error * 1e6
            );
            assert!(
                error < previous_error * 0.30,
                "Expected second-order spatial convergence"
            );
            assert!(
                (envelope / exact - 1.).abs() < 1e-4,
                "Incorrect collision damping"
            );
            assert!(
                (frequency / (discrete_omega / (2. * PI)) - 1.).abs() < 2e-7,
                "Incorrect temporal phase evolution"
            );
            previous_error = error;
        }
    }
    #[test]
    fn horizontal_source_does_not_compress_unmagnetized_1d_plasma() {
        let mut s = uniform(64);
        s.source.amplitude = 1e-3;
        s.source.angle = PI / 2.;
        for _ in 0..200 {
            s.step();
        }
        assert!(s.y[..64].iter().all(|v| v.abs() < 1e-18));
        assert!(s.y[64..129].iter().any(|v| v.abs() > 0.));
    }
    #[test]
    fn hermite_preserves_linear_gradient() {
        let z = [0., 10.];
        let n = [1e11, 2e11];
        let d = [1e10, 1e10];
        assert!((density_at(3., &z, &n, &d) / 1.3e11 - 1.).abs() < 1e-14);
    }
    #[test]
    fn magnetized_unforced_energy_is_conserved() {
        let mut s = uniform(64);
        s.g.bx = 3e6;
        s.g.bz = 4e6;
        for j in 1..64 {
            s.y[64 + j] = (PI * j as f64 / 64.).sin();
        }
        let e = s.energy();
        for _ in 0..1000 {
            s.step();
        }
        assert!((s.energy() / e - 1.).abs() < 2e-5);
    }
    #[test]
    fn gradient_equilibrium_and_energy_balance() {
        let n = 128;
        let ne = (0..n)
            .map(|j| 2e11 * (1. + 0.3 * (PI * (j as f64 + 0.5) / n as f64).cos()))
            .collect();
        let mut s = Simulation::new(n, 0., 2., ne, 300., 35e-6, 0.8, 0.);
        s.g.sponge.fill(0.);
        for j in 1..n {
            s.y[3 * n + 2 + j] = (PI * j as f64 / n as f64).sin();
        }
        let e = s.energy();
        for _ in 0..2000 {
            s.step();
        }
        assert!((s.energy() / e - 1.).abs() < 1e-5);
        assert!(s.y[..n].iter().sum::<f64>().abs() < 1e-16);
    }
    #[test]
    fn doubling_field_doubles_density_and_time_refinement_converges() {
        let run = |amp: f64, refine: usize| {
            let mut s = uniform(64);
            s.source.amplitude = amp;
            s.dt /= refine as f64;
            for _ in 0..500 * refine {
                s.step();
            }
            s.y
        };
        let a = run(1e-3, 1);
        let b = run(2e-3, 1);
        let fine = run(1e-3, 2);
        let norm = a[..64].iter().map(|x| x * x).sum::<f64>().sqrt();
        for i in 0..a.len() {
            assert!((2. * a[i] - b[i]).abs() < 1e-12);
        }
        let error = a[..64]
            .iter()
            .zip(&fine)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
            / norm;
        assert!(error < 1e-4, "{error}");
    }
}

pub mod ten_moment;
