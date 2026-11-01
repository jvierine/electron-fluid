//! Nonlinear electron moments, independent of mesh and Maxwell implementation.
//!
//! State is [n, n*u_x, n*u_y, n*u_z, M_xx, M_xy, M_xz, M_yy, M_yz, M_zz],
//! where M_ij = n*u_i*u_j + P_ij/m. The caller may use consistent normalized
//! units. Gaussian (zero central third moment) closure: heat flux is zero.
//! Collision source is number-conserving BGK relaxation toward a stationary
//! Maxwellian bath, NOT a Coulomb/Fokker--Planck collision operator.
pub type Moments = [f64; 10];
pub const PAIRS: [(usize, usize); 6] = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)];
pub fn matrix(q: &Moments) -> [[f64; 3]; 3] {
    let mut m = [[0.; 3]; 3];
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        m[i][j] = q[4 + k];
        m[j][i] = q[4 + k];
    }
    m
}
pub fn pressure(q: &Moments) -> [[f64; 3]; 3] {
    let mut p = matrix(q);
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] -= q[1 + i] * q[1 + j] / q[0];
        }
    }
    p
}
pub fn admissible(q: &Moments) -> bool {
    if !q.iter().all(|v| v.is_finite()) || q[0] <= 0. {
        return false;
    }
    let p = pressure(q);
    p[0][0] > 0.
        && p[0][0] * p[1][1] - p[0][1] * p[0][1] > 0.
        && p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[1][2])
            - p[0][1] * (p[0][1] * p[2][2] - p[1][2] * p[0][2])
            + p[0][2] * (p[0][1] * p[1][2] - p[1][1] * p[0][2])
            > 0.
}
/// Conservative flux along direction a. No small-amplitude expansion.
pub fn flux(q: &Moments, a: usize) -> Moments {
    let n = q[0];
    let u = [q[1] / n, q[2] / n, q[3] / n];
    let m = matrix(q);
    let mut f = [0.; 10];
    f[0] = q[1 + a];
    for i in 0..3 {
        f[1 + i] = m[i][a];
    }
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        f[4 + k] = u[a] * m[i][j] + u[i] * m[j][a] + u[j] * m[i][a] - 2. * n * u[i] * u[j] * u[a];
    }
    f
}
/// Acceleration is the electric force per unit electron mass; gyro is e B/m
/// for electrons, so magnetic acceleration is gyro cross u. nu is BGK rate;
/// theta is the stationary bath variance k_B T_b/m. Retains full wave B.
pub fn source(q: &Moments, acceleration: [f64; 3], gyro: [f64; 3], nu: f64, theta: f64) -> Moments {
    let m = matrix(q);
    let mut s = [0.; 10];
    let l = [
        [0., -gyro[2], gyro[1]],
        [gyro[2], 0., -gyro[0]],
        [-gyro[1], gyro[0], 0.],
    ];
    for i in 0..3 {
        s[1 + i] = q[0] * acceleration[i] - nu * q[1 + i];
        for j in 0..3 {
            s[1 + i] += l[i][j] * q[1 + j];
        }
    }
    for (k, &(i, j)) in PAIRS.iter().enumerate() {
        s[4 + k] = acceleration[i] * q[1 + j] + acceleration[j] * q[1 + i]
            - nu * (m[i][j] - if i == j { q[0] * theta } else { 0. });
        for a in 0..3 {
            s[4 + k] += l[i][a] * m[a][j] + l[j][a] * m[i][a];
        }
    }
    s
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gaussian_flux_and_pressure() {
        let q = [2., 2., 4., 6., 4., 4., 6., 12., 12., 24.];
        assert!(admissible(&q));
        let p = pressure(&q);
        assert_eq!(p, [[2., 0., 0.], [0., 4., 0.], [0., 0., 6.]]);
        let f = flux(&q, 2);
        assert_eq!(f[0], 6.);
        assert_eq!(f[4], 12.);
        assert_eq!(f[9], 108.);
    }
    #[test]
    fn force_work_and_gyro_energy() {
        let q = [2., 2., 4., 6., 4., 4., 6., 12., 12., 24.];
        let a = [0.3, -0.2, 0.7];
        let s = source(&q, a, [0.4, 0.5, 0.6], 0., 1.);
        assert!(
            (0.5 * (s[4] + s[7] + s[9]) - (a[0] * q[1] + a[1] * q[2] + a[2] * q[3])).abs() < 1e-13
        );
    }
    #[test]
    fn bgk_bath_has_correct_moments() {
        let q = [2., 2., 4., 6., 4., 4., 6., 12., 12., 24.];
        let s = source(&q, [0.; 3], [0.; 3], 0.2, 1.);
        assert_eq!(s[0], 0.);
        assert_eq!(s[1], -0.4);
        assert!((0.5 * (s[4] + s[7] + s[9]) + 0.1 * (40. - 6.)).abs() < 1e-14);
    }
}
