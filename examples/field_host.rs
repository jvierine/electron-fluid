//! Minimal external host: uniform Maxwell-Ampere plasma oscillation.
//! Run: cargo run --example field_host
use arecibo_fluid::coupling::{FieldCoupledFluid,YeeFluidCoupler};
use arecibo_fluid::maxwell_coupling::CouplingGeometry;
fn main()->Result<(),Box<dyn std::error::Error>> {
    let n=8;let widths=vec![1.;n];let wp=vec![vec![0.1;n];3];
    let zero=vec![vec![0.;n];3];let rho=vec![vec![1.;n];3];let buffer=vec![0.;n];
    let geometry=CouplingGeometry{widths:&widths,kx:0.,gyro:[0.;3],plasma_frequency:&wp,
        collision_rate:&zero,sqrt_density:&rho,thermal_speed:0.01,density_buffer:&buffer};
    let dt=0.01;let mut fluid=YeeFluidCoupler::new(geometry,dt,
        vec![vec![[0.;2];n];3],vec![[0.;2];n])?;
    let mut electric=vec![vec![[0.;2];n];3];electric[0].fill([1.,0.]);
    for _ in 0..1000 {
        fluid.advance(&electric)?;
        for q in 0..3 {for i in 0..n {for r in 0..2 {
            electric[q][i][r]-=dt*wp[q][i]*fluid.current()[q][i][r];
        }}}
    }
    println!("E_x after 1000 coupled steps: {} V/m",electric[0][0][0]);Ok(())
}
