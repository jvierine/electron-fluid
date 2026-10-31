use arecibo_fluid::coupling::{FieldCoupledFluid,YeeFluidCoupler};
use arecibo_fluid::maxwell_coupling::CouplingGeometry;

#[test]
fn independent_ampere_host_reproduces_plasma_oscillation() {
    let n=8;let dt=0.01;let wp=0.2_f64;
    let widths=vec![0.1;n];let plasma=vec![vec![wp;n];3];
    let zero=vec![vec![0.;n];3];let rho=vec![vec![1.;n];3];let buffer=vec![0.;n];
    let g=CouplingGeometry{widths:&widths,kx:0.,gyro:[0.;3],plasma_frequency:&plasma,
        collision_rate:&zero,sqrt_density:&rho,thermal_speed:0.01,density_buffer:&buffer};
    let theta=2.*(wp*dt/2.).asin();let mut j=vec![vec![[0.;2];n];3];
    j[0].fill([-(theta/2.).sin(),0.]);
    let mut fluid=YeeFluidCoupler::new(g,dt,j,vec![[0.;2];n]).unwrap();
    let mut e=vec![vec![[0.;2];n];3];e[0].fill([1.,0.]);
    for step in 1..10001 {
        fluid.advance(&e).unwrap();
        // An independent uniform-field Ampere host, with no FDTD dependency.
        for i in 0..n { e[0][i][0]-=dt*wp*fluid.current()[0][i][0]; }
        assert!((e[0][0][0]-(step as f64*theta).cos()).abs()<2e-12);
    }
    assert!(fluid.density().iter().flatten().all(|v|*v==0.));
}
#[test]
fn nonuniform_periodic_continuity_conserves_total_charge() {
    let n=64;let dt=0.001;
    let widths:Vec<_>=(0..n).map(|i|0.1+0.01*(i as f64).sin()).collect();
    let plasma=vec![vec![0.2;n];3];let zero=vec![vec![0.;n];3];
    let rho:Vec<Vec<f64>>=(0..3).map(|q|(0..n).map(|i|1.+0.1*(i as f64+q as f64/2.).sin()).collect()).collect();
    let buffer=vec![0.;n];let g=CouplingGeometry{widths:&widths,kx:0.,gyro:[0.;3],
        plasma_frequency:&plasma,collision_rate:&zero,sqrt_density:&rho,
        thermal_speed:0.01,density_buffer:&buffer};
    let mut j=vec![vec![[0.;2];n];3];
    for i in 0..n { j[2][i]=[(i as f64).cos(),(i as f64).sin()]; }
    let mut f=YeeFluidCoupler::new(g,dt,j,vec![[0.;2];n]).unwrap();
    let e=vec![vec![[0.;2];n];3];for _ in 0..100 { f.advance(&e).unwrap(); }
    for q in 0..2 {
        let integral:f64=(0..n).map(|i| f.density()[i][q]*rho[0][i]*(widths[i]+widths[(i+n-1)%n])/2.).sum();
        assert!(integral.abs()<1e-15,"charge integral {integral}");
    }
    assert!(f.density().iter().any(|v|v[0].abs()>1e-4));
    let saved=f.density().to_vec();
    let mut bad=e.clone();bad[0][0][0]=f64::NAN;
    assert!(f.advance(&bad).is_err());assert_eq!(saved,f.density());
}
#[test]
fn uniform_drag_matches_exact_trapezoidal_decay() {
    let n=5;let dt=0.01;let nu=0.3;
    let widths=vec![1.;n];let plasma=vec![vec![0.2;n];3];let rates=vec![vec![nu;n];3];
    let rho=vec![vec![1.;n];3];let buffer=vec![0.;n];
    let g=CouplingGeometry{widths:&widths,kx:0.,gyro:[0.;3],plasma_frequency:&plasma,
        collision_rate:&rates,sqrt_density:&rho,thermal_speed:0.01,density_buffer:&buffer};
    let mut f=YeeFluidCoupler::new(g,dt,vec![vec![[1.,0.];n];3],vec![[0.;2];n]).unwrap();
    for _ in 0..1000 { f.advance(&vec![vec![[0.;2];n];3]).unwrap(); }
    let exact=((1.-nu*dt/2.)/(1.+nu*dt/2.)).powi(1000);
    assert!((f.current()[0][0][0]-exact).abs()<1e-14);
}
