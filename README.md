# Electron fluid

Rust solver for linear electron momentum and continuity, including pressure,
magnetization and collisional drag. It runs natively or as WebAssembly.

Two interfaces are available:

- `Simulation`: standalone vertical fluid with an electrostatic Poisson field.
- `FieldCoupledFluid` / `YeeFluidCoupler`: accepts an external model's electric
  field and returns current and density on a nonuniform Yee grid. It has no
  Maxwell-library dependency and performs no additional Poisson solve.

```sh
cargo test --release
cargo run --release --example field_host
cargo doc --no-deps
```

The crate is named `electron-fluid`; the Rust library name remains
`arecibo_fluid` for compatibility. See [the coupling contract](docs/coupling.md)
for units, staggering, boundaries and integration order. A complete coupled
application is at https://github.com/jvierine/arecibo-hf-conversion.

The linear closure includes electron-ion and electron-neutral drag. It does
not include Landau damping, mobile ions, trapping or nonlinear instability.
Tests cover analytic dispersion, ringdown, charge and energy accounting,
and an independent field host. Numerical validation does not validate the
missing kinetic physics. Earlier experiment-specific notes are retained in
[the development record](docs/arecibo_history.md).
