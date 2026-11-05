# Changelog

## 0.1.0 — 2026-11-05

- Extract equilibrium, PEST and Boozer reconstruction and plotting from cuMES
  into the `equilibrium_plot` library and `equilibrium-plot` command.
- Retain Matplotlib figures, coil overlays, field lines, figure settings, and
  the optional PyVista backend. Keep the existing cuMES launcher compatible.
- Read native and Boozer results through `equilibrium-io>=0.3.0`. Supported
  files use cuMES schemas; external VMEC wout files are not supported.
