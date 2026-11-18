# Equilibrium plotting

Create equilibrium and magnetic-coordinate figures from saved cuMES and Boozer
results without building or importing the solver. `equilibrium-io` owns the
native and Boozer file readers. This package owns reconstruction, coordinate
grids, coil overlays and rendering.

```sh
python -m pip install "equilibrium-io @ git+https://github.com/12ff54e/equilibrium-io.git@v0.3.0"
python -m pip install .
equilibrium-plot --state equilibrium.bin --output-dir figures
equilibrium-plot --boozer-state boozer.bin --coordinate-system boozer \
  --coordinate-only --output-dir figures
```

`python -m equilibrium_plot` provides the same CLI. Install `.[hdf5]` for HDF5
input and `.[pyvista]` for the optional VTK renderer, selected by
`--backend pyvista`. Matplotlib is the default renderer. Native NetCDF input
uses SciPy, already required for interpolation. Files use cuMES schemas;
external VMEC wout files are not supported.

Native results produce perspective, top, combined and R-Z slice figures.
`--field-lines` adds LCFS field lines; `--coils PATH` adds coil geometry.
`--out BASE` prefixes output names, while `--output-dir DIR` writes the
standard names directly into that directory. PEST reconstruction needs native
state; Boozer plots use saved Boozer v7/v8 data.

For custom figures, pass an explicit immutable configuration:

```python
from dataclasses import replace
from equilibrium_plot.application import main
from equilibrium_plot.defaults import FIGURE_PARAMETERS

settings = replace(FIGURE_PARAMETERS, save_dpi=200, field_colormap="viridis")
main(settings, ["--state", "equilibrium.bin", "--output-dir", "figures"])
```

The public `FigureConfig`, `CoordinateGrid`, `make_boozer_grid`,
`make_pest_grid`, and figure functions are available from `equilibrium_plot`.
The plotting adapters translate reader errors into CLI diagnostics. The reader
API in `equilibrium_io` raises ordinary Python exceptions.

Run the analytic geometry, coordinate mapping and rendering tests against an
installed wheel:

```sh
python -m pip install build
python -m build --wheel
python -m pip install dist/*.whl
python -m unittest discover -s tests
```

This initial release keeps the cuMES numerical reconstruction and figure
defaults. It requires `equilibrium-io>=0.3.0`. When using local wheels, install
the I/O wheel before the plotting wheel.
