"""Construct magnetic-coordinate plot grids from equilibrium-io data.

The mixed (s, theta_b, alpha) grid preserves the source field-period angle.
The displacement nu is in physical radians: alpha_b = alpha + nfp * nu.
PEST plots use native geometry and physical lambda.
"""

from dataclasses import dataclass

import numpy as np

from equilibrium_io.boozer import (
    ASYMMETRIC_FAMILY_NAMES, ASYMMETRIC_FOURIER_CONVENTION, BoozerData,
    COORDINATE_CONVENTION, FAMILY_NAMES, FOURIER_CONVENTION, MAGICS,
    MAGNETIC_CONVENTION, PERIOD, PROFILE_CONVENTION, PROFILE_UNITS,
    SCHEMAS, TRACING_CONVENTION, load_boozer as _load_boozer,
)


def load_boozer(path):
    """Keep command-line error handling at the plotting boundary."""
    try:
        return _load_boozer(path)
    except (ValueError, OSError) as error:
        raise SystemExit(str(error)) from error


@dataclass(frozen=True)
class CoordinateGrid:
    """Geometry and |B| on a regular magnetic-coordinate angular grid."""

    coordinate: str
    angle_label: str
    toroidal_label: str
    s: np.ndarray
    theta: np.ndarray
    zeta: np.ndarray
    r: np.ndarray
    z: np.ndarray
    b: np.ndarray


def _synthesize_mixed_fields(data, theta, alpha):
    """Synthesize R, Z, and physical-angle nu on a uniform mixed grid."""
    mt = data.mode_m[:, None] * theta[None, :]
    nz = data.mode_n[:, None] * alpha[None, :]
    cos_mt, sin_mt = np.cos(mt), np.sin(mt)
    cos_nz, sin_nz = np.cos(nz), np.sin(nz)

    def combine(first, first_t, first_z, second, second_t, second_z):
        return (np.einsum("sm,mz,mt->szt", data.families[first], first_z,
                          first_t, optimize=True) +
                np.einsum("sm,mz,mt->szt", data.families[second], second_z,
                          second_t, optimize=True))

    r = combine("rmncc", cos_mt, cos_nz, "rmnss", sin_mt, sin_nz)
    z = combine("zmnsc", sin_mt, cos_nz, "zmncs", cos_mt, sin_nz)
    nu = combine("numnsc", sin_mt, cos_nz, "numncs", cos_mt, sin_nz)
    if "rmnsc" in data.families:
        r += combine("rmnsc", sin_mt, cos_nz, "rmncs", cos_mt, sin_nz)
        z += combine("zmncc", cos_mt, cos_nz, "zmnss", sin_mt, sin_nz)
        nu += combine("numncc", cos_mt, cos_nz, "numnss", sin_mt, sin_nz)
    return r, z, nu


def _resample_uniform_axis(values, target_count, axis):
    """Linearly resample one uniformly periodic array axis."""
    source_count = values.shape[axis]
    if source_count == target_count:
        return values.copy()
    position = np.arange(target_count) * source_count / target_count
    lower = np.floor(position).astype(int) % source_count
    upper = (lower + 1) % source_count
    fraction = position - np.floor(position)
    shape = [1] * values.ndim
    shape[axis] = target_count
    fraction = fraction.reshape(shape)
    return ((1.0 - fraction) * np.take(values, lower, axis=axis) +
            fraction * np.take(values, upper, axis=axis))


def _resample_field(data, ntheta, nzeta):
    field = _resample_uniform_axis(data.b, ntheta, axis=-1)
    return _resample_uniform_axis(field, nzeta, axis=-2)


def _periodic_interp(target, source, values):
    source = np.mod(source, PERIOD)
    order = np.argsort(source)
    source = source[order]
    values = np.asarray(values)[order]
    if np.min(np.diff(source)) <= 64.0 * np.finfo(float).eps:
        raise ValueError("magnetic poloidal map contains duplicate angles")
    return np.interp(target, source, values, period=PERIOD)


def _remap_mixed_to_uniform_boozer(
        source_alpha, target_alpha_b, nfp, nu, *fields):
    """Invert ``alpha_b=alpha+nfp*nu`` at fixed ``(s, theta_b)``."""
    source_alpha = np.asarray(source_alpha, dtype=float)
    target_alpha_b = np.asarray(target_alpha_b, dtype=float)
    nu = np.asarray(nu, dtype=float)
    fields = tuple(np.asarray(field, dtype=float) for field in fields)
    if nfp < 1 or source_alpha.ndim != 1 or target_alpha_b.ndim != 1 or \
            source_alpha.size < 1 or target_alpha_b.size < 1 or nu.ndim != 3 or \
            nu.shape[1] != source_alpha.size or any(
                field.shape != nu.shape for field in fields):
        raise ValueError(
            "mixed-to-Boozer remap arrays have inconsistent extents")
    if source_alpha.size == 1:
        return [np.repeat(field, target_alpha_b.size, axis=1)
                for field in fields]

    outputs = [np.empty((nu.shape[0], target_alpha_b.size, nu.shape[2]))
               for _ in fields]
    for surface in range(nu.shape[0]):
        for poloidal in range(nu.shape[2]):
            mapped = source_alpha + nfp * nu[surface, :, poloidal]
            closed = np.unwrap(np.append(mapped, mapped[0] + PERIOD))
            if np.any(np.diff(closed) <= 0.0) or not np.all(
                    np.isfinite(closed)):
                raise ValueError(
                    "source-alpha to Boozer-alpha map is not orientation "
                    "preserving")
            for output, field in zip(outputs, fields):
                output[surface, :, poloidal] = _periodic_interp(
                    target_alpha_b, mapped, field[surface, :, poloidal])
    return outputs


def make_boozer_grid(data, ntheta=None, nzeta=None):
    """Return geometry and |B| on a grid uniform in both Boozer angles."""
    ntheta = data.ntheta if ntheta is None else int(ntheta)
    nzeta = data.nzeta if nzeta is None else int(nzeta)
    if ntheta < 4 or nzeta < 1:
        raise ValueError("coordinate-grid angular dimensions are invalid")
    theta = PERIOD * np.arange(ntheta) / ntheta
    alpha = PERIOD * np.arange(nzeta) / nzeta
    r_mixed, z_mixed, nu = _synthesize_mixed_fields(data, theta, alpha)
    b_mixed = _resample_field(data, ntheta, nzeta)
    r, z, b = _remap_mixed_to_uniform_boozer(
        alpha, alpha, data.nfp, nu, r_mixed, z_mixed, b_mixed)
    return CoordinateGrid("Boozer", r"$\theta_b$", r"$\alpha_b$",
                          data.s.copy(), theta, alpha, r, z, b)


def make_pest_grid(s, theta, zeta, r, z, b, lambda_physical):
    """Remap native-state fields from theta onto a regular PEST grid.

    ``lambda_physical`` is the normalized displacement reconstructed from the
    native cuMES lambda coefficients, and ``theta_p = theta + lambda``. Every
    real-space input has shape ``[nonaxis_surface, zeta, theta]``.
    """
    s = np.asarray(s, dtype=float)
    theta = np.asarray(theta, dtype=float)
    zeta = np.asarray(zeta, dtype=float)
    arrays = [np.asarray(value, dtype=float)
              for value in (r, z, b, lambda_physical)]
    r, z, b, lambda_physical = arrays
    expected = (s.size, zeta.size, theta.size)
    if theta.ndim != 1 or zeta.ndim != 1 or s.ndim != 1 or \
            theta.size < 4 or zeta.size < 1 or any(
                value.shape != expected for value in arrays):
        raise ValueError("native PEST plotting arrays have inconsistent extents")
    if any(not np.all(np.isfinite(value))
           for value in [s, theta, zeta] + arrays):
        raise ValueError("native PEST plotting arrays contain non-finite data")
    if np.any(b <= 0.0) or np.any(np.diff(s) <= 0.0):
        raise ValueError("native PEST plotting grid is nonphysical")

    r_p = np.empty_like(r)
    z_p = np.empty_like(z)
    b_p = np.empty_like(b)
    target = theta
    for surface in range(s.size):
        for toroidal in range(zeta.size):
            theta_p = target + lambda_physical[surface, toroidal]
            closed = np.unwrap(np.append(theta_p, theta_p[0] + PERIOD))
            if np.any(np.diff(closed) <= 0.0):
                raise ValueError(
                    "native theta-to-PEST map is not orientation preserving")
            r_p[surface, toroidal] = _periodic_interp(
                target, theta_p, r[surface, toroidal])
            z_p[surface, toroidal] = _periodic_interp(
                target, theta_p, z[surface, toroidal])
            b_p[surface, toroidal] = _periodic_interp(
                target, theta_p, b[surface, toroidal])
    return CoordinateGrid("PEST", r"$\theta_p$", r"$\alpha$", s.copy(),
                          target, zeta, r_p, z_p, b_p)


def interpolate_zeta(values, zeta):
    """Periodically interpolate an ``[..., zeta, theta]`` array."""
    nzeta = values.shape[-2]
    grid_index = (zeta % PERIOD) * nzeta / PERIOD
    lower = int(np.floor(grid_index)) % nzeta
    fraction = grid_index - np.floor(grid_index)
    upper = (lower + 1) % nzeta
    return (1.0 - fraction) * values[..., lower, :] + \
        fraction * values[..., upper, :]
