#!/usr/bin/env python3
"""
Fill a circular on-sky region with smaller circles and plot them
in a correct TAN WCS projection.

Now with:
- configurable spacing (circles can overlap strongly),
- allowing circles to spill slightly outside the big circle,
- fill_circle_with_circles returns a pandas DataFrame.

Requires: numpy, astropy, matplotlib, pandas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.wcs import WCS


# ----------------------------------------------------------------------
# Geometry: fill big circle with small circles in tangent plane
# ----------------------------------------------------------------------
def fill_circle_with_circles(
    ra0_deg,
    dec0_deg,
    big_radius_deg,
    small_radius_deg,
    spacing_factor=0.7,
    spill_factor=0.3,
    require_full_inside=False,
):
    """
    Generate RA, Dec positions for small circles that hex-pack a larger circle.

    Returns a pandas DataFrame with columns:
        RAdeg, DEdeg, x_offset_deg, y_offset_deg
    """
    centre = SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg, frame="icrs")
    off_frame = SkyOffsetFrame(origin=centre)

    R = float(big_radius_deg)
    r = float(small_radius_deg)

    if r <= 0 or R <= 0:
        raise ValueError("Both big_radius_deg and small_radius_deg must be > 0.")
    if r > R:
        raise ValueError("small_radius_deg must be <= big_radius_deg.")

    spacing_factor = float(spacing_factor)
    if not (0 < spacing_factor <= 1.0):
        raise ValueError("spacing_factor must be in (0, 1].")

    spill_factor = float(spill_factor)
    if spill_factor < 0:
        raise ValueError("spill_factor must be >= 0.")

    # Hexagonal spacing
    dx = 2.0 * r * spacing_factor
    dy = np.sqrt(3.0) * r * spacing_factor

    # Cover slightly more than the circle to allow spill/overlap
    y_vals = np.arange(-R - r, R + r + dy, dy)
    x_vals = np.arange(-R - r, R + r + dx, dx)

    if require_full_inside:
        max_r2 = (R - r) ** 2
    else:
        max_r2 = (R + spill_factor * r) ** 2

    xs = []
    ys = []

    for j, y in enumerate(y_vals):
        row_shift = 0.0 if (j % 2 == 0) else dx / 2.0
        for x0 in x_vals:
            x = x0 + row_shift
            if x * x + y * y <= max_r2:
                xs.append(x)
                ys.append(y)

    xs = np.asarray(xs)
    ys = np.asarray(ys)

    # Convert offsets back to sky coordinates
    offs = SkyCoord(
        lon=xs * u.deg,
        lat=ys * u.deg,
        frame=off_frame
    )
    sky = offs.icrs

    # Build DataFrame
    df = pd.DataFrame({
        "GAIA_RA": sky.ra.deg,
        "GAIA_DEC": sky.dec.deg,
        "x_offset_deg": xs,
        "y_offset_deg": ys,
    })

    return df

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, SkyOffsetFrame

def three_circles_triangular(
    ra0_deg,
    dec0_deg,
    small_radius_deg,
    orientation_deg=0.0,
    ra_offset=None,
    dec_offset=None,
):
    """
    Generate RA, Dec positions for three circle centres in an equilateral-triangle
    configuration such that the input (ra0_deg, dec0_deg) lies on the circumference
    of all three circles (i.e. it is a common intersection point).

    Returns a pandas DataFrame with columns:
        GAIA_RA, GAIA_DEC, x_offset_deg, y_offset_deg

    Notes
    -----
    - The three circle centres are placed on a ring of radius = small_radius_deg
      around the centre point, at position angles separated by 120 degrees.
    - In the local tangent-plane (SkyOffsetFrame), x_offset_deg ~ "east" and
      y_offset_deg ~ "north" offsets in degrees.
    """
    
    if ra_offset:
        ra0_deg = ra0_deg + ra_offset
    elif dec_offset:
        dec0_deg = dec0_deg + dec_offset
    
    centre = SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg, frame="icrs")
    off_frame = SkyOffsetFrame(origin=centre)

    r = float(small_radius_deg)
    if r <= 0:
        raise ValueError("small_radius_deg must be > 0.")

    theta0 = np.deg2rad(float(orientation_deg))
    thetas = theta0 + np.deg2rad(np.array([0.0, 120.0, 240.0]))

    xs = r * np.cos(thetas)
    ys = r * np.sin(thetas)
    offs = SkyCoord(lon=xs * u.deg, lat=ys * u.deg, frame=off_frame)
    sky = offs.icrs

    df = pd.DataFrame({
        "GAIA_RA": sky.ra.deg,
        "GAIA_DEC": sky.dec.deg,
        "x_offset_deg": xs,
        "y_offset_deg": ys,
    })
    return df



# ----------------------------------------------------------------------
# WCS creation
# ----------------------------------------------------------------------
def make_tan_wcs(ra0_deg,
                 dec0_deg,
                 fov_deg,
                 pixscale_arcmin=2.0):
    """
    Create a TAN WCS centered on (ra0, dec0) with a given field-of-view.
    """
    pixscale_deg = pixscale_arcmin / 60.0
    nx = int(np.ceil(fov_deg / pixscale_deg))
    ny = int(np.ceil(fov_deg / pixscale_deg))

    w = WCS(naxis=2)
    w.wcs.crval = [ra0_deg, dec0_deg]
    w.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cd = np.array([[-pixscale_deg, 0.0],
                         [0.0,            pixscale_deg]])

    return w, nx, ny


# ----------------------------------------------------------------------
# Plotting in WCS using DataFrame
# ----------------------------------------------------------------------
def plot_circles_wcs(
    df_circles,
    ra0_deg,
    dec0_deg,
    big_radius_deg,
    small_radius_deg,
    pixscale_arcmin=2.0,
    min_fov_factor=2.1,
    min_fov_deg=1.0,
    circle_npts=180,
    draw_small_circles=True,
):
    """
    Plot the circles returned by fill_circle_with_circles() in a TAN WCS.
    """
    ra_small = df_circles["GAIA_RA"].to_numpy()
    dec_small = df_circles["GAIA_DEC"].to_numpy()

    center_sc = SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg, frame="icrs")

    effective_R = big_radius_deg
    auto_fov = max(min_fov_factor * effective_R, float(min_fov_deg))

    w, nx, ny = make_tan_wcs(
        ra0_deg, dec0_deg, auto_fov, pixscale_arcmin=pixscale_arcmin
    )

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection=w)

    # Draw big circle
    theta = np.linspace(0, 2 * np.pi, 720)
    big_gc = center_sc.directional_offset_by(theta * u.rad,
                                             big_radius_deg * u.deg)
    ax.plot(big_gc.ra.deg, big_gc.dec.deg,
            transform=ax.get_transform('world'),
            color='k', lw=1.5)

    # Plot centres
    ax.scatter(
        ra_small, dec_small,
        transform=ax.get_transform('world'),
        s=4, alpha=0.6, color="tab:blue"
    )

    # Draw small circle outlines
    if draw_small_circles:
        phi = np.linspace(0, 2 * np.pi, circle_npts)
        for ra_c, dec_c in zip(ra_small, dec_small):
            c_sc = SkyCoord(ra_c * u.deg, dec_c * u.deg, frame="icrs")
            small_gc = c_sc.directional_offset_by(phi * u.rad,
                                                  small_radius_deg * u.deg)
            ax.plot(small_gc.ra.deg,
                    small_gc.dec.deg,
                    transform=ax.get_transform('world'),
                    lw=0.4, color="gray", alpha=0.7)

    # WCS limits
    half = auto_fov / 2.0
    ra_min, ra_max = ra0_deg - half, ra0_deg + half
    dec_min, dec_max = dec0_deg - half, dec0_deg + half
    box = SkyCoord([ra_min, ra_max, ra_min, ra_max] * u.deg,
                   [dec_min, dec_max, dec_max, dec_min] * u.deg,
                   frame="icrs")
    xp, yp = w.world_to_pixel(box)
    #ax.set_xlim(np.min(xp), np.max(xp))
    #ax.set_ylim(np.min(yp), np.max(yp))

    ax.set_aspect(1.0)
    ax.set_xlabel("RA (ICRS)")
    ax.set_ylabel("Dec (ICRS)")
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()

    return fig, ax, w


# ----------------------------------------------------------------------
# Example usage
# ----------------------------------------------------------------------
if __name__ == "__main__":
    ra0 = 120.0
    dec0 = 30.0
    big_R = 1.0
    small_r = 0.1

    df = fill_circle_with_circles(
        ra0,
        dec0,
        big_R,
        small_r,
        spacing_factor=0.7,
        spill_factor=0.3
    )

    plot_circles_wcs(
        df,
        ra0,
        dec0,
        big_R,
        small_r
    )
