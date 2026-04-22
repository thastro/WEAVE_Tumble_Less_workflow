import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np


def plot_coords(df, ra_col='RAdeg', de_col='DEdeg', figsize=(10, 6),
                projection='aitoff', coordinates='radec', hue=None,
                cmap='viridis', markersize=0.7, ax=None, color='k', alpha=1.0):

    if ax is None:
        ax = plt.gca()

    sky_coords = SkyCoord(ra=df[ra_col]*u.deg, dec=df[de_col]*u.deg, frame='icrs', unit=u.deg)
    ax.grid(True)

    if projection is not None:

        if coordinates == 'galactic':
            sky_coords = sky_coords.galactic
            sky_coord_x = sky_coords.l.wrap_at('180d').radian
            sky_coord_y = sky_coords.b.radian
            ax.set_xlabel('gal. longitude')
            ax.set_ylabel('gal. latitude')

        elif coordinates == 'radec':
            sky_coord_x = sky_coords.ra.wrap_at('180d').radian
            sky_coord_y = sky_coords.dec.radian
            ax.set_xlabel('RA')
            ax.set_ylabel('DEC')

        if hue is None:
            ax.scatter(sky_coord_x, sky_coord_y, marker='.', s=markersize, color=color,
                       alpha=alpha)

            return ax

        else:
            plot = ax.scatter(sky_coord_x, sky_coord_y, marker='.', s=markersize, c=df[hue],
                              cmap=cmap, alpha=alpha)

            return plot
    else:
        sky_coord_x = sky_coords.ra
        sky_coord_y = sky_coords.dec

        if hue is None:
            ax.scatter(sky_coord_x, sky_coord_y, marker='.', s=markersize, color=color,
                       alpha=alpha)

            return ax

        else:
            plot = ax.scatter(sky_coord_x, sky_coord_y, marker='.', s=markersize, c=df[hue],
                              cmap=cmap, alpha=alpha)

        return plot


def plot_twilight_times(twilight_data, color='darkslategray', ax=None):

    if ax is None:
        ax = plt.gca()

    for _, twilights in twilight_data.iterrows():
        if twilights['morning'][0] and twilights['morning'][1]:
            # Convert astropy Time objects to datetime if necessary
            start_time = twilights['morning'][0].datetime
            end_time = twilights['morning'][1].datetime

            # Plot morning twilight as a shaded region
            ax.axvspan(start_time, end_time, color=color, alpha=0.3)

        if twilights['evening'][0] and twilights['evening'][1]:
            start_time = twilights['evening'][0].datetime
            end_time = twilights['evening'][1].datetime

            # Plot evening twilight as a shaded region
            ax.axvspan(start_time, end_time, color=color, alpha=0.3)

    return ax


def plot_altitudes(altitudes, times, alt_lim=30.0, title='month',
                   linewidth=1.0, mean=False, ax=None):

    if ax is None:
        ax = plt.gca()

    ax.axhline(alt_lim, linewidth=0.5, alpha=0.8, ls='--', color='k')

    ax.set_xlabel(rf'{title}')
    ax.set_ylabel('Altitude deg')
    ax.tick_params(axis='x', labelrotation=45)
    ax.set_ylim(0, 90)

    if mean:
        ax.plot(times.datetime, altitudes[1], color='tab:green', linewidth=linewidth,
                label='mean-alt')
        ax.legend()
        return ax
    # ax.plot(times.datetime, altitudes, color='tab:blue', linewidth=linewidth, alpha=0.5)
    for i in range(0, len(altitudes)):
        ax.plot(times.datetime, altitudes[i], color='tab:blue', linewidth=linewidth, alpha=0.5)

    return ax


def wcs_plot(df, ra_col='RAdeg', de_col='DEdeg',
             color='grey', markersize=0.7, alpha=1.0,
             marker='o', fill=False, kepler_fov=False,
             label=None, c=None, cmap='viridis', facecolor=None,
             edgecolor=None, lw=None, zorder=1, ax=None):

    from scipy.spatial import ConvexHull
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    if ax is None:
        ax = plt.gca()

    # Get coordinate types from axis WCS
    wcs = ax.wcs
    ctype1 = wcs.wcs.ctype[0]
    ctype2 = wcs.wcs.ctype[1]

    # Check if coordinate types are Galactic
    is_galactic = ctype1.startswith('GLON') and ctype2.startswith('GLAT')

    if kepler_fov:
        label = 'Kepler FOV'
        for i in np.arange(0, len(df), 4):
            subset = df.iloc[i:i+4][[ra_col, de_col]]
            points = subset.values

            # Transform points if necessary
            if is_galactic:
                skycoord = SkyCoord(ra=points[:, 0]*u.degree, dec=points[:, 1]*u.degree,
                                    frame='icrs')
                points = np.column_stack((skycoord.galactic.l.degree,
                                          skycoord.galactic.b.degree))

            hull = ConvexHull(points)
            vertices = points[hull.vertices]

            if fill:
                # Fill the entire convex hull area
                ax.fill(vertices[:, 0], vertices[:, 1], color=color, alpha=alpha,
                        transform=ax.get_transform('world'), label=label, zorder=zorder)
            else:
                # Plot only the outline of the convex hull
                ax.plot(np.append(vertices[:, 0], vertices[0, 0]),
                        np.append(vertices[:, 1], vertices[0, 1]),
                        color=color, alpha=alpha,
                        transform=ax.get_transform('world'),
                        label=label, zorder=zorder)

            label = None

        return ax

    # Transform input coordinates if necessary
    if is_galactic:
        skycoord = SkyCoord(ra=df[ra_col].values*u.degree, dec=df[de_col].values*u.degree,
                            frame='icrs')
        x_data = skycoord.galactic.l.degree
        y_data = skycoord.galactic.b.degree
    else:
        x_data = df[ra_col].values
        y_data = df[de_col].values

    if c is None:
        plot = ax.scatter(x_data, y_data, transform=ax.get_transform('world'),
                          color=color, s=markersize, alpha=alpha,
                          marker=marker, label=label, facecolor=facecolor, edgecolor=edgecolor,
                          lw=lw)
    else:
        plot = ax.scatter(x_data, y_data, transform=ax.get_transform('world'),
                          s=markersize, alpha=alpha,
                          marker=marker, label=label, c=c,
                          cmap=cmap, facecolor=facecolor, edgecolor=edgecolor,
                          lw=lw)
    return plot


def create_wcs_ax(fig, figsize=(10, 10), coords=('RA---', 'DEC--'), projection='AIT',
                  center=(270, 45), grid=True, nrows=1, ncols=1, index=1):
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    def _create_header():
        header = fits.Header()
        header['CTYPE1'] = coords[0] + projection
        header['CTYPE2'] = coords[1] + projection

        # Determine if the coordinate system is Galactic
        is_galactic = coords[0].startswith('GLON') and coords[1].startswith('GLAT')

        # Transform center coordinates if necessary
        if is_galactic:
            # Center is given in RA-Dec, need to convert to Galactic
            skycoord = SkyCoord(ra=center[0]*u.degree, dec=center[1]*u.degree, frame='icrs')
            header['CRVAL1'] = skycoord.galactic.l.degree
            header['CRVAL2'] = skycoord.galactic.b.degree
        else:
            # Center is in RA-Dec
            header['CRVAL1'] = center[0]
            header['CRVAL2'] = center[1]

        header['CRPIX1'] = 0
        header['CRPIX2'] = 0
        header['CDELT1'] = 0.5  # deg/pixel
        header['CDELT2'] = 0.5   # deg/pixel
        return header

    header = _create_header()
    wcs = WCS(header)
    ax = fig.add_subplot(nrows, ncols, index, projection=wcs)

    if grid:
        ax.grid(color='gray', ls='-')

    return ax


def plot_WTL_wcs_circles(coverage_df,
                         ra_col='ra', de_col='dec',
                         ax=None, label=None, alpha=0.3,
                         facecolor='none', edgecolor='k', scaling_factor=1.0):
    from astropy.visualization.wcsaxes import SphericalCircle
    from astropy.coordinates import SkyCoord

    if ax is None:
        ax = plt.gca()

    if label is not None:
        label = label

    for bid in coverage_df['batchId'].unique():
        batch = coverage_df[coverage_df['batchId'] == bid]
        central_stars = batch[batch['guide'] == 1]

        for _, c in central_stars.iterrows():
            # Create a SkyCoord object with the correct frame
            center_coord = SkyCoord(ra=c[ra_col]*u.deg, dec=c[de_col]*u.deg, frame='icrs')

            # Use the appropriate transform for the coordinate frame
            circle = SphericalCircle((center_coord.ra, center_coord.dec), 1.0*u.deg,
                                     edgecolor=edgecolor, facecolor=facecolor, alpha=alpha,
                                     label=None, transform=ax.get_transform('icrs'), lw=2)
            ax.add_patch(circle)

    circle_label_size = calculate_scaling_factor(ax, scaling_factor)

    ax.scatter([], [], s=circle_label_size, label=label,
               edgecolor=edgecolor, alpha=0.5, facecolor=facecolor)


def calculate_scaling_factor(ax, radius):
    fig = ax.figure
    renderer = fig.canvas.get_renderer()
    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Calculate the width and height of the axis in pixels
    bbox = ax.get_window_extent(renderer)
    width, height = bbox.width, bbox.height

    # Calculate the data range
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]

    # Calculate the number of pixels per data unit
    x_pixels_per_unit = width / x_range
    y_pixels_per_unit = height / y_range

    # Use the average scaling factor for x and y directions
    scaling_factor = (x_pixels_per_unit + y_pixels_per_unit) / 2

    # Calculate scatter size in points^2, considering the radius in data units
    scatter_size = (radius * scaling_factor) ** 2

    return scatter_size


def plot_wtl_field(wtl_df, wtl_df2=None, coords='galactic', weave_fov=False,
                   scaling_factor=1.0, fontsize=16, show_coverage=False, star_color='k',
                   plato_color='forestgreen', host_color='crimson', candidate_color='darkblue',
                   ax1_title=None, ax2_title=None, kepler_fov=False,
                   ra_col='ra', de_col='dec'):

    from astropy.visualization.wcsaxes import add_scalebar
    from matplotlib.font_manager import FontProperties

    lop1 = wtl_df[~wtl_df['sourceFlag'].isin([8, 9])]['LOP'].iloc[0]
    lop_field1 = wtl_df[wtl_df['LOP'] == lop1]
    circle_label = '$\pi$ deg$^2$ WEAVE FOV'

    if show_coverage:
        weave_fov = False
        circle_label = 'complete guide fibre coverage'

    figsize = (10, 10)
    ncols = 1
    if wtl_df2 is not None:
        ncols = 2
        lop2 = wtl_df2[~wtl_df2['sourceFlag'].isin([8, 9])]['LOP'].iloc[0]
        lop_field2 = wtl_df2[wtl_df2['LOP'] == lop2]
        center2 = (93.5, -43)
        figsize = (20, 10)
    if lop1 == 'N':
        center = (277, 53)
    if lop1 == 'S':
        center = (93.5, -43)
    if coords == 'galactic':
        coords = ('GLON-', 'GLAT-')
        x_label = 'Galactic Longitude (deg)'
        y_label = 'Galactic Latitude (deg)'
    if coords == 'radec':
        coords = ('RA---', 'DEC--')
        x_label = 'RA (hms)'
        y_label = 'DEC (deg)'

    fig = plt.figure(figsize=figsize)

    ax1 = create_wcs_ax(fig=fig, center=center, coords=coords,
                        projection='AIT', nrows=1, ncols=ncols, index=1)

    wcs_plot(wtl_df[wtl_df['sourceFlag'] != 9], ra_col=ra_col, de_col=de_col,
             ax=ax1, alpha=1, color=star_color)
    ax1.scatter([], [], color=star_color, s=40, label='Additional GDR3 stars')

    wcs_plot(lop_field1, ra_col=ra_col, de_col=de_col,
             ax=ax1, color=plato_color, alpha=0.7, label=None, markersize=3)
    ax1.scatter([], [], color=plato_color, label='PLATO LOP', s=40)

    if weave_fov:
        plot_WTL_wcs_circles(wtl_df, ra_col=ra_col, de_col=de_col,
                             ax=ax1, label=circle_label, alpha=0.5,
                             scaling_factor=scaling_factor)

    if show_coverage:
        full_guides = wtl_df[wtl_df['batchScore'].isin([3, 33, 333])]
        plot_WTL_wcs_circles(full_guides, ra_col=ra_col, de_col=de_col,
                             ax=ax1, label=circle_label, alpha=0.3, facecolor='lightgreen',
                             scaling_factor=scaling_factor)

    wcs_plot(wtl_df[wtl_df['sourceFlag'].isin([4, 5])], ra_col=ra_col, de_col=de_col, ax=ax1,
             alpha=1, label='known hosts', marker='o', markersize=80, facecolor='none',
             edgecolor=host_color, lw=2)
    wcs_plot(wtl_df[wtl_df['sourceFlag'].isin([7])], ra_col=ra_col, de_col=de_col, ax=ax1,
             alpha=1, label='candidate hosts', marker='D', markersize=65, facecolor='none',
             edgecolor=candidate_color, lw=2)

    fontprops = FontProperties(size=16)
    add_scalebar(ax1, 2*u.deg, label="2 deg", color="k", size_vertical=0.3,
                 fontproperties=fontprops)

    ax1.set_aspect('equal', adjustable='datalim')
    ax1.set_title(ax1_title, fontsize=fontsize)
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.legend(loc='upper left')
    ax1.set_xlabel(x_label, fontsize=fontsize)
    ax1.set_ylabel(y_label, fontsize=fontsize)

    if wtl_df2 is not None:
        ax2 = create_wcs_ax(fig=fig, center=center2, coords=coords,
                            projection='AIT', nrows=1, ncols=ncols, index=2)
        wcs_plot(wtl_df2[wtl_df2['sourceFlag'] != 9], ra_col=ra_col, de_col=de_col,
                 ax=ax2, alpha=1, color=star_color)
        wcs_plot(lop_field2, ra_col=ra_col, de_col=de_col, ax=ax2,
                 color=plato_color, alpha=0.7, label='PLATO LOP', markersize=3)
        wcs_plot(wtl_df2[wtl_df2['sourceFlag'].isin([4, 5])], ra_col=ra_col, de_col=de_col,
                 ax=ax2, alpha=1, label='Known hosts', marker='o', markersize=80,
                 facecolor='none', edgecolor=host_color, lw=2)
        wcs_plot(wtl_df2[wtl_df2['sourceFlag'].isin([7])], ra_col=ra_col, de_col=de_col,
                 ax=ax2, alpha=1, label='Candidate hosts', marker='D', markersize=65,
                 facecolor='none',
                 edgecolor=candidate_color, lw=2)

        if weave_fov:
            plot_WTL_wcs_circles(wtl_df2, ra_col=ra_col, de_col=de_col,
                                 ax=ax2, label=None, alpha=0.5)

        if show_coverage:
            full_guides2 = wtl_df2[wtl_df2['batchScore'].isin([3, 33, 333])]
            plot_WTL_wcs_circles(full_guides2, ra_col=ra_col, de_col=de_col,
                                 ax=ax2, label=None, alpha=0.3, facecolor='lightgreen')

        wcs_plot(wtl_df2[wtl_df2['sourceFlag'].isin([4, 5])], ra_col=ra_col, de_col=de_col,
                 ax=ax2, alpha=1, label='Known hosts', marker='o', markersize=80,
                 facecolor='none', edgecolor=host_color, lw=2)
        wcs_plot(wtl_df2[wtl_df2['sourceFlag'].isin([7])], ra_col=ra_col, de_col=de_col,
                 ax=ax2, alpha=1, label='Candidate hosts', marker='D', markersize=65,
                 facecolor='none',
                 edgecolor=candidate_color, lw=2)

        add_scalebar(ax2, 2*u.deg, label="2 deg", color="k", size_vertical=0.3,
                     fontproperties=fontprops)

        ax2.set_aspect('equal', adjustable='datalim')
        ax2.coords['ra'].tick_params(labelsize=18)
        ax2.coords['dec'].tick_params(labelsize=18)
        ax2.set_title(ax1_title, fontsize=fontsize)
        ax2.set_xlabel(x_label, fontsize=fontsize)
        ax2.set_ylabel(' ', fontsize=fontsize)

        return ax1, ax2

    return ax1


def plot_wtl_field_general(
    wtl_df,
    wtl_df2=None,
    coords='galactic',
    weave_fov=False,
    scaling_factor=1.0,
    fontsize=16,
    show_coverage=False,
    star_color='k',
    host_color='crimson',
    host_size=10, 
    candidate_color='darkblue',
    obi_color='forestgreen',
    obi_label='Object of interest',
    obi_flags=None,
    host_flags=(4, 5),
    candidate_flags=(7,),
    ax1_title=None,
    ax2_title=None,
    ra_col='ra',
    de_col='dec',
    center=None,):
    """
    Generalised WTL field plotter.

    Parameters
    ----------
    wtl_df : DataFrame
        Main field catalogue; all sources are plotted.
    wtl_df2 : DataFrame or None, optional
        Second field catalogue to plot side-by-side.
    coords : {'galactic', 'radec'}
        Coordinate system for the WCS projection.
    weave_fov : bool
        If True, plot WEAVE FOV circles for all pointings in each field.
    scaling_factor : float
        Scaling factor passed to plot_WTL_wcs_circles.
    fontsize : int
        Base font size for labels and titles.
    show_coverage : bool
        If True, plot only fields with full guide-fibre coverage
        (batchScore in [3, 33, 333]) and colour them differently.
        This also overrides weave_fov=False.
    star_color : str
        Color for plotting all stars.
    host_color : str
        Edge color for host markers.
    candidate_color : str
        Edge color for candidate host markers.
    host_flags : tuple or list
        sourceFlag values to treat as "known hosts".
    candidate_flags : tuple or list
        sourceFlag values to treat as "candidate hosts".
    ax1_title : str or None
        Title for first axis.
    ax2_title : str or None
        Title for second axis (if None, uses ax1_title).
    ra_col, de_col : str
        Column names for longitude and latitude in degrees.
    center : None, tuple, or tuple of tuples
        If None, use median positions of the data as in the original version.
        If a single tuple (ra, dec) is provided, use it as the centre for the first
        (or only) plot. If a second field is plotted and `center` is a tuple of two
        tuples, e.g. ((ra1, dec1), (ra2, dec2)), use them for the first and second
        plots respectively.
    """

    from astropy.visualization.wcsaxes import add_scalebar
    from matplotlib.font_manager import FontProperties
    import astropy.units as u

    # Label for FOV circles
    circle_label = r'$\pi$ deg$^2$ WEAVE FOV'
    if show_coverage:
        weave_fov = False  # coverage plot replaces generic FOV circles
        circle_label = 'complete guide fibre coverage'

    # Figure layout
    figsize = (10, 10)
    ncols = 1
    if wtl_df2 is not None:
        ncols = 2
        figsize = (20, 10)

    # WCS coordinate system
    if coords == 'galactic':
        wcs_coords = ('GLON-', 'GLAT-')
        x_label = 'Galactic Longitude (deg)'
        y_label = 'Galactic Latitude (deg)'
    elif coords == 'radec':
        wcs_coords = ('RA---', 'DEC--')
        x_label = 'RA (hms)'
        y_label = 'DEC (deg)'
    else:
        raise ValueError("coords must be 'galactic' or 'radec'")

    # Field centres
    if center is None:
        # Original behaviour: use median position of sources
        center1 = (wtl_df[ra_col].median(), wtl_df[de_col].median())
        center2 = None
        if wtl_df2 is not None:
            center2 = (wtl_df2[ra_col].median(), wtl_df2[de_col].median())
    else:
        # centre explicitly provided
        if wtl_df2 is None:
            # Single panel
            if isinstance(center[0], (tuple, list)):
                # User accidentally passed a tuple of tuples; take the first
                center1 = tuple(center[0])
            else:
                center1 = tuple(center)
            center2 = None
        else:
            # Two panels
            # If a single (ra, dec) is given, use for both; if two, use separately
            if len(center) == 2 and not isinstance(center[0], (tuple, list)):
                center1 = tuple(center)
                center2 = tuple(center)
            else:
                center1 = tuple(center[0])
                center2 = tuple(center[1])

    fig = plt.figure(figsize=figsize)

    # ------------------------------------------------------------------
    # First field
    # ------------------------------------------------------------------
    ax1 = create_wcs_ax(
        fig=fig,
        center=center1,
        coords=wcs_coords,
        projection='AIT',
        nrows=1,
        ncols=ncols,
        index=1,
    )

    # All stars
    wcs_plot(
        wtl_df[~wtl_df['sourceFlag'].isin([-1, 9])],
        ra_col=ra_col,
        de_col=de_col,
        ax=ax1,
        alpha=1.0,
        color=star_color,
    )
    # Dummy handle for legend
    ax1.scatter([], [], color=star_color, s=40, label='All stars')

    # WEAVE FOV circles
    if weave_fov:
        plot_WTL_wcs_circles(
            wtl_df,
            ra_col=ra_col,
            de_col=de_col,
            ax=ax1,
            label=circle_label,
            alpha=0.5,
            scaling_factor=scaling_factor,
        )

    # Guide-fibre coverage circles (subset by batchScore)
    if show_coverage and 'batchScore' in wtl_df.columns:
        full_guides = wtl_df[wtl_df['batchScore'].isin([3, 33, 333])]
        if not full_guides.empty:
            plot_WTL_wcs_circles(
                full_guides,
                ra_col=ra_col,
                de_col=de_col,
                ax=ax1,
                label=circle_label,
                alpha=0.3,
                facecolor='lightgreen',
                scaling_factor=scaling_factor,
            )

    # Known hosts
    if host_flags is not None and len(host_flags) > 0:
        hosts1 = wtl_df[wtl_df['sourceFlag'].isin(host_flags)]
        if not hosts1.empty:
            wcs_plot(
                hosts1,
                ra_col=ra_col,
                de_col=de_col,
                ax=ax1,
                alpha=1.0,
                label='Known hosts',
                marker='o',
                markersize=host_size,
                facecolor='none',
                edgecolor=host_color,
                lw=2,
            )

    # Candidate hosts
    if candidate_flags is not None and len(candidate_flags) > 0:
        cand1 = wtl_df[wtl_df['sourceFlag'].isin(candidate_flags)]
        if not cand1.empty:
            wcs_plot(
                cand1,
                ra_col=ra_col,
                de_col=de_col,
                ax=ax1,
                alpha=1.0,
                label='Candidate hosts',
                marker='D',
                markersize=host_size,
                facecolor='none',
                edgecolor=candidate_color,
                lw=2,
            )

    # Objects of interest
    if obi_flags is not None and len(obi_flags) > 0:
        obi1 = wtl_df[wtl_df['sourceFlag'].isin(obi_flags)]
        if not obi1.empty:
            wcs_plot(
                obi1,
                ra_col=ra_col,
                de_col=de_col,
                ax=ax1,
                alpha=1.0,
                label=obi_label,
                marker='^',
                markersize=host_size,
                facecolor='none',
                edgecolor=obi_color,
                lw=2,
            )

    # Scalebar and cosmetics
    fontprops = FontProperties(size=16)
    add_scalebar(
        ax1,
        2 * u.deg,
        label="2 deg",
        color="k",
        size_vertical=0.3,
        fontproperties=fontprops,
    )

    ax1.set_aspect('equal', adjustable='datalim')
    if ax1_title is not None:
        ax1.set_title(ax1_title, fontsize=fontsize)
    ax1.tick_params(axis='both', which='major', labelsize=18)
    ax1.legend(loc='upper left', markerscale=1.2)
    ax1.set_xlabel(x_label, fontsize=fontsize)
    ax1.set_ylabel(y_label, fontsize=fontsize)

    # ------------------------------------------------------------------
    # Second field (optional)
    # ------------------------------------------------------------------
    if wtl_df2 is not None:
        ax2 = create_wcs_ax(
            fig=fig,
            center=center2,
            coords=wcs_coords,
            projection='AIT',
            nrows=1,
            ncols=ncols,
            index=2,
        )

        # All stars
        wcs_plot(
            wtl_df2,
            ra_col=ra_col,
            de_col=de_col,
            ax=ax2,
            alpha=1.0,
            color=star_color,
        )
        ax2.scatter([], [], color=star_color, s=40, label='All stars')

        # WEAVE FOV circles
        if weave_fov:
            plot_WTL_wcs_circles(
                wtl_df2,
                ra_col=ra_col,
                de_col=de_col,
                ax=ax2,
                label=None,
                alpha=0.5,
                scaling_factor=scaling_factor,
            )

        # Guide-fibre coverage
        if show_coverage and 'batchScore' in wtl_df2.columns:
            full_guides2 = wtl_df2[wtl_df2['batchScore'].isin([3, 33, 333])]
            if not full_guides2.empty:
                plot_WTL_wcs_circles(
                    full_guides2,
                    ra_col=ra_col,
                    de_col=de_col,
                    ax=ax2,
                    label=None,
                    alpha=0.3,
                    facecolor='lightgreen',
                    scaling_factor=scaling_factor,
                )

        # Known hosts
        if host_flags is not None and len(host_flags) > 0:
            hosts2 = wtl_df2[wtl_df2['sourceFlag'].isin(host_flags)]
            if not hosts2.empty:
                wcs_plot(
                    hosts2,
                    ra_col=ra_col,
                    de_col=de_col,
                    ax=ax2,
                    alpha=1.0,
                    label='Known hosts',
                    marker='o',
                    markersize=80,
                    facecolor='none',
                    edgecolor=host_color,
                    lw=2,
                )

        # Candidate hosts
        if candidate_flags is not None and len(candidate_flags) > 0:
            cand2 = wtl_df2[wtl_df2['sourceFlag'].isin(candidate_flags)]
            if not cand2.empty:
                wcs_plot(
                    cand2,
                    ra_col=ra_col,
                    de_col=de_col,
                    ax=ax2,
                    alpha=1.0,
                    label='Candidate hosts',
                    marker='D',
                    markersize=65,
                    facecolor='none',
                    edgecolor=candidate_color,
                    lw=2,
                )

        add_scalebar(
            ax2,
            2 * u.deg,
            label="2 deg",
            color="k",
            size_vertical=0.3,
            fontproperties=fontprops,
        )

        ax2.set_aspect('equal', adjustable='datalim')
        if ax2_title is not None:
            ax2.set_title(ax2_title, fontsize=fontsize)
        elif ax1_title is not None:
            ax2.set_title(ax1_title, fontsize=fontsize)

        ax2.tick_params(axis='both', which='major', labelsize=18)
        ax2.set_xlabel(x_label, fontsize=fontsize)
        ax2.set_ylabel(' ', fontsize=fontsize)

        return ax1, ax2

    return ax1



def plot_tl_fov_coverage(field_stars_df, center=(0, 0), coords='radec',
                         outpath=None, figsize=(8, 8), alpha=1, mag_col=None,
                         host_flag_labels=None, plot_fov=True, show_central_star=False,
                         star_color='k', star_markersize=5):
    """
    Parameters
    ----------
    field_stars_df : DataFrame
    host_flag_labels : dict or None
        Optional mapping like {"Known planet hosts": 1, "Candidate planet hosts": 2}.
        Values must not be -1 or 9.
    """
    plot_df = field_stars_df.copy()

    # ------------------------------------------------------------------
    # Optional host categories: validate and adjust selection
    # ------------------------------------------------------------------
    if host_flag_labels is not None:
        # assert that forbidden sourceFlag values are not used
        flags = set(host_flag_labels.values())
        invalid = {f for f in flags if f in (-1, 9, 0)}
        assert not invalid, (
            f"sourceFlag values {sorted(invalid)} are not allowed in host_flag_labels "
            "(cannot use -1, 9 or 0)."
        )

        # Plot 0 (science), 9 (central guide) and all host flags
        allowed_flags = {0, 9} | flags
        plot_stars = plot_df[plot_df['sourceFlag'].isin(allowed_flags)]
    else:
        # Original behaviour: only 0 and 9 are considered
        plot_stars = plot_df[plot_df['sourceFlag'].isin([0, 9])]

    fig = plt.figure(figsize=figsize)

    if coords == 'radec':
        ax = create_wcs_ax(fig=fig, center=center, coords=('RA---', 'DEC--'))
        coord_x = 'GAIA_RA'
        coord_y = 'GAIA_DEC'
    elif coords == 'galactic':
        ax = create_wcs_ax(fig=fig, center=center, coords=('GLON-', 'GLAT-'))
        coord_x = 'GAIA_GAL_LONG'
        coord_y = 'GAIA_GAL_LAT'
    else:
        raise ValueError("coords must be 'radec' or 'galactic'")

    # Define marker styles for host categories (no '*')
    host_markers = ['d', 's', '^', 'v', 'P', 'X', 'h', 'H', 'o']
    if host_flag_labels is not None and len(host_flag_labels) > len(host_markers):
        raise ValueError(
            f"Not enough marker styles for {len(host_flag_labels)} host categories; "
            f"max supported is {len(host_markers)}."
        )

    regular_markersize = star_markersize
    host_markersize = 4 * regular_markersize

    for bid in plot_df['batchId'].unique():
        plot_ss = plot_stars[(plot_stars['batchId'] == bid)]
        if plot_ss.empty:
            continue

        # ------------------------------------------------------------------
        # Science targets (sourceFlag == 0), excluding host flags if provided
        # ------------------------------------------------------------------
        if host_flag_labels is None:
            sci_mask = ~plot_ss['sourceFlag'].isin([-1, 9])
        else:
            host_flags_set = set(host_flag_labels.values())
            sci_mask = (~plot_ss['sourceFlag'].isin([-1, 9]))

        sci_ss = plot_ss[sci_mask]
        if not sci_ss.empty:
            wcs_plot(
                sci_ss,
                ra_col=coord_x, de_col=coord_y,
                kepler_fov=False, alpha=alpha, fill=True, zorder=0,
                color=star_color, marker='*', markersize=regular_markersize, ax=ax,
                label='Science targets',
            )

        # ------------------------------------------------------------------
        # Central guide star(s): sourceFlag == 9
        # ------------------------------------------------------------------
        if show_central_star:
            guide_ss = plot_ss[plot_ss['sourceFlag'] == 9]
            if not guide_ss.empty:
                wcs_plot(
                    guide_ss,
                    ra_col=coord_x, de_col=coord_y,
                    kepler_fov=False, alpha=alpha, fill=True, zorder=0,
                    color='k', marker='D', markersize=15, ax=ax,
                    label='Central guide star',
                )

        # ------------------------------------------------------------------
        # Host categories, if provided
        # ------------------------------------------------------------------
        if host_flag_labels is not None:
            for (label, flag), marker in zip(host_flag_labels.items(), host_markers):
                cat_ss = plot_ss[plot_ss['sourceFlag'] == flag]
                if cat_ss.empty:
                    continue
                wcs_plot(
                    cat_ss,
                    ra_col=coord_x, de_col=coord_y,
                    kepler_fov=False, alpha=1.0, fill=True, zorder=1,
                    facecolor='none',edgecolor='k', marker=marker,
                    markersize=host_markersize, ax=ax, lw=1.5,
                    label=label,
                )

    # Plot WEAVE FOV
    if plot_fov:
        plot_WTL_wcs_circles(plot_df, ra_col=coord_x, de_col=coord_y, alpha=0.5)

    # Remove duplicated legend items
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), markerscale=1.2, loc='upper left')
    ax.axis('equal')

    if mag_col is not None:
        mag_type = mag_col.split('_')[-1]
        mag_ss = field_stars_df[~field_stars_df['sourceFlag'].isin([-1, 9])]
        ax.set_title(
            f'On-sky TL fields | '
            f'{mag_ss[mag_col].min():.2f} < {mag_type}'
            f' < {mag_ss[mag_col].max():.2f} mag'
        )
    else:
        ax.set_title('On-sky TL fields')

    if coords == 'radec':
        ax.set_xlabel('RA (hms)')
        ax.set_ylabel('Dec (deg)')
    elif coords == 'galactic':
        ax.set_xlabel('LON (deg)')
        ax.set_ylabel('LAT (deg)')

    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=150)

    plt.show()

    
    
def plot_input_fits_cat(cat_df,
                        center=(0, 0),
                        figsize=(12, 6),
                        coords='radec',
                        outpath=None,
                        mag_col='GAIA_MAG_G',
                        alpha=0.8):

    fig = plt.figure(figsize=figsize)

    # choose columns + WCS coord types
    if coords == 'radec':
        coord_x, coord_y = 'GAIA_RA', 'GAIA_DEC'
        wcs_coords = ('RA---', 'DEC--')
    elif coords == 'galactic':
        coord_x, coord_y = 'GAIA_GAL_LONG', 'GAIA_GAL_LAT'
        wcs_coords = ('GLON-', 'GLAT-')
    else:
        raise ValueError("coords must be 'radec' or 'galactic'")

    mag_type = mag_col.split('_')[-1]

    # left panel: plain scatter
    ax0 = fig.add_subplot(1, 2, 1)
    ax0.scatter(cat_df[coord_x],
                cat_df[coord_y],
                marker='*',
                s=10,
                alpha=alpha,
                color='crimson',
                label=(
                    f'Science targets; n = {cat_df.shape[0]},\n'
                    fr'{cat_df[mag_col].min():.2f} < {mag_type} < {cat_df[mag_col].max():.2f}')
                )

    ax0.axis('equal')
    ax0.grid(True)
    ax0.set_title('Deprojected')
    ax0.legend(loc='upper left', markerscale=2.2)

    # right panel: WCS axis (placed in slot 2 of the 1×2 grid)
    ax1 = create_wcs_ax(
        fig=fig,
        coords=wcs_coords,
        projection='AIT',
        center=center,          # use the function's center
        grid=True,
        nrows=1, ncols=2, index=2
    )

    wcs_plot(cat_df, ra_col=coord_x, de_col=coord_y,
             kepler_fov=False, alpha=alpha, fill=True, zorder=0,
             color='crimson', marker='*', markersize=10, ax=ax1)

    ax1.axis('equal')
    ax1.set_title('On-sky')

    if coords == 'radec':
        ax0.set_xlabel('RA (deg)')
        ax0.set_ylabel('Dec (deg)')
        ax1.set_xlabel('RA (hms)')
        ax1.set_ylabel('Dec (deg)')

    elif coords == 'galactic':
        ax0.set_xlabel('LON (deg)')
        ax0.set_ylabel('LAT (deg)')
        ax1.set_xlabel('LON (deg)')
        ax1.set_ylabel('LAT (deg)')

    plt.tight_layout()

    if outpath:
        plt.savefig(outpath, dpi=150)
    plt.show()
    plt.close('all')


if __name__ == '__main__':
    pass
