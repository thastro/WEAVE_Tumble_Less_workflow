import logging
import numpy as np
import pandas as pd
from scipy.spatial import distance
from astropy.coordinates import SkyCoord
import astropy.units as u





def assign_stars_to_circles(catalogue_df, circles,
                            ra_col='GAIA_RA',
                            de_col='GAIA_DEC',
                            src_id_col='CNAME',
                            max_passes=None,
                            outpath='./resources/WG/'):
    """
    Assign science targets and central guide stars to a list of circle centers,
    and build a guide-to-field lookup table.

    For each circle center (RA, Dec), the function:
      1) assigns all science targets within 0.98° to that circle (`fieldId`),
      2) selects the closest guide star within 0.98° as the circle's central guide
         (`centId`), and
      3) records, for every guide star inside any circle, the set of fields it
         belongs to and which field (if any) it serves as the central guide.

    It returns an augmented catalogue containing science targets and the identified
    central guide stars (with `guide=1` and correct `fieldId`), plus a lookup table
    mapping each guide to the fields it appears in and its central field.

    Args:
        catalogue_df (pd.DataFrame): Input catalogue with at least:
                                     - `sourceFlag` (9 = guide star; anything else = science),
                                     - RA/Dec columns (see `ra_col`, `de_col`),
                                     - an ID column used for guides (see `src_id_col`).
        circles (list[tuple[float, float]]): Circle centers as [(ra_deg, dec_deg), ...].
        ra_col (str, optional): Column name for Right Ascension (degrees).
                                Defaults to 'GAIA_RA'.
        de_col (str, optional): Column name for Declination (degrees). Defaults to 'GAIA_DEC'.
        src_id_col (str, optional): Identifier column used for guide-star IDs
                                    (e.g., 'CNAME' or 'GAIA_ID'). Defaults to 'CNAME'.
        max_passes (int or None, optional): Unused placeholder for API compatibility.
                                            Defaults to None.
        outpath (str or None, optional): If not None, write the guide lookup table to
                                         `<outpath>/guide_look_up.parquet`.
                                          Defaults to './resources/WG/'.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]:
            - pd.DataFrame: `catalogue_df_out` — concatenation of:
                * all science targets with assigned `fieldId` where applicable (`guide=0`), and
                * the *central* guide stars for each circle with `fieldId` and `guide=1`.
            - pd.DataFrame: `guide_look_up_df` — one row per guide with columns:
                * `src_id_col`: the guide identifier,
                * `lutFieldId`: underscore-joined list of all fieldIds containing this guide,
                * `centId`: the fieldId where the guide serves as the central guide
                            (NaN if none).

    Notes:
        - Membership and central-guide assignment use a hard-coded angular threshold of 0.98°.
          The `radius` parameter is currently ignored.
        - Only *central* guide stars are appended to the returned catalogue; all guide stars
          (central and non-central) appear in the lookup table.
        - If `outpath` is provided, the lookup table is saved as Parquet
          ('guide_look_up.parquet') without the DataFrame index.
    """

    science_targets = catalogue_df[catalogue_df['sourceFlag'] != 9].copy()
    science_targets['guide'] = 0

    guide_stars = catalogue_df[catalogue_df['sourceFlag'] == 9].copy()
    guide_stars['guide'] = 0
    guide_stars['fieldId'] = pd.NA

    science_coords = SkyCoord(ra=science_targets[ra_col].values * u.deg,
                              dec=science_targets[de_col].values * u.deg)
    guide_coords = SkyCoord(ra=guide_stars[ra_col].values * u.deg,
                            dec=guide_stars[de_col].values * u.deg)

    guide_look_up_dict = {}  # GAIA_ID -> (set of fieldIds, centId if any)
    cent_ids = []            # List of (central_id, field_number)

    for field_number, (ra_center, dec_center) in enumerate(circles, start=1):
        center_coord = SkyCoord(ra=ra_center * u.deg, dec=dec_center * u.deg)

        # Assign science targets. Circle radius set t0 0.98 deg, ensuring
        # that targets are not too close to the field plate edge.
        separations = science_coords.separation(center_coord)
        in_circle = separations < 0.98 * u.deg
        science_targets.loc[in_circle, 'fieldId'] = field_number

        # Guide star selection
        guide_separations = guide_coords.separation(center_coord)
        in_circle_guide = guide_separations < 0.98 * u.deg

        if in_circle_guide.any():
            idx_closest = np.argmin(guide_separations[in_circle_guide])
            central_idx = guide_stars.index[in_circle_guide][idx_closest]
            central_id = guide_stars.at[central_idx, src_id_col]
            cent_ids.append((central_id, field_number))
        else:
            logging.warning(f'NOOOOOOooooo guide for field {field_number}. '
                            'Make sure you added the guide catalogue correctly!')
            continue

        # Add all guide stars to lookup dict
        for idx in guide_stars.index[in_circle_guide]:
            gaia_id = guide_stars.at[idx, src_id_col]
            if gaia_id not in guide_look_up_dict:
                guide_look_up_dict[gaia_id] = (set(), None)
            guide_look_up_dict[gaia_id][0].add(field_number)

    # Now assign centId based on cent_ids list (overwrites previous guesses)
    for central_id, field_number in cent_ids:
        if central_id in guide_look_up_dict:
            guide_look_up_dict[central_id] = (
                guide_look_up_dict[central_id][0],
                field_number
            )
        else:
            # In case somehow central star was not previously added
            guide_look_up_dict[central_id] = ({field_number}, field_number)

    # Create guide lookup table
    guide_look_up_df = pd.DataFrame([
        {
            src_id_col: gaia_id,
            'lutFieldId': '_'.join(map(str, sorted(field_ids))),
            'centId': cent_id if cent_id is not None else np.nan
        }
        for gaia_id, (field_ids, cent_id) in guide_look_up_dict.items()
    ])

    # Recover central guide stars only and assign correct fieldId
    cent_df = pd.DataFrame(cent_ids, columns=[src_id_col, 'fieldId'])
    cent_df['guide'] = 1

    central_guide_stars = cent_df.merge(
        guide_stars, on=src_id_col, how='left', suffixes=('', '_drop'))
    central_guide_stars['fieldId'] = cent_df['fieldId']
    central_guide_stars['guide'] = 1

    central_guide_stars = central_guide_stars.drop(
        columns=[c for c in central_guide_stars if c.endswith('_drop')])

    # Overwrite with correct fieldId and guide flag
    central_guide_stars['fieldId'] = cent_df['fieldId']
    central_guide_stars['guide'] = 1

    # Final catalogue
    catalogue_df = pd.concat([science_targets, central_guide_stars], ignore_index=True)

    if outpath is not None:
        guide_look_up_df.to_parquet(outpath + 'guide_look_up.parquet', index=False)

    return catalogue_df, guide_look_up_df


def assign_batch_id_optimized(circle_coords_df,
                              fields_selection_df,
                              batch_size=3,
                              merge_threshold=3,
                              ra_col='RAdeg',
                              de_col='DEdeg'):
    """
    Forces as many 3-field batches as possible.
    Any leftover 1- or 2-field batches are merged where possible.
    """

    coords = circle_coords_df[[ra_col + '_c', de_col + '_c']].values
    field_ids = circle_coords_df['fieldId_c'].values
    # Calculate pairwise distances
    from scipy.spatial.distance import pdist, squareform
    dist_matrix = squareform(pdist(coords))

    # Create list of unassigned field indices
    remaining = set(range(len(coords)))
    batches = []
    batch_id = 1

    while len(remaining) >= 3:
        # Try forming a triangle (3 fields) with minimal max separation
        triplet = None
        min_max_dist = np.inf

        for i in remaining:
            for j in remaining:
                if j <= i:
                    continue
                for k in remaining:
                    if k <= j or k == i:
                        continue
                    dists = [dist_matrix[i, j], dist_matrix[i, k], dist_matrix[j, k]]
                    max_d = max(dists)
                    if max_d < min_max_dist:
                        min_max_dist = max_d
                        triplet = (i, j, k)

        if triplet is not None:
            batches.append((batch_id, [field_ids[i] for i in triplet]))
            remaining -= set(triplet)
            batch_id += 1
        else:
            break

    # Handle remaining fields (1 or 2)
    leftovers = list(remaining)
    if leftovers:
        batches.append((batch_id, [field_ids[i] for i in leftovers]))
        batch_id += 1

    # Map back to DataFrames
    field_to_batch = {fid: b_id for b_id, fids in batches for fid in fids}
    circle_coords_df['batchId'] = circle_coords_df['fieldId_c'].map(field_to_batch)
    fields_selection_df['batchId'] = fields_selection_df['fieldId'].map(field_to_batch)

    # Clean up
    circle_coords_df['batchId'] = circle_coords_df['batchId'].astype('Int64')
    fields_selection_df['batchId'] = fields_selection_df['batchId'].astype('Int64')

    return fields_selection_df.sort_values(by=['batchId', ra_col]), batch_id - 1


def assign_batch_id_optimized_fast(circle_coords_df,
                                   fields_selection_df,
                                   batch_size=3,
                                   merge_threshold=3,
                                   ra_col='RAdeg',
                                   de_col='DEdeg'):

    circle_coords_df = circle_coords_df.copy()
    fields_selection_df = fields_selection_df.copy()

    coords = circle_coords_df[[f"{ra_col}_c", f"{de_col}_c"]].values
    field_ids = circle_coords_df['fieldId_c'].values
    n_fields = len(coords)

    # Compute pairwise distances (condensed and full matrix)
    dist_condensed = distance.pdist(coords)
    dist_matrix = distance.squareform(dist_condensed)

    # Build sorted list of all unique field pairs by distance
    pairs = [
        (i, j, dist_matrix[i, j])
        for i in range(n_fields)
        for j in range(i + 1, n_fields)
    ]
    pairs.sort(key=lambda x: x[2])  # sort by distance

    assigned = set()
    batches = []
    batch_id = 1

    # Greedily assign triplets
    for i, j, _ in pairs:
        if i in assigned or j in assigned:
            continue
        # Find third candidate closest to both i and j
        candidates = [
            k for k in range(n_fields)
            if k not in assigned and k != i and k != j
        ]
        best_k = None
        best_max_sep = np.inf
        for k in candidates:
            dists = [dist_matrix[i, k], dist_matrix[j, k], dist_matrix[i, j]]
            max_sep = max(dists)
            if max_sep < best_max_sep:
                best_max_sep = max_sep
                best_k = k
        if best_k is not None:
            batches.append((batch_id, [field_ids[i], field_ids[j], field_ids[best_k]]))
            assigned.update([i, j, best_k])
            batch_id += 1

    # Handle leftovers (0, 1, or 2)
    unassigned = [idx for idx in range(n_fields) if idx not in assigned]
    if unassigned:
        batches.append((batch_id, [field_ids[idx] for idx in unassigned]))
        batch_id += 1

    # Map to batchId
    field_to_batch = {fid: b_id for b_id, fids in batches for fid in fids}
    circle_coords_df['batchId'] = circle_coords_df['fieldId_c'].map(field_to_batch)
    fields_selection_df['batchId'] = fields_selection_df['fieldId'].map(field_to_batch)

    return fields_selection_df.sort_values(by=['batchId', ra_col]), batch_id - 1


def get_science_covering_guide_star_circles_fixed_grid(df,
                                                       ra_col='GAIA_RA',
                                                       dec_col='GAIA_DEC',
                                                       mag_col='GAIA_MAG_G',
                                                       flag_col='sourceFlag',
                                                       central_star_mag_lim=16.0,
                                                       spacing=1.5,
                                                       match_radius=1.0,
                                                       circles_df=None,
                                                       filter_circles_to_science_bbox=True):
    """
    Determine circle center coordinates (RA, Dec) anchored on guide stars to cover
    science targets. If a set of precomputed circles is provided, their centers will
    be snapped to the nearest guide stars; otherwise, a fixed grid of candidate
    centers will be generated within the science bounding box and then snapped.

    Args:
        df (pd.DataFrame): Input catalogue containing RA, Dec, and a sourceFlag column,
                        where sourceFlag==9 identifies guide stars.
        ra_col (str): Name of the column containing Right Ascension in degrees.
        dec_col (str): Name of the column containing Declination in degrees.
        flag_col (str): Name of the column indicating source type
                        (9 = guide star, -1 = sky position).
        spacing (float): Grid spacing in degrees, only used if no `circles_df` is given.
        match_radius (float): Maximum angular distance (in degrees) for snapping grid
                            or input circle centers to a guide star. This value also
                            defines the science coverage radius of each circle.
        circles_df (pd.DataFrame or None): Optional DataFrame with precomputed circle
                                        centers to be aligned with guide stars.
        filter_circles_to_science_bbox (bool): If True and `circles_df` is provided,
                                            restricts candidate circles to those
                                            lying within the science bounding box,
                                            preserving consistency with the default
                                            grid-based logic.

    Returns:
        list[tuple[float, float]]: A list of (RA, Dec) coordinates representing the
                                final, guide-star-centered circle positions.
    """

    from collections import defaultdict

    # Split guide and science stars
    guides = df[df[flag_col] == 9].copy()
    science = df[df[flag_col] != 9].copy()
    guides = guides[guides[mag_col] <= central_star_mag_lim]

    assert not guides.empty, "No guide stars available. Check input DataFrame!"
    assert not science.empty, "No science targets available. Check input DataFrame!"

    # Science bounding box defining the grid area
    ra_min, ra_max = science[ra_col].min(), science[ra_col].max()
    dec_min, dec_max = science[dec_col].min(), science[dec_col].max()

    # Candidate circle centres (grid OR provided circle centres)
    if circles_df is None:
        # Generate grid over the science bbox
        ra_grid = np.arange(ra_min, ra_max + spacing, spacing)
        dec_grid = np.arange(dec_min, dec_max + spacing, spacing)
        # assert if grid is empty (==0)

        assert ra_grid.size > 0 and dec_grid.size > 0, "Grid generation failed!"

        ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
        cand_ra = ra_mesh.ravel()
        cand_dec = dec_mesh.ravel()

    else:
        # Use provided circle centres
        if filter_circles_to_science_bbox:
            mask = ((circles_df[ra_col].between(ra_min, ra_max, inclusive='both')) &
                    (circles_df[dec_col].between(dec_min, dec_max, inclusive='both')))

            sub = circles_df.loc[mask]
        else:
            sub = circles_df

        assert not sub.empty, "No circles remain after filtering to science bbox!"

        cand_ra = sub[ra_col].to_numpy(dtype=float)
        cand_dec = sub[dec_col].to_numpy(dtype=float)

    # Build candidate coords
    grid_points = SkyCoord(ra=cand_ra * u.deg, dec=cand_dec * u.deg)

    # Guide star coordinates
    gs_coords = SkyCoord(ra=guides[ra_col].values * u.deg,
                         dec=guides[dec_col].values * u.deg)

    # Snap candidate points to nearest guide star (within radius)
    idx_nearest, sep2d, _ = grid_points.match_to_catalog_sky(gs_coords)
    mask_within_radius = sep2d < match_radius * u.deg
    assert np.any(mask_within_radius), ("No grid/circle centers within match_radius "
                                        "of any guide star!")

    snapped_guide_coords = gs_coords[idx_nearest[mask_within_radius]]

    # De-duplicate (RA, Dec) pairs (round to micro-degree to avoid float fuzz)
    rounded_ra = np.round(snapped_guide_coords.ra.deg, 6)
    rounded_dec = np.round(snapped_guide_coords.dec.deg, 6)
    coord_pairs = np.rec.fromarrays([rounded_ra, rounded_dec], names='ra,dec')
    _, unique_idx = np.unique(coord_pairs, return_index=True)
    unique_circles = snapped_guide_coords[unique_idx]

    # Science coordinates
    sc_coords = SkyCoord(ra=science[ra_col].values * u.deg,
                         dec=science[dec_col].values * u.deg)

    # Keep only circles that cover at least one science target
    circle_map = defaultdict(set)
    for i, center in enumerate(unique_circles):
        sep = center.separation(sc_coords)
        science_in_circle = np.where(sep < match_radius * u.deg)[0]
        if len(science_in_circle) > 0:
            circle_map[i] = set(science_in_circle)

    assert circle_map, "No circles cover any science targets!"

    # Final selected circles
    final_circle_idxs = list(circle_map.keys())
    final_centers = unique_circles[final_circle_idxs]

    # Return as list of (RA, Dec) tuples
    return list(zip(final_centers.ra.deg, final_centers.dec.deg))


def get_science_covering_guide_star_circles_greedy(df,
                                                   ra_col='GAIA_RA',
                                                   dec_col='GAIA_DEC',
                                                   flag_col='sourceFlag',
                                                   cname_col='CNAME',
                                                   gaia_id_col='GAIA_ID',
                                                   mag_col='GAIA_MAG_G',
                                                   central_star_mag_lim=16.0,
                                                   spacing=1.5,
                                                   circle_radius=1.0,
                                                   match_radius=1.0,
                                                   move_radius=None,
                                                   min_science_per_circle=1,
                                                   random_state=42):
    """
    Compute an optimized set of guide-star-centered circles that cover science targets.
    The algorithm seeds a grid in a local tangent plane over the science footprint,
    snaps candidate centers to the nearest guide within `match_radius`, augments with
    additional guides to cover all coverable targets, solves a greedy set-cover to
    minimize the circle count, and optionally performs a local refinement step that
    swaps centers to nearby guides within `move_radius`. Targets that are not coverable
    by any guide within `circle_radius` are returned separately rather than causing the
    routine to fail.

    Args:
        df (pd.DataFrame): Input catalogue containing RA/Dec, a source type flag, IDs,
            and magnitudes. Guide stars are indicated by `flag_col == 9`.
        ra_col (str, optional): Column name for Right Ascension in degrees.
            Defaults to 'GAIA_RA'.
        dec_col (str, optional): Column name for Declination in degrees.
            Defaults to 'GAIA_DEC'.
        flag_col (str, optional): Column indicating source type
            (9 = guide star, -1 = sky position). Defaults to 'sourceFlag'.
        cname_col (str, optional): Column with the target name/identifier used in
            reports of uncovered targets. Defaults to 'CNAME'.
        gaia_id_col (str, optional): Column with the Gaia source identifier used in
            reports of uncovered targets. Defaults to 'GAIA_ID'.
        mag_col (str, optional): Column with guide-star magnitude used to filter
            central guide candidates. Defaults to 'GAIA_MAG_G'.
        central_star_mag_lim (float, optional): Maximum allowed magnitude for guide
            stars that may serve as circle centers. Defaults to 16.0.
        spacing (float, optional): Grid spacing (degrees) in the tangent plane used to
            seed candidate centers. Defaults to 1.5.
        circle_radius (float, optional): Science coverage radius (degrees) for each
            circle. Defaults to 1.0.
        match_radius (float, optional): Maximum distance (degrees) to snap a seeded grid
            node to a guide star. Also used when scoring candidates during seeding.
            Defaults to 1.0.
        move_radius (float, optional): Maximum distance (degrees) a chosen center may be
            swapped/nudged to a nearby guide during local refinement. If None,
            `match_radius` is used. Defaults to None.
        min_science_per_circle (int, optional): Minimum number of science targets that a
            candidate circle must cover to be retained. Defaults to 1.
        random_state (int, optional): Seed for the internal random number generator used
            for shuffling/tie-breaking in the greedy set-cover. Defaults to 42.

    Returns:
        list[tuple[float, float]]: Final (RA_deg, Dec_deg) circle centers anchored on
            guide stars after greedy selection and local refinement.
        pd.DataFrame: DataFrame of science targets that are not coverable by any guide
            within `circle_radius`. Contains columns `CNAME`, `GAIA_ID`, `RA`, `Dec`.
            Empty if all coverable targets are covered.

    Notes:
        - Coverage and geometric operations are performed in a local tangent-plane
          (astropy `SkyOffsetFrame`) centered on the median science position.
        - When no viable candidates remain (`min_science_per_circle` filter) or no
          coverable science is found, the function returns an empty result.
    """

    from astropy.coordinates import SkyOffsetFrame

    rng = np.random.default_rng(random_state)
    if move_radius is None:
        move_radius = match_radius

    guides = df[df[flag_col] == 9].copy()
    science = df[df[flag_col] != 9].copy()

    guides = guides[guides[mag_col] <= central_star_mag_lim]

    if guides.empty or science.empty:
        return [], pd.DataFrame(columns=['CNAME', 'GAIA_ID', 'RA', 'Dec'])

    # Tangent plane around science median
    sci_coords = SkyCoord(science[ra_col].to_numpy() * u.deg,
                          science[dec_col].to_numpy() * u.deg)
    origin = SkyCoord(np.median(sci_coords.ra.deg) * u.deg,
                      np.median(sci_coords.dec.deg) * u.deg)
    off = SkyOffsetFrame(origin=origin)
    sci_off = sci_coords.transform_to(off)
    sci_xy = np.column_stack([sci_off.lon.deg, sci_off.lat.deg])  # (Nsci,2)

    gs_coords = SkyCoord(guides[ra_col].to_numpy() * u.deg,
                         guides[dec_col].to_numpy() * u.deg)
    gs_off = gs_coords.transform_to(off)
    gs_xy = np.column_stack([gs_off.lon.deg, gs_off.lat.deg])     # (Ngs,2)

    # Science footprint bounds (convex-hull + buffer; bbox fallback)
    use_poly = False
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(sci_xy)
        hull_pts = sci_xy[hull.vertices]
        try:
            from shapely.geometry import Polygon, Point
            poly = Polygon(hull_pts).buffer(circle_radius)
            minx, miny, maxx, maxy = poly.bounds
            use_poly = True
        except Exception:
            use_poly = False
    except Exception:
        use_poly = False

    if not use_poly:
        xmin, ymin = sci_xy.min(axis=0) - circle_radius
        xmax, ymax = sci_xy.max(axis=0) + circle_radius
        minx, miny, maxx, maxy = xmin, ymin, xmax, ymax

    # --- Grid inside buffered footprint ---
    xs = np.arange(minx, maxx + spacing*0.5, spacing)
    ys = np.arange(miny, maxy + spacing*0.5, spacing)
    gx, gy = np.meshgrid(xs, ys)
    grid_xy = np.column_stack([gx.ravel(), gy.ravel()])

    if use_poly:
        from shapely.geometry import Point
        inside = np.array([poly.contains(Point(x, y)) for x, y in grid_xy])
        grid_xy = grid_xy[inside]

    if grid_xy.size == 0:
        logging.warning("No grid points inside buffered science footprint.")
        return [], pd.DataFrame(columns=['CNAME', 'GAIA_ID', 'RA', 'Dec'])

    # Helper funtions
    def coverage_vector_for_center_xy(center_xy):
        dx = sci_xy[:, 0] - center_xy[0]
        dy = sci_xy[:, 1] - center_xy[1]
        return (dx*dx + dy*dy) <= (circle_radius**2)

    def coverage_matrix_for_coords(centers_coords):
        cen_off = centers_coords.transform_to(off)
        cen_xy = np.column_stack([cen_off.lon.deg, cen_off.lat.deg])
        dx = sci_xy[:, 0][:, None] - cen_xy[:, 0][None, :]
        dy = sci_xy[:, 1][:, None] - cen_xy[:, 1][None, :]
        return (dx*dx + dy*dy) <= (circle_radius**2)

    # Snapping: for each grid node, choose guide within match_radius
    # that covers MOST science
    # Dist from grid -> guides
    xg = grid_xy[:, 0][:, None]            # (Ng,1)
    yg = grid_xy[:, 1][:, None]
    xs_g = gs_xy[:, 0][None, :]            # (1,Ngs)
    ys_g = gs_xy[:, 1][None, :]

    d2_g = (xg - xs_g)**2 + (yg - ys_g)**2
    within = d2_g <= (match_radius**2)

    cand_idxs = []
    for i in range(grid_xy.shape[0]):
        gids = np.where(within[i])[0]
        if gids.size == 0:
            continue
        # Pick guide that maximizes science coverage
        best_g = None
        best_score = -1
        for gidx in gids:
            vec = coverage_vector_for_center_xy(gs_xy[gidx])
            score = int(vec.sum())
            if score > best_score:
                best_score = score
                best_g = gidx
        cand_idxs.append(best_g)

    if len(cand_idxs) == 0:
        # No snapped candidates from grid
        cand_idxs = np.array([], dtype=int)
    else:
        cand_idxs = np.unique(np.array(cand_idxs, dtype=int))

    cand_coords = gs_coords[cand_idxs] if cand_idxs.size else SkyCoord([], [], unit='deg')

    # --- Feasibility: which science are coverable at all? (any guide within circle_radius) ---
    dx_all = sci_xy[:, 0][:, None] - gs_xy[:, 0][None, :]
    dy_all = sci_xy[:, 1][:, None] - gs_xy[:, 1][None, :]
    d2_all = dx_all*dx_all + dy_all*dy_all
    coverable_mask = (d2_all <= (circle_radius**2)).any(axis=1)
    uncoverable_mask = ~coverable_mask

    if np.any(uncoverable_mask):
        # Prepare warning table
        cn = science[cname_col] if cname_col in science.columns else pd.Series(['<NA>']*np.sum(uncoverable_mask))
        gid = science[gaia_id_col] if gaia_id_col in science.columns else pd.Series(['<NA>']*np.sum(uncoverable_mask))
        uncv = pd.DataFrame({
            'CNAME': np.asarray(cn)[uncoverable_mask],
            'GAIA_ID': np.asarray(gid)[uncoverable_mask],
            'RA': sci_coords.ra.deg[uncoverable_mask],
            'Dec': sci_coords.dec.deg[uncoverable_mask],
        })
    else:
        uncv = pd.DataFrame(columns=['CNAME', 'GAIA_ID', 'RA', 'Dec'])

    # Augment: if snapped candidates don't cover all coverable science, add all nearby guides
    cov = coverage_matrix_for_coords(cand_coords) if cand_idxs.size else np.zeros((sci_xy.shape[0], 0), bool)
    not_yet_cov = coverable_mask & (~cov.any(axis=1) if cov.shape[1] > 0 else np.ones(cov.shape[0] if cov.size else sci_xy.shape[0], bool))
    if np.any(not_yet_cov):
        add_gidx = np.where((d2_all[not_yet_cov, :] <= (circle_radius**2)).any(axis=0))[0]
        cand_idxs = np.unique(np.concatenate([cand_idxs, add_gidx])) if cand_idxs.size else add_gidx
        cand_coords = gs_coords[cand_idxs]
        cov = coverage_matrix_for_coords(cand_coords)

    # Filter weak circles
    keep_c = cov.sum(axis=0) >= max(1, min_science_per_circle)
    if not np.any(keep_c):
        print(uncv)
        return []
    cand_idxs = cand_idxs[keep_c]
    cand_coords = cand_coords[keep_c]
    cov = cov[:, keep_c]

    # --- Set cover over coverable subset ---
    cov_eff = cov[coverable_mask, :]  # (Neff, Ncand)
    if cov_eff.shape[0] == 0:
        print(uncv)
        return []

    uncovered = np.ones(cov_eff.shape[0], dtype=bool)
    order = np.arange(cov_eff.shape[1])
    rng.shuffle(order)

    chosen_rel = []
    while uncovered.any():
        gains = (cov_eff[uncovered, :]).sum(axis=0)
        j = order[np.argmax(gains[order])]
        if gains[j] == 0:
            break
        chosen_rel.append(j)
        newly = cov_eff[:, j] & uncovered
        uncovered[newly] = False

    # Redundancy prune (currently disabled). First tests show little to no effect.
    # Decision to be made after more thorough evaluation.
    chosen_rel = list(chosen_rel)

    def prune(choices):
        changed = True
        while changed:
            changed = False
            for k in list(choices):
                trial = [t for t in choices if t != k]
                if len(trial) == 0:
                    continue
                if cov_eff[:, trial].any(axis=1).all():
                    choices.remove(k)
                    changed = True
        return choices

    # chosen_rel = prune(chosen_rel)

    # Local refinement (nudge to nearby guides to remove another circle)
    chosen = np.array(chosen_rel, dtype=int)

    # Precompute chosen center XY for neighborhood queries
    chosen_xy = np.column_stack([gs_coords[cand_idxs[chosen]].transform_to(off).lon.deg,
                                 gs_coords[cand_idxs[chosen]].transform_to(off).lat.deg])

    # Iteratively try improving the current guide-star selection by testing
    # local substitutions: for each chosen guide, look for nearby alternatives
    # (within move_radius) that maintain full science-target coverage but allow
    # one of the current guides to be pruned. If a beneficial swap is found,
    # update the chosen set and repeat until no improvement remains.
    improved = True
    while improved and len(chosen) > 1:
        improved = False
        # current coverage matrix of chosen
        for col_idx in range(len(chosen)):
            # Nearby alternative guides within move_radius of this chosen center
            cx, cy = chosen_xy[col_idx]
            d2_to_all_guides = (gs_xy[:,0] - cx)**2 + (gs_xy[:,1] - cy)**2
            alt_idxs = np.where(d2_to_all_guides <= (move_radius**2))[0]
            if alt_idxs.size == 0:
                continue

            # Try swapping this column with each alternative; then prune and see if count drops
            for g_alt in alt_idxs:
                alt_vec_full = ((sci_xy[:,0] - gs_xy[g_alt,0])**2 +
                                (sci_xy[:,1] - gs_xy[g_alt,1])**2) <= (circle_radius**2)
                alt_vec = alt_vec_full[coverable_mask]  # restrict to coverable set

                # Build trial matrix: replace column col_idx with alt_vec
                trial_cols = cov_eff[:, chosen].copy()
                trial_cols[:, col_idx] = alt_vec

                # If replacement loses coverage, skip
                if not trial_cols.any(axis=1).all():
                    continue

                # Trial prune
                active = list(range(trial_cols.shape[1]))
                changed = True
                while changed:
                    changed = False
                    for k in list(active):
                        trial = [t for t in active if t != k]
                        if len(trial) == 0:
                            continue
                        if trial_cols[:, trial].any(axis=1).all():
                            active.remove(k)
                            changed = True

                if len(active) < len(chosen):
                    if g_alt not in cand_idxs:
                        cand_idxs = np.append(cand_idxs, g_alt)
                        # Extend cov_eff with this new column
                        new_col = alt_vec[:, None]
                        cov_eff = np.hstack([cov_eff, new_col])
                        # Map trial column index (col_idx) to this new absolute column (last)
                        chosen[col_idx] = cov_eff.shape[1] - 1
                    else:
                        # Set chosen[col_idx] to absolute index of existing column for g_alt
                        # Find its absolute column in cov (cand_idxs alignment)
                        abs_idx = np.where(cand_idxs == g_alt)[0][0]
                        # Map to cov_eff column indices (same order as cand_idxs after filters)
                        chosen[col_idx] = abs_idx

                    # Now rebuild 'chosen' from 'active' relative positions
                    # 'trial_cols' corresponds to cov_eff[:, chosen] ordering; map back:
                    chosen = chosen[active]
                    # Recompute chosen_xy after change
                    chosen_xy = np.column_stack([
                        gs_coords[cand_idxs[chosen]].transform_to(off).lon.deg,
                        gs_coords[cand_idxs[chosen]].transform_to(off).lat.deg
                    ])
                    improved = True
                    break  # Break alt loop; restart outer while
            if improved:
                break
        # End for col_idx
    # End while improved
    if len(chosen):
        final_coords = gs_coords[cand_idxs[chosen]]
    else:
        final_coords = SkyCoord([], [], unit='deg')

    centers = list(zip(final_coords.ra.deg, final_coords.dec.deg))

    return centers


def prune_circles(catalogue_df, circles,
                  ra_col="GAIA_RA",
                  dec_col="GAIA_DEC",
                  flag_col="sourceFlag",
                  ref_circle_idx=None,
                  ref_ra=None, ref_dec=None,
                  circle_radius_deg=0.98,
                  chunk_size_targets=100_000):
    """
    Parameters
    ----------
    catalogue_df : pd.DataFrame
        Must contain RA/Dec and sourceFlag.
    circles : list[tuple]
        Either [(ra, dec), ...] or [(id, ra, dec), ...] in degrees.
    ref_circle_idx : int, optional
        Index into `circles` to use as the reference circle.
    ref_ra, ref_dec : float, optional
        Reference circle centre if you don't pass ref_circle_idx.
    circle_radius_deg : float
        Circle radius (assumed same for all).
    chunk_size_targets : int
        Chunk size for separation computations.

    Returns
    -------
    kept_circle_indices : list[int]
        Indices in the ORIGINAL `circles` list that remain after pruning.
    removed_circle_indices : list[int]
        Indices in the ORIGINAL `circles` list that were removed as redundant.
    kept_circles : list[tuple]
        The corresponding circle tuples from the ORIGINAL `circles` list.
    coverage_lut : pd.DataFrame
        Long-form table: columns ['circle_idx','target_index'] mapping each kept/intersecting
        circle to the science targets (row index into `science_df`) it covers.
    science_df : pd.DataFrame
        Subset of `catalogue_df` where sourceFlag != 9, index reset.
    intersecting_circle_indices : list[int]
        All intersecting circle indices BEFORE pruning.
    circle_ids : list
        List of IDs if provided in the input tuples; otherwise equals kept_circle_indices.
    """

    ids = []
    ras = []
    decs = []

    for tup in circles:
        if len(tup) == 2:
            ra, dec = tup
            ids.append(None)
            ras.append(float(ra))
            decs.append(float(dec))
        elif len(tup) == 3:
            cid, ra, dec = tup
            ids.append(cid)
            ras.append(float(ra))
            decs.append(float(dec))
        else:
            raise ValueError("Each circle tuple must be (ra, dec) or (id, ra, dec).")

    ras = np.asarray(ras, dtype=float)
    decs = np.asarray(decs, dtype=float)

    # Reference circle selection
    if (ref_circle_idx is None) == (ref_ra is None or ref_dec is None):
        raise ValueError("Provide either ref_circle_idx OR (ref_ra and ref_dec), "
                         "but not both.")

    if ref_circle_idx is not None:
        if not (0 <= ref_circle_idx < len(circles)):
            raise IndexError("ref_circle_idx out of range.")

        ref_ra = ras[ref_circle_idx]
        ref_dec = decs[ref_circle_idx]

    # Science targets only
    science_df = catalogue_df[catalogue_df[flag_col] != 9].reset_index(drop=True)
    sci_coords = SkyCoord(science_df[ra_col].to_numpy()*u.deg,
                          science_df[dec_col].to_numpy()*u.deg, frame="icrs")

    # Intersecting circles
    ref_center = SkyCoord(ref_ra*u.deg, ref_dec*u.deg, frame="icrs")
    all_centers = SkyCoord(ras*u.deg, decs*u.deg, frame="icrs")
    sep_to_ref = all_centers.separation(ref_center).deg
    intersect_mask = sep_to_ref <= (360.0 + 1e-12)
    intersecting_idx = np.nonzero(intersect_mask)[0].tolist()

    # Coverage computation
    R = circle_radius_deg
    coverage = {}           # idx -> set(target_index)
    coverage_rows = []      # (circle_idx, target_index)

    N_sci = len(science_df)
    target_chunks = [slice(i, min(i+chunk_size_targets, N_sci))
                     for i in range(0, N_sci, chunk_size_targets)]

    inter_ras = ras[intersecting_idx]
    inter_decs = decs[intersecting_idx]
    inter_centers = SkyCoord(inter_ras*u.deg, inter_decs*u.deg, frame="icrs")

    for j, cidx in enumerate(intersecting_idx):
        c_center = inter_centers[j:j+1]
        covered_idx = []
        for sl in target_chunks:
            seps = sci_coords[sl].separation(c_center).deg
            inside = np.where(seps <= R + 1e-12)[0]
            if inside.size:
                covered_idx.extend((inside + sl.start).tolist())
        s = set(covered_idx)
        coverage[cidx] = s
        coverage_rows.extend([(cidx, k) for k in s])

    # Redundancy pruning
    remaining = set(intersecting_idx)
    removed = []

    # Remove smallest-coverage circles first
    ordered = sorted(remaining, key=lambda c: len(coverage[c]))
    for c in ordered:
        if c not in remaining:
            continue
        others = remaining - {c}
        union_others = set().union(*(coverage[o] for o in others)) if others else set()
        if coverage[c].issubset(union_others):
            remaining.remove(c)
            removed.append(c)

    # Removed_circle_indices = removed[:]
    kept_circles = [circles[i] for i in sorted(remaining)]

    return kept_circles


def define_fields(catalogue_df,
                  fixed_grid=False,
                  circles_df=None,
                  prune_fields=False,
                  batch_size=3,
                  max_passes=5,
                  merge_threshold=3,
                  radius=0.98,
                  central_star_mag_lim=16.,
                  grid_offset=1.0,
                  ra_col='GAIA_RA',
                  de_col='GAIA_DEC',
                  mag_col='GAIA_MAG_G',
                  src_id_col='GAIA_ID',
                  verbose=False,
                  lut_outpath='./resources/WG/'):
    """
    Main function to place circles and assign stars to them in a star catalogue.
    This function places circles either around bright stars or using optimized placement,
    assigns stars to those circles, removes fields without science targets, and
    organizes fields into batches.

    Returns:
        tuple: (stars_in_circles, guide_look_up_df)
    """

    if circles_df is None and fixed_grid:
        print('Using fixed grid approach for field definition.')
        circles = get_science_covering_guide_star_circles_fixed_grid(
            catalogue_df,
            ra_col='GAIA_RA',
            dec_col='GAIA_DEC',
            flag_col='sourceFlag',
            mag_col='GAIA_MAG_G',
            central_star_mag_lim=central_star_mag_lim,
            spacing=grid_offset,
            match_radius=1.0)

    elif circles_df is not None and (fixed_grid or not fixed_grid):
        print('Using provided circles dataframe for field definition.')
        circles = get_science_covering_guide_star_circles_fixed_grid(
            catalogue_df,
            ra_col='GAIA_RA',
            dec_col='GAIA_DEC',
            flag_col='sourceFlag',
            mag_col='GAIA_MAG_G',
            central_star_mag_lim=central_star_mag_lim,
            spacing=grid_offset,
            match_radius=1.0,
            circles_df=circles_df,
            filter_circles_to_science_bbox=False)

    elif circles_df is None and not fixed_grid:
        print('Using automatic (greedy) approach for field definition.')
        circles = get_science_covering_guide_star_circles_greedy(
            catalogue_df,
            ra_col='GAIA_RA',
            dec_col='GAIA_DEC',
            flag_col='sourceFlag',
            mag_col='GAIA_MAG_G',
            central_star_mag_lim=central_star_mag_lim,
            spacing=grid_offset,
            match_radius=1.0)

    # Remove redundant fields (no science targets or covered by other circles)
    if len(circles) == 0:
        logging.error("No fields were placed. Check input catalogue and parameters.")

    if prune_fields:
        pruned_circles = prune_circles(catalogue_df, circles,
                                       ra_col="GAIA_RA", dec_col="GAIA_DEC",
                                       flag_col="sourceFlag", ref_circle_idx=0,
                                       circle_radius_deg=radius)

        circles = pruned_circles

    # Assign stars to circles
    catalogue_df, guide_look_up_df = assign_stars_to_circles(catalogue_df, circles,
                                                             max_passes=max_passes,
                                                             ra_col=ra_col,
                                                             de_col=de_col,
                                                             src_id_col=src_id_col,
                                                             outpath=lut_outpath)

    circles_with_ids = [{ra_col+'_c': c[0], de_col+'_c': c[1], 'fieldId_c': i+1}
                        for i, c in enumerate(circles)]
    circles_df = pd.DataFrame(circles_with_ids)
    if verbose:
        print('----------------------------------------------')
        print(f"Total placed circles: {circles_df.shape[0]}")

    stars_in_circles = catalogue_df[catalogue_df['fieldId'].notnull()]
    stars_in_circles = stars_in_circles.astype({'fieldId': 'int32'})
    if verbose:
        print('----------------------------------------------')
        print('Field defnition completed!')
        count_sc = stars_in_circles[~stars_in_circles['sourceFlag'].isin([9])].shape[0]
        count_gs = stars_in_circles[stars_in_circles['sourceFlag'].isin([9])].shape[0]
        print('Nr. of targets: \n'
              f'Science: {count_sc}\n'
              f'Guide: {count_gs}')

    # Standard batch assignment
    if (batch_size < 3) and (batch_size > 0):
        stars_in_circles, num_batches = assign_batch_id_optimized(
            circles_df, stars_in_circles, batch_size=batch_size,
            merge_threshold=merge_threshold, ra_col=ra_col, de_col=de_col)

    elif batch_size == 3:
        stars_in_circles, num_batches = assign_batch_id_optimized_fast(
            circles_df, stars_in_circles, batch_size=batch_size,
            merge_threshold=merge_threshold, ra_col=ra_col, de_col=de_col)
    else:
        num_batches = -1
    if verbose:
        print('----------------------------------------------')
        print('Batch ID assignment completed!')
        print(f'Nr. of field batches: {num_batches}')
        print('----------------------------------------------')
    stars_in_circles = stars_in_circles.astype({src_id_col: 'str'})

    return stars_in_circles, guide_look_up_df
