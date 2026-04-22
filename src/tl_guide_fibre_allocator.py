import logging

import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs

from tl_protofielding import assign_targprog


def check_guide_fibre_status(fibre_status_xml,
                             fibre_ids_to_check=None):
    """
    Checks if any fibres in the list have a status other than 'enabled' in a given
    fibre status XML file.

    Args:
        fibre_status_xml (str): Path to the XML file.
        fibre_ids_to_check (list of int): List of fibre IDs to check, if None, defaults to
                                          WEAVE guide fibre ids

    Returns:
        dict: {'A': [...], 'B': [...]} with fibre IDs of status not 'enabled'.
    """

    if fibre_ids_to_check is None:
        fibre_ids_to_check = [249, 123, 1005, 879, 753, 627, 501, 375]

    tree = ET.parse(fibre_status_xml)
    root = tree.getroot()

    result = {'A': [], 'B': []}
    fibre_ids_set = set(map(str, fibre_ids_to_check))

    for plate_tag, key in [('PLATE_A', 'A'), ('PLATE_B', 'B')]:
        plate = root.find(plate_tag)
        for fibre in plate.findall('fibre'):
            fid = fibre.attrib['id']
            if fid in fibre_ids_set and fibre.attrib['status'] != 'enabled':
                result[key].append(int(fid))

    return result


def zpn_projection(batch,
                   p_coeff=[203.00898850, 0., 3.03143062, 0., -0.08049567, 0., 0.],
                   ra_col='RAdeg',
                   de_col='DEdeg'):

    # c1 = 205.0 gives marginally better results
    # Define Polynomial distortion coefficients
    dist_coeff = [(2, i + 1, coeff) for i, coeff in enumerate(p_coeff)]

    fields = []
    for field_id in batch['fieldId'].unique():
        # Set field
        field = batch[batch['fieldId'] == field_id]

        # Get central star corrdinates
        central_star = field[field['guide'] == 1]
        central_ra, central_de = central_star[ra_col], central_star[de_col]

        # Calculate vectors relative to the central star for merging onto plate
        field['Vx'] = field[ra_col] - central_ra.values[0]
        field['Vy'] = field[de_col] - central_de.values[0]

        # Convert to world coordinates
        world_coords = np.column_stack((field[ra_col], field[de_col]))
        world_coords = SkyCoord(ra=world_coords[:, 0], dec=world_coords[:, 1], unit='deg')

        # Setup ZPN projection in astropy WCS
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0, 0]
        w.wcs.crval = [central_ra.values[0], central_de.values[0]]
        w.wcs.ctype = ['RA---ZPN', 'DEC--ZPN']
        w.wcs.set_pv(dist_coeff)

        # Calculate ZPN projection
        plate_coords = w.world_to_pixel(world_coords)

        # Stack plate coordinates and flip Y axis to account for telescope mirroring effect
        plate_coords = np.column_stack((plate_coords[0], -plate_coords[1]))

        # Set plate coordinates for each star in the field
        field['X'] = plate_coords[:, 0]
        field['Y'] = plate_coords[:, 1]

        # Append field to fields list and create the projected batch DataFrame
        fields.append(field)
        batch_projected = pd.concat(fields, ignore_index=True)

        # import matplotlib.pyplot as plt
        # from matplotlib.patches import Circle

        # fig, ax = plt.subplots(figsize=(5, 5))

        # circle_plate = Circle((0, 0), 205.0, edgecolor='tab:blue',
        #                      facecolor='none', alpha=0.5)

        # ax.add_patch(circle_plate)

        # ax.scatter(batch_projected[batch_projected['sourceFlag'] == 9]['X'],
        #          batch_projected[batch_projected['sourceFlag'] == 9]['Y'], s=5)
        # for _, f in batch_projected.iterrows():
        #   ax.text(f['X'], f['Y'], f['fieldId'])
        # ax.set_xlim(-210, 210)
        # ax.set_ylim(-210, 210)
        # ax.axis('equal')
        # plt.grid(True)
        # batch_id = batch['batchId'].iloc[0]
        # plt.gca().invert_xaxis()
        # plt.gca().invert_yaxis()
        # plt.title(f'Batch {batch_id} - ZPN Projection')
        # plt.tight_layout()
        # plt.savefig('./temp_plots/zpn_projection.png', dpi=150)
        # plt.show()
        # plt.close('all')

    return batch_projected


def assign_fibre_ids(df,
                     fibre_ids,
                     mag_col='GAIA_MAG_G',
                     src_id='GAIA_ID'):
    """
    Allocate guide-fibre IDs to reachable guide stars, aiming to ensure that each field in a
    batch has at least one central guide and two additional guides (when possible).
    Stars competing for the same fibre are resolved by brightness (the smallest value in
    `mag_col` gets selected). The alglorithm runs until a stable configuration is reached or
    a maximum number of iterations is exceeded. A stable configuration is one where no further
    changes occur in the allocation of fibre IDs. To enable this the allocation workflow
    maintains a live pool of free fibre IDs and, on each pass,
    reconciles the current frame against the initial assignment.
    When multiple candidates contend for an already-used fibre, the brightest star keeps
    (or receives) that fibre. After a stable configuration is reached, one of the remaining
    free fibres is assigned to the central guide star (with exception of fibreid 1005, this
    fibre is not able to reach the centre). Fobre 1005 is therefore prioritized, when
    selecting off-center fibres.

    Args:
        df (pd.DataFrame): Input table of candidate stars. Must contain:
            - a unique identifier column named by `src_id`,
            - a magnitude column named by `mag_col` (smaller = brighter),
            - a writable/nullable `fibreId` column (created if missing),
            - any additional columns your calling code requires for per-field logic
              (e.g., `fieldId`, `guide`), if used upstream.
        fibre_ids (np.ndarray): Array of available guide-fibre IDs that can be assigned.
        mag_col (str, optional): Column with magnitudes used for brightness ranking.
            Defaults to 'gaiaV'.
        src_id (str, optional): Column name holding the star identifier used to match and
            assign fibres. Defaults to 'GAIA_ID'.

    Returns:
        pd.DataFrame: A copy of `df` with the `fibreId` column populated for stars that
        received an allocation. Previously assigned fibres are preserved when possible,
        and conflicts are resolved in favor of the brightest star.

    Notes:
        - The allocation pool excludes any IDs already present in `df['fibreId']`. If the
          special constant `AVOID_CENTRAL_FIBRE_ID` has been assigned anywhere, that ID
          is permanently removed from the free pool and will not be reused.
        - Conflict resolution compares the current working frame with an `initial_assignment`
          snapshot to stabilize assignments across passes. Ensure your calling code defines
          `initial_assignment = df.copy()` (or equivalent) before iterative updates, and
          maintains a consistent `free_fibre_ids` pool across calls if you allocate in
          multiple stages.
        - Brightness is defined by ascending values in `mag_col` (i.e., lower is brighter).
        - The function assumes that input rows are already filtered to “reachable” stars for
          the guide fibres; reachability checks are not performed here.

    Raises:
        KeyError: If required columns (e.g., `src_id`, `mag_col`) are missing from `df`.
        ValueError: If `fibre_ids` is empty or contains non-integer values.
    """

    def _update_free_fibre_ids(df, fibre_ids):
        # Return fibre IDs that are not (yet) allocated
        alloc_fib_ids = df['fibreId'].dropna().astype(int, errors='ignore')
        free = fibre_ids[~np.isin(fibre_ids, alloc_fib_ids)]
        # If 1005 is assigned anywhere, never consider it free again
        if (alloc_fib_ids == _AVOID_CENTRAL_FIBRE_ID).any():
            free = free[free != _AVOID_CENTRAL_FIBRE_ID]
        return free

    def _update_fibres(frame):
        # Compare the current frame with the initial assignment to find common stars
        compare_frame = initial_assignment[initial_assignment['fibreId'].isin(free_fibre_ids)]
        common_stars = frame.merge(compare_frame[[src_id, 'fibreId']], on=src_id,
                                   how='inner', suffixes=('', '_from_compare'))

        for fibre_id in common_stars['fibreId_from_compare'].unique():
            stars_with_same_fibre = common_stars[common_stars['fibreId_from_compare'] ==
                                                 fibre_id]

            # Find the brightest star among the stars with the same fibreId
            brightest_star_index = stars_with_same_fibre[mag_col].idxmin()
            fibre_id_to_assign = stars_with_same_fibre.loc[brightest_star_index,
                                                           'fibreId_from_compare']
            # Assign the fibre ID to the brightest star in the frame
            frame.loc[frame[src_id] ==
                      stars_with_same_fibre.loc[brightest_star_index,
                                                src_id], 'fibreId'] = fibre_id_to_assign

        return frame

    def _assign_last_fibre(df):
        """
        Assigns the last fibre ID to the brightest star in the DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing star data with fibre IDs.

        Returns:
            pd.DataFrame: DataFrame with the last fibre ID assigned to the brightest star.
        """
        # Group fields and find brightest star index in each group
        grouped = df.groupby('fieldId')
        max_indices_per_group = grouped.apply(lambda x: x[mag_col].idxmin() if
                                              len(x) > 1 else None).dropna()
        # Find the index of the star with the highest magnitude among the max indices
        if not max_indices_per_group.empty:
            index_to_drop = df.loc[max_indices_per_group][mag_col].idxmax()
            df_dropped = df.drop(index_to_drop)
        else:
            df_dropped = df.copy()

        return df_dropped

    def _drop_and_update(df):
        """
        Drops rows with NaN fibre IDs and updates the fibre IDs. Incldes the logic to avoid
        dropping fibre ID 1005 once it is assigned. This is due to the fact that fibre 1005
        cannot reach the centre of the field plate, so it should never be assigned as the
        central fibre.

        Args:
            df (pd.DataFrame): DataFrame containing star data with fibre IDs.

        Returns:
            tuple: A tuple containing:
                - df (pd.DataFrame): DataFrame with dropped rows and updated fibre IDs.
                - free_fibre_ids (np.array): Updated array of free fibre IDs.
        """

        df_3_stars = df.dropna(subset=['fibreId'])
        valid_fields = df_3_stars.groupby('fieldId').filter(lambda x: len(x) > 2)

        to_drop = []
        for _, g in valid_fields.groupby('fieldId'):
            # do not drop the star holding 1005 (keeps 1005 from becoming free)
            cand = g[g['fibreId'] != _AVOID_CENTRAL_FIBRE_ID]
            if cand.empty:
                continue  # nothing to drop without touching 1005
            to_drop.append(cand[mag_col].idxmax())  # drop faintest among non-1005

        if to_drop:
            df = df.drop(to_drop)

        free_fibre_ids = _update_free_fibre_ids(df, fibre_ids)
        return df, free_fibre_ids

    # ######################### Main Logic #########################

    # Fibre 1005 cannot reach the centre of the field plate, the variable below ensures that.
    # THIS MUST NOT BE CHANGED!!!
    _AVOID_CENTRAL_FIBRE_ID = 1005

    # Initial assignment of stars with relevant columns
    # before the per-field loop
    df = df.sort_values([mag_col, src_id], kind="mergesort")

    initial_assignment = df[[src_id, 'fieldId', 'fibreId', mag_col,
                             'batchId', 'gateDist', 'sourceFlag']]

    field_count = df['fieldId'].value_counts()

    # Iterate over fields with at least 1 star, sorted by number of stars
    for field_id in field_count[field_count >= 1].sort_values().index:
        # Get DataFrame for the current field and sort by magnitude
        field_df = df[df['fieldId'] == field_id].sort_values([mag_col, src_id],
                                                             kind="mergesort")

        # Handle different stars with the same fibreId by ordering them by brightness
        for index, star in field_df.iterrows():
            bright_dupl = (
                df[(df['fibreId'] == star['fibreId']) & (df['fieldId'] != star['fieldId'])]
                .sort_values([mag_col, src_id], kind="mergesort")
                .groupby('fieldId', sort=False)
                .first())

            if not bright_dupl.empty:
                brightest_in_other_field = bright_dupl.iloc[0]
                other_field_id = bright_dupl.index[0]
                brightest_star = df[(df['fieldId'] == other_field_id) &
                                    (df['fibreId'] != star['fibreId'])][mag_col].min()

                # If the brightest duplicate is brighter than the brightest star in the
                # other field, handle the conflict
                if brightest_in_other_field[mag_col] < brightest_star:
                    fibre_id_dupl = brightest_in_other_field['fibreId']
                    if df.loc[index, 'fibreId'] == fibre_id_dupl:
                        try:
                            df.drop(index, inplace=True)
                            df.loc[df['fibreId'] == star['fibreId'],
                                   'fibreId'] = star['fibreId']

                        except KeyError:
                            continue

                    continue

            # Set fibre ID to NaN for non-brightest duplicates
            df.loc[(df['fibreId'] == star['fibreId']) &
                   (df.index != index), 'fibreId'] = pd.NA

    # Drop duplicate stars based on source ID, keeping the first occurrence
    # (arbitrary choice, keeping 'last' also works)
    if df['fibreId'].isna().any():
        df = df.drop_duplicates(subset=[src_id, 'GAIA_DR'], keep='first')

    # First initial update of fibre IDs acts as intial assignment
    free_fibre_ids = _update_free_fibre_ids(df, fibre_ids)

    if free_fibre_ids.size == 0:
        df, free_fibre_ids = _drop_and_update(df)

    # Stability loop to ensure all stars are assigned a fibre ID
    iters = 0
    iters_max = 20
    stable = False

    while not stable:
        # Save the current state of free fibre IDs
        previous_free_fibre_ids = free_fibre_ids.copy()
        # Update the DataFrame with the current free fibre IDs
        df = _update_fibres(df)
        # Drop and update the DataFrame, and get the new set of free fibre IDs
        df, free_fibre_ids = _drop_and_update(df)
        # Check if the free fibre IDs have stabilized (no changes)
        stable = np.array_equal(free_fibre_ids, previous_free_fibre_ids)
        iters += 1

        # Break out of the loop when it is cycling between two ids continuously.
        # It is very likely that the allocation was successful. But a warning is outputted to
        # advise the user to check the fibre allocation output (verbose=True).
        if iters >= iters_max:
            logging.warning(f"Max iterations reached ({iters_max}) "
                            f"for batch {int(df['batchId'].unique()[0])}. "
                            f"User is advised to visually confirm fibre allocation. "
                            f"Set verbose=True to see more details.")
            break

    # Housekeeping: drop NaN and duplicates (Safety first!)
    df = df.dropna(subset=['fibreId'])
    df = df.sort_values(by=mag_col, ascending=True).drop_duplicates(subset='fibreId',
                                                                    keep='first')
    if df.empty:
        df = None

    # Assign last (central) guide fibre
    if free_fibre_ids.size == 0 and df is not None and not df.empty:
        df = _assign_last_fibre(df)
        # Update ids to account for last fibre assignment
        free_fibre_ids = _update_free_fibre_ids(df, fibre_ids)
        # central_fibre_id = free_fibre_ids[-1]
        central_fibre_id = next((fid for fid in reversed(free_fibre_ids) if fid != 1005), None)

        if central_fibre_id is None:
            logging.warning(f"No available fibre ID found other than 1005!")

    else:
        central_fibre_id = next((fid for fid in reversed(free_fibre_ids) if fid != 1005), None)
        if central_fibre_id is None:
            print(f"No available central fibre found other than fibreid=1005 "
                  f"in batch {int(df['batchId'].unique()[0])}. "
                  "Recommend setting maximize_batches=True")
            central_fibre_id = 1005

    return df, central_fibre_id


def assign_guide_fibres(merged_df,
                        max_gate_angle,
                        hemisphere,
                        plate,
                        guide_mag_limit=None,
                        mag_col='gaiaV',
                        ra_col='ra',
                        de_col='dec',
                        src_id='sourceId',
                        guide_reach_factor=12.0,
                        broken_fibres=None):

    """
    Assign WEAVE guide-fibre IDs to guide stars field-by-field using the instrument
    gate geometry, plate, hemisphere, and a broken-fibre mask. Only guide stars
    (`sourceFlag == 9`) are considered. For each field, the routine:
      1) identifies the central guide (row(s) with `guide == 1`) and switches to the
         superposition-frame coordinates (`Vx`, `Vy`) for geometry,
      2) determines the target gate positions/angles for the active plate/hemisphere,
         masking out any broken fibres,
      3) filters candidate guides (`guide == 0`) to those inside each gate's triangular
         reach region (limited by `max_gate_angle` and `guide_reach_factor`), and
      4) calls `assign_fibre_ids(...)` to allocate specific `fibreId`s to the valid
         candidates, resolving conflicts by brightness (`mag_col`).

    The output is the per-batch DataFrame with guide assignments in `fibreId`: central
    guides keep their role and are added back; non-central assigned guides are tagged
    with `guide == 2`. If a field has no valid assignments, it is skipped.

    Args:
        merged_df (pd.DataFrame): Input table containing **only the batch under
            consideration**. Must include at least:
            - `fieldId` (int): field identifier,
            - `sourceFlag` (int): 9 for guide stars,
            - `guide` (int): 1 for the central guide of the field, 0 for other candidates,
            - `X`, `Y` (float): star coordinates in the superposition frame,
            - `Vx`, `Vy` (float): central guide coordinates in the same frame,
            - columns named by `mag_col`, `ra_col`, `de_col`, `src_id` (see below).
        max_gate_angle (float): Half-opening angle (degrees) of a gate's reach cone.
        hemisphere (str): Observation hemisphere, 'N' (north) or 'S' (south). Affects
            fibre-ID ordering/rotation.
        plate (str): Plate identifier, 'A' or 'B'. Used with `broken_fibres`.
        guide_mag_limit (float, optional): Maximum magnitude for guide candidates. If
            provided, fainter stars are excluded before allocation. Defaults to None.
        mag_col (str, optional): Magnitude column used for ranking (lower = brighter).
            Defaults to 'gaiaV'.
        ra_col (str, optional): Right Ascension column (degrees). Defaults to 'ra'.
        de_col (str, optional): Declination column (degrees). Defaults to 'dec'.
        src_id (str, optional): Source identifier column used when allocating fibres.
            Defaults to 'sourceId'.
        guide_reach_factor (float, optional): Multiplier used to offset the gate-cone
            edges to reduce fibre crossings (acts on cone length correction).
            Defaults to 12.0 (do not change unless you understand the geometry).
        broken_fibres (dict | None, optional): Mapping like
            `{ 'A': [ids...], 'B': [ids...] }` listing broken fibre IDs per plate.
            Broken fibres are masked out and never allocated. Defaults to None
            (treated as no broken fibres).

    Returns:
        pd.DataFrame | None: The batch DataFrame with assigned guide fibres:
            - `fibreId` set on selected non-central guides,
            - central guides re-added with their `fibreId` (returned by `assign_fibre_ids`),
            - non-central, fibre-assigned guides marked with `guide == 2`,
            - columns `Vx`, `Vy` dropped from the output.
            Returns `None` if no valid fields produced assignments.

    Notes:
        - The function **filters to `sourceFlag == 9` internally**; science targets are
          not considered here.
        - Fibre IDs are hard-coded for WEAVE and rotated between hemispheres; a small
          fixed offset (`fibre_pa_offset`) is applied to target angles.
        - Reach regions are modeled as triangles (gate vertex plus two rays at
          `±max_gate_angle`), and candidate guides are selected if their (`X`, `Y`)
          lies inside the triangle.
        - Final per-fibre selection among candidates is brightness-based and delegated
          to `assign_fibre_ids(...)`, which also preserves/stabilizes allocations.

    Raises:
        KeyError: If required columns are missing (e.g., `fieldId`, `sourceFlag`,
          `guide`, `X`, `Y`, `Vx`, `Vy`, `mag_col`, `src_id`).
        ValueError: If `hemisphere` not in {'N','S'} or `plate` not in {'A','B'}.
    """

    def _is_point_in_triangle(pt, A, B, C):
        """
        Determine if a point is inside a triangle formed by points A, B, and C.

        Args:
            pt (tuple): Point to check, given as (x, y).
            A (tuple): First vertex of the triangle, given as (x, y).
            B (tuple): Second vertex of the triangle, given as (x, y).
            C (tuple): Third vertex of the triangle, given as (x, y).

        Returns:
            bool: True if the point is inside the triangle or on its edge, False otherwise.
        """

        # Creating numpy arrays for points and vectors
        A, B, C, pt = np.array(A), np.array(B), np.array(C), np.array(pt)
        AB, BC, CA = B - A, C - B, A - C
        AP, BP, CP = pt - A, pt - B, pt - C

        # Computing cross products using numpy
        cross_AB_AP = np.cross(AB, AP)
        cross_BC_BP = np.cross(BC, BP)
        cross_CA_CP = np.cross(CA, CP)

        # Checking if the point is on the same side of each edge or on the edge
        return (cross_AB_AP >= 0 and cross_BC_BP >= 0 and cross_CA_CP >= 0) or \
            (cross_AB_AP <= 0 and cross_BC_BP <= 0 and cross_CA_CP <= 0)

    def _euclidean_distance(x1, y1, x2, y2):
        """
        Calculate euclidian distance between two points (x1, y1) & (x2, y2).

        Args:
            x1 (float): X coordinate of the first point.
            y1 (float): Y coordinate of the first point.
            x2 (float): X coordinate of the second point.
            y2 (float): Y coordinate of the second point.

        Returns:
            float: Euclidean distance between the two points.
        """

        return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    def _calculate_angle(A, B):
        """
        Calculate the angle in degrees between two vectors A and B using the dot product.

        Args:
            A (np.ndarray): First vector.
            B (np.ndarray): Second vector.

        Returns:
            float: Angle in degrees between the two vectors.
        """
        dot_product = np.dot(A, B)
        norm_a = np.linalg.norm(A)
        norm_b = np.linalg.norm(B)
        cos_theta = dot_product / (norm_a * norm_b)

        # Ensuring the value is within the valid range for arccos, due to floating point errors
        cos_theta = min(1, max(cos_theta, -1))
        angle = np.arccos(cos_theta)

        # Convert from radians to degrees
        angle_deg = np.degrees(angle)

        return angle_deg

    def _filter_single_stars(df):
        """
        Filters out fields that only have one star.

        Args:
            df (pd.DataFrame): DataFrame containing star data with 'fieldId'
                               and 'fibreId' columns.

        Returns:
            pd.DataFrame: DataFrame with fields that have more than one star.
        """

        single_star_counts = df.groupby('fieldId').size()
        ss_field_ids = single_star_counts[single_star_counts == 1].index.to_list()
        single_stars_df = df[df['fieldId'].isin(ss_field_ids)]

        ss_fibre_ids = single_stars_df['fibreId'].unique()
        mask = ~df['fibreId'].isin(ss_fibre_ids) | df['fieldId'].isin(ss_field_ids)

        return df[mask]

    # ######################### Main Logic #########################

    batch_frames = []
    central_stars = []
    fibre_pa_offset = 0.5442  # Offset of the fibre tier relative to the unit circle angles

    # Hardcoded guide fibre ids for WEAVE (DO NOT CHANGE!!!)
    fibre_ids_n = np.array([249, 123, 1005, 879, 753, 627, 501, 375])
    fibre_ids_s = np.roll(fibre_ids_n, -4)  # rolled by 4, accounting for negative Declination
    angles = np.array([0, 45, 90, 135, 180, 225, 270, 315])

    # Check if broken fibres are provided, if not, set to empty dict.
    # Also take hemisphere into account.
    if hemisphere == 'N':
        broken_fibre_mask = ~np.isin(fibre_ids_n, broken_fibres.get(plate, []))
        fibre_ids = fibre_ids_n[broken_fibre_mask]
        target_angles_deg = angles + fibre_pa_offset
        target_angles_deg = target_angles_deg[broken_fibre_mask]

    elif hemisphere == 'S':
        broken_fibre_mask = ~np.isin(fibre_ids_s, broken_fibres.get(plate, []))
        fibre_ids = fibre_ids_s[broken_fibre_mask]
        target_angles_deg = angles + fibre_pa_offset
        target_angles_deg = target_angles_deg[broken_fibre_mask]

    # import matplotlib.pyplot as plt
    # from matplotlib.patches import Circle
    # fig, ax = plt.subplots(figsize=(5, 5))
    # field_radius = 205.0
    # circle = Circle((0, 0), field_radius, fill=False)

    # Assure that only guide stars are used for guide fibre assignment
    merged_df = merged_df[merged_df['sourceFlag'] == 9]

    for field_id in merged_df['fieldId'].unique():

        field_stars = merged_df[(merged_df['fieldId'] == field_id) & (merged_df['guide'] == 0)]
        central_star = merged_df[(merged_df['fieldId'] == field_id) &
                                 (merged_df['guide'] == 1)]

        # Switch to superposition field coordinates
        central_x = central_star['Vx'].values[0]
        central_y = central_star['Vy'].values[0]

        if central_star.empty:
            # Skip this field if there is no central star
            logging.warning('No central star found for field '
                            f'{field_id}. Skipping field.')
            continue

        central_stars.append(central_star)

        # Filtering other stars based on angular distance
        max_dist_to_center = field_stars.apply(lambda row: _euclidean_distance(
                                               central_x, central_y, row['X'],
                                               row['Y']), axis=1).max()

        min_dist_from_center = max_dist_to_center / 8

        valid_stars = field_stars[field_stars.apply(lambda row: _euclidean_distance(
                                                    central_x, central_y, row['X'],
                                                    row['Y']) > min_dist_from_center,
                                                    axis=1)]

        # ax.scatter(field_stars['X'], field_stars['Y'], s=2, alpha=1.0)

        # Reconstruct field/fibre gate geometry check if stars are within the reachable
        # region of the guide fibres. Then assign the fibre IDs to the valid stars.

        gate_dist = 273.76237
        field_frames = []

        for angle, fibre_id in zip(target_angles_deg, fibre_ids):
            gate_x = central_x + np.cos(np.deg2rad(angle)) * gate_dist
            gate_y = central_y + np.sin(np.deg2rad(angle)) * gate_dist
            mga_left = np.deg2rad(180 + angle + max_gate_angle)
            mga_right = np.deg2rad(180 + angle - max_gate_angle)

            # Ensures minimal crossing of fibres by shifting the cone border away from centre
            correction = ((gate_dist/np.cos(np.deg2rad(max_gate_angle))) -
                          gate_dist) * -guide_reach_factor

            # Endpoints for cone borders
            l_end_ra = gate_x + (gate_dist+correction) * np.cos(mga_left * u.rad)
            l_end_de = gate_y + (gate_dist+correction) * np.sin(mga_left * u.rad)
            r_end_ra = gate_x + (gate_dist+correction) * np.cos(mga_right * u.rad)
            r_end_de = gate_y + (gate_dist+correction) * np.sin(mga_right * u.rad)

            gate_coords = (gate_x, gate_y)
            left_coords = (l_end_ra.value, l_end_de.value)
            right_coords = (r_end_ra.value, r_end_de.value)

            # ax.scatter(gate_x, gate_y, s=10)
            # ax.text(gate_x, gate_y, rf'{fibre_id}')
            # ax.plot([gate_x, l_end_ra.value], [gate_y, l_end_de.value], 'k--')
            # ax.plot([gate_x, r_end_ra.value], [gate_y, r_end_de.value], 'k--')
            # print('TEEEEEST')
            # print(gate_coords)
            # print(left_coords)
            # print(right_coord
            # print('-----------------')

            mask = valid_stars.apply(lambda row: _is_point_in_triangle((row['X'],
                                                                       row['Y']),
                                                                       gate_coords,
                                                                       left_coords,
                                                                       right_coords),
                                     axis=1).squeeze()

            if mask.shape != ():
                valid_stars_filt = valid_stars[mask]
                if guide_mag_limit is not None:
                    try:
                        valid_stars_filt = valid_stars_filt[valid_stars_filt[mag_col] <=
                                                            guide_mag_limit]
                    except KeyError:
                        logging.warning(
                            f"No guides in fibre reach at Batch-Id:"
                            f"{merged_df['batchId'].iloc[0]}")
                        continue

                field_distances = valid_stars_filt.apply(lambda row: _euclidean_distance(
                                                         gate_coords[0], gate_coords[1],
                                                         row['X'], row['Y']), axis=1)

                gate_coords_np = np.array([gate_coords[0], gate_coords[1]])

                # print(gate_coords_np, fibre_id)
                field_angles = valid_stars_filt.apply(lambda row:
                                                      _calculate_angle(
                                                        np.array([row['X'], row['Y']]),
                                                        gate_coords_np), axis=1)

                valid_stars_filt['gateDist'] = field_distances
                valid_stars_filt['gateAngle'] = field_angles
                valid_stars_filt['fibreId'] = fibre_id

                field_frames.append(valid_stars_filt.drop_duplicates(subset=['fibreId'],
                                                                     keep='first'))
            else:
                continue

        if len(field_frames) > 0:
            combined_field_frame = pd.concat(field_frames).sort_values(by=mag_col,
                                                                       ascending=False)
            # print(combined_field_frame[[src_id, 'fieldId', 'fibreId', mag_col,
            #                            'batchId']])

        else:
            continue
        batch_frames.append(combined_field_frame)

    # ax.add_patch(circle)
    # ax.grid(True)
    # ax.axis('equal')
    # plt.title(rf'{merged_df["batchId"].iloc[0]} - {field_id}')
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.tight_layout()
    # plt.savefig('./temp_plots/guide_fibre_assignment.png', dpi=150)
    # plt.close('all')

    if len(batch_frames) > 0:
        # Create batch frame from individual fields
        combined_batch_df = pd.concat(batch_frames)

        # combined_batch_df = _filter_single_stars(combined_batch_df)

        # print('HERE0', combined_batch_df[[src_id, 'batchId', 'fieldId', 'fibreId', 'guide',
        # 'sourceFlag']].head())
        # print(combined_batch_df[[src_id, 'batchId', 'fieldId', 'fibreId',
        # 'guide', 'sourceFlag']])
        combined_batch_df, central_fibre = assign_fibre_ids(combined_batch_df, fibre_ids,
                                                            mag_col=mag_col, src_id=src_id)

        for star in central_stars:
            star['fibreId'] = central_fibre

        central_stars_df = pd.concat(central_stars)

        drop_columns = ['Vx', 'Vy']
        central_stars_df = central_stars_df.drop(columns=drop_columns[2:])

        if combined_batch_df is not None:
            combined_batch_df = combined_batch_df.drop(columns=drop_columns)
            combined_batch_df['guide'] = 2
            combined_batch_df = pd.concat([combined_batch_df, central_stars_df],
                                          ignore_index=True)

            return combined_batch_df

        else:
            print("combined_batch_df is None. Continuing with next batch.")
            pass
    else:
        print("No valid fields found for guide fibre assignment. Continuing with next batch.")
        pass


def maximize_final_batches(df,
                           master_df,
                           max_field_distance=4.0,
                           ra_col='ra',
                           de_col='dec',
                           src_id='sourceId',
                           mag_col='gaiaV',
                           guide_mag_limit=15.0,
                           max_gate_angle=13.8,
                           hemisphere='N',
                           plate='A',
                           guide_reach_factor=12.0,
                           broken_fibres=None):
    """
    Maximizes the final batches by trying to find as many valid 3-field configurations
    as possible. If this is not possible anymore then it tries to find pairs of fields
    and finally single fields. The function will return a DataFrame with the assigned
    fields, guide stars and their scores. Batch IDs, scores and field scores are dynamically
    updated during the process.

    Args:
        df (pd.DataFrame): DataFrame containing the fields to be assigned to batches.
        master_df (pd.DataFrame): Initial DataFrame containing all fields and guide stars.
        max_field_distance (float): Maximum distance in degrees between fields.
                                    Defaults to 4.0.
        ra_col (str): RA column in the dataframe. Defaults to 'ra'.
        de_col (str): DEC column in the dataframe. Defaults to 'dec'.
        src_id (str): Source ID column in the dataframe. Defaults to 'sourceId'.
        mag_col (str): Magnitude column in the dataframe. Defaults to 'gaiaV'.
        guide_mag_limit (float): Maximum magnitude (Gaia G mag) of the guide stars.
                                 Defaults to 15.0.
        max_gate_angle (float): Max reachable angle of the guide fibres in degrees.
                                Defaults to 13.8. Best not to change!!!!
        hemisphere (str): Observational hemisphere, either 'north' or 'south'.
                          Defaults to 'north'.
        plate (str): Plate identifier, either 'A' or 'B'. Defaults to 'A'.
        guide_reach_factor (float): Factor defining the region reacbale for the guide fibres.
                                    This avoids crossing of the fibres. Defaults to 12.0.
                                    Best not to change!!!!
        broken_fibres (dict, optional): Dictionary containing broken fibre IDs for each plate.
                                        This is created using the check_fibre_status function.
                                        If not provided, it will be set to None, resulting in
                                        an all guide-fibres-enabled case.

    Raises:
        ValueError: If the number of fields in a batch is higher than 3.
        ValueError: If no available fibre ID is found other than 1005.

    Returns:
        pd.DataFrame: DataFrame containing the assigned fields with their batch ID,
                      guide stars, and scores.
    """

    from itertools import combinations

    # ############################ Helper functions ############################

    def _haversine(ra1, dec1, ra2, dec2):

        """
        Calculates the Haversine distance in degrees between two points
        on the celestial sphere.

        Args:
            ra1 (float): Right Ascension of the first point in degrees.
            dec1 (float): Declination of the first point in degrees.
            ra2 (float): Right Ascension of the second point in degrees.
            dec2 (float): Declination of the second point in degrees.

        Returns:
            float: Haversine distance in degrees between the two points.
        """

        ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
        d_ra = ra2 - ra1
        d_dec = dec2 - dec1
        a = np.sin(d_dec / 2)**2 + np.cos(dec1) * np.cos(dec2) * np.sin(d_ra / 2)**2

        return 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)) * (180/np.pi)

    def _within_distance(fields, max_dist):

        """
        Checks if all fields are within a specified maximum distance from each other.

        Args:
            fields (list): List of fields, each containing 'ra' and 'dec' keys.
            max_dist (float): Maximum distance in degrees to check between fields.

        Returns:
            bool: True if all fields are within the maximum distance, False otherwise.
        """

        coords = SkyCoord(ra=[f[ra_col] for f in fields] * u.deg,
                          dec=[f[de_col] for f in fields] * u.deg)
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                if coords[i].separation(coords[j]).deg > max_dist:
                    return False

        return True

    def _assign_batch_and_guides(field_ids, next_batch_id):

        """
        Assigns a batch ID and guide stars to the fields specified by field_ids.

        Args:
            field_ids (list): List of field IDs to be assigned to the batch.
            next_batch_id (int): The next batch ID to be assigned.

        Returns:
            pd.DataFrame: DataFrame containing the assigned fields with their batch ID
                    and guide stars.
        """

        subset = df[df['fieldId'].isin(field_ids)]
        guides = master_df[(master_df['fieldId'].isin(field_ids)) &
                           (~master_df['guide'].isin([1, 2])) &
                           (master_df['sourceFlag'] == 9)]
        combined = pd.concat([subset, guides], ignore_index=True)
        combined = combined.copy()
        combined['batchId'] = next_batch_id

        try:
            projected = zpn_projection(combined, ra_col=ra_col, de_col=de_col)
            guide_star_df = assign_guide_fibres(projected,
                                                hemisphere=hemisphere,
                                                plate=plate,
                                                guide_reach_factor=guide_reach_factor,
                                                guide_mag_limit=guide_mag_limit,
                                                max_gate_angle=max_gate_angle,
                                                src_id=src_id,
                                                mag_col=mag_col,
                                                ra_col=ra_col,
                                                de_col=de_col,
                                                broken_fibres=broken_fibres)
            if guide_star_df is not None:
                guide_star_df = calculate_field_score(guide_star_df)
                if (guide_star_df['fieldScore'] >= 3).all():
                    guide_star_df = calculate_batch_score(guide_star_df)
                    return guide_star_df

        except Exception:
            logging.warning(f"Failed to assign batch and guides for fields: {field_ids}")
            return None

        # Explicit fallback in case guide_star_df is None or scores are too low
        return None

    def _try_pairings(remaining_df,
                      start_batch_id):
        """
        Attempts to pair fields from the remaining DataFrame into batches of two.

        Args:
            remaining_df (pd.DataFrame): DataFrame containing fields to be paired.
            start_batch_id (int): The starting batch ID for the pairings.

        Returns:
            tuple: A tuple containing:
                - results (list): List of DataFrames with paired fields.
                - used (set): Set of field IDs that have been used in pairings.
                - next_batch_id (int): The next batch ID to be used.
        """

        central = remaining_df[remaining_df['guide'] == 1].drop_duplicates('fieldId')
        ra_idx = central.columns.get_loc(ra_col)
        de_idx = central.columns.get_loc(de_col)

        used = set()
        results = []
        batch_id = start_batch_id
        for f1, f2 in combinations(central.itertuples(index=False), 2):
            if f1.fieldId in used or f2.fieldId in used:
                continue
            if max_field_distance is not None:
                sep = _haversine(f1[ra_idx], f1[de_idx], f2[ra_idx], f2[de_idx])
                if sep > max_field_distance:
                    continue
            res = _assign_batch_and_guides([f1.fieldId, f2.fieldId], batch_id)
            if res is not None:
                used.update([f1.fieldId, f2.fieldId])
                results.append(res)
                batch_id += 1

        return results, used, batch_id

    def _final_cleanup(df_combined):
        """
        Cleans up the combined DataFrame by removing invalid pairs and updating scores.
        Used after all pairings are done.

        Args:
            df_combined (pd.DataFrame): DataFrame containing the combined pairs data.

        Returns:
            pd.DataFrame: Cleaned DataFrame with the final updated scores and valid pairs.
        """

        invalid_batch_ids, overlap_ids = _get_invalid_pair_batch_ids(df_combined)

        df_combined = _switch_central_fibres(df_combined, overlap_ids, plate=plate,
                                             broken_fibres=broken_fibres)
        df_combined = _redo_fibre_assignment(df_combined, master_df, invalid_batch_ids)
        invalid_offset_ids = _get_invalid_offsets(df_combined, ra_col, de_col,
                                                  max_sep=max_field_distance)
        df_combined = _split_invalid_batches(df_combined, invalid_offset_ids)
        df_combined = _update_batch_score(df_combined)
        
        for col in ["guide", "batchScore", "fieldScore"]:
            # turn '' / ' ' / None / 'nan' into proper missing
            df_combined[col] = (
                df_combined[col]
                .replace(r"^\s*$", pd.NA, regex=True)      # empty/whitespace strings -> NA
                .replace({"nan": pd.NA, "None": pd.NA})    # common string junk -> NA
            )
            # coerce anything else non-numeric to reminder missing
            df_combined[col] = pd.to_numeric(df_combined[col], errors="coerce")

        df_combined = df_combined.astype({'guide': 'Int64',
                                          'batchScore': 'Int64',
                                          'fieldScore': 'Int64'})

        return df_combined

    def _update_batch_score(df):
        """
        Updates the batch score based on the field scores.

        Args:
            df (pd.DataFrame): DataFrame containing the pairs data.

        Raises:
            ValueError: If the number of fields in a batch is higher than 3.

        Returns:
            pd.DataFrame: DataFrame with updated 'fieldScore' and 'batchScore' columns.
        """

        batch_groups = df.groupby('batchId')
        for batchId, group in batch_groups:
            count = group['fieldId'].nunique()
            if (group['fieldScore'] >= 3).all():
                if count == 1:
                    df.loc[df['batchId'] == batchId, 'batchScore'] = 3
                elif count == 2:
                    df.loc[df['batchId'] == batchId, 'batchScore'] = 33
                elif count == 3:
                    df.loc[df['batchId'] == batchId, 'batchScore'] = 333
                else:
                    raise ValueError(f"Number of fields in batch higher than 3! "
                                     f"(BatchId={batchId}, Fields={count})")

        return df.astype({'fieldScore': int, 'batchScore': int})

    def _get_invalid_pair_batch_ids(df):
        """
        Identifies invalid batches based on fibre ID conditions.

        Args:
            df (pd.DataFrame): DataFrame containing the pairs data with columns:

        Returns:
            tuple: A tuple containing:
                   - invalid_batches (list): List of batch IDs that do not meet the conditions.
                   - overlap_batches (list): List of batch IDs that have special
                                             overlap conditions.
        """

        invalid_batches = []  # Batches that fail any condition
        overlap_batches = []  # Batches where special overlap conditions are satisfied

        for batch in df[df['guide'].isin([1])]['batchId'].unique():
            df_batch = df[df['batchId'] == batch]
            fields = df_batch['fieldId'].unique()

            if len(fields) < 3:
                continue

            field1 = (
                df_batch[
                    (df_batch['fieldId'] == fields[0]) &
                    (df_batch['guide'].isin([1, 2]))
                ][['batchId', 'guide', 'fieldId', 'fibreId']]
                .reset_index(drop=True)
            )

            field2 = (
                df_batch[
                    (df_batch['fieldId'] == fields[1]) &
                    (df_batch['guide'].isin([1, 2]))
                ][['batchId', 'fieldId', 'guide', 'fibreId']]
                .reset_index(drop=True)
            )

            # Check that central guide fibreId is the same
            guide1_valid = (
                field1[field1['guide'] == 1]['fibreId'].reset_index(drop=True) ==
                field2[field2['guide'] == 1]['fibreId'].reset_index(drop=True)
            ).all()

            # Check that off-centre guide fibreIds do not overlap
            guide2_valid = (
                set(field1[field1['guide'] == 2]['fibreId']) &
                set(field2[field2['guide'] == 2]['fibreId']) == set()
            )

            # Check for overlap between guide == 1 and guide == 2
            guide1_fibreIds = (
                set(field1[field1['guide'] == 1]['fibreId']) |
                set(field2[field2['guide'] == 1]['fibreId'])
            )
            guide2_fibreIds = (
                set(field1[field1['guide'] == 2]['fibreId']) |
                set(field2[field2['guide'] == 2]['fibreId'])
            )
            overlap = guide1_fibreIds & guide2_fibreIds

            # Ensure guide == 2 fibreIds do not overlap
            no_guide2_overlap = (
                set(field1[field1['guide'] == 2]['fibreId']) &
                set(field2[field2['guide'] == 2]['fibreId']) == set()
            )

            # Add to overlap_batches if special conditions are satisfied
            if (overlap or not guide1_valid) and no_guide2_overlap:
                overlap_batches.append(batch)
            elif not (guide1_valid and guide2_valid):
                invalid_batches.append(batch)

        return invalid_batches, overlap_batches

    def _switch_central_fibres(pairs_df,
                               overlap_batch_ids,
                               plate='A',
                               broken_fibres=None):
        """
        Reassigns central guide fibre IDs for overlapping batches.

        Args:
            pairs_df (pd.DataFrame): DataFrame containing the pairs data.
            overlap_batch_ids (list): List of batch IDs where fibre IDs need to be switched.
            plate (str): Plate identifier ('A' or 'B'). Default is 'A'.
            broken_fibres (dict): Dictionary containing broken fibre IDs for each plate. If
                                  None, uses default fibre IDs, all fibres are treated as
                                  "enabled"
        Raises:
            ValueError: If an invalid plate value is provided.

        Returns:
            _type_: _description_
        """

        fibre_ids = np.array([249, 123, 1005, 879, 753, 627, 501, 375])

        if broken_fibres is None:
            broken_fibres = {'A': [], 'B': []}

        if plate == 'A':
            # Define available fibre IDs for Plate A
            broken_fibre_mask = ~np.isin(fibre_ids, broken_fibres.get(plate, []))
            fibre_ids = fibre_ids[broken_fibre_mask]
        elif plate == 'B':
            # Define available fibre IDs for Plate B
            broken_fibre_mask = ~np.isin(fibre_ids, broken_fibres.get(plate, []))
            fibre_ids = fibre_ids[broken_fibre_mask]

        else:
            raise ValueError("Invalid plate value. Use 'A' or 'B'.")

        # Iterate through overlap_batch_ids to modify pairs_df
        for batch in overlap_batch_ids:
            # Extract the batch data
            batch_df = pairs_df[pairs_df['batchId'] == batch].copy()

            # Get guide == 2 fibre IDs to avoid conflicts
            guide2_fibre_ids = batch_df[batch_df['guide'] == 2]['fibreId'].unique()

            # Find a replacement fibre ID for guide == 1
            replacement_found = False
            for fibre_id in fibre_ids:
                if fibre_id not in guide2_fibre_ids and fibre_id != 1005:
                    replacement_fibre_id = fibre_id
                    replacement_found = True
                    break

            if not replacement_found:
                # Raise a warning then continue
                print(f"WARNING: No valid replacement fibre ID found for batch {batch}.")
                continue

            # Update the central guide (guide == 1) fibre ID directly in pairs_df
            pairs_df.loc[(pairs_df['batchId'] == batch) &
                         (pairs_df['guide'] == 1), 'fibreId'] = replacement_fibre_id

        return pairs_df

    def _split_invalid_batches(pairs_df,
                               invalid_batch_ids):
        """
        Splits invalid batches into two separate batches, each containing one field.
        Ensures the new batchId is unique.

        Parameters:
        pairs_df (pd.DataFrame): The DataFrame containing the pairs data.
        invalid_batch_ids (list): List of batchIds to be split.

        Returns:
        pd.DataFrame: Updated pairs_df with split batches.
        """

        # Get the maximum existing batchId to ensure uniqueness for new batchIds
        max_batch_id = pairs_df['batchId'].max()

        # Iterate over each invalid batchId
        for batch_id in invalid_batch_ids:
            # Extract the invalid batch
            batch_df = pairs_df[pairs_df['batchId'] == batch_id]

            # Get unique fields in the batch
            field_ids = batch_df['fieldId'].unique()

            if len(field_ids) != 2:
                raise ValueError(f"Batch {batch_id} does not have exactly"
                                 "two fields to split.")

            # Assign one field to keep the old batchId
            pairs_df.loc[
                (pairs_df['batchId'] == batch_id) &
                (pairs_df['fieldId'] == field_ids[0]),
                'batchId'] = batch_id

            # Assign a new batchId to the other field
            new_batch_id = max_batch_id + 1
            pairs_df.loc[
                (pairs_df['batchId'] == batch_id) &
                (pairs_df['fieldId'] == field_ids[1]),
                'batchId'] = new_batch_id

            # Increment the max_batch_id for the next new batch
            max_batch_id += 1

        return pairs_df

    def _reset_fibre_ids(pairs_df,
                         invalid_batch_ids):
        """
        Resets the fibreIds for fields in invalid batches.

        Parameters:
        pairs_df (pd.DataFrame): The DataFrame containing the pairs data.
        invalid_batch_ids (list): List of batchIds to be reset.

        Returns:
        pd.DataFrame: Updated pairs_df with reset fibreIds.
        """

        # Iterate over each invalid batchId
        for batch_id in invalid_batch_ids:
            # Extract the invalid batch
            batch_df = pairs_df[pairs_df['batchId'] == batch_id]

            # Reset the fibreId for each field in the batch
            for field_id in batch_df['fieldId'].unique():
                pairs_df.loc[(pairs_df['batchId'] == batch_id) &
                             (pairs_df['fieldId'] == field_id), 'fibreId'] = pd.NA #np.nan
                pairs_df.loc[(pairs_df['batchId'] == batch_id) &
                             (pairs_df['fieldId'] == field_id) &
                             (pairs_df['guide'] == 2), 'guide'] = pd.NA #np.nan

        return pairs_df

    def _redo_fibre_assignment(pairs_df,
                               master_df,
                               invalid_batch_ids):
        """
        Redoes the fibre assignment for fields in invalid batches.

        Parameters:
        pairs_df (pd.DataFrame): The DataFrame containing the pairs data.
        invalid_batch_ids (list): List of batchIds to be reassigned.

        Returns:
        pd.DataFrame: Updated pairs_df with new fibre assignments.
        """

        guide_dfs = []
        field_ids_to_remove = []

        print('Recalculating fibre assignment for invalid batches...')
        for batch_id in invalid_batch_ids:
            # Extract the invalid batch
            batch_no_guide = pairs_df[pairs_df['batchId'] == batch_id]
            guides = master_df[
                (master_df['fieldId'].isin(batch_no_guide['fieldId'].unique())) &
                (~master_df['guide'].isin([1, 2])) &
                (master_df['sourceFlag'] == 9)
            ]

            field_ids_to_remove.extend(guides['fieldId'].unique())
            # print(batch_no_guide["fieldId"].unique())
            # print(guides[['batchId', 'fieldId', 'guide', 'sourceFlag']])
            batch = pd.concat([batch_no_guide, guides], ignore_index=True).astype({src_id:
                                                                                  'int64'})
            batch.loc[:, 'batchId'] = batch_id

            merged_fields = zpn_projection(batch=batch, ra_col=ra_col, de_col=de_col)
            # print(merged_fields[merged_fields['sourceFlag'] != 9][['batchId', "fieldId",
            # "guide", 'sourceFlag']])
            guide_star_df = assign_guide_fibres(merged_fields,
                                                hemisphere=hemisphere,
                                                plate=plate,
                                                guide_reach_factor=guide_reach_factor,
                                                guide_mag_limit=guide_mag_limit,
                                                max_gate_angle=max_gate_angle,
                                                src_id=src_id,
                                                mag_col=mag_col,
                                                ra_col=ra_col,
                                                de_col=de_col,
                                                broken_fibres=broken_fibres)
            guide_star_df = calculate_field_score(guide_star_df)
            guide_star_df = calculate_batch_score(guide_star_df)
            batch_score = guide_star_df['batchScore'].unique()[0]
            if guide_star_df is not None and (batch_score == 3 or batch_score == 33):
                guide_dfs.append(guide_star_df)
                # print(guide_star_df[['batchId', "fieldId", "guide", "fibreId",
                #                     'fieldScore', 'batchScore']])
            elif guide_star_df is not None and batch_score == 32:
                guide_star_df = guide_star_df[guide_star_df['fieldScore'].isin([3])]
                guide_star_df = calculate_batch_score(guide_star_df)
                guide_dfs.append(guide_star_df)
                # print(guide_star_df[['batchId', "fieldId", "guide", "fibreId",
                #                    'fieldScore', 'batchScore']])
            elif guide_star_df is None:
                logging(f"Second guide star assignment failed at Batch-Id: {batch_id}")
        # print(pairs_df[(pairs_df['guide'].isin([1,2]))][['batchId', "fieldId", "guide",
        # "fibreId", 'fieldScore', 'batchScore', src_id]])
        # combined_guide_frame = pd.concat(guide_dfs).astype({src_id: 'str'})
        if not guide_dfs:
            print("Warning: No valid batches after reassignment.")
            return pairs_df

        combined_guide_frame = pd.concat(guide_dfs).astype({src_id: 'str'})

        pairs_df = pairs_df[~pairs_df['fieldId'].isin(field_ids_to_remove)]
        pairs_df = pd.concat([pairs_df, combined_guide_frame], ignore_index=True)
        num_redid_batches = combined_guide_frame['batchId'].unique().shape[0]
        print("Successfully recalculated fibre assignment for \n"
              f"{num_redid_batches} invalid batches.")

        return pairs_df.dropna(subset=[ra_col, de_col])

    def _get_invalid_offsets(df,
                             ra_col,
                             de_col,
                             num_fields=2,
                             max_sep=4):
        """Identifies batches with invalid offsets based on the separation between fields.

        Args:
            df (pd.DataFrame): DataFrame containing the fields and their attributes.
            ra_col (pd.Series or str): Column name or Series containing the right ascension
                                       values.
            de_col (pd.Series or str): Column name or Series containing the declination values.
            num_fields (int): Number of fields expected in each batch. Defaults to 2.
            max_sep (float): Maximum allowed separation in degrees between fields in a batch.

        Returns:
            _type_: _description_
        """

        print('Checking for invalid offsets...')
        df_guides = df[df['guide'].isin([1])].copy()
        grouped = df_guides.groupby('batchId')
        invalid_batch_ids = []
        for batch_id, fields_df in grouped:
            if len(fields_df) == num_fields:
                field1 = fields_df.iloc[0]
                field2 = fields_df.iloc[1]

                ra1, dec1 = field1[ra_col], field1[de_col]
                ra2, dec2 = field2[ra_col], field2[de_col]

                coord1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg, frame='icrs')
                coord2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg, frame='icrs')

                separation_deg = coord1.separation(coord2).deg

                if separation_deg > max_sep:
                    invalid_batch_ids.append(batch_id)
        print(f'Found {len(invalid_batch_ids)} batches with invalid offsets:')
        print(invalid_batch_ids)
        return invalid_batch_ids

    pairs_df = df[df['fieldScore'] > 2].astype({src_id: 'str'})

    count = pairs_df.groupby('batchId')['fieldId'].nunique()
    # print('COUNT',  count)
    relevant_batch_ids = count[count.isin([1])].index
    # print('rel: ',  relevant_batch_ids)

    # Check if there's only one relevant batch ID, in which case skip pairing
    if len(relevant_batch_ids) == 1:
        updated_df = _update_batch_score(pairs_df)

        # Count the number of batches with 1 and 2 fields
        field_counts = updated_df.groupby('batchId')['fieldId'].nunique()
        single_field_batches = (field_counts == 1).sum()
        two_field_batches = (field_counts == 2).sum()
        three_field_batches = (field_counts == 3).sum()

        print(updated_df['batchId'].nunique(), 'batches have been created!')
        print(single_field_batches, 'batches have 1 field.')
        print(two_field_batches, 'batches have 2 fields.')
        print(three_field_batches, 'batches have 3 fields.')

        return updated_df.astype({'guide': 'int', 'batchScore': 'int', 'fieldScore': 'int'})

    df = df[df['fieldScore'] > 2].astype({src_id: 'str'})
    unassigned = df.groupby('batchId').filter(lambda g: g['fieldId'].nunique() == 1)
    central = unassigned[unassigned['guide'] == 1].drop_duplicates('fieldId')

    triplet_results = []
    used_fields = set()
    next_batch_id = df['batchId'].max() + 1

    # Try to form valid triplets
    for triplet in combinations(central.itertuples(index=False), 3):
        f1, f2, f3 = triplet
        if any(f.fieldId in used_fields for f in [f1, f2, f3]):
            continue
        if max_field_distance is not None and not _within_distance([f1._asdict(),
                                                                    f2._asdict(),
                                                                    f3._asdict()],
                                                                   max_field_distance):
            continue

        result = _assign_batch_and_guides([f1.fieldId, f2.fieldId, f3.fieldId], next_batch_id)

        if result is not None:
            triplet_results.append(result)
            used_fields.update([f1.fieldId, f2.fieldId, f3.fieldId])
            next_batch_id += 1

    # Try to form valid pairs with remaining fields
    remaining = central[~central['fieldId'].isin(used_fields)]
    pair_results, used_pair_fields, next_batch_id = _try_pairings(
        df[df['fieldId'].isin(remaining['fieldId'])],
        next_batch_id
    )

    # print('PAIRS', pair_results, used_pair_fields)
    # Try to upgrade pairs to triplets using leftover singles

    all_used = used_fields.union(used_pair_fields)
    final_remaining = central[~central['fieldId'].isin(all_used)]
    upgraded_pairs = []
    singles_to_remove = set()
    for pair_df in pair_results:
        batch_id = pair_df['batchId'].iloc[0]
        pair_fields = pair_df[pair_df['guide'] == 1]['fieldId'].unique()
        for _, single_row in final_remaining.iterrows():
            if single_row['fieldId'] in singles_to_remove:
                continue
            all_fields = list(pair_fields) + [single_row['fieldId']]
            if max_field_distance is not None:
                if not _within_distance([df[df['fieldId'] == fid].iloc[0].to_dict()
                                         for fid in all_fields], max_field_distance):
                    continue
            result = _assign_batch_and_guides(all_fields, batch_id)
            if result is not None:
                upgraded_pairs.append(result)
                singles_to_remove.add(single_row['fieldId'])
                break  # only one single can be added to each pair

    # Update results: remove upgraded pairs and used singles
    pair_results = [
        p for p in pair_results
        if p['batchId'].iloc[0] not in [
            up['batchId'].iloc[0] for up in upgraded_pairs
            ]
        ]

    final_remaining = final_remaining[~final_remaining['fieldId'].isin(singles_to_remove)]

    # Try to form 1+1+1 triplets from remaining singles
    new_triplets = []
    used_in_new_triplets = set()
    single_central = final_remaining[~final_remaining['fieldId'].isin(singles_to_remove)]
    for f1, f2, f3 in combinations(single_central.itertuples(index=False), 3):
        if any(f.fieldId in used_in_new_triplets for f in [f1, f2, f3]):
            continue
        if max_field_distance is not None and not _within_distance([f1._asdict(), f2._asdict(),
                                                                    f3._asdict()],
                                                                   max_field_distance):
            continue

        result = _assign_batch_and_guides([f1.fieldId, f2.fieldId, f3.fieldId], next_batch_id)

        if result is not None:
            new_triplets.append(result)
            used_in_new_triplets.update([f1.fieldId, f2.fieldId, f3.fieldId])
            next_batch_id += 1

    final_remaining = final_remaining[~final_remaining['fieldId'].isin(used_in_new_triplets)]

    # Try to form pairs from any leftover singles
    pair_leftovers, used_leftovers, next_batch_id = _try_pairings(
        df[df['fieldId'].isin(final_remaining['fieldId'])], next_batch_id)
    final_remaining = final_remaining[~final_remaining['fieldId'].isin(used_leftovers)]

    # Step 4: Assign any truly remaining singles to their own batches
    all_used = used_fields.union(used_pair_fields)
    final_remaining = central[~central['fieldId'].isin(all_used)]
    single_df = df[df['fieldId'].isin(final_remaining['fieldId'])].copy()
    unique_ids = single_df['fieldId'].unique()
    single_df = single_df.copy()
    single_df['batchId'] = single_df.groupby('fieldId').ngroup() + next_batch_id
    next_batch_id += len(unique_ids)

    df_combined = pd.concat([df[~df['fieldId'].isin(central['fieldId'])],
                             single_df] + triplet_results + pair_results, ignore_index=True)
    df_combined = df_combined.dropna(subset=[ra_col, de_col])

    # Final checks and fixes
    df_combined = _final_cleanup(df_combined)

    print(df_combined['batchId'].nunique(), 'batches created:')
    field_counts = df_combined.groupby('batchId')['fieldId'].nunique().value_counts()
    print(field_counts.get(1, 0), 'single-field batches')
    print(field_counts.get(2, 0), 'two-field batches')
    print(field_counts.get(3, 0), 'three-field batches')

    return df_combined


def calculate_batch_score(df):
    """Calculate the batch score based on the number of fields with successfully allocated
       guide stars/fibres.

    Args:
        df (pd.DataFrame): DataFrame containing the fields and their guide star assignments.

    Returns:
        pd.DataFrame: DataFrame with an additional 'batchScore' column indicating the score
    """

    df = df.sort_values(by='batchId')
    df['batchScore'] = np.nan
    workon_df = df.dropna(subset=['fibreId'])
    workon_df = df[df['sourceFlag'] == 9]

    count_vals = []
    for bid in workon_df['batchId'].unique():
        counts = workon_df[workon_df['batchId'] ==
                           bid]['fieldId'].value_counts().sort_values(ascending=False)
        count_str = ''.join(counts.astype(str).values)
        count_vals.append(int(count_str) if count_str.isdigit() else 0)

    for bid, v in zip(workon_df['batchId'].unique(), count_vals):
        df.loc[df['batchId'] == bid, 'batchScore'] = int(v)
    return df


def calculate_field_score(df):
    """Calculate the field score based on the number of fibres assigned to each field.

    Args:
        df (pd.DataFrame): DataFrame containing the fields and their fibre assignments.

    Returns:
        pd.DataFrame: DataFrame with an additional 'fieldScore' column indicating the score
    """

    df = df.sort_values(by='batchId')
    workon_df = df.dropna(subset=['fibreId'])
    workon_df = df[df['sourceFlag'] == 9]

    counts = workon_df['fieldId'].value_counts()
    indices, values = counts.index, counts.values

    for i, v in zip(indices, values):
        df.loc[df['fieldId'] == i, 'fieldScore'] = v

    return df


def sample_sky_fibre_coordinates(
    batch_df, guide_df,
    num_samples=20,
    radius_mm=5.0,
    plate_coord_range=(-190.0, 190.0),
    seed=None,
    # WCS / column settings to mirror zpn_projection
    p_coeff=[203.00898850, 0., 3.03143062, 0., -0.08049567, 0., 0.],
    ra_col='GAIA_RA',
    de_col='GAIA_DEC',
    field_col='fieldId',
    guide_col='guide',
    # columns we populate explicitly and should NOT auto-fill with ''
    skip_cols=('sourceFlag', 'TARGUSE',
               'TARGCLASS', 'X', 'Y',
               'batchScore', 'fieldScore',
               'fieldId', 'batchId', 'fibreId')):

    """
    For each field in df (grouped by field_col), sample 'num_samples' sky-fibre plate
    coordinates (X, Y) within plate_coord_range, ensuring each is >= radius_mm away from
    *all* other points in that same field (existing science/guides and newly sampled points).
    Convert sampled (X, Y) to (RA, Dec) via inverse ZPN using the field's central (guide==1)
    star as CRVAL. Append the new rows and return df.

    Notes:
      - Uses the same ZPN setup as zpn_projection and accounts for the Y flip:
        world->pixel used (x, y) then plate Y = -y, so inverse uses (x, -plateY).
      - Fills any columns present in df but missing in the new rows with '' (except skip_cols).
      - Sets TARGUSE='S', TARGCLASS='SKY', sourceFlag=-1 for the new sky fibres.
    """

    from astropy import wcs

    rng = np.random.default_rng(seed)

    dist_coeff = [(2, i + 1, coeff) for i, coeff in enumerate(p_coeff)]

    all_new_rows = []

    # Include guides so spacing accounts for them too
    batch_df = pd.concat([batch_df, guide_df], ignore_index=True)

    # Work per field for correct WCS & spacing
    for fid, fdf in batch_df.groupby(field_col, sort=False):
        # existing plate coords for spacing checks (clean NaN/inf)
        existing = fdf[['X', 'Y']].to_numpy(dtype=float, copy=True)
        if existing.size:
            existing = existing[np.isfinite(existing).all(axis=1)]
        else:
            existing = np.empty((0, 2), dtype=float)

        # Identify central guide star to define CRVAL
        central = fdf[fdf[guide_col] == 1]
        if central.empty:
            raise ValueError(f"No central guide star (guide==1) found for fieldId={fid}")
        crval_ra = float(central.iloc[0][ra_col])
        crval_de = float(central.iloc[0][de_col])

        # Build WCS for this field (RA---ZPN/DEC--ZPN + PV distortion)
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0, 0]
        w.wcs.crval = [crval_ra, crval_de]
        w.wcs.ctype = ['RA---ZPN', 'DEC--ZPN']
        w.wcs.set_pv(dist_coeff)

        samples = []
        need = int(num_samples)
        max_iters = 20

        for _ in range(max_iters):
            if need <= 0:
                break

            # Draw more candidates to keep iterations low
            n_try = max(need * 4, 16)
            cand_x = rng.uniform(plate_coord_range[0], plate_coord_range[1], size=n_try)
            cand_y = rng.uniform(plate_coord_range[0], plate_coord_range[1], size=n_try)
            candidates = np.column_stack((cand_x, cand_y))

            # First filter by distance to EXISTING points in this field
            if existing.shape[0] > 0:
                diff = candidates[:, None, :] - existing[None, :, :]
                min_dist_to_existing = np.sqrt(np.sum(diff**2, axis=2)).min(axis=1)
                candidates = candidates[min_dist_to_existing >= radius_mm]

            # Greedy accept to enforce pairwise spacing among newly added points
            accepted_now = []
            for cx, cy in candidates:
                if need <= 0:
                    break
                if accepted_now:
                    v = np.array([cx, cy])
                    diff2 = np.array(accepted_now) - v
                    # distance to the closest newly accepted point
                    if np.sqrt((diff2**2).sum(axis=1)).min() < radius_mm:
                        continue
                accepted_now.append([cx, cy])
                need -= 1

            if accepted_now:
                chosen = np.asarray(accepted_now, dtype=float)
                samples.extend(chosen.tolist())
                # Update existing so next iteration sees these too
                existing = np.vstack((existing, chosen))

        if need > 0:
            logging.warning(f"[fieldId={fid}] Could only find {num_samples - need} "
                            f"valid sky coordinates after {max_iters} iterations.\n"
                            "Please check if sky target distribution is sufficiently uniform!")

        samples = np.array(samples, dtype=float)  # (N,2) with columns X,Y (plate)
        # Invert the Y flip used in zpn_projection: pixel_y = -plate_y
        pixel_x = samples[:, 0]
        pixel_y = -samples[:, 1]

        # pixel -> world (RA/Dec in degrees)
        sky = w.pixel_to_world(pixel_x, pixel_y)
        ra_vals = np.asarray(sky.ra.to(u.deg).value, dtype=float)
        de_vals = np.asarray(sky.dec.to(u.deg).value, dtype=float)

        # Build frame for this field's new sky fibres
        new_rows = pd.DataFrame({
            'X': samples[:, 0],
            'Y': samples[:, 1],
            ra_col: ra_vals,
            de_col: de_vals,
            field_col: fid,
            'TARGUSE': 'S',
            'TARGCLASS': 'SKY',
            'sourceFlag': -1,
            'batchId': fdf['batchId'].iloc[0],
            'fieldId': fid,
            'fibreId': '',  # no fibreId for sky targets
        })

        all_new_rows.append(new_rows)

    # Combine all fields' new rows
    new_rows = pd.concat(all_new_rows, ignore_index=True)

    # Fill columns that exist in df but not in new_rows with ''
    skip_set = set(skip_cols)
    missing_cols = [col for col in batch_df.columns
                    if col not in new_rows.columns and col not in skip_set]
    for col in missing_cols:
        new_rows[col] = ''

    # Ensure same column order as df (and include any extra new columns at the end)
    shared_cols = [c for c in batch_df.columns if c in new_rows.columns]
    other_cols = [c for c in new_rows.columns if c not in batch_df.columns]
    new_rows = new_rows[shared_cols + other_cols]

    # Keep original science + any prior sky (-1), then append the new sky rows
    # df_keep = batch_df[batch_df['sourceFlag'].isin([0, -1])]
    df_keep = batch_df[batch_df['sourceFlag'] != 9]

    return pd.concat([df_keep, new_rows], ignore_index=True)


def decode_byte_columns(df: pd.DataFrame, to_string_dtype: bool = True) -> pd.DataFrame:
    obj_cols = [c for c in df.columns if df[c].dtype == 'object']
    for c in obj_cols:
        s = df[c]
        if s.map(lambda v: isinstance(v, (bytes, np.bytes_))).any():
            df[c] = s.map(
                lambda v: (
                    "" if isinstance(v, (bytes, np.bytes_)) and len(v) == 0
                    else (v.decode("utf-8", "ignore") if isinstance(v, (bytes, np.bytes_))
                          else v)
                )
            )
    if to_string_dtype:
        # only cast columns that are all strings after decode
        for c in obj_cols:
            if df[c].map(lambda v: (v is None) or isinstance(v, str)).all():
                df[c] = df[c].astype("string")
    return df


def assign_sky_attributes(df):

    assert (df['TARGUSE'] == 'S').any(), "No sky targets found in the DataFrame."

    df.loc[df['sourceFlag'] == -1, 'TARGPRIO'] = 1.0

    return df


def split_batches_at_max_sep(
    df,
    max_sep_deg=10.0,
    batch_col="batchId",
    field_col="fieldId",
    guide_col="guide",
    ra_col="GAIA_RA",
    dec_col="GAIA_DEC",
    targprog_col="TARGPROG",
    batchscore_col="batchScore",
    keep_original_batch_for_largest_group=True,
    return_mapping=False):

    """
    Split batches so that within each resulting batch, all guide-field pairs satisfy:
        separation <= max_sep_deg

    Supported: batches with 1–3 guide fields (typical WTL case).

    Updates:
      - batchId: new ids are float but integer-valued (e.g. 75.0), unique vs existing ids
      - TARGPROG: always "{batchIdInt}_{fieldIdInt}" (e.g. "75_11")
      - batchScore: depends on number of fields in batch after split:
            1 -> 3, 2 -> 33, 3 -> 333 (generalizes to n -> int("3"*n))

    Only these columns are modified; all rows remain in their original fieldId.
    """

    if max_sep_deg <= 0:
        raise ValueError("max_sep_deg must be > 0.")

    required = {batch_col, field_col, guide_col, ra_col, dec_col, targprog_col, batchscore_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"df is missing required columns: {sorted(missing)}")

    df_out = df.copy()

    # ---- treat 75 and 75.0 as the same used integer id ----
    def _to_int_or_none(x):
        if x is None:
            return None
        try:
            if isinstance(x, float) and np.isnan(x):
                return None
        except Exception:
            pass
        if isinstance(x, (bool, np.bool_)):
            return None
        try:
            xf = float(x)
        except Exception:
            return None
        if not np.isfinite(xf):
            return None
        if abs(xf - round(xf)) < 1e-9:
            return int(round(xf))
        return None

    used_int_ids = set()
    for v in df_out[batch_col].tolist():
        iv = _to_int_or_none(v)
        if iv is not None:
            used_int_ids.add(iv)

    next_int = (max(used_int_ids) + 1) if used_int_ids else 1

    def _new_batch_id_float() -> float:
        nonlocal next_int
        while next_int in used_int_ids:
            next_int += 1
        bid = next_int
        used_int_ids.add(bid)
        next_int += 1
        return float(bid)

    def _score_from_n(n: int) -> int:
        if n <= 0:
            return np.nan
        return int("3" * int(n))

    # ---- guide rows as centres; ensure 1 per (batch, field) ----
    df_guides = df_out[df_out[guide_col] == 1].copy()
    df_guides = df_guides.drop_duplicates(subset=[batch_col, field_col], keep="first")

    mapping_rows = []

    for batch_id, g in df_guides.groupby(batch_col, sort=False):
        n = len(g)
        if n <= 1:
            continue
        if n > 3:
            raise ValueError(
                f"Batch {batch_id!r} has {n} guide fields; this function currently supports up to 3."
            )

        fields = g[field_col].tolist()
        ra_vals = g[ra_col].to_numpy(dtype=float)
        dec_vals = g[dec_col].to_numpy(dtype=float)
        coords = SkyCoord(ra=ra_vals * u.deg, dec=dec_vals * u.deg, frame="icrs")

        # Pairwise separations
        # For n=2: sep01
        # For n=3: sep01, sep02, sep12
        def sep(i, j):
            return coords[i].separation(coords[j]).deg

        if n == 2:
            s01 = sep(0, 1)
            if s01 <= max_sep_deg:
                continue  # keep as is
            # split into two singles: keep original batch for one field, new batch for the other
            if keep_original_batch_for_largest_group:
                # both size 1; deterministic: keep original for smallest fieldId string
                keep_i = 0 if str(fields[0]) <= str(fields[1]) else 1
            else:
                keep_i = 0
            move_i = 1 - keep_i
            new_id = _new_batch_id_float()
            moved_field = fields[move_i]

            mask = (df_out[batch_col] == batch_id) & (df_out[field_col] == moved_field)
            df_out.loc[mask, batch_col] = new_id

            mapping_rows.append(
                dict(
                    old_batchId=float(batch_id) if _to_int_or_none(batch_id) is not None else batch_id,
                    new_batchId=new_id,
                    fieldId=moved_field,
                    component_index=1,
                    n_fields_in_component=1,
                )
            )
            continue

        # n == 3
        s01 = sep(0, 1)
        s02 = sep(0, 2)
        s12 = sep(1, 2)

        # If all pairs are within threshold, keep triplet
        if max(s01, s02, s12) <= max_sep_deg:
            continue

        # Determine which pairs are allowed
        pairs = [
            ((0, 1), s01),
            ((0, 2), s02),
            ((1, 2), s12),
        ]
        ok_pairs = [(ij, d) for ij, d in pairs if d <= max_sep_deg]

        if len(ok_pairs) == 0:
            # no pair qualifies => split into three singles
            # keep original batch for one single; move the other two to new batchIds
            order = sorted([0, 1, 2], key=lambda i: str(fields[i]))
            keep_i = order[0] if keep_original_batch_for_largest_group else 0
            for i in [0, 1, 2]:
                if i == keep_i:
                    continue
                new_id = _new_batch_id_float()
                moved_field = fields[i]
                mask = (df_out[batch_col] == batch_id) & (df_out[field_col] == moved_field)
                df_out.loc[mask, batch_col] = new_id
                mapping_rows.append(
                    dict(
                        old_batchId=float(batch_id) if _to_int_or_none(batch_id) is not None else batch_id,
                        new_batchId=new_id,
                        fieldId=moved_field,
                        component_index=i,
                        n_fields_in_component=1,
                    )
                )
            continue

        # At least one pair qualifies => create one double (best/closest ok pair) + one single
        (i, j), _dmin = min(ok_pairs, key=lambda t: t[1])
        k = ({0, 1, 2} - {i, j}).pop()

        # Decide which group keeps original batchId
        # groups: double={i,j} (size2), single={k} (size1)
        if keep_original_batch_for_largest_group:
            keep_double = True  # size2 wins
        else:
            # deterministic: keep group whose fieldId-list is lexicographically smallest
            double_key = [str(fields[x]) for x in sorted([i, j])]
            single_key = [str(fields[k])]
            keep_double = double_key <= single_key

        if keep_double:
            # move single to a new batch
            new_id = _new_batch_id_float()
            moved_field = fields[k]
            mask = (df_out[batch_col] == batch_id) & (df_out[field_col] == moved_field)
            df_out.loc[mask, batch_col] = new_id
            mapping_rows.append(
                dict(
                    old_batchId=float(batch_id) if _to_int_or_none(batch_id) is not None else batch_id,
                    new_batchId=new_id,
                    fieldId=moved_field,
                    component_index=1,
                    n_fields_in_component=1,
                )
            )
        else:
            # keep single, move the double to a new batch
            new_id = _new_batch_id_float()
            moved_fields = [fields[i], fields[j]]
            mask = (df_out[batch_col] == batch_id) & (df_out[field_col].isin(moved_fields))
            df_out.loc[mask, batch_col] = new_id
            for moved_field in moved_fields:
                mapping_rows.append(
                    dict(
                        old_batchId=float(batch_id) if _to_int_or_none(batch_id) is not None else batch_id,
                        new_batchId=new_id,
                        fieldId=moved_field,
                        component_index=0,
                        n_fields_in_component=2,
                    )
                )

    # ---- update TARGPROG for all rows ----
    batch_int = pd.to_numeric(df_out[batch_col], errors="coerce").round().astype("Int64")
    field_int = pd.to_numeric(df_out[field_col], errors="coerce").round().astype("Int64")

    if batch_int.isna().any() or field_int.isna().any():
        bad = df_out[batch_int.isna() | field_int.isna()][[batch_col, field_col, targprog_col]].head(10)
        raise ValueError(
            "Cannot update TARGPROG because some batchId/fieldId values are not integer-valued or are NaN. "
            f"Examples:\n{bad}"
        )

    df_out[targprog_col] = batch_int.astype(str) + "_" + field_int.astype(str)

    # ---- update batchScore based on number of guide fields per batch ----
    guides_after = df_out[df_out[guide_col] == 1].drop_duplicates(subset=[batch_col, field_col], keep="first")
    n_fields_per_batch = guides_after.groupby(batch_col)[field_col].nunique()
    # batch_score_map = n_fields_per_batch.map(_score_from_n).to_dict()
    # df_out[batchscore_col] = df_out[batch_col].map(batch_score_map)

    if return_mapping:
        mapping_df = pd.DataFrame(mapping_rows)
        if not mapping_df.empty:
            mapping_df["old_batchId"] = pd.to_numeric(mapping_df["old_batchId"], errors="coerce").astype(float)
            mapping_df["new_batchId"] = pd.to_numeric(mapping_df["new_batchId"], errors="coerce").astype(float)
        return df_out, mapping_df

    return df_out


def run_guide_fibre_allocation(df_field_selection,
                               wg_comb_df,
                               lookup_df,
                               hemisphere,
                               use_degrees=True,
                               max_gate_angle=12.5,
                               plate='A',
                               split_batch_at_max_sep=True,
                               maximize_batches=False,
                               max_field_distance=10.0,
                               guide_mag_limit=None,
                               incl_guide_layer=False,
                               perfect_fields_only=False,
                               verbose=False,
                               mag_col='GAIA_MAG_G',
                               ra_col='GAIA_RA',
                               de_col='GAIA_DEC',
                               src_id_col='CNAME',
                               guide_reach_factor=18.0,
                               broken_fibres=None,
                               add_sky_targets=False,
                               max_num_sky_targets=30,
                               outpath='./'):

    """
    Constructs and optionally saves an XML configuration based on various DataFrames and field
    parameters.

    This function creates an XML structure for observation configuration using data from
    multiple DataFrames and specific field selection criteria. It also handles saving the XML
    if a filename is provided.

    Args:
        df_field_selection (DataFrame): DataFrame with field selection data.
        df_circle_coordinates (DataFrame): DataFrame with circle coordinates.
        field_id (int): The ID of the field to configure.
        observation_df (DataFrame): DataFrame with observation details.
        configure_df (DataFrame): DataFrame with configuration details.
        surveys_df (DataFrame): DataFrame with survey details.
        use_degrees (bool):
        mixed_mode (bool):
        mIFU_loc (str):
        filename (str, optional): File name to save the XML. If None, the XML is not saved.
        Defaults to None.

    Returns:
        None
    """

    def _match_field_ids_w_lookup_table(lut_val, field_ids):
        parts = lut_val.split('_')
        return any(part in field_ids for part in parts)

    def _explode_batch_by_lookup(batch, lookup_df, batch_field_ids,
                                 src_id_col,
                                 field_id_col="fieldId",
                                 lut_field_col="lutFieldId"):

        """
        Explode rows in `batch` where guide is NaN into one row per field membership
        (from lookup_df.lutFieldId), but only for memberships present in `batch_field_ids`.
        """

        # Prepare types
        batch_field_ids_set = set(map(str, batch_field_ids))

        # Split batch into (A) rows to explode (guide is NaN) and (B) the rest
        mask_to_explode = batch["guide"].isna()
        batch_rest = batch.loc[~mask_to_explode].copy()
        batch_to_explode = batch.loc[mask_to_explode].copy()

        if batch_to_explode.empty:
            return batch.copy()  # nothing to do

        # Map src_id -> lutFieldId and split into rows
        # Merge only the keys to avoid duplicating full rows early
        keys = batch_to_explode[[src_id_col]].drop_duplicates()
        lut = keys.merge(
            lookup_df[[src_id_col, lut_field_col]],
            on=src_id_col,
            how="left",
        )

        # Drop rows without lut info and split "1_2_3" -> ["1","2","3"], then explode
        lut = lut.dropna(subset=[lut_field_col]).copy()
        if lut.empty:
            # No lookups found; return original batch unchanged
            return batch.copy()

        lut["__fields_list__"] = lut[lut_field_col].astype(str).str.split("_")
        lut = lut.explode("__fields_list__", ignore_index=True)

        # Keep only memberships present in current iteration's batch_field_ids
        lut = lut[lut["__fields_list__"].isin(batch_field_ids_set)].copy()
        if lut.empty:
            # No memberships intersect with current batch fields
            return pd.concat([batch_rest, batch_to_explode], ignore_index=True)

        # Prepare exploded mapping: (src_id, fieldId_mapped)
        lut.rename(columns={"__fields_list__": "__mapped_field__"}, inplace=True)

        # Merge mapping back to rows-to-explode -> duplicates rows per membership
        exploded = batch_to_explode.merge(lut[[src_id_col, "__mapped_field__"]],
                                          on=src_id_col,
                                          how="inner")

        # Set fieldId to the mapped value; cast to match original dtype if numeric
        # Try to cast to integer where possible, but be safe with np.integer
        mapped = pd.to_numeric(exploded["__mapped_field__"], errors="ignore")
        if np.issubdtype(
            getattr(
                batch[field_id_col].dtype, "type",
                type(batch[field_id_col].dtype)), np.integer):

            exploded[field_id_col] = pd.to_numeric(mapped, errors="coerce").astype("Int64")
        else:
            exploded[field_id_col] = mapped.astype(str)

        exploded = exploded.drop(columns=["__mapped_field__", lut_field_col], errors="ignore")

        # Concatenate back: rest + exploded
        out = pd.concat([batch_rest, exploded], ignore_index=True).astype({'fieldId': float})
        out.loc[out['guide'].isna(), 'guide'] = 0

        return out

    print(f'Allocation in progress for plate {plate} - {hemisphere}')
    print('------------------------------------------------------'
          '------------------')

    guide_dfs = []
    science_dfs = []
    batch_ids = np.sort(df_field_selection['batchId'].unique())

    # Loop through field batches and assign guide stars and fibres
    for batch_id in batch_ids:
        batch = df_field_selection[df_field_selection['batchId'] == batch_id]
        batch['fibreId'] = pd.NA #np.nan

        batch_field_ids = set(batch['fieldId'].astype(str).unique())

        matched_rows = lookup_df[lookup_df['lutFieldId'].apply(
            _match_field_ids_w_lookup_table, field_ids=batch_field_ids)]

        central_guide_ids = batch[batch['guide'] == 1][src_id_col].unique()
        src_ids_from_lut = matched_rows[src_id_col].unique()
        src_ids_from_lut = np.setdiff1d(src_ids_from_lut, central_guide_ids)

        if src_ids_from_lut.shape[0] > 6000:
            rng = np.random.default_rng(seed=42)  # fixed seed for reproducibility
            src_ids_from_lut = rng.choice(src_ids_from_lut, 6000, replace=False)

        # Now get all the guide stars for each field and add them to the batch
        guide_star_rows = wg_comb_df[wg_comb_df[src_id_col].isin(src_ids_from_lut)]
        if guide_star_rows.empty:
            logging.warning(f'No guide stars found for Batch-Id: {batch_id} '
                            'Skipping this batch.')
            continue

        # Add the queried guide stars to the current batch.
        batch = pd.concat([batch, guide_star_rows], ignore_index=True).astype(
                                                                           {src_id_col: 'str'})

        # Set batchId for the added guide stars that have no batchId set
        batch.loc[batch['batchId'].isna(), 'batchId'] = batch['batchId'].unique()[0]

        # Explode the batch to account for guide stars that are part of multiple fields.
        batch = _explode_batch_by_lookup(batch=batch, lookup_df=lookup_df,
                                         batch_field_ids=batch_field_ids,
                                         src_id_col=src_id_col,
                                         field_id_col="fieldId",
                                         lut_field_col="lutFieldId")

        # Merge fields and apply ZPN projection (onto field plate).
        merged_fields = zpn_projection(batch=batch, ra_col=ra_col, de_col=de_col)
        science_df = merged_fields[merged_fields['sourceFlag'] != 9].copy()

        # Allocate guide fibres to accomodate the superposition of fields.
        guide_star_df = assign_guide_fibres(merged_fields,
                                            max_gate_angle,
                                            hemisphere=hemisphere,
                                            guide_mag_limit=guide_mag_limit,
                                            plate=plate,
                                            ra_col=ra_col,
                                            de_col=de_col,
                                            mag_col=mag_col,
                                            src_id=src_id_col,
                                            guide_reach_factor=guide_reach_factor,
                                            broken_fibres=broken_fibres)

        # if True gives only fields with 2 secondary and 1 central guide star
        if perfect_fields_only and guide_star_df is not None:
            guide_counts = guide_star_df[guide_star_df['guide'] == 2].groupby('fieldId').size()
            valid_fields = guide_counts[guide_counts >= 2].index
            guide_star_df = guide_star_df[(guide_star_df['fieldId'].isin(valid_fields))]

        add_sky = False

        if add_sky_targets:
            n_fields = batch["fieldId"].nunique()
            n_sci = batch.loc[~batch["sourceFlag"].isin([-1, 9])].shape[0]

            # Handle single field batches
            if n_fields < 2:
                if 100 <= n_sci < 300:
                    num_samples = int(n_sci * 0.25)
                if 50 <= n_sci < 100:
                    num_samples = int(n_sci * 0.35)                    
                elif 30 <= n_sci < 50:
                    num_samples = max_num_sky_targets
                elif 300 <= n_sci <= 900:
                    num_samples = max_num_sky_targets
                elif n_sci > 900:
                    logging.warning(
                        "More than 900 science targets! "
                        "Adjust target count to enable sufficient sky fibre coverage!"
                    )
                    num_samples = 10
                else:
                    num_samples = 15
            # Add up to max_num_sky_targets if more than 2 fields
            else:
                num_samples = max_num_sky_targets
            """
            # For batches with only one field:
            if (batch['fieldId'].unique().shape[0] < 2
                    and 50 <= batch[~batch['sourceFlag'].isin([-1, 9])].shape[0] < 300):
                # Add approx. 25% sky targets if less than 300 science targets
                num_samples = int(batch[~batch['sourceFlag'].isin([-1, 9])].shape[0] * 0.25)
            elif (batch['fieldId'].unique().shape[0] < 2
                    and 30 <= batch[~batch['sourceFlag'].isin([-1, 9])].shape[0] < 50):
                # Add approx. 25% sky targets if less than 300 science targets
                num_samples = max_num_sky_targets
            elif ((batch['fieldId'].unique().shape[0] < 2)
                    and (batch[~batch['sourceFlag'].isin([-1, 9])].shape[0]) => 300):
                # Add up to 80 sky targets if more than 300 science targets
                num_samples = 80

            elif ((batch['fieldId'].unique().shape[0] < 2)
                    and (batch[~batch['sourceFlag'].isin([-1, 9])].shape[0]) < 30):
                # Add at least 10 sky targets if less than 30 science targets
                num_samples = 10

            # For batches with 2 or more fields:
            elif batch['fieldId'].unique().shape[0] >= 2:
                # Add up to max_num_sky_targets if more than 2 fields
                num_samples = max_num_sky_targets
            """
            # Sample sky fibre coordinates and add to science_df
            science_df = sample_sky_fibre_coordinates(
                science_df, guide_star_df, num_samples=num_samples, radius_mm=10.0,
                plate_coord_range=(-190.0, 190.0), seed=42)

            # Set flag to indicate that sky targets were added
            add_sky = True

        if science_df is not None:
            science_dfs.append(science_df)

        # Create a DataFrame of guide stars for each batch
        if guide_star_df is not None:
            guide_dfs.append(guide_star_df)

            # Verbose output
            if verbose:
                print('Batch ID: ', batch_id)
                print('Guide fibre assignment: ')
                print('------------------------------------------------------'
                      '------------------')
                print(guide_star_df[[src_id_col, 'fieldId', 'guide',
                                     'fibreId', 'sourceFlag', mag_col]
                                    ].sort_values(by=['fieldId', 'guide']))
                print('------------------------------------------------------'
                      '------------------')
                print('Nr. of fields in batch', batch['fieldId'].unique().shape[0])
                print('Nr. of guide stars in batch', guide_star_df.shape[0])
                print('Nr. of non-guide targets in batch',
                      batch[batch['sourceFlag'] != 9].shape[0])
                counts = (
                    batch[~batch['sourceFlag'].isin([-1, 9])]     # exclude guides
                    .groupby('fieldId')
                    .size()
                    )

                if add_sky:
                    n_sky = (science_df['sourceFlag'] == -1).sum()
                    n_fields = science_df['fieldId'].nunique()
                    if n_fields == 0:
                        print("No fields found; cannot compute sky targets per field!")
                    else:
                        print(f"Added {int(n_sky / n_fields)} sky targets per field!")

                print('------------------------------------------------------'
                      '------------------')
                if add_sky:
                    print("Nr. of non-guide targets per field (excluding added sky targets):")
                elif not add_sky:
                    print("Nr. of non-guide targets per field:")
                print(counts)
                print('------------------------------------------------------'
                      '------------------')
                print('------------------------------------------------------'
                      '------------------')

        elif guide_star_df is None:
            print('No guides found at Batch-Id:', batch_id)
            continue

        if perfect_fields_only:
            batch = batch[(batch['fieldId'].isin(guide_star_df['fieldId'].unique())) |
                          (~batch['fibreId'].isna())]

    combined_guide_frame = pd.concat(guide_dfs).astype({src_id_col: 'str',
                                                        'fibreId': 'Int64',
                                                        'fieldId': 'Int64',
                                                        'batchId': 'Int64'})

    combined_science_frame = pd.concat(science_dfs).astype({src_id_col: 'str'})

    # df_field_selection.loc[:, 'fibreId'] = np.nan
    # df_field_selection.loc[:, 'gateDist'] = np.nan
    # df_field_selection.loc[:, 'gateAngle'] = np.nan
    # df_field_selection.loc[:, 'X'] = np.nan
    # df_field_selection.loc[:, 'Y'] = np.nan
    # print(merged_fields[merged_fields['sourceFlag'] != 9][['X', 'Y']])
    # print(merged_fields[merged_fields['sourceFlag']!=9].shape)

    field_selection_with_guides = pd.concat([combined_science_frame,
                                             combined_guide_frame],
                                            ignore_index=True)

    # print(field_selection_with_guides.shape)
    # If True output does not contain stars from the guide cataloge except alreadey
    # assigned guidestars
    if not incl_guide_layer:
        field_selection_with_guides = field_selection_with_guides[
                                        ~((field_selection_with_guides['sourceFlag'] == 9) &
                                          (field_selection_with_guides['guide'].isin([0, -1]))
                                          )]

    # Calculate field and batch scores
    field_selection_with_guides = calculate_field_score(field_selection_with_guides)
    field_selection_with_guides = calculate_batch_score(field_selection_with_guides)

    # Drop fields that have no guide stars
    field_selection_with_guides = field_selection_with_guides[
                                  ~field_selection_with_guides['fieldScore'].isna()]

    if perfect_fields_only:
        field_selection_with_guides = field_selection_with_guides[
                                      field_selection_with_guides['fieldScore'] >= 2]   


    if maximize_batches:
        print('Maximizing fields per batch...')
        field_selection_with_guides = maximize_final_batches(
            field_selection_with_guides,
            df_field_selection,
            ra_col=ra_col,
            de_col=de_col,
            src_id=src_id_col,
            mag_col=mag_col,
            plate=plate,
            guide_mag_limit=guide_mag_limit,
            hemisphere=hemisphere,
            max_gate_angle=max_gate_angle,
            guide_reach_factor=guide_reach_factor,
            max_field_distance=max_field_distance,
            broken_fibres=broken_fibres,
        )

        field_selection_with_guides = calculate_field_score(field_selection_with_guides)
        field_selection_with_guides = calculate_batch_score(field_selection_with_guides)

        
    if split_batch_at_max_sep:
        print(f'Splitting batches with fields with a separation higher than {max_field_distance} deg...')
        print(type(field_selection_with_guides))
        field_selection_with_guides = split_batches_at_max_sep(field_selection_with_guides,
                                                               max_sep_deg=max_field_distance,
                                                               return_mapping=False)
        
        field_selection_with_guides = calculate_field_score(field_selection_with_guides)
        field_selection_with_guides = calculate_batch_score(field_selection_with_guides)            
    

    # Check for insufficient guide fibre allocations
    for batch_id, group in field_selection_with_guides.groupby("batchId"):
        bad_fields = group.loc[group["fieldScore"] < 3, "fieldId"].unique()
        if len(bad_fields):
            logging.warning(
                f"Field(s) {', '.join(map(str, bad_fields))} in batch {batch_id} "
                f"have insufficient guide fibre coverage. \n"
                "The affected field cannot be observed! \n"
                f"Please rerun the field definition routine!!! \n"
                "TIP: Slightly adjusting the grid_offset value during field definition might "
                "do the trick."
            )

    print('Field selection complete!')
    print('Number of batches:', len(field_selection_with_guides['batchId'].unique()))
    print('Number of fields:', len(field_selection_with_guides['fieldId'].unique()))
    print('------------------------------------------------------'
          '------------------')

    field_selection_with_guides = field_selection_with_guides.astype({src_id_col: 'str'})

    # print(field_selection_with_guides[
    #             field_selection_with_guides['sourceFlag'] == -1][src_id])

    # Calculations for final print outs:
    field_num = df_field_selection[df_field_selection['sourceFlag'] != 9].shape[0]
    cov_num = field_selection_with_guides[
        ~field_selection_with_guides['sourceFlag'].isin([9, -1])].shape[0]

    batch_scores = np.sort(field_selection_with_guides['batchScore']
                           .dropna().astype(int).unique())

    batch_counts = (field_selection_with_guides
                    .groupby('batchScore')['batchId'].nunique().sort_index())

    print('Number of batches per batch score:')
    for score, count in batch_counts.items():
        print(f'batchScore {int(score)} : {count} batches')

    print('------------------------------------------------------'
          '------------------')
    print('Nr. science targets before allocation:', field_num)
    print('Nr. science targets after allocation:', cov_num)
    print(f'Science target coverage factor: {cov_num/field_num} '
          '(nr. initial objects/nr. objects in\n'
          'allocated fields)')
    print('------------------------------------------------------'
          '------------------')
    print('Unique batch scores: ', batch_scores)
    print(('Batch score gives the nr. of successfully allocated guide fibres/field:\n'
           '333 = 3 fields with 3 good guide allocations --> best score,\n'
           '111 = 3 fields with one allocation each --> worst score\n'
           'FIELDS WITH A SCORE < 3 CANNOT BE OBSERVED!'
           ))

    print('------------------------------------------------------'
          '------------------')
    print('------------------------------------------------------'
          '------------------')

    # Replace byte strings (b'') in the dataframe with regular strings
    # field_selection_with_guides = decode_byte_columns(field_selection_with_guides)

    # Save the DataFrame to a CSV file
    field_selection_with_guides.to_csv(
        f'{outpath}tl_coverage_plate_{plate}_{hemisphere}.csv', index=False)

    # Assign TARGPROG values to enable correct xml creation
    field_selection_with_guides = assign_targprog(field_selection_with_guides)
    field_selection_with_guides = assign_sky_attributes(field_selection_with_guides)

    return field_selection_with_guides


###############################################################################################
###############################################################################################
###############################################################################################


if __name__ == "__main__":

    # Example DataFrame for <observation> element
    observation_df = pd.DataFrame([{
        'num_sky_fibres': '50',
        'max_calibration': '20',
        'maximum_gate_angle': '10'
    }])

    # Example DataFrame for <configure> element
    configure_df = pd.DataFrame([{
            'num_sky_fibres': '100',
            'max_calibration': '30',
            'maximum_gate_angle': '15',
            'plate': 'PLATE_B'
        }])

    # Example DataFrame for <surveys> element
    surveys_df = pd.DataFrame([
            {'name': 'survey_1', 'max_fibres': '700'}
        ])

    month = 'Aug'
    field = 'lopn'
    hemisphere = 'N'
    plate = 'PLATE_A'

    weave_df = pd.DataFrame([{'datamver': "8.00", 'author': 'th721@ast.cam.ac.uk'}])

    observation_df = pd.DataFrame([{'name': 'WTLS_test',
                                    'program': 'WTLS',
                                    'casuid': '-1',
                                    'chained': 'False',
                                    'coordinate_system': 'ICRS',
                                    'linkedgroup': '%%%',
                                    'ob_class': 'science',
                                    'ob_priority': '1.0',
                                    'obsgroup_validity': '',
                                    'obstemp': '%%%%%',
                                    'progtemp': '9%%%%',
                                    'tac_alloc': '%%%',
                                    'tac_id': '%%%',
                                    'trimester': '%%%'}])

    # Example DataFrame for <configure> element
    configure_df = pd.DataFrame([{
        'num_sky_fibres': '10',
        'max_calibration': '30',
        'maximum_gate_angle': '13.8',
        'plate': f'{plate}'
    }])

    # Example DataFrame for <surveys> element
    surveys_df = pd.DataFrame([
        {'name': 'WTLS_test_run', 'max_fibres': '700'}
    ])

    exposure_df = pd.DataFrame([{'type': 'science',
                                 'exp_time': '30',
                                 'arm': 'both',
                                 'cal_lamp': 'None',
                                 'cal_lamp_filter_A': 'WCALU_clear1',
                                 'cal_lamp_filter_B': 'EMPTY',
                                 'ff_ilu_1': '0',
                                 'ff_ilu_2': '0',
                                 'ff_ilu_3': '0',
                                 'speed': 'slow'}])

    calib_df = pd.DataFrame({
                'arm': ['red', 'blue', 'red', 'blue', 'red', 'blue'],
                'cal_lamp': ['QTH', 'QTH', 'ThArCr', 'ThArCr', 'ThArCr', 'ThArCr'],
                'cal_lamp_filter_A': ['WCALU_clear1', 'WCALU_clear1', 'WCALU_clear1',
                                      'WCALU_clear1', 'WCALU_clear1', 'WCALU_clear1'],
                'cal_lamp_filter_B': ['EMPTY', 'EMPTY', 'EMPTY', 'EMPTY', 'EMPTY', 'EMPTY'],
                'exp_time': [3, 9, 10, 120, 10, 120],
                'ff_ilu_1': [0, 0, 0, 0, 0, 0],
                'ff_ilu_2': [0, 0, 0, 0, 0, 0],
                'ff_ilu_3': [0, 0, 0, 0, 0, 0],
                'order': [1, 1, 2, 2, 4, 4],
                'speed': ['fast', 'fast', 'fast', 'fast', 'fast', 'fast'],
                'type': ['fibre_lampflat', 'fibre_lampflat', 'arc', 'arc', 'arc', 'arc'],
    })
