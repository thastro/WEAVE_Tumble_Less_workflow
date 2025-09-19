import logging
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from scipy.spatial import cKDTree
import healpy as hp
from astropy.coordinates import SkyCoord
import astropy.units as u


def add_healpix_id(df, ra_col='ra', de_col='dec', hp_col='healpixId', nside=32, nested=True):
    """
    Adds a new column to the dataframe containing the HEALPix ID based on RA and Dec
    coordinates.

    Parameters:
        df: pandas DataFrame containing the RA and Dec coordinates
        ra_col: name of the column containing the RA coordinates (default 'RA')
        dec_col: name of the column containing the Dec coordinates (default 'Dec')
        nside: NSIDE parameter for HEALPix (default 32)
        nested: whether to use nested ordering for HEALPix (default True)

    Returns:
     df: pandas DataFrame with a new column 'healpix_id' containing the HEALPix IDs
    """
    df = df.copy()
    # Convert RA and Dec to theta and phi
    theta = np.radians(90.0 - df[de_col].values)  # colatitude
    phi = np.radians(df[ra_col].values)            # longitude

    # Compute HEALPix IDs
    healpix_id = hp.ang2pix(nside, theta, phi, nest=nested)

    # Add HEALPix IDs to the dataframe
    df.loc[:, hp_col] = np.int64(healpix_id)

    return df


def remove_single_star_fields(df, df_other):
    # Identify 'fieldId's with only one star
    # fields_with_plato_targets = df[df['sourceFlag'] != 9]
    counts = df.groupby('fieldId').size()
    single_star_fields = counts[counts == 1].index.tolist()

    # Remove entries with these 'fieldId's from the original DataFrame
    filtered_df = df[df['fieldId'].isin(single_star_fields) == False]

    # Remove entries with these 'fieldId's from the other DataFrame
    filtered_df_other = df_other[df_other['fieldId_c'].isin(single_star_fields) == False]

    return filtered_df, filtered_df_other


def place_optimized_circles(catalogue_df,
                            bright_stars_df=None,
                            central_star_mag_lim=None,
                            offset=1.0,
                            n_clusters=None,
                            ra_col='RAdeg',
                            de_col='DEdeg',
                            mag_col='gaiaV'):

    """
    Place circles either around bright stars or find optimal circle placement based on a
    catalogue.

    If a DataFrame of bright stars is provided, places circles around these stars. If not
    provided, it invokes `find_optimal_circle_centers` to determine the optimal placement of
    circles.

    Args:
        catalogue_df (DataFrame): A DataFrame containing the catalogue of stars.
        bright_stars_df (DataFrame, optional): A DataFrame containing bright stars. Defaults
        to None.

    Returns:
        list: A list of tuples representing the circle centers (RA deg, DEC deg).
    """

    if bright_stars_df is not None:
        return [(row[ra_col], row[de_col]) for index, row in bright_stars_df.iterrows()]
    else:
        if n_clusters is None:
            print('heuristic field selection...')
            return find_optimal_circle_centers_spherical(catalogue_df, offset=offset,
                                               central_star_mag_lim=central_star_mag_lim,
                                               ra_col=ra_col, de_col=de_col, mag_col=mag_col)
        print('k-means field selection...')
        return find_optimal_circle_centers_kmeans(catalogue_df, offset=offset,
                                                  central_star_mag_lim=central_star_mag_lim,
                                                  n_clusters=n_clusters,
                                                  ra_col=ra_col, de_col=de_col,
                                                  mag_col=mag_col)


def find_optimal_circle_centers_kmeans(catalogue_df,
                                       central_star_mag_lim=None,
                                       n_clusters=10,
                                       offset=1.0, mag_col='gaiaV',
                                       ra_col='RAdeg', de_col='DEdeg'):
    """
    Find an efficient placement of circle centers to cover stars in a catalogue,
    balancing between the brightness of stars and spatial coverage.

    Args:
        catalogue_df (DataFrame): DataFrame containing the star catalogue.
        central_star_mag_lim (float, optional): Magnitude limit for the central star.
                                                Defaults to None.
        n_clusters (int): Number of clusters to form for initial grouping of stars.
                          Defaults to 10.
        offset (float): Minimum distance between circle centers, acting as the circle's
                        radius. Defaults to 1.0.

    Returns:
        list: A list of tuples representing the efficient circle centers (RAdeg, DEdeg).
    """

    # Filter stars brighter than the magnitude limit if specified
    if central_star_mag_lim is not None:
        # Define which stars can be used as central guide stars.
        catalogue_df = catalogue_df[(catalogue_df[mag_col] <= central_star_mag_lim) &
                                    (catalogue_df['sourceFlag'] != 5) &
                                    (catalogue_df['sourceFlag'] != 7) &
                                    (catalogue_df['sourceFlag'] != 4)
                                    (catalogue_df['sourceFlag'] == 9)]

    # Cluster stars based on their spatial positions
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    positions = catalogue_df[[ra_col, de_col]]
    kmeans.fit(positions)
    catalogue_df['cluster'] = kmeans.labels_

    circles = []

    # Process each cluster
    for cluster in range(n_clusters):
        cluster_df = catalogue_df[catalogue_df['cluster'] == cluster]
        # Sort stars in the cluster by brightness
        cluster_df.sort_values(by=mag_col, ascending=True, inplace=True)

        while not cluster_df.empty:
            brightest_star = cluster_df.iloc[0]
            circle_center = (brightest_star[ra_col], brightest_star[de_col])

            # Check for overlap with existing circles
            if all(np.sqrt((circle_center[0] - x[0])**2
                           + (circle_center[1] - x[1])**2) >= offset for x in circles):

                circles.append(circle_center)
                # Remove stars within the circle's area, including other clusters
                catalogue_df = catalogue_df[np.sqrt((catalogue_df[ra_col] -
                                            circle_center[0])**2 +
                                            (catalogue_df[de_col] -
                                            circle_center[1])**2) > offset]
                # Update the cluster_df for the next iteration
                cluster_df = catalogue_df[catalogue_df['cluster'] == cluster]

    return circles


def find_optimal_circle_centers_o(catalogue_df,
                                central_star_mag_lim=None,
                                offset=1.0, mag_col='gaiaV',
                                ra_col='RAdeg', de_col='DEdeg'):
    """
    Find the optimal placement of circle centers to cover the most stars in a catalogue,
    with each circle centered around the brightest star in the remaining subset.

    This function iterates over the star catalogue, placing circles centered on the brightest
    stars while ensuring no overlaps between circles.

    Args:
        catalogue_df (DataFrame): A DataFrame containing the catalogue of stars.
        central_star_mag_lim (float, optional): The magnitude limit for the central star.
                                                Stars brighter than this limit will be used
                                                as central guide stars. Defaults to None.
    Returns:
        list: A list of tuples representing the optimal circle centers (RAdeg, DEdeg).
    """

    circles = []
    remaining_stars = catalogue_df[(catalogue_df['sourceFlag'] != 5) &
                                   (catalogue_df['sourceFlag'] != 7) &
                                   (catalogue_df['sourceFlag'] == 9)]

    # Sort by brightness brightest first
    remaining_stars.sort_values(by=mag_col, ascending=True, inplace=True)

    if central_star_mag_lim is not None:

        while remaining_stars[mag_col].min() <= central_star_mag_lim:
            # Always choose the brightest star in the remaining subset
            brightest_star = remaining_stars.iloc[0]
            circle_center = (brightest_star[ra_col], brightest_star[de_col])

            # Check if this circle overlaps with any existing circle
            if all(np.sqrt((circle_center[0] - x[0])**2 + (circle_center[1] - x[1])**2) >=
                   offset for x in circles):
                circles.append(circle_center)

                # Remove stars covered by the new circle
                remaining_stars = remaining_stars[np.sqrt((remaining_stars[ra_col] -
                                                           circle_center[0])**2 +
                                                          (remaining_stars[de_col] -
                                                           circle_center[1])**2) > offset]
    else:

        while not remaining_stars.empty:
            # Always choose the brightest star in the remaining subset
            brightest_star = remaining_stars.iloc[0]
            circle_center = (brightest_star[ra_col], brightest_star[de_col])

            # Check if this circle overlaps with any existing circle
            if all(np.sqrt((circle_center[0] - x[0])**2 + (circle_center[1] - x[1])**2) >=
                   offset for x in circles):
                circles.append(circle_center)

                # Remove stars covered by the new circle
                remaining_stars = remaining_stars[np.sqrt((remaining_stars[ra_col] -
                                                           circle_center[0])**2 +
                                                          (remaining_stars[de_col] -
                                                           circle_center[1])**2) > offset]

    return circles


def find_optimal_circle_centers(catalogue_df, central_star_mag_lim=None, offset=1.0,
                                mag_col='gaiaV', ra_col='RAdeg', de_col='DEdeg'):
    """
    Find the optimal placement of circle centers to cover the most stars in a catalogue,
    with each circle centered around the brightest star in the remaining subset.

    Args:
        catalogue_df (DataFrame): A DataFrame containing the catalogue of stars.
        central_star_mag_lim (float, optional): The magnitude limit for the central star.
                                                Stars brighter than this limit will be used
                                                as central guide stars. Defaults to None.
    Returns:
        list: A list of tuples representing the optimal circle centers (RAdeg, DEdeg).
    """

    circles = []
    remaining_stars = catalogue_df[(catalogue_df['sourceFlag'] != 5) &
                                   (catalogue_df['sourceFlag'] != 7) &
                                   (catalogue_df['sourceFlag'] != 4) &
                                   (catalogue_df['sourceFlag'] == 9)]

    # Sort by brightness brightest first
    remaining_stars.sort_values(by=mag_col, ascending=True, inplace=True)

    # Create a KD-tree for efficient spatial queries
    star_positions = remaining_stars[[ra_col, de_col]].values
    tree = cKDTree(star_positions)

    if central_star_mag_lim is not None:
        while not remaining_stars.empty and remaining_stars[
                                            mag_col].iloc[0] <= central_star_mag_lim:
            # Always choose the brightest star in the remaining subset
            brightest_star = remaining_stars.iloc[0]
            circle_center = (brightest_star[ra_col], brightest_star[de_col])

            # Check if this circle overlaps with any existing circle
            if all(np.linalg.norm(np.array(circle_center) - np.array(x)) >= offset
                   for x in circles):

                circles.append(circle_center)

                # Remove stars covered by the new circle using KD-tree
                indices_to_remove = tree.query_ball_point(circle_center, offset)
                remaining_stars = remaining_stars.drop(
                                  remaining_stars.index[indices_to_remove]
                                  )
                star_positions = remaining_stars[[ra_col, de_col]].values
                tree = cKDTree(star_positions)
    else:
        while not remaining_stars.empty:
            # Always choose the brightest star in the remaining subset
            brightest_star = remaining_stars.iloc[0]
            circle_center = (brightest_star[ra_col], brightest_star[de_col])

            # Check if this circle overlaps with any existing circle
            if all(np.linalg.norm(np.array(circle_center) - np.array(x)) >= offset
                   for x in circles):
                circles.append(circle_center)

                # Remove stars covered by the new circle using KD-tree
                indices_to_remove = tree.query_ball_point(circle_center, offset)
                remaining_stars = remaining_stars.drop(remaining_stars.index[indices_to_remove])
                star_positions = remaining_stars[[ra_col, de_col]].values
                tree = cKDTree(star_positions)

    return circles


def find_optimal_circle_centers_spherical(catalogue_df, central_star_mag_lim=None, offset=1.8,
                                          mag_col='gaiaV', ra_col='RAdeg', de_col='DEdeg'):
    """
    Find optimal circle centers using true spherical angular distances.

    Args:
        catalogue_df (DataFrame): Catalogue of stars.
        central_star_mag_lim (float, optional): Magnitude limit for central star.
        offset (float): Radius of each circle in degrees.
        mag_col (str): Magnitude column name.
        ra_col (str): RA column name (degrees).
        de_col (str): Dec column name (degrees).

    Returns:
        list of tuples: List of (RA, Dec) circle centers.
    """
    circles = []

    # Filter stars
    remaining_stars = catalogue_df[
        (catalogue_df['sourceFlag'] == 9)
    ].copy()
    remaining_stars.sort_values(by=mag_col, inplace=True)

    star_coords = SkyCoord(ra=remaining_stars[ra_col].values * u.deg,
                           dec=remaining_stars[de_col].values * u.deg)

    while len(star_coords) > 0:
        # Take brightest remaining star
        center_ra = star_coords[0].ra.deg
        center_dec = star_coords[0].dec.deg
        center = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg)

        # Check overlap with previous circles
        if circles:
            previous_centers = SkyCoord([ra * u.deg for ra, _ in circles],
                                        [dec * u.deg for _, dec in circles])
            sep_to_others = center.separation(previous_centers)
            if np.any(sep_to_others < offset * u.deg):
                star_coords = star_coords[1:]  # Skip this star
                continue

        # Add this as a new circle center
        circles.append((center_ra, center_dec))

        # Remove all stars within offset radius of this center
        separations = center.separation(star_coords)
        star_coords = star_coords[separations >= offset * u.deg]

    return circles



def find_optimal_circle_centers_spherical_n(catalogue_df, offset=1.8, central_star_mag_lim=None,
                                               mag_col='gaiaV',
                                               ra_col='RAdeg',
                                               de_col='DEdeg'):
    """
    Optimized version: Place non-overlapping fields centered on guide stars to
    maximize coverage of science targets.

    Much faster via precomputed spatial indexing.

    Returns:
        List of (RA, Dec) tuples for selected field centers.
    """
    guide_df = catalogue_df[catalogue_df['sourceFlag'] == 9].copy()
    science_df = catalogue_df[catalogue_df['sourceFlag'] != 9].copy()

    if guide_df.empty or science_df.empty:
        return []

    guide_coords = SkyCoord(ra=guide_df[ra_col].values * u.deg,
                            dec=guide_df[de_col].values * u.deg)
    science_coords = SkyCoord(ra=science_df[ra_col].values * u.deg,
                              dec=science_df[de_col].values * u.deg)

    # Step 1: Precompute which science targets are within offset of which guide stars
    idx_center, idx_target, sep2d, _ = guide_coords.search_around_sky(guide_coords, offset * u.deg)

    # Map: guide index → list of science indices it can cover
    from collections import defaultdict
    guide_to_science = defaultdict(set)
    for i_center, i_target in zip(idx_center, idx_target):
        guide_to_science[i_target].add(i_center)

    # Step 2: Track coverage
    covered = np.zeros(len(science_df), dtype=bool)
    placed_centers = []
    blocked_guides = np.zeros(len(guide_df), dtype=bool)  # avoid placing too close to existing centers

    guide_tree = guide_coords  # For checking field overlap

    while True:
        best_guide = None
        best_new_targets = set()
        best_count = 0

        for gid, science_ids in guide_to_science.items():
            if blocked_guides[gid]:
                continue

            # Remove already-covered science targets
            uncovered = [i for i in science_ids if not covered[i]]
            if len(uncovered) == 0:
                continue

            if len(uncovered) > best_count:
                best_count = len(uncovered)
                best_guide = gid
                best_new_targets = uncovered

        if best_guide is None:
            break

        # Place field at best guide star
        ra = guide_coords[best_guide].ra.deg
        dec = guide_coords[best_guide].dec.deg
        placed_centers.append((ra, dec))

        # Mark science targets as covered
        covered[best_new_targets] = True

        # Mark nearby guide stars as blocked (to avoid overlap)
        new_center = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
        sep_to_guides = new_center.separation(guide_coords)
        blocked_guides |= sep_to_guides < offset * u.deg

    print(f"[FAST] Placed {len(placed_centers)} fields covering {covered.sum()} / {len(science_df)} science targets.")
    return placed_centers


def assign_stars_to_circles(catalogue_df, circles,
                            ra_col='GAIA_RA', de_col='GAIA_DEC',
                            src_id_col='CNAME', max_passes=None, radius=1.0,
                            outpath='./resources/WG/'):

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

        # Assign science targets
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


    # Merge with original guide star info to get full row
    #central_guide_stars = pd.merge(cent_df, guide_stars, on=src_id_col, how='left')
    #print(central_guide_stars.columns)

    # Overwrite with correct fieldId and guide flag
    central_guide_stars['fieldId'] = cent_df['fieldId']
    central_guide_stars['guide'] = 1


    # Final catalogue
    catalogue_df = pd.concat([science_targets, central_guide_stars], ignore_index=True)
    #print('----------------------------------------------')
    #print(f"Nr. of central guide stars added: {len(cent_df)}")

    if outpath is not None:
        guide_look_up_df.to_parquet(outpath + 'guide_look_up.parquet', index=False)

    return catalogue_df, guide_look_up_df



def assign_batch_id_optimized_o(circle_coords_df,
                              fields_selection_df,
                              batch_size=5,
                              merge_threshold=3,
                              ra_col='RAdeg',
                              de_col='DEdeg'):

    # Improved initial estimate of clusters
    #n_clusters = max(1, int(np.ceil(len(circle_coords_df) / (batch_size * 0.9))))
    preferred_batch_size = 3
    n_clusters = max(1, int(np.ceil(len(circle_coords_df) / preferred_batch_size)))

    while True:
        # Fit K-Means
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(circle_coords_df[
                                                                   [ra_col+'_c',
                                                                    de_col+'_c']
                                                                   ])
        circle_coords_df['batchId'] = kmeans.labels_

        # Check the size of the largest batch
        largest_batch_size = circle_coords_df['batchId'].value_counts().max()

        # Break loop if largest batch size is within limit, adjust otherwise
        if largest_batch_size <= batch_size:
            break
        else:
            n_clusters_adjustment = max(1, int((largest_batch_size - batch_size) / batch_size))
            n_clusters += n_clusters_adjustment

    # Aggressive batch merging
    centroids = kmeans.cluster_centers_
    for batch_id in range(n_clusters):
        batch_mask = circle_coords_df['batchId'] == batch_id
        #if batch_mask.sum() < merge_threshold:
        if batch_mask.sum() != preferred_batch_size:            

            small_batch_centroid = centroids[batch_id]
            distances = np.linalg.norm(centroids - small_batch_centroid, axis=1)
            distances[batch_id] = np.inf  # Ignore distance to itself

            # Merge strategy: Consider more than nearest batch for potential merge
            sorted_distances_indices = np.argsort(distances)
            for nearest_batch_id in sorted_distances_indices:
                if (circle_coords_df[
                    circle_coords_df['batchId'] ==
                    nearest_batch_id].shape[0] +
                    batch_mask.sum() <= batch_size + 1):

                    circle_coords_df.loc[batch_mask, 'batchId'] = nearest_batch_id
                    break

    # Post-merge optimization: Repeat merging for any remaining small batches

    # Remap batch IDs to be continuous
    unique_batch_ids = sorted(circle_coords_df['batchId'].unique())
    batch_id_mapping = {old_id: new_id + 1 for new_id, old_id in enumerate(unique_batch_ids)}
    circle_coords_df['batchId'] = circle_coords_df['batchId'].map(batch_id_mapping)

    # Mapping batches to fields_selection_df
    field_to_batch_map = dict(zip(circle_coords_df['fieldId_c'], circle_coords_df['batchId']))
    fields_selection_df['batchId'] = fields_selection_df['fieldId'].map(field_to_batch_map)

    return fields_selection_df.sort_values(by=['batchId', ra_col]), len(unique_batch_ids)



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
    import numpy as np
    from scipy.spatial import distance

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


def assign_exoplanet_batch_id(fields_df,
                              circle_coord_df,
                              max_exo_fields=5,
                              exo_radius=5.0,
                              candidates_only=False,
                              ra_col='RAdeg',
                              de_col='DEdeg'):

    def _create_even_sub_batches(batch_fields, max_exo_fields):
        sub_batches = []
        total_fields = len(batch_fields)

        if total_fields <= max_exo_fields:
            return [batch_fields]

        optimal_sub_batch_size = max(max_exo_fields, min(max_exo_fields, total_fields // 2 +
                                     total_fields % 2))

        start_index = 0
        while start_index < total_fields:
            end_index = start_index + optimal_sub_batch_size
            if end_index > total_fields or total_fields - start_index <= max_exo_fields:
                end_index = total_fields
            sub_batches.append(batch_fields[start_index:end_index])
            start_index = end_index
        return sub_batches

    exo_field_ids = fields_df[(fields_df['sourceFlag'] == 5) |
                              (fields_df['sourceFlag'] == 7) |
                              (fields_df['sourceFlag'] == 4)]['fieldId'].unique()

    if candidates_only:
        exo_field_ids = fields_df[fields_df['sourceFlag'] == 7]['fieldId'].unique()

    batch_id = 1
    all_batch_ids = []
    fields_collection = []

    for exo_field_id in exo_field_ids:
        exo_field = fields_df[fields_df['fieldId'] == exo_field_id]
        exo_field_guide = exo_field[exo_field['sourceFlag'].isin([4, 5, 7])]
        exo_field_coords = (exo_field_guide[ra_col].values[0], exo_field_guide[de_col].values[0])

        batch_fields = []
        for field_id in fields_df['fieldId'].unique():
            field = fields_df[fields_df['fieldId'] == field_id]
            field_guide = field[field['guide'] == 1]
            if not field_guide.empty:
                field_coords = (field_guide[ra_col].values[0], field_guide[de_col].values[0])
                distance = np.sqrt((exo_field_coords[0] - field_coords[0])**2 +
                                   (exo_field_coords[1] - field_coords[1])**2)
                if distance <= exo_radius:
                    batch_fields.append((field_id, distance))
        fields_collection.append(batch_fields)

    fields_collection = [[x[0] for x in sorted(lst, key=lambda y: y[1])]
                         for lst in fields_collection]

    result_lists = []
    seen = set()
    for lst in fields_collection:
        new_list = [x for x in lst if not (x in seen or seen.add(x))]
        if len(new_list) == 1 and result_lists:
            result_lists[-1].extend(new_list)
        else:
            result_lists.append(new_list)
    fields_collection = result_lists

    for batch_fields in fields_collection:
        if batch_id in all_batch_ids:
            fields_df['batchId'] += fields_df['batchId'] >= batch_id
            circle_coord_df['batchId'] += circle_coord_df['batchId'] >= batch_id
            all_batch_ids = [x + 1 if x >= batch_id else x for x in all_batch_ids]

        sub_batch_fields = _create_even_sub_batches(batch_fields, max_exo_fields)
        for sub_batch in sub_batch_fields:
            fields_df.loc[fields_df['fieldId'].isin(sub_batch), 'batchId'] = batch_id
            circle_coord_df.loc[circle_coord_df['fieldId_c'].isin(sub_batch),
                                'batchId'] = batch_id
            all_batch_ids.append(batch_id)
            batch_id += 1

    for _, circle in circle_coord_df.iterrows():
        if ~np.isnan(circle['batchId']):
            ra_circle, de_circle = circle[ra_col + '_c'], circle[de_col+'_c']
            in_circle = np.sqrt((fields_df[ra_col] - ra_circle)**2 +
                                (fields_df[de_col] - de_circle)**2) < 0.98

            for i, row in fields_df[in_circle].iterrows():
                star_coords = (row[ra_col], row[de_col])
                distances = circle_coord_df.apply(
                    lambda x: np.sqrt((star_coords[0] - x[ra_col + '_c'])**2 +
                                      (star_coords[1] - x[de_col + '_c'])**2), axis=1)
                closest_circle = circle_coord_df.loc[distances.idxmin()]
                fields_df.at[i, 'batchId'] = closest_circle['batchId']
                fields_df.at[i, 'fieldId'] = closest_circle['fieldId_c']
        else:
            continue

    fields_df['batchId'] = fields_df['batchId'].astype('Int64')
    fields_df = fields_df.dropna(subset=['batchId'])
    circle_coord_df['batchId'] = circle_coord_df['batchId'].astype('Int64')
    circle_coord_df = circle_coord_df.dropna(subset=['batchId'])

    return fields_df, circle_coord_df, len(fields_df['batchId'].unique())


def assign_batch_id(circle_coords_df,
                    fields_selection_df,
                    batch_size=5,
                    merge_threshold=3,
                    ra_col='RAdeg',
                    de_col='DEdeg'):

    """
    Assigns batch IDs to circles using K-Means clustering, ensuring that no batch is larger
    than `batch_size` and merging smaller batches (less than `merge_threshold`) with
    adjacent batches, allowing merged batches to have up to `batch_size + 1` circles.

    Args:
        circle_coords_df (pd.DataFrame): DataFrame with coordinates of circles.
                                         Must include 'ra_c' and 'dec_c' columns.
        fields_selection_df (pd.DataFrame): DataFrame for mapping batch IDs.
                                            Must include 'fieldId' column.
        batch_size (int): Maximum number of circles per batch. Defaults to 5.
        merge_threshold (int): Batches smaller than this will be merged. Defaults to 3.

    Returns:
        pd.DataFrame: Updated `fields_selection_df` with 'batchId', sorted by 'batchId' and
                              'RAdeg'.
    """

    # Estimate initial number of clusters
    n_clusters = max(1, np.ceil(len(circle_coords_df) / batch_size).astype(int))

    while True:
        # Fit K-Means
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(circle_coords_df[[ra_col,
                                                                                     de_col]
                                                                                    ])
        circle_coords_df['batchId'] = kmeans.labels_

        # Check the size of the largest batch
        largest_batch_size = circle_coords_df['batchId'].value_counts().max()

        # Break loop if largest batch size is within limit
        if largest_batch_size <= batch_size:
            break

        # Increase number of clusters for next iteration
        n_clusters += 1

    # Merge small batches
    centroids = kmeans.cluster_centers_
    for batch_id in range(n_clusters):
        # Identify small batches
        batch_mask = circle_coords_df['batchId'] == batch_id
        #if batch_mask.sum() < merge_threshold:
        if batch_mask.sum() != preferred_batch_size:            # Calculate centroid of the small batch
            small_batch_centroid = centroids[batch_id]

            # Find nearest batch centroid
            distances = np.linalg.norm(centroids - small_batch_centroid, axis=1)
            distances[batch_id] = np.inf  # Ignore distance to itself
            nearest_batch_id = np.argmin(distances)

            # Merge if combined size is <= batch_size + 1
            nearest_batch_mask = circle_coords_df['batchId'] == nearest_batch_id
            if nearest_batch_mask.sum() + batch_mask.sum() <= batch_size + 1:
                circle_coords_df.loc[batch_mask, 'batchId'] = nearest_batch_id

    # Remap batch IDs to be continuous from 1 to n
    unique_batch_ids = sorted(circle_coords_df['batchId'].unique())
    batch_id_mapping = {old_id: new_id + 1 for new_id, old_id in enumerate(unique_batch_ids)}
    circle_coords_df['batchId'] = circle_coords_df['batchId'].map(batch_id_mapping)

    # Mapping batches to fields_selection_df based on fieldId
    field_to_batch_map = dict(zip(circle_coords_df['fieldId_c'], circle_coords_df['batchId']))
    fields_selection_df['batchId'] = fields_selection_df['fieldId'].map(field_to_batch_map)

    final_dataset = fields_selection_df.sort_values(by=['batchId', ra_col])

    return final_dataset, len(unique_batch_ids)


def get_science_covering_guide_star_circles_fixed(
    df,
    ra_col='GAIA_RA',
    dec_col='GAIA_DEC',
    mag_col='GAIA_MAG_G',
    flag_col='sourceFlag',
    central_star_mag_lim=16.0,
    spacing=1.5,
    match_radius=1.0,
    circles_df=None,
    filter_circles_to_science_bbox=True
):
    """
    Returns guide-star-centered circle coordinates (RA, Dec) that cover science targets.
    If `circles_df` is provided, snaps those circle centres to the nearest guide stars.
    Otherwise generates an evenly spaced grid inside the science bounding box and snaps to guides.

    Parameters
    ----------
    df : DataFrame
        Input with RA, Dec, and sourceFlag columns (9 = guide star).
    ra_col, dec_col : str
        Column names for RA/Dec (deg) in `df`.
    flag_col : str
        Column name for the sourceFlag (9 = guide star).
    spacing : float
        Grid spacing in degrees (used only when `circles_df` is None).
    match_radius : float
        Max separation (deg) for snapping grid/circle points to a guide star, and
        also used as the science coverage radius for each circle.
    circles_df : DataFrame or None
        Optional DataFrame of precomputed circle centres to be snapped to guides.
    filter_circles_to_science_bbox : bool
        If True and `circles_df` is provided, only use circles within the science
        bounding box (inclusive). Keeps behaviour close to the original grid logic.

    Returns
    -------
    List[Tuple[float, float]]
        List of (RA, Dec) tuples for valid circle centers.
    """

    from collections import defaultdict

    # Split guide and science stars
    guides = df[df[flag_col] == 9].copy()
    science = df[df[flag_col] != 9].copy()
    guides = guides[guides[mag_col] <= central_star_mag_lim]

    if science.empty or guides.empty:
        return []  # nothing to do

    # Science bounding box
    ra_min, ra_max = science[ra_col].min(), science[ra_col].max()
    dec_min, dec_max = science[dec_col].min(), science[dec_col].max()

    # --- Candidate circle centres (grid OR provided circles) ---
    if circles_df is None:
        # Generate grid over the science bbox
        ra_grid = np.arange(ra_min, ra_max + spacing, spacing)
        dec_grid = np.arange(dec_min, dec_max + spacing, spacing)
        if ra_grid.size == 0 or dec_grid.size == 0:
            return []
        ra_mesh, dec_mesh = np.meshgrid(ra_grid, dec_grid)
        cand_ra = ra_mesh.ravel()
        cand_dec = dec_mesh.ravel()
    else:
        # Use provided circle centres
        if filter_circles_to_science_bbox:
            mask = (
                (circles_df[ra_col].between(ra_min, ra_max, inclusive='both')) &
                (circles_df[dec_col].between(dec_min, dec_max, inclusive='both'))
            )
            sub = circles_df.loc[mask]
        else:
            sub = circles_df
        if sub.empty:
            return []
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
    if not np.any(mask_within_radius):
        return []

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

    if not circle_map:
        return []

    final_circle_idxs = list(circle_map.keys())
    final_centers = unique_circles[final_circle_idxs]

    # Return as list of (RA, Dec) tuples
    return list(zip(final_centers.ra.deg, final_centers.dec.deg))


def get_science_covering_guide_star_circles(
    df,
    ra_col='GAIA_RA',
    dec_col='GAIA_DEC',
    flag_col='sourceFlag',
    cname_col='CNAME',
    gaia_id_col='GAIA_ID',
    mag_col='GAIA_MAG_G',
    central_star_mag_lim=16.0,  # max mag for central guide star
    spacing=1.5,               # grid spacing (deg) in tangent plane
    circle_radius=1.0,         # circle radius (deg)
    match_radius=1.0,          # grid node can snap to any guide within this (deg)
    move_radius=None,          # how far a chosen center may "nudge" to another guide (deg); 
    min_science_per_circle=1,  # drop very weak circles
    random_state=42
):
    """
    Science-footprint grid -> smart snap to guides -> augment for full cover ->
    greedy set cover -> redundancy prune -> local refinement (swap to nearby guides)
    to remove unnecessary circles. Non-failing: warns & returns uncoverable targets.

    Returns
    -------
    centers : List[(float,float)]
        Selected (RA_deg, Dec_deg) of guide-centred circles after refinement.
    uncovered_df : pd.DataFrame
        Uncoverable science targets (no guide within circle_radius) with CNAME/GAIA_ID/RA/Dec.
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

    # --- Tangent plane around science median ---
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

    # --- Helpers ---
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

    # --- Smart snapping: for each grid node, choose guide within match_radius that covers MOST science ---
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
        uncv = pd.DataFrame(columns=['CNAME','GAIA_ID','RA','Dec'])

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

    # Redundancy prune
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

    chosen_rel = prune(chosen_rel)

    # Local refinement (nudge to nearby guides to remove another circle)
    chosen = np.array(chosen_rel, dtype=int)

    # Precompute chosen center XY for neighborhood queries
    chosen_xy = np.column_stack([
        gs_coords[cand_idxs[chosen]].transform_to(off).lon.deg,
        gs_coords[cand_idxs[chosen]].transform_to(off).lat.deg
    ])

    improved = True
    while improved and len(chosen) > 1:
        improved = False
        # current coverage matrix of chosen
        cov_ch = cov_eff[:, chosen]  # (Neff, Nchosen)
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
                    # Accept improvement: update chosen to reflect alt swap + pruning
                    # Rebuild chosen indices from active columns relative to the trial matrix
                    # First, update the corresponding chosen index to the new guide candidate
                    # We need to ensure the alt guide is in cand_idxs; if not, append it.
                    if g_alt not in cand_idxs:
                        cand_idxs = np.append(cand_idxs, g_alt)
                        # extend cov_eff with this new column
                        new_col = alt_vec[:, None]
                        cov_eff = np.hstack([cov_eff, new_col])
                        # map trial column index (col_idx) to this new absolute column (last)
                        chosen[col_idx] = cov_eff.shape[1] - 1
                    else:
                        # set chosen[col_idx] to absolute index of existing column for g_alt
                        # find its absolute column in cov (cand_idxs alignment)
                        abs_idx = np.where(cand_idxs == g_alt)[0][0]
                        # map to cov_eff's column indices (same order as cand_idxs after filters)
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
                    break  # break alt loop; restart outer while
            if improved:
                break
        # end for col_idx
    # end while improved

    final_coords = gs_coords[cand_idxs[chosen]] if len(chosen) else SkyCoord([], [], unit='deg')
    centers = list(zip(final_coords.ra.deg, final_coords.dec.deg))
    return centers


def define_fields(catalogue_df,
                  fixed_grid=False,
                  circles_df=None,
                  batch_size=3,
                  max_passes=5,
                  merge_threshold=3,
                  radius=0.98,
                  central_star_mag_lim=16.,
                  grid_offset=1.0,
                  n_clusters=None,
                  exo_subset=False,
                  exo_radius=2.5,
                  max_exo_fields=4,
                  candidates_only=False,
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
        tuple: (stars_in_circles, circles_df)
    """

    if circles_df is None and fixed_grid:
        print('Using fixed grid approach for field definition.')
        circles = get_science_covering_guide_star_circles_fixed(
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
        circles = get_science_covering_guide_star_circles_fixed(
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
        print('Using automatic approach for field definition.')
        circles = get_science_covering_guide_star_circles(
            catalogue_df,
            ra_col='GAIA_RA',
            dec_col='GAIA_DEC',
            flag_col='sourceFlag',
            mag_col='GAIA_MAG_G',
            central_star_mag_lim=central_star_mag_lim,
            spacing=grid_offset,
            match_radius=1.0)

    # Assign stars to circles
    catalogue_df, guide_look_up_df = assign_stars_to_circles(catalogue_df, circles,
                                                             max_passes=max_passes,
                                                             radius=radius,
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
        print('Nr. of targets: \n'
           f'Science: {stars_in_circles[~stars_in_circles["sourceFlag"].isin([9])].shape[0]}\n'
           f'Guide: {stars_in_circles[stars_in_circles["sourceFlag"].isin([9])].shape[0]}')

    # Optionally assign batches to exoplanet candidate fields
    if exo_subset:
        (stars_in_circles,
         circles_df,
         num_batches) = assign_exoplanet_batch_id(stars_in_circles,
                                                  circle_coord_df=circles_df,
                                                  exo_radius=exo_radius,
                                                  max_exo_fields=max_exo_fields,
                                                  candidates_only=candidates_only,
                                                  ra_col=ra_col, de_col=de_col)
                                                  
        print(f'Nr. of field batches: {num_batches}')
        print('Nr. of targets in exoplanet batches:', stars_in_circles.shape)
        return stars_in_circles, circles_df

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
