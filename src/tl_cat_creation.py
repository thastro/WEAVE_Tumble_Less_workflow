import pandas as pd
import numpy as np
import os
from astropy.io import fits
import healpy as hp
from astropy.table import Table
from workflow.utils.get_resources import get_guide_cat


def query_weave_guide_stars(healpix_ids,
                            out_dir='./resources/WG',
                            split_nr=5,
                            mag_low=14.5,
                            mag_high=16,
                            cache=False,
                            neighbour_hp=True):
    """
    Query guide stars from the WEAVE guide star catalogue and save them in either one or
    multiple parquet files, depending on the `split_nr` parameter. The orginal fits files
    are removed after reading the data, unless `cache` is set to True. Also filters the
    stars based on the `mag_low` and `mag_high` parameters. The mag_low and mag_high should
    not be more than 2 magnitudes apart, due to autoguider requirements.

    Args:
        healpix_ids (list): List of HEALPix IDs to query the guide stars for.
        out_dir (str): Directory where the output parquet files will be saved.
        split_nr (int, optional): Desired nr. of fits files to be combined in the same
                                  parquet file. Defaults to 5.
        mag_low (float): Minimum magnitude of the guide stars to be included.
                                   Defaults to 14.5.
        mag_high (float): Maximum magnitude of the guide stars to be included.
                                  Defaults to 16.
        cache (bool): If True, the original fits files will not be removed after
                                reading. Defaults to False.
        neighbour_hp (bool): If True, the function will also query the neighbouring
                             HEALPix indices for guide stars. Defaults to True.
    """

    nside = 32
    nest = True

    if neighbour_hp:
        expanded_ids = set(healpix_ids)
        for hp_id in healpix_ids:
            neighbours = hp.get_all_neighbours(nside, hp_id, nest=nest)
            valid_neighbours = [n for n in neighbours if n >= 0]  # exclude -1 (invalid)
            expanded_ids.update(valid_neighbours)
        healpix_ids = sorted(expanded_ids)

    dfs = []
    hpids = []
    batch = 0

    for healpix_id in healpix_ids:
        try:
            file_path = get_guide_cat(healpix_index=healpix_id, directory=out_dir)
        except SystemExit:
            print(f"SystemExit triggered (likely HTTP 404) for HEALPix ID {healpix_id}, skipping.")
            continue
        except Exception as e:
            print(f"Unexpected error retrieving guide stars for HEALPix ID {healpix_id}: {e}")
            continue

            continue
        if file_path is None:
            continue
        try:
            with fits.open(file_path) as hdul:
                data = hdul[1].data
                df = Table(data).to_pandas()
                df = df[((df['GAIA_MAG_G'] >= mag_low) &
                         (df['GAIA_MAG_G'] <= mag_high)) &
                        (df['TARGPROG'].isin(['ANY']))]
                dfs.append(df)
            hpid = file_path.split('/')[-1].replace('.fits', '').split('_')[1]
            hpids.append(hpid)
        finally:
            if os.path.exists(file_path) and not cache:
                os.remove(file_path)

        if split_nr is not None and len(dfs) == split_nr:
            combined = pd.concat(dfs, ignore_index=True)
            parquet_name = 'WG_all.parquet' if split_nr is None else 'WG_' + '_'.join(hpids) + '.parquet'
            parquet_path = os.path.join(out_dir, parquet_name)
            combined.to_parquet(parquet_path, index=False)
            dfs.clear()
            hpids.clear()
            batch += 1

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        parquet_name = 'WG_all.parquet' if split_nr is None else 'WG_' + '_'.join(hpids) + '.parquet'
        parquet_path = os.path.join(out_dir, parquet_name)
        combined.to_parquet(parquet_path, index=False)


def decode_byte_columns(df: pd.DataFrame, to_string_dtype: bool = True) -> pd.DataFrame:
    obj_cols = [c for c in df.columns if df[c].dtype == 'object']
    for c in obj_cols:
        s = df[c]
        if s.map(lambda v: isinstance(v, (bytes, np.bytes_))).any():
            df[c] = s.map(
                lambda v: (
                    "" if isinstance(v, (bytes, np.bytes_)) and len(v) == 0
                    else (v.decode("utf-8", "ignore") if isinstance(v, (bytes, np.bytes_)) else v)
                )
            )
    if to_string_dtype:
        # only cast columns that are all strings after decode
        for c in obj_cols:
            if df[c].map(lambda v: (v is None) or isinstance(v, str)).all():
                df[c] = df[c].astype("string")
    return df


def combine_wg_with_input_cat(input_df, wg_df, id_col='CNAME'):

    """
    Combine the input DataFrame with the WEAVE guide star DataFrame, ensuring that only
    unique entries are retained based on the specified ID column. The function will
    append new rows from the WEAVE guide star DataFrame that do not exist in the input
    DataFrame, while preserving the columns that match the input DataFrame.

    Args:
        input_df (pd.DataFrame): The input DataFrame to which new rows will be added.
        wg_df (pd.DataFrame): The WEAVE guide star DataFrame containing guide stars from
                              previously queried helapix regions.
        id_col (str, optional): The column name used to identify unique entries.
                                Defaults to 'CNAME'.

    Returns:
        pd.DataFrame: A new DataFrame that combines the input DataFrame with new rows
                      from the WEAVE guide star DataFrame, ensuring uniqueness based on
                      the specified ID column.
    """
    input_df = input_df.copy()
    input_df.loc[:, 'sourceFlag'] = 0

    new_rows = wg_df[~wg_df[id_col].isin(input_df[id_col])].copy()

    # Select only columns that match input_df
    common_cols = [col for col in input_df.columns if col in new_rows.columns]
    new_rows_subset = new_rows[common_cols].copy()
    new_rows_subset.loc[:, 'sourceFlag'] = 9

    concat_df = pd.concat([input_df, new_rows_subset], ignore_index=True)

    del input_df, wg_df, new_rows, new_rows_subset

    return concat_df


if __name__ == "__main__":
    # Do Nothing
    pass
