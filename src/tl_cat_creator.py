import os
import logging

import pandas as pd
import numpy as np

import healpy as hp
from astropy.io import fits
from astropy.table import Table


def fill_wasp_fits_template(catalogue_path,
                            fits_template_path,
                            output_path='./',
                            mapping={},
                            header_updates=None,
                            defaults=None,
                            allow_passthrough_same_name=True,
                            extname=None,
                            fill_missing_with_tnull=True,
                            chunk_read_kwargs=None,
                            dataframe_postprocess=None,
                            passthrough_static=None,
                            checksum=True,
                            overwrite=False):
    """
    Write a WASP-compliant FITS catalogue from a user-provided table (CSV or Parquet)
    by mapping its columns onto a survey-specific WASP Catalogue Template. Columns
    with identical names are auto-mapped; differing names must be provided via
    `mapping`. Template columns missing from the source can be filled from `defaults`
    or, when available, from the template's per-column TNULLn values (recommended for
    strict WASP compatibility).

    Args:
        catalogue_path (str): Path to the input catalogue (.csv or .parquet).
        fits_template_path (str): Path to the WASP FITS template
                                  (Primary HDU + binary table in HDU 1).
        output_path (str): Destination path (including filename) for the output FITS.
        mapping (dict): Mapping from source columns to template columns.
            - Basic: `{"dr3_source_id": "GAIA_ID", "ra": "GAIA_RA", ...}
            - Advanced (with transform): {"gmag": ("GAIA_MAG_G", lambda s: s.astype(float))}
        header_updates (dict, optional): Key-value pairs merged into the
                                         primary header (HDU 0).
        defaults (dict, optional): Fallback values for template columns when source data
            or a mapping is not provided.
        allow_passthrough_same_name (bool, optional): If True, copy source columns whose
            names match the template exactly. Defaults to True.
        extname (str, optional): If provided, set as `EXTNAME`
                                 on HDU 1 (e.g., "WS2025XX-001").
        fill_missing_with_tnull (bool, optional): If True, fill missing values using the
            template column's TNULLn where defined. Defaults to True.
        chunk_read_kwargs (dict, optional): Extra keyword arguments forwarded to `pd.read_csv`
            (useful for large CSVs and chunked reads).
        dataframe_postprocess (callable, optional): Hook applied after loading, before mapping:
            `df = dataframe_postprocess(df)`.
        passthrough_static (dict, optional): Constant values to set on specific
                                             template columns (e.g., {"APS_FLAG": ""}).
        checksum (bool, optional): If True, write FITS with checksums.
                                   Defaults to True.
        overwrite (bool, optional): If True, allow overwriting an existing `output_path`.
                                    Defaults to False.

    Returns:
        None: The function writes the FITS file to output_path.

    Notes:
        - Auto-mapping applies only where source and template column names are identical;
          all other mappings must be specified in mapping.
        - When both defaults and TNULLn exist, the fill order is typically: mapped data
          → defaults → TNULLn (if fill_missing_with_tnull=True).
        - Use `dataframe_postprocess for last-mile cleaning (dtype fixes, unit conversions,
          renaming before mapping, etc.).
    """

    # Safety: do not overwrite the template by mistake
    if os.path.abspath(output_path) == os.path.abspath(fits_template_path):
        raise ValueError("Output path must not be the same as the FITS template path.")

    # Load catalogue
    ext_in = os.path.splitext(catalogue_path)[1].lower()
    if ext_in in (".parquet", ".pq"):
        df = pd.read_parquet(catalogue_path)
    elif ext_in in (".csv", ".txt"):
        if chunk_read_kwargs:
            chunks = pd.read_csv(catalogue_path, **chunk_read_kwargs)
            if hasattr(chunks, "__iter__"):
                df = pd.concat(chunks, ignore_index=True)
            else:
                df = chunks
        else:
            df = pd.read_csv(catalogue_path)
    else:
        raise ValueError(f"Unsupported catalogue extension: {ext_in}")

    if dataframe_postprocess is not None:
        df = dataframe_postprocess(df)

    # Open FITS template HDU
    with fits.open(fits_template_path, mode="readonly") as hdul:
        prim_hdr = hdul[0].header.copy()
        tab_hdr = hdul[1].header.copy()

        # Detach columns from file: rebuild ColDefs from name/format only
        tmpl_cols = fits.ColDefs([
            fits.Column(name=c.name, format=c.format)
            for c in hdul[1].columns
        ])
        tmpl_col_names = list(tmpl_cols.names)

    # Primary header updates (on copy)
    if header_updates:
        for k, v in header_updates.items():
            prim_hdr[k] = v
    prim_hdr["EXTEND"] = True

    defaults = defaults.copy() if defaults else {}
    passthrough_static = passthrough_static.copy() if passthrough_static else {}

    # Build the output data frame (template columns only, in order)
    out_df = pd.DataFrame(index=np.arange(len(df)), columns=tmpl_col_names, dtype="object")

    def _fits_tform_to_dtype(tform):
        if tform is None:
            return (None, None)
        tform = tform.strip().upper()
        if tform.endswith("A"):
            try:
                n = int(tform[:-1])
            except ValueError:
                n = None
            return (np.bytes_, n)
        return {
            "E": (np.float32, None),
            "D": (np.float64, None),
            "I": (np.int16,  None),
            "J": (np.int32,  None),
            "K": (np.int64,  None),
        }.get(tform, (None, None))

    # Normalize mapping: target -> (source, transform_or_None)
    normalized_map = {}
    for src, tgt in mapping.items():
        if isinstance(tgt, tuple) and len(tgt) == 2:
            target, transform = tgt
        else:
            target, transform = tgt, None
        normalized_map[target] = (src, transform)

    # Fill out_df following priority:
    # 1) explicit mapping, 2) passthrough same-name, 3) static passthrough 4) defaults
    for tgt_col in tmpl_col_names:
        filled = False

        if tgt_col in normalized_map:
            src_col, transform = normalized_map[tgt_col]
            if src_col in df.columns:
                series = df[src_col]
                if transform is not None:
                    series = transform(series)
                out_df[tgt_col] = series
                filled = True
            else:
                logging.warning(f"Source column '{src_col}' for target ' "
                                "{tgt_col}' not found in catalogue.")

        if not filled and allow_passthrough_same_name and tgt_col in df.columns:
            out_df[tgt_col] = df[tgt_col]
            filled = True

        if not filled and tgt_col in passthrough_static:
            out_df[tgt_col] = passthrough_static[tgt_col]
            filled = True

        if not filled and tgt_col in defaults:
            out_df[tgt_col] = defaults[tgt_col]

    # Enforce template types, TNULLs, and ASCII
    prim_hdu = fits.PrimaryHDU(header=prim_hdr)
    bin_hdu = fits.BinTableHDU.from_columns(tmpl_cols, nrows=len(out_df))

    if extname:
        if "EXTNAME" in bin_hdu.header:
            del bin_hdu.header["EXTNAME"]
        bin_hdu.header["EXTNAME"] = extname

    # We need the *template* table header (tab_hdr) for TNULLn, TTYPEn, etc.
    for idx, col in enumerate(tmpl_cols, start=1):
        tgt = col.name
        tform = col.format
        np_dtype, str_len = _fits_tform_to_dtype(tform)
        tnull_key = f"TNULL{idx}"
        tnull_val = tab_hdr.get(tnull_key, None) if fill_missing_with_tnull else None

        if tgt not in out_df.columns:
            logging.warning(f"Target FITS column '{tgt}' "
                            "missing in output frame; filling with defaults/TNULL.")
            if np_dtype is np.bytes_:
                bin_hdu.data[tgt] = np.array(["" for _ in range(len(out_df))],
                                             dtype=f"S{str_len if str_len else 1}")
            else:
                if tnull_val is not None:
                    bin_hdu.data[tgt] = np.array([tnull_val] * len(out_df))
                else:
                    fill = 0
                    if np_dtype in (np.float32, np.float64):
                        fill = np.nan
                    bin_hdu.data[tgt] = np.array([fill] * len(out_df))
            continue

        series = out_df[tgt]

        # Fixed-width ASCII strings
        if np_dtype is np.bytes_:
            ser = series.fillna("").astype(str).values.astype(f"S{str_len if str_len else 1}")
            bin_hdu.data[tgt] = ser
        else:
            # Numerics
            arr = series.values
            if tnull_val is not None:
                mask = pd.isna(arr)
                if mask.any():
                    if np_dtype in (np.int16, np.int32, np.int64):
                        tmp = np.array(arr, dtype="object")
                        tmp[mask] = tnull_val
                        try:
                            arr = tmp.astype(np_dtype)
                        except Exception:
                            arr = np.array([tnull_val if pd.isna(x) else x for x in series],
                                           dtype=np.int64)
                    else:
                        arr = arr.astype(np_dtype, copy=False)
                else:
                    arr = arr.astype(np_dtype, copy=False)
            else:
                try:
                    arr = arr.astype(np_dtype, copy=False) if np_dtype is not None else arr
                except Exception:
                    if np_dtype in (np.float32, np.float64):
                        arr = pd.to_numeric(series, errors="coerce").astype(np_dtype)
                    elif np_dtype in (np.int16, np.int32, np.int64):
                        coerced = pd.to_numeric(series, errors="coerce")
                        arr = coerced.fillna(0).astype(np_dtype)
                    else:
                        pass

            bin_hdu.data[tgt] = arr

    # Preserve the original table header exactly
    bin_hdu.header.clear()
    for card in tab_hdr.cards:
        bin_hdu.header.append(card, end=True)

    # Write output FITS
    hdul_out = fits.HDUList([prim_hdu, bin_hdu])
    hdul_out.writeto(output_path, overwrite=overwrite, checksum=checksum)
    return hdul_out


def add_healpix_id(df,
                   ra_col='ra',
                   de_col='dec',
                   hp_col='healpixId',
                   nside=32,
                   nested=True):
    """
    Adds a new column to the dataframe containing the HEALPix ID based on RA and Dec
    coordinates.

    Parameters:
        df: pandas DataFrame containing the RA and Dec coordinates
        ra_col: name of the column containing the RA coordinates (default 'RA')
        dec_col: name of the column containing the Dec coordinates (default 'Dec')
        hp_col: name of the new column to store HEALPix IDs (default 'healpixId')
        nside: NSIDE parameter for HEALPix (default 32)
        nested: whether to use nested ordering for HEALPix (default True)

    Returns:
     df: The inputted pandas DataFrame with a new column 'healpix_id' containing
         the HEALPix IDs tailored to the WEAVE guide star cataogue HP grid.
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


def query_weave_guide_stars(healpix_ids,
                            out_dir='./resources/WG',
                            base_name='WG_',
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
        base_name (str): Base name for the output parquet files.
        split_nr (int, optional): Desired nr. of fits files to be combined in the same
                                  parquet file. Defaults to 5.
        mag_low (float): Minimum magnitude (Gaia G-band) of the guide stars to be included.
                         Defaults to 14.5.
        mag_high (float): Maximum magnitude (Gaia G-band) of the guide stars to be included.
                          Defaults to 16.
        cache (bool): If True, the downloaded fits files will not be removed after
                      reading. Defaults to False.
        neighbour_hp (bool): If True, the function will also query the neighbouring
                             HEALPix indices for guide stars. This allows for better
                             guide fibre allocatcion later .Defaults to True.
    """

    from workflow.utils.get_resources import get_guide_cat

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
            print(f"SystemExit triggered (likely HTTP 404) for HEALPix ID {healpix_id}, "
                  "skipping.")
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
            if split_nr is None:
                parquet_name = f"{base_name}_all.parquet"
            else:
                parquet_name = f"{base_name}_{'_'.join(hpids)}.parquet"

            parquet_path = os.path.join(out_dir, parquet_name)
            combined.to_parquet(parquet_path, index=False)
            dfs.clear()
            hpids.clear()
            batch += 1

    if dfs:
        combined = pd.concat(dfs, ignore_index=True)
        if split_nr is None:
            suffix = "all"
        else:
            suffix = "_".join(hpids)

        parquet_name = f"{base_name}_{suffix}.parquet"
        parquet_path = os.path.join(out_dir, parquet_name)
        combined.to_parquet(parquet_path, index=False)


def combine_wg_with_input_cat(input_df, wg_df, src_id_col='CNAME'):

    """
    Combine the input DataFrame with the WEAVE guide star DataFrame, ensuring that only
    unique entries are retained based on the specified ID column. The function will
    append new rows from the WEAVE guide star DataFrame that do not exist in the input
    DataFrame, while preserving the columns that match the input DataFrame.

    Args:
        input_df (pd.DataFrame): The input DataFrame to which new rows will be added.
        wg_df (pd.DataFrame): The WEAVE guide star DataFrame containing guide stars from
                              previously queried helapix regions.
        id_col (str, optional): The name of the column containing the unique target identifier.
                                Defaults to 'CNAME'.

    Returns:
        pd.DataFrame: A new DataFrame that combines the input DataFrame with new rows
                      from the WEAVE guide star DataFrame, ensuring uniqueness based on
                      the specified ID column.
    """

    input_df = input_df.copy()
    input_df.loc[:, 'sourceFlag'] = 0

    new_rows = wg_df[~wg_df[src_id_col].isin(input_df[src_id_col])].copy()

    # Select only columns that match input_df
    common_cols = [col for col in input_df.columns if col in new_rows.columns]
    new_rows_subset = new_rows[common_cols].copy()
    new_rows_subset.loc[:, 'sourceFlag'] = 9

    concat_df = pd.concat([input_df, new_rows_subset], ignore_index=True)

    del input_df, wg_df, new_rows, new_rows_subset

    return concat_df


if __name__ == "__main__":

    # Example Usage
    # Load an example input catalogue as a pandas DataFrame.
    with fits.open("./WS2025XX-XXX.fits") as hdu:
        cat_df = Table(hdu[1].data).to_pandas()

    cat_df = add_healpix_id(cat_df, ra_col='GAIA_RA', de_col='GAIA_DEC',
                            hp_col='healpixId', nside=32, nested=True)

    healpix_ids = cat_df['healpixId'].unique()

    # Query guide stars. In real use case neighbour_hp should be set to True,
    # here set to False to limit the number of queried stars.
    query_weave_guide_stars(healpix_ids, split_nr=None,
                            out_dir='./resources/WG',
                            base_name='WG_',
                            mag_low=14.5,
                            mag_high=15,
                            cache=False,
                            neighbour_hp=False)

    wg_df = pd.read_parquet('./resources/WG/WG_all.parquet')
    cat_comb_df = combine_wg_with_input_cat(cat_df, wg_df, src_id_col='CNAME')
    print(cat_comb_df)
