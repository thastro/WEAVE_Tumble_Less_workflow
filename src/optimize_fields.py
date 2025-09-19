import os
import contextlib
import numpy as np
import pandas as pd
from collections import namedtuple


Result = namedtuple("Result", ["offset", "n_science", "n_fields",
                    "coverage_frac", "field_stars", "lut_df"])


def evaluate_offset(offset,
                    cat_comb_df,
                    fs,
                    bright_stars_df=None,
                    max_passes=5,
                    batch_size=3,
                    merge_threshold=3,
                    radius=0.98,
                    central_star_mag_lim=18.,
                    n_clusters=None,
                    exo_subset=False,
                    exo_radius=2.5,
                    max_exo_fields=3,
                    candidates_only=False,
                    ra_col='GAIA_RA', de_col='GAIA_DEC',
                    mag_col='GAIA_MAG_G', src_id_col='CNAME',
                    verbose=True,
                    fixed_grid=False):
    """Run define_fields once and compute metrics for this offset."""

    if fixed_grid:
        if verbose:
            print(f"Using fixed grid with offset={offset}")
        fixed_grid = True
    else:
        if verbose:
            print(f"Using adaptive grid with offset={offset}")
        fixed_grid = False

    # ---- suppress prints from define_fields ----
    with open(os.devnull, 'w') as devnull, \
         contextlib.redirect_stdout(devnull), \
         contextlib.redirect_stderr(devnull):
        field_stars, lut_df = fs.define_fields(
            cat_comb_df,
            max_passes=max_passes,
            batch_size=batch_size,
            merge_threshold=merge_threshold,
            radius=radius,
            central_star_mag_lim=central_star_mag_lim,
            grid_offset=offset,
            n_clusters=n_clusters,
            exo_subset=exo_subset,
            exo_radius=exo_radius,
            max_exo_fields=max_exo_fields,
            candidates_only=candidates_only,
            ra_col=ra_col, de_col=de_col,
            mag_col=mag_col, src_id_col=src_id_col,
            verbose=False,                 # if the function supports it
            fixed_grid=fixed_grid
        )

    # Identify science vs field rows (assuming sourceFlag==9 denotes field/guide definitions)
    n_science = field_stars[~field_stars["sourceFlag"].isin([9])].shape[0]
    n_fields = field_stars[field_stars["sourceFlag"].isin([9])].shape[0]

    total_candidates = cat_comb_df[cat_comb_df['sourceFlag'] != 9].shape[0]
    coverage_frac = n_science / max(1, total_candidates)

    return Result(offset, n_science, n_fields, coverage_frac, field_stars, lut_df)


def sweep_offsets(cat_comb_df,
                  fs,
                  offsets=None,
                  **kwargs):
    """
    Evaluate a set of offsets and return a DataFrame of metrics
    plus the best result according to lexicographic objective:
      (coverage_frac desc, n_fields asc, offset asc).
    """
    if offsets is None:
        # Reasonable sweep for grid/snapping-style offsets. Adjust as needed.
        # Use a finer grid if results are sensitive.
        # offsets = np.round(np.linspace(0.30, 0.90, 25), 3)
        offsets = np.round(np.arange(0.30, 1.30 + 0.01, 0.01), 2)

    results = []
    best = None

    for off in offsets:
        r = evaluate_offset(off, cat_comb_df, fs, **kwargs)
        results.append(r)

        # Lexicographic ranking:
        if best is None:
            best = r
        else:
            # Prefer higher coverage
            if r.coverage_frac > best.coverage_frac + 1e-9:
                best = r
            elif abs(r.coverage_frac - best.coverage_frac) <= 1e-9:
                # Among equal coverage, prefer fewer fields
                if r.n_fields < best.n_fields:
                    best = r
                elif r.n_fields == best.n_fields:
                    # Prefer smaller offset for stability/reproducibility
                    if r.offset < best.offset:
                        best = r

    # Pack results into a table for inspection
    df = pd.DataFrame([{
        "offset": r.offset,
        "n_science": r.n_science,
        "n_fields": r.n_fields,
        "coverage_frac": r.coverage_frac,
    } for r in results]).sort_values(["coverage_frac", "n_fields", "offset"],
                                     ascending=[False, True, True],
                                     ignore_index=True)

    # Attach the full best artifacts so you can use them immediately
    best_payload = {
        "offset": best.offset,
        "field_stars": best.field_stars,
        "lut_df": best.lut_df,
        "metrics": {
            "n_science": best.n_science,
            "n_fields": best.n_fields,
            "coverage_frac": best.coverage_frac
        },
        "sweep_table": df
    }
    return best_payload


if __name__ == "__main__":
    pass
