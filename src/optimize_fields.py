import os
import contextlib
import numpy as np
import pandas as pd
from collections import namedtuple




def evaluate_offset(offset,
                    cat_comb_df,
                    fs,
                    bright_stars_df=None,
                    max_passes=5,
                    batch_size=3,
                    merge_threshold=3,
                    radius=0.98,
                    central_star_mag_lim=18.,
                    ra_col='GAIA_RA', de_col='GAIA_DEC',
                    mag_col='GAIA_MAG_G', src_id_col='CNAME',
                    verbose=False,
                    fixed_grid=False,
                    prune_fields=False):
    """Run define_fields once and compute metrics for this offset."""
    Result = namedtuple("Result", ["offset", "n_science", "n_fields",
                        "coverage_frac", "field_stars", "lut_df"])

    if fixed_grid:
        if verbose:
            print(f"Using fixed grid with offset={offset}")
        fixed_grid = True
    else:
        if verbose:
            print(f"Using adaptive grid with offset={offset}")
        fixed_grid = False

    # suppress prints from define_fields
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
            ra_col=ra_col, de_col=de_col,
            mag_col=mag_col, src_id_col=src_id_col,
            verbose=False,                 # if the function supports it
            fixed_grid=fixed_grid,
            prune_fields=prune_fields
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


def _run_define_fields(cat_comb_df, fs, circles_df=None, suppress_io=True, **kwargs):
    """
    Helper: run fs.define_fields with optional circles_df passed through kwargs.
    Returns (field_stars, lut_df).
    """
    # Always forward circles_df if provided (this is how we control which fieldIds are enabled)
    if circles_df is not None:
        kwargs = dict(kwargs)
        kwargs["circles_df"] = circles_df
        # If your pipeline expects this flag to bypass bbox trimming, keep it here:
        #kwargs.setdefault("filter_circles_to_science_bbox", False)

    # Default to quiet runs unless caller overrides
    kwargs.setdefault("verbose", False)

    if suppress_io:
        with open(os.devnull, 'w') as devnull, \
             contextlib.redirect_stdout(devnull), \
             contextlib.redirect_stderr(devnull):
            field_stars, lut_df = fs.define_fields(cat_comb_df, **kwargs)
    else:
        field_stars, lut_df = fs.define_fields(cat_comb_df, **kwargs)

    return field_stars, lut_df


def evaluate_optimal_circles(cat_comb_df, circles_df, fs, **kwargs):
    """
    Evaluate fs.define_fields() for a given circles_df subset and return coverage metrics.
    """
    Result = namedtuple("Result", ["circles_df", "n_science", "n_fields",
                               "coverage_frac", "field_stars", "lut_df"])

    field_stars, lut_df = _run_define_fields(cat_comb_df, fs, circles_df=circles_df, **kwargs)

    # Identify science vs. field rows (sourceFlag == 9 denotes guide/field rows)
    n_science = field_stars[~field_stars["sourceFlag"].isin([9])].shape[0]
    n_fields = field_stars[field_stars["sourceFlag"].isin([9])].shape[0]

    total_candidates = cat_comb_df[cat_comb_df['sourceFlag'] != 9].shape[0]
    coverage_frac = n_science / max(1, total_candidates)

    return Result(circles_df, n_science, n_fields, coverage_frac, field_stars, lut_df)


def sweep_optimal_circles(cat_comb_df,
                          fs,
                          circles_df=None,
                          guide_flag_col="guide",
                          ra_col="GAIA_RA",
                          dec_col="GAIA_DEC",
                          field_id_col="fieldId",
                          metric="n_science",
                          tol=0,
                          randomize_order=False,
                          **kwargs):
    """
    Greedy backward elimination of fieldIds:
    - Start with all provided fieldIds (from circles_df or from cat_comb_df[guide==1]).
    - At each pass, try removing one fieldId at a time.
    - Permanently drop a fieldId if removal keeps the *maximum* science coverage.
    - Stop when no more fieldIds can be removed without reducing coverage.

    Parameters
    ----------
    cat_comb_df : pd.DataFrame
    fs          : object with define_fields(...) -> (field_stars, lut_df)
    circles_df  : pd.DataFrame or None
        If None, it is built from cat_comb_df rows with guide_flag_col==1,
        using [ra_col, dec_col, field_id_col] columns.
    metric      : {"n_science","coverage_frac"}
        Metric used to decide if removal preserves maximum coverage.
    tol         : int or float
        Allowed tolerance for metric when comparing to baseline (e.g., 0 means must be >= exact).
    randomize_order : bool
        If True, shuffle field removal trial order each pass.

    Returns
    -------
    final_result : Result
        Result for the final kept set of fieldIds.
    removed_field_ids : list
        List of fieldIds that were removed by the sweep.
    attempts_log : pd.DataFrame
        Log of all attempted removals with their resulting metrics and decisions.
    baseline_result : Result
        Result for the initial full set before pruning.
    """
    # --- Build circles_df if not provided ---
    if circles_df is None:
        guides = cat_comb_df[cat_comb_df['guide'] == 1][[ra_col, dec_col, field_id_col]].copy()
        if guides.empty:
            raise ValueError("No guides found to build circles_df (check guide_flag_col and column names).")
        # Minimal circles_df expected by your pipeline:
        circles_df = guides.rename(columns={ra_col: "GAIA_RA", dec_col: "GAIA_DEC", field_id_col: "fieldId"}).copy()

    # --- Establish the baseline with all fields present ---
    baseline_result = evaluate_optimal_circles(cat_comb_df, circles_df, fs, **kwargs)
    baseline_metric = baseline_result.n_science if metric == "n_science" else baseline_result.coverage_frac

    # We will attempt to remove fields while keeping >= baseline_metric - tol
    current_circles = circles_df.copy()
    all_field_ids = list(pd.unique(current_circles["fieldId"]))
    kept_field_ids = set(all_field_ids)
    removed_field_ids = []
    attempts_rows = []

    # Greedy backward elimination
    improved = True
    while improved and len(kept_field_ids) > 0:
        improved = False

        trial_order = list(kept_field_ids)
        if randomize_order:
            rng = np.random.default_rng(42)
            rng.shuffle(trial_order)

        for fid in trial_order:
            # Trial: remove fid
            trial_circles = current_circles[current_circles["fieldId"] != fid]
            trial_res = evaluate_optimal_circles(cat_comb_df, trial_circles, fs, **kwargs)
            trial_metric = trial_res.n_science if metric == "n_science" else trial_res.coverage_frac

            # Decide if we can drop fid without losing target coverage
            keep_ok = (trial_metric >= baseline_metric - tol)

            attempts_rows.append({
                "removed_try": fid,
                "trial_n_science": trial_res.n_science,
                "trial_n_fields": trial_res.n_fields,
                "trial_coverage_frac": trial_res.coverage_frac,
                "decision_keep_removed": bool(keep_ok),
            })

            if keep_ok:
                # Commit the removal
                current_circles = trial_circles
                kept_field_ids.remove(fid)
                removed_field_ids.append(fid)
                improved = True
                # Continue with this new current set (no need to re-check previously tried this pass)
        # end for fid
    # end while improved

    # Final result after all removable fields are dropped
    final_result = evaluate_optimal_circles(cat_comb_df, current_circles, fs, **kwargs)

    attempts_log = pd.DataFrame(attempts_rows, columns=[
        "removed_try", "trial_n_science", "trial_n_fields",
        "trial_coverage_frac", "decision_keep_removed"
    ])

    return final_result, removed_field_ids, attempts_log, baseline_result

def sweep_optimal_circles(df, cat_comb_df, fs, **kwargs):

    guides = df[df['guide'] == 1].copy()

    field_ids = guides['fieldId'].unique()
    #field_ids = [17]
    ids_to_remove = []
    for fid in field_ids:
        guide_subset = guides[~guides['fieldId'].isin([fid])]

        with open(os.devnull, 'w') as devnull, \
            contextlib.redirect_stdout(devnull), \
            contextlib.redirect_stderr(devnull):

            field_stars, _ = fs.define_fields(cat_comb_df, grid_offset=1.615,
                                              central_star_mag_lim=16.,
                                              src_id_col='GAIA_ID',
                                              verbose=False,
                                              fixed_grid=True,
                                              circles_df=guide_subset)

        n_science = field_stars[~field_stars["guide"].isin([1,2])].shape[0]
        # print(fid, n_science)
        n_fields = field_stars[field_stars["guide"].isin([1])].shape[0]
        total_candidates = cat_comb_df[~cat_comb_df['sourceFlag'].isin([9, -1])].shape[0]
        coverage_frac = n_science / max(1, total_candidates)

        # if the id removal leads to the same number of science targets as the
        # total_candidates, we can add it to the ids_to_remove list
        if n_science == total_candidates:
            ids_to_remove.append(fid)
        for fid in ids_to_remove:
            guides = guides[~guides['fieldId'].isin([fid])]

    print(f"Removed {len(ids_to_remove)} fields without affecting the science coverage.")

    return ids_to_remove


def sweep_optimal_circles(df, cat_comb_df, fs, coverage_tol=0, **kwargs):
    """
    Iterative backward elimination:
      - Start with all fieldIds found in df[df['guide']==1].
      - At each step, try removing ONE fieldId.
      - Commit a removal only if coverage (n_science) stays >= current coverage - coverage_tol.
      - After committing, recompute the baseline using the new kept set, then restart.
    Returns the list of removed fieldIds.
    """

    # --- helper: run define_fields quietly on a given guides subset ---
    def run_with_guides(guide_subset):
        with open(os.devnull, 'w') as devnull, \
             contextlib.redirect_stdout(devnull), \
             contextlib.redirect_stderr(devnull):
            field_stars, _ = fs.define_fields(cat_comb_df, circles_df=guide_subset, **kwargs)
        # coverage metric consistent with your code
        n_science = field_stars[~field_stars["guide"].isin([1, 2])].shape[0]
        return n_science

    # Build guides table (don’t mutate the original df)
    guides_all = df[df['guide'] == 1].copy()
    # Ensure one row per fieldId (keep first RA/Dec; adjust if you store multiple)
    if not {'GAIA_RA','GAIA_DEC','fieldId'}.issubset(guides_all.columns):
        raise ValueError("Expected columns GAIA_RA, GAIA_DEC, fieldId in df where guide == 1.")
    guides_all = guides_all[['GAIA_RA','GAIA_DEC','fieldId']].drop_duplicates('fieldId')

    kept_ids = list(guides_all['fieldId'].unique())
    removed_ids = []

    # Baseline coverage with all kept_ids
    current_guides = guides_all[guides_all['fieldId'].isin(kept_ids)]
    current_cov = run_with_guides(current_guides)

    changed = True
    while changed and kept_ids:
        changed = False

        # Iterate over a snapshot of current kept_ids (don’t mutate while iterating)
        for fid in list(kept_ids):
            trial_ids = [x for x in kept_ids if x != fid]
            trial_guides = guides_all[guides_all['fieldId'].isin(trial_ids)]
            trial_cov = run_with_guides(trial_guides)

            # Commit removal only if coverage doesn’t drop (within tolerance)
            if trial_cov >= (current_cov - coverage_tol):
                kept_ids.remove(fid)
                removed_ids.append(fid)
                current_cov = trial_cov  # update baseline to the new state
                changed = True
                break  # restart loop from the new state (crucial for interaction effects)

    print(f"Removed {len(removed_ids)} fields without affecting coverage (tol={coverage_tol}).")
    return removed_ids


if __name__ == "__main__":
    pass
