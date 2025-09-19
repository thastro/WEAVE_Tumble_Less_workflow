from .check_populated_ifu_fits_cat import check_populated_ifu_fits_cat
import logging
try:
    from .add_photometry_from_ps import add_photometry_from_ps
    from .add_photometry_from_sdss import add_photometry_from_sdss
    from .mark_star_fibres import mark_star_fibres
    from .populate_column_based_on_targid import populate_column_based_on_targid
except:
    logging.warning("Cannot import stag6 modules")
