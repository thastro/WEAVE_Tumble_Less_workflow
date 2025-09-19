from . import utils

from . import stage1
from . import stage2
from . import stage3
# Stage 4 is configure
from . import stage5
from . import stage6
from . import stage7
# Stage 8 is WASP
from . import stage9
from . import stage10


from .stage1 import create_ifu_driver_cat as stage1_main
from .stage1 import check_ifu_driver_cat as stage1_check
from .stage1 import plot_ifu_driver_cat as stage1_plot
from .stage2 import create_xml_files as stage2_main
from .stage3 import add_guide_and_calib_stars as stage3_main
# Stage 4 is configure
from .stage5 import create_ifu_fits_cat as stage5_main
import logging
try:
    from .stage6 import check_populated_ifu_fits_cat as stage6_check
    from .stage6 import add_photometry_from_ps as stage6_ps
    from .stage6 import add_photometry_from_sdss as stage6_sdss
    from .stage6 import mark_star_fibres as stage6_star
    from .stage6 import populate_column_based_on_targid as stage6_pop_on_targid
except:
    logging.warning("Cannot import stage6 modules")
from .stage7 import create_combo_fits_cat as stage7_main
# Stage 8 is WASP
from .stage9 import fill_xmls_with_fits_info as stage9_main
from .stage10 import plot_data_from_xml as stage10_plot

