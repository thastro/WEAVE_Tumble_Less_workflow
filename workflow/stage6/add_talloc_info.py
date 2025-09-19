import os
import numpy as np
import argparse
from astropy.io import fits
#==============================================================================
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Add TALLOC code to header')
    parser.add_argument('fits_catalogue',  help='name of the FITS catalogue')
    parser.add_argument('talloc_code',     help='TALLOC code')
    args = parser.parse_args()

    fits.setval(args.fits_catalogue, 'TACALLOC', value=args.talloc_code, ext=0)

    exit()
