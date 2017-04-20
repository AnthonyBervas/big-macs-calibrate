#-*- coding: utf-8 -*-s
import numpy as np
from clusters import data
from astropy.io import fits
from astropy.io.fits import getdata
from astropy.table import Table
from argparse import ArgumentParser

description = """Convert hdf5 file to fits file for BIGMACS."""
prog = "hdf5_to_fits_BIGMACS.py"
usage = """%s [options] input""" % prog

parser = ArgumentParser(prog=prog, usage=usage, description=description)
parser.add_argument('input', help='hdf5 file to convert in fits file')

args = parser.parse_args()

# Read hdf5 input file
d = data.read_hdf5(args.input)

# Get number of filters
nfilters = len(d['deepCoadd_meas'].group_by('id').groups[0])
print("We have data for %d filters in this dataset" % nfilters)

cat = 'deepCoadd_meas'
oid = 'id'

print("Number of %s sources before cuts : %d" % (cat, len(d[cat]) / nfilters))

data.correct_for_extinction(d[cat], d['extinction'], ifilt='i_old')

# Get data only for stars
mag_type = 'modelfit_CModel_mag_extcorr'

filt = d['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0
filt &= d['deepCoadd_meas']['detect_isPrimary'] == 1
filt &= d[cat]['modelfit_CModel_flag'] == 0
filt &= d[cat]['modelfit_CModel_flux'] > 0

filt &= (d[cat]['modelfit_CModel_flux'] / d[cat]
         ['modelfit_CModel_fluxSigma']) > 10

filtS = d['deepCoadd_meas']['base_ClassificationExtendedness_value'] < 0.5

s = d[cat][filt & filtS].group_by(oid)
f = (s.groups.indices[1:] - s.groups.indices[:-1]) == nfilters
star = s.groups[f]

filters = set(star['filter'])

mags = {f: star[star['filter'] == f][mag_type] for f in filters}

mags_err = {f: star[star['filter'] == f]['modelfit_CModel_magSigma'] for f in filters}

# Get RA & DEC tables
X_WORLD = star['coord_ra_deg'][star['filter'] == 'r']
Y_WORLD = star['coord_dec_deg'][star['filter'] == 'r']

# Create fits file with magnitudes and coordinates
col1 = fits.Column(name='X_WORLD', format='E', array=X_WORLD)
col2 = fits.Column(name='Y_WORLD', format='E', array=Y_WORLD)
cols = fits.ColDefs([col1, col2])

for f in filters:
	col = fits.Column(name='mag_%s'%f, format='E', array=mags[f])
	col_err = fits.Column(name='mag_%s_err'%f, format='E', array=mags_err[f])
	cols.add_col(col)
	cols.add_col(col_err)

tbhdu = fits.BinTableHDU.from_columns(cols)

# Header of fits file
prihdr = fits.Header()
prihdr['COMMENT'] = "RA,DEC / Stars Magnitude / Stars Magnitude Errors."
prihdu = fits.PrimaryHDU(header=prihdr)

# Header + First extension
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('mag.fits')

# Check fits file content
hdu_list = fits.open('mag.fits')
hdu_list.info()

# Check first extension
g = getdata('mag.fits', 1)
t = Table(g)
print t