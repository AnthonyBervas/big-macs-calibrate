#-*- coding: utf-8 -*-
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
parser.add_argument('hold', help='Select the hold filter during BIGMACS ZPs fit')

args = parser.parse_args()

# Read hdf5 input file
d = data.read_hdf5(args.input)

# Get number of filters
nfilters = len(d['deepCoadd_meas'].group_by('id').groups[0])
print("We have data for %d filters in this dataset" % nfilters)

cat = 'deepCoadd_meas'
oid = 'id'

print("Number of %s sources before cuts : %d" % (cat, len(d[cat]) / nfilters))

# Correction of magnitudes for extinction
data.correct_for_extinction(d[cat], d['extinction'], ifilt='i_new')

# Get data only for stars et select the extinction for the hold filter.
mag_type_hold = 'modelfit_CModel_mag_extcorr'
mag_type = 'modelfit_CModel_mag'

filt = d['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0
filt &= d['deepCoadd_meas']['detect_isPrimary'] == 1
filt &= d[cat]['modelfit_CModel_flag'] == 0
filt &= d[cat]['modelfit_CModel_flux'] > 0

# Mask to remove bad magnitudes
# filt &= (d[cat]['modelfit_CModel_flux'] / d[cat]['modelfit_CModel_fluxSigma']) > 10

# Mask to select stars
filtS = d['deepCoadd_meas']['base_ClassificationExtendedness_value'] < 0.5

s = d[cat][filt & filtS].group_by(oid)
f = (s.groups.indices[1:] - s.groups.indices[:-1]) == nfilters
star = s.groups[f]

# Set bad magnitudes for big macs
for n in range(len(star)):
	if (star['modelfit_CModel_flux'][n] / star['modelfit_CModel_fluxSigma'][n]) <= 10:
		star['modelfit_CModel_mag'][n] = 99

star_p = star

# Classify magnitudes by filters
filters = set(star_p['filter'])

mags = {f: star_p[star_p['filter'] == f][mag_type] for f in filters}
mags_err = {f: star_p[star_p['filter'] == f]['modelfit_CModel_magSigma'] for f in filters}
mags_hold = {f: star_p[star_p['filter'] == f][mag_type_hold] for f in filters}
mags_err_hold = {f: star_p[star_p['filter'] == f]['modelfit_CModel_magSigma'] for f in filters}

# Get RA & DEC tables
X_WORLD = star_p['coord_ra_deg'][star_p['filter'] == args.hold]
Y_WORLD = star_p['coord_dec_deg'][star_p['filter'] == args.hold]

# Create fits file with magnitudes and coordinates
col1 = fits.Column(name='X_WORLD', format='E', array=X_WORLD)
col2 = fits.Column(name='Y_WORLD', format='E', array=Y_WORLD)
cols = fits.ColDefs([col1, col2])

for f in filters:
	if f != args.hold:
		col = fits.Column(name='mag_%s'%f, format='E', array=mags[f])
		col_err = fits.Column(name='mag_%s_err'%f, format='E', array=mags_err[f])
		cols.add_col(col)
		cols.add_col(col_err)
	else:
		col3 = fits.Column(name='mag_r', format='E', array=mags_hold[args.hold])
		col2_err = fits.Column(name='mag_r_err', format='E', array=mags_err_hold[args.hold])
		cols.add_col(col3)
		cols.add_col(col2_err)

tbhdu = fits.BinTableHDU.from_columns(cols)

# Header of fits file
prihdr = fits.Header()
prihdr['COMMENT'] = ""
prihdu = fits.PrimaryHDU(header=prihdr)

# Write fits file
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto('mag.fits')


# Check fits file content
hdu_list = fits.open('mag.fits')
hdu_list.info()

# Check first extension
# g = getdata('mag.fits', 1)
# t = Table(g)
# print t.more()
