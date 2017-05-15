#-*- coding: utf-8 -*-

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
parser.add_argument('output', help='Name of the output file')
parser.add_argument('--extinction', help='Output of clusters_extinction (hdf5 file)')
parser.add_argument('--mag', type=str, help='Stars magnitude column', default='modelfit_CModel_mag')
parser.add_argument('--cut', type=float, help='Select the value of the magnitude cut ((g-i) > value)')

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
if args.extinction is not None:
    data.correct_for_extinction(d[cat], te=d['extinction'], ifilt='i_new', mag=args.mag)
    mag = args.mag + '_extcorr'
else:
    mag = args.mag

mag_err = args.mag + 'Sigma'

# Filters to select stars
filt = d['deepCoadd_meas']['base_ClassificationExtendedness_flag'] == 0
filt &= d['deepCoadd_meas']['detect_isPrimary'] == 1
filt &= d[cat]['modelfit_CModel_flag'] == 0
filt &= d[cat]['modelfit_CModel_flux'] > 0

# Mask to remove bad magnitudes
filt &= (d[cat]['modelfit_CModel_flux'] / d[cat]['modelfit_CModel_fluxSigma']) > 10

# Mask to select stars
filtS = d['deepCoadd_meas']['base_ClassificationExtendedness_value'] < 0.5

s = d[cat][filt & filtS].group_by(oid)
f = (s.groups.indices[1:] - s.groups.indices[:-1]) == nfilters
star = s.groups[f]

# Set bad magnitudes for big macs
#for n in range(len(star)):
#	if (star['modelfit_CModel_flux'][n] / star['modelfit_CModel_fluxSigma'][n]) <= 10:
#		star['modelfit_CModel_mag_extcorr'][n] = 99

# Magnitude cut
def get_mag(g, f):
    return g[mag][g['filter'] == f]
def set_mag(g, f, nv):
    g[mag][g['filter'] == f] = nv

if args.cut is not None:
        print len(star.groups)
        for group in star.groups:
            if (get_mag(group, 'g') - get_mag(group, 'i')) > args.cut:
                set_mag(group, 'u', 99)

# Classify magnitudes by filters
filters = set(star['filter'])

mags = {f: star[star['filter'] == f][mag] for f in filters}
mags_err = {f: star[star['filter'] == f][mag_err] for f in filters}

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
prihdr['COMMENT'] = ""
prihdu = fits.PrimaryHDU(header=prihdr)

# Write fits file
thdulist = fits.HDUList([prihdu, tbhdu])
thdulist.writeto(args.output + '.fits')

# Check fits file content
hdu_list = fits.open(args.output + '.fits')
hdu_list.info()

# Check first extension
g = getdata(args.output + '.fits', 1)
t = Table(g)
print t
