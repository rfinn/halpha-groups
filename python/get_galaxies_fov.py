#!/usr/bin/env python

"""
GOAL 
* get a list of all galaxies that fall within the UAT coadds

PROCEDURE
make a table that's row-matched to AGC with columns
agc number
in FOV
in filter
image name
image2 - if galaxy is in more than one coadd
instrument
halpha filter
r-band filter

"""

import build_web_common as buildweb

import glob
import os
import numpy as np

from astropy.io import fits
from astropy.table import Table, Column
from astropy import units as u


import argparse


##############################################################
# dictionary of Halpha filters 
##############################################################
lmin={'ha4':6573., 'ha8':6606., 'ha12':6650., 'ha16':6682., 'INT197':6540.5}
lmax={'ha4':6669., 'ha8':6703., 'ha12':6747., 'ha16':6779., 'INT197':6615.5}



def parse_uat_coadd_name(imagename):
    '''
    example name: UAT-144.582+17.103-HDI-20160414-NRGb049-h01-R.fits
    '''

    t = imagename.split('-')
    #print("number of fields from splitting imagename = ",len(t))
    if len(t) == 7:
        instrument = t[2]
    else: # negative declination
        instrument = t[3]

    thisfilter = t[-1].split('.fits')[0]

    return instrument, thisfilter


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description ='Get galaxies in FOV of the UAT coadded images.  Expecting the agc catalog to be in {homedir}/research/AGC/agcnorthminus1.full200617.fits', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--coaddir', dest = 'coaddir', default = '.', help = 'directory for coadds.  default is current directory')
    parser.add_argument('--outdir', dest = 'outdir', default =  '.', help = 'output directory for table with galaxies in FOV.  default is current directory')
    #parser.add_argument('--agcpath', dest = 'agcpath', default = None, help = 'photometric catalog to use for bootrapping photometry')    
    args = parser.parse_args()

    homedir = os.getenv("HOME")
    
    #############################################################
    # read in agc
    #############################################################
    agcfile = os.path.join(homedir,'research/AGC','agcnorthminus1.full200617.fits')
    agc = Table.read(agcfile)

    # get super velocity for agc, which prioritizes vopt if available
    agcsuperv = agc['vopt']

    novoptflag = agcsuperv == 0
    agcsuperv[novoptflag] = agc['v21'][novoptflag]

    #############################################################
    # create an output table with same length as agc
    #############################################################

    infov = Column(np.zeros(len(agc),'bool'), name='infov')
    infilter = Column(np.zeros(len(agc),'bool'), name='infilter')
    imagename = Column(np.zeros(len(agc),dtype='S53'), name='imagename')
    imagename2 = Column(np.zeros_like(imagename), name='imagename2')
    instrument = Column(np.zeros(len(agc),dtype='S3'), name='instrument')
    hafilter = Column(np.zeros(len(agc),dtype='S4'), name='hafilter')
    rfilter = Column(np.zeros(len(agc),dtype='S4'), name='rfilter')
    superv = Column(agcsuperv*u.km/u.s, name='supervel')#,units=u.km/u.s)    
    groupcat = Table([agc['AGCnr'],agc['radeg'],agc['decdeg'], superv, infov, infilter, instrument, rfilter, hafilter,imagename, imagename2])

    #############################################################
    # get list of R-band images
    #############################################################
    Rfilelist = glob.glob(os.path.join(args.coaddir,'UAT*R.fits'))
    rfilelist = glob.glob(os.path.join(args.coaddir,'UAT*r.fits'))                             
    allfilelist = Rfilelist + rfilelist
    
    #############################################################
    # loop over images
    #############################################################
    for fname in allfilelist:
        print(f"Working on {fname}")
    
        #############################################################
        # get galaxies in footprint
        #############################################################
        imx, imy, keepflag = buildweb.get_galaxies_fov(fname, agc['radeg'], agc['decdeg'])

        #############################################################
        # set in fov flag
        #############################################################
        groupcat['infov'][keepflag] = True



        #############################################################
        # save parent image name, halpha filter
        #############################################################

        groupcat['imagename'][keepflag] = os.path.basename(fname)

        # TODO - if imagename is already populated, save as imagename2
        # but do we also then need to save instrument2, rfilter2, hafilter2, in case they are different?
        #############################################################
        # save parent image name, halpha filter
        #############################################################
        
        instrument, rfilter = parse_uat_coadd_name(fname)
        #print(rfilter, instrument)
        groupcat['rfilter'][keepflag] = rfilter
        groupcat['instrument'][keepflag] = instrument

    
        rheader = fits.getheader(fname)
        try:
            haimage = rheader['HAIMAGE']
        except KeyError:
            print("WARNING: no HAIMAGE in header of ", fname)
            print("continuing to next image...")
            continue
                      
        instrument, hafilter = parse_uat_coadd_name(rheader['HAIMAGE'])
        groupcat['hafilter'][keepflag] = hafilter

        #############################################################
        # check redshifts against filter throughput
        #############################################################
        zmax=((lmax[hafilter])/6563.)-1
        zmin=((lmin[hafilter])/6563.)-1
        vflag = (agcsuperv/3.e5 > zmin) & (agcsuperv/3.e5 < zmax)

        #############################################################
        # set in filter flag
        #############################################################
        groupcat['infilter'][keepflag & vflag] = True
    
    # write out file
    groupcat.write('uat_groups_halpha_sample.fits', format='fits', overwrite=True)
