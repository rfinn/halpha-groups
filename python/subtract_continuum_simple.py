#!/usr/bin/env python


'''

GOAL:

To create a continuum-subtracted image given
(1) an R-band and Halpha image, and
(2) the scale factor to apply to the R-band image.

PROCEDURE:

- scale R band by a factor
- subtract from Halpha
- When running in interactive mode, you can iterate until the you are happy
- save resultant continuum-subtracted image if desired
  
EXAMPLES:

(1) to run within ipython type:

%run ~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale .055 --mosaic

The mosaic flag increases the size of the images that are displayed so you can see the results better.

(2) to run from the command line type:

~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale .055

(3) to run on cutout images rather than mosaics:

 ~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364

INPUT/OUPUT:

REQUIRED MODULES:

astropy
argparse
matplotlib
scipy

EXTRA NOTES:

WRITTEN BY:
   Rose Finn

'''
from astropy.io import fits
import argparse
from argparse import RawDescriptionHelpFormatter
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import sys


def plot_results(ha, r, cs, vmin = .5, vmax = 98):
        
    v1,v2=scoreatpercentile(ha,[vmin,vmax])#.5,99
    
    fig = plt.figure(1,figsize=figure_size)
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    #Halpha plus continuum
    plt.subplot(1,3,1)
    plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
    plt.title('Halpha + cont')
    #plt.gca().set_yticks(())
    plt.xlabel('NSA ID '+id,fontsize=14)
    #R
    plt.subplot(1,3,2)
    plt.imshow(r,cmap='gray_r',vmin=v1/scale,vmax=v2/scale,origin='lower')
    plt.title('R')
    plt.gca().set_yticks(())
    #Continuum subtracted image
    plt.subplot(1,3,3)
    v3,v4=scoreatpercentile(cs,[vmin,vmax])
    plt.imshow(cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
    plt.gca().set_yticks(())
    plt.title('contsub, scale = %4.4f'%(scale))
    plt.show(block=False)
    return fig

def write_cs_image(cs,ha_header,outimage,fig=None):
    """
    PARAMS:
    - cs: continuum subtracted image
    - ha_header: header from the Halpha image, will be used to save the CS image
    - outimage: name for the CS image
    
    OPTIONAL:
    - fig: figure instance reference, use this to save a png version of the figure.
    """
    newfile = fits.PrimaryHDU()
    newfile.data = cs
    newfile.header = ha_header
    fits.writeto(outimage, newfile.data, header = newfile.header, overwrite=True)
    output = outimage.split('.fits')
    if fig is not None:
        fig.savefig(output[0]+'.png')


parser = argparse.ArgumentParser(description ='This program subtracts scaled R-band image from Halpha using a constant scale factor.\n \nTo subtract mosaics:\n~/github/halpha-groups/python/subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale 0.0445 --mosaic \n\nTo subtract cutouts:\n~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364', formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('--pointing',dest = 'pointing', default = None, help = 'Cluster prefix of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
parser.add_argument('--id',dest = 'id', default = None, help = 'NSAID of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image.  Use this if you are subtracting mosaic images rather than cutouts.')
parser.add_argument('--ha', dest = 'ha', default = None, help = 'Halpha image.  Use this if you are subtracting mosaic images rather than cutouts.')
parser.add_argument('--scale', dest = 'scale', default = 0.0445, help = 'factor to scale R-band image by before subtracting from Halpha image')
parser.add_argument('--zp', dest = 'zp', default = False, action = 'store_true', help = 'Use the PHOTZP header keywords to determine the scale factor for continuum subtraction.')
parser.add_argument('--interactive', dest = 'interactive', default = False, action = 'store_true', help = 'Set this to run the continuum subtraction interactively. Default is False.')
parser.add_argument('--mosaic', dest = 'mos', default = False, action = 'store_true', help = 'set this if subtracting mosaic images rather than cutouts.  It will make the figure bigger so you can see the mosaics better. Default is False.')

args = parser.parse_args()

figure_size=(10,4)
if args.mos:
    figure_size=(14,6)
    
if args.id:
    rimage = args.pointing+'-'+args.id+'-R.fits'
    haimage = args.pointing+'-'+args.id+'-Ha.fits'
    id = args.id
else:
    rimage = args.r
    haimage = args.ha
    t = args.r.split('-')
    id = t[0]
r,r_header = fits.getdata(rimage,header=True)
ha,ha_header = fits.getdata(haimage,header=True)


if args.zp:
    outname_suffix = 'CS_ZP.fits'
else:
    outname_suffix = 'CS.fits'

# define the output continuum-subtracted image name
outimage = haimage.replace('.fits','-CS-ZP.fits')


if args.zp:
    # get PHOTZP from image headers
    try:
        RZP = r_header['PHOTZP']
    except KeyError:
        print("Error: no PHOTZP found in r-band image header")
        print("make sure you ran getzp.py on the r-band image.  For example:")
        print("\t python ~/github/HalphaImaging/python3/getzp.py  --image {rimage} --instrument h --filter R")
        sys.exit()
    try:
        HZP = ha_header['PHOTZP']
    except KeyError:
        print("Error: no PHOTZP found in Halpha image header")
        print("make sure you ran getzp.py on the Halpha image.  For example:")
        print("\t python ~/github/HalphaImaging/python3/getzp.py  --image {rimage} --instrument h --filter ha")
        sys.exit()

    # calculate the scale ratio for the r-band image from the difference in PHOTZP
    # m2 - m1 = 2.5 log10 (f1/f2)
    # let 1 = Halpha, 2 = r, then f1/f2 is the scale factor such that
    # CS = Halpha - scale * R

    scale = 10.**((HZP-RZP)/2.5)
else:
    scale = float(args.scale)

# get the CS image
cs = ha - scale*r

if not args.interactive:
    # plot the results
    fig = plot_results(ha, r, cs)
    
    # write out results
    write_cs_image(cs,ha_header,outimage,fig=fig)    

else:  
    vmin = .5
    vmax = 98
    adjust_scale = True
    while adjust_scale:
        plt.close('all')
        fig = plot_results(ha, r, cs, vmin = vmin, vmax = vmax)
        #cs = ha  - scale*r
        t=raw_input('enter new value for scale factor;\n   a to adjust contrast percentages(e.g. [.5,98], default is [50,99]). Enter min max.;\n   w to write output and quit;\n   q to quit without saving\n') #TASH after lunch, work on writing code for a
        try:
            scale=float(t)
        except ValueError:
            if t.find('a') > -1:
                t=raw_input('enter new value for contrast scale;\n   example: .5 99\n' )
                try:
                    vmin,vmax = t.split()
                    vmin = float(vmin)
                    vmax = float(vmax)
                except:
                    print("Didn't understand input entry, try again")
                    continue
                adjust_scale = False
                if t.find('q') > -1:
                    break
       
                if t.find('w') > -1:
                    write_cs_image(cs,ha_header,outimage,fig=fig)    
