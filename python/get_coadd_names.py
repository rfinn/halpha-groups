#!/usr/bin/env python
"""
run from directory that contains the coadds

this will create the list of coadded images that I can then use 
as input for other files


"""
import glob
import argparse
import os

# get list of r-band coadded images
parser = argparse.ArgumentParser(description ='Get the list of coadded images.  Run from the coadd directory.')
args = parser.parse_args()



a = glob.glob('UAT*HDI*-R.fits')
b = glob.glob('UAT*MOS*-R.fits')

#d = glob.glob('UAT*BOK*-r.fits')
#######################################################
# updating to use the shifted r-band images
#######################################################    
rfiles = a + b 

print(f"number of targets = {len(rfiles)}")

# write out as a csv file
outfile = open('uat-coadds.csv','w')
outfile4 = open('uat-coadds-test.csv','w')
outfile2 = open('uat-coadds-fullpath.txt','w')
outfile3 = open('uat-coadds-fullpath-test.txt','w')
coadd_dir = os.getcwd()
for i in range(len(rfiles)):
    #basname = rfiles[i].replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    outfile.write(f"{rfiles[i]}\n")
    outfile2.write(f"{coadd_dir}/{rfiles[i]}\n")
    if i < 2:
        outfile3.write(f"{coadd_dir}/{rfiles[i]}\n")
        outfile4.write(f"{rfiles[i]}\n")
outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()



######################################################
# REPEAT FOR THE HALPHA IMAGES
######################################################
a = glob.glob('UAT*ha4.fits')
b = glob.glob('UAT*ha8.fits')
#b = glob.glob('UAT*HDI*-r.fits')
c = glob.glob('UAT*ha12.fits')
d = glob.glob('UAT*ha16.fits')

hfiles = a + b + c + d

hfiles.sort()
print(f"number of targets = {len(hfiles)}")

# write out as a csv file
outfile = open('uat-coadds-ha.csv','w')
outfile4 = open('uat-coadds-ha-test.csv','w')
outfile2 = open('uat-coadds-ha-fullpath.txt','w')
outfile3 = open('uat-coadds-ha-fullpath-test.txt','w')
for i in range(len(hfiles)):
    #basname = rfiles[i].replace("-r-shifted.fits","").replace("-r.fits","").replace("-R.fits","")
    outfile.write(f"{hfiles[i]}\n")
    outfile2.write(f"{coadd_dir}/{hfiles[i]}\n")
    if i < 2:
        outfile3.write(f"{coadd_dir}/{hfiles[i]}\n")
        outfile4.write(f"{hfiles[i]}\n")
outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()
