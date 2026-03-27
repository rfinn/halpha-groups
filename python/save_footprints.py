#!/usr/bin/env python3

from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table


def find_fits_files():
    """Find matching FITS files in the current directory."""
    files = set()
    files.update(Path(".").glob("UAT*R.fits"))
    files.update(Path(".").glob("UAT*r.fits"))
    return sorted(files)


def get_footprint_from_header(filename):
    """Return footprint as (4,2) array [RA, Dec]."""
    with fits.open(filename) as hdul:
        header = hdul[0].header
        wcs = WCS(header)
        footprint = wcs.calc_footprint()

    return footprint  # shape (4,2)


def main():
    files = find_fits_files()

    if not files:
        print("No matching FITS files found.")
        return

    rows = []

    for f in files:
        try:
            fp = get_footprint_from_header(f)

            # Ensure shape is (4,2)
            fp = np.asarray(fp).reshape(4, 2)

            row = {
                "image": f.name,

                # --- scalar columns ---
                "ra1": fp[0, 0], "dec1": fp[0, 1],
                "ra2": fp[1, 0], "dec2": fp[1, 1],
                "ra3": fp[2, 0], "dec3": fp[2, 1],
                "ra4": fp[3, 0], "dec4": fp[3, 1],

                # --- array column ---
                "footprint": fp  # shape (4,2)
            }

            rows.append(row)

        except Exception as e:
            print(f"Skipping {f.name}: {e}")

    if not rows:
        print("No valid footprints were generated.")
        return

    table = Table(rows=rows)

    # Optional: explicitly define dtype for clarity
    table["footprint"] = np.array(table["footprint"], dtype=float)

    # --- write outputs ---
    table.write("image_footprints.fits", overwrite=True)
    table.write("image_footprints.ecsv", overwrite=True)

    print(f"Wrote {len(table)} rows to:")
    print("  image_footprints.fits")
    print("  image_footprints.ecsv")


if __name__ == "__main__":
    main()
