import os, re, sys
import matplotlib.pyplot as plt
import numpy as np
import avro.schema
from avro.datafile import DataFileReader, DataFileWriter
from avro.io import DatumReader, DatumWriter
from astropy.io import fits

# projection transformation order
order = 2
deg2rad = np.pi / 180.
rad2deg = 180. / np.pi

# tangent projection coordinates
def xieta(x, y, PV): # all in degrees
  r = np.sqrt(x**2 + y**2)
  xicomp = PV[0, 0] + PV[0, 1] * x + PV[0, 2] * y + PV[0, 3] * r + PV[0, 4] * x**2 + PV[0, 5] * x * y + PV[0, 6] * y**2 + PV[0, 7] * x**3 + PV[0, 8] * x**2 * y + PV[0, 9] * x * y**2 + PV[0, 10] * y**3
  etacomp = PV[1, 0] + PV[1, 1] * y + PV[1, 2] * x + PV[1, 3] * r + PV[1, 4] * y**2 + PV[1, 5] * y * x + PV[1, 6] * x**2 + PV[1, 7] * y**3 + PV[1, 8] * y**2 * x + PV[1, 9] * y * x**2 + PV[1, 10] * x**3
  return (xicomp, etacomp)

# RA DEC given pixel coordinates
def RADEC(i, j, CD11, CD12, CD21, CD22, CRPIX1, CRPIX2, CRVAL1, CRVAL2, PV):

  # i, j to x, y
  x = CD11 * (i - CRPIX1) + CD12 * (j - CRPIX2) # deg 
  y = CD21 * (i - CRPIX1) + CD22 * (j - CRPIX2) # deg

  # x, y to xi, eta
  (xi, eta) = xieta(x, y, PV)

  # xi, eta to RA, DEC
  num1 = (xi * deg2rad) / np.cos(CRVAL2 * deg2rad) # rad
  den1 = 1. - (eta * deg2rad) * np.tan(CRVAL2 * deg2rad) # rad
  alphap = np.arctan2(num1, den1) # rad
  RA  = CRVAL1 + alphap * rad2deg # deg
  num2 = (eta * deg2rad + np.tan(CRVAL2 * deg2rad)) * np.cos(alphap) # rad
  DEC = np.arctan2(num2, den1) * rad2deg # deg

  return (RA / 15., DEC) # hr deg

# apply previous transformation
def applytransformation(order, x1, y1, sol):

    # this is slow, but I prefer fewer bugs than speed at the moment...

    x1t = sol[0] + sol[2] * x1 + sol[3] * y1
    y1t = sol[1] + sol[4] * x1 + sol[5] * y1
    if order > 1:
        x1t = x1t + sol[6] * x1 * x1 + sol[7] * x1 * y1 + sol[8] * y1 * y1
        y1t = y1t + sol[9] * x1 * x1 + sol[10] * x1 * y1 + sol[11] * y1 * y1
    if order > 2:
        x1t = x1t + sol[12] * x1 * x1 * x1 + sol[13] * x1 * x1 * y1 + sol[14] * x1 * y1 * y1 + sol[15] * y1 * y1 * y1
        y1t = y1t + sol[16] * x1 * x1 * x1 + sol[17] * x1 * x1 * y1 + sol[18] * x1 * y1 * y1 + sol[19] * y1 * y1 * y1

    return (x1t, y1t)


SHAREDDIR = "/home/apps/astro/data/SHARED"
ARCHIVEDIR = "/home/apps/astro/data/ARCHIVE"

# check fields

fields = os.listdir(SHAREDDIR)
timeseries = []
idsource = 0
npsf = 21
npsf2 = npsf**2

filtername = {}
airmass = {}
exptime = {}

# astrometric variables
NAXIS = np.zeros((2, 2))
CD = np.zeros((2, 2, 2))
nPV1 = 2
nPV2 = 11
PV = np.zeros((2, nPV1, nPV2))
CRVAL = np.zeros((2, 2))
CRPIX = np.zeros((2, 2))

schema = avro.schema.parse(open("hits.avsc", "rb").read())
writer = DataFileWriter(open("hits-demo.avro", "wb"), DatumWriter(), schema)

for field in fields:
    print(field)

    ccds = os.listdir("%s/%s" % (SHAREDDIR, field))
    
    for ccd in ccds:
        
        print(ccd)

        candidatesdir = "%s/%s/%s/CANDIDATES" % (SHAREDDIR, field, ccd)
        calibrationsdir = "%s/%s/%s/CALIBRATIONS" % (SHAREDDIR, field, ccd)

        if os.path.exists(candidatesdir):

            cands = os.listdir(candidatesdir)

            for cand in cands:
                
                if cand[:4] == "time":
                    
                    print(cand)

                    timeseries.append("%s/%s" % (candidatesdir, cand))
                    idsource += 1
                    data = np.load("%s/%s" % (candidatesdir, cand))
                    i, j, ncoincidence, npossible = data[:4]
                    ncoincidence = int(ncoincidence) + 1
                    MJDs     = data[4: 4 + ncoincidence]
                    fluxes   = data[4 + ncoincidence: 4 + 2 * ncoincidence]
                    e_fluxes = data[4 + 2 * ncoincidence: 4 + 3 * ncoincidence]
                    mags     = data[4 + 3 * ncoincidence: 4 + 4 * ncoincidence]
                    e1_mags  = data[4 + 4 * ncoincidence: 4 + 5 * ncoincidence]
                    e2_mags  = data[4 + 5 * ncoincidence: 4 + 6 * ncoincidence]
                    probs    = data[4 + 6 * ncoincidence: 4 + 7 * ncoincidence]
                    labels   = data[4 + 7 * ncoincidence: 4 + 8 * ncoincidence]
                    ipixs    = data[4 + 8 * ncoincidence: 4 + 9 * ncoincidence]
                    jpixs    = data[4 + 9 * ncoincidence: 4 + 10 * ncoincidence]
                    ipixs = np.array(ipixs, dtype = int)
                    jpixs = np.array(jpixs, dtype = int)
                    
                    MJDref = ("%14.8f" % MJDs[0]).rstrip("0")
                    if not MJDref in filtername.keys():
                        fitsname = "%s/%s/%s/%s/%s_%s_%s_image_crblaster.fits" % (ARCHIVEDIR, field, MJDref, ccd, field, ccd, MJDref)
                        header = fits.open(fitsname)[0].header
                        filtername[MJDref] = header["FILTER"][0]
                        airmass[MJDref] = float(header["AIRMASS"])
                        exptime[MJDref] = float(header["EXPTIME"])                                
                                
                    for icand in range(len(MJDs))[1:]:
                        
                        sci, ref = re.findall("(\d+).*-(\d+).*", labels[icand])[0]
                        sci = int(sci)
                        ref = int(ref)
                        candfile = "%s/cand_%s_%s_%s_grid%02i_lanczos2.npy" % (candidatesdir, field, ccd, labels[icand], ref)
                        psffile = "%s/psf_%s_%s_%s_grid%02i_lanczos2.npy" % (calibrationsdir, field, ccd, labels[icand], ref)

                        if not os.path.exists(candfile) or not os.path.exists(psffile):
                            print("WARNING: %s or %s does not exist" % candfile, psffile) 
                            continue

                        # get science fits header data
                        MJDsci = ("%14.8f" % float(MJDs[icand])).rstrip("0")
                        if not MJDsci in airmass.keys():
                            fitsname = "%s/%s/%s/%s/%s_%s_%s_image_crblaster.fits" % (ARCHIVEDIR, field, MJDsci, ccd, field, ccd, MJDsci)
                            header = fits.open(fitsname)[0].header
                            airmass[MJDsci] = float(header["AIRMASS"])
                            exptime[MJDsci] = float(header["EXPTIME"])

                        # load candidate information
                        canddata = np.load(candfile)
                        psf = np.array(np.load(psffile).flatten(), dtype = float) # need to do this conversion for AVRO to ingest
                        mask = (canddata[:, 0] ==  ipixs[icand]) & (canddata[:, 1] == jpixs[icand])
                        imSNR = canddata[mask, -4 * npsf2: -3 * npsf2][0]
                        im1   = canddata[mask, -3 * npsf2: -2 * npsf2][0]
                        im2   = canddata[mask, -2 * npsf2: -1 * npsf2][0]
                        imt   = canddata[mask, -npsf2:][0]
                        

                        # get astrometric solution
                        matchRADEC = np.load("%s/matchRADEC_%s_%s_%02i-%02i.npy" % (calibrationsdir, field, ccd, sci, ref))
                        (afluxADUB, e_afluxADUB, rmsdeg, CRVAL[0, 0], CRVAL[0, 1], CRPIX[0, 0], CRPIX[0, 1], CD[0, 0, 0], CD[0, 0, 1], CD[0, 1, 0], CD[0, 1, 1], nPV1, nPV2, ordersol) = matchRADEC[0:14]
                        # unpack sol_astrometry_RADEC terms
                        nend = 20
                        if ordersol == 2:
                            nend = 26
                        elif ordersol == 3:
                            nend = 34
                        sol_astrometry_RADEC = matchRADEC[14: nend]
                        # unpack PV terms
                        PV[0] = matchRADEC[nend: nend + int(nPV1 * nPV2)].reshape((int(nPV1), int(nPV2)))

                        # apply astrometric solution
                        (RA, DEC) = RADEC(jpixs[icand], ipixs[icand], CD[0, 0, 0], CD[0, 0, 1], CD[0, 1, 0], CD[0, 1, 1], CRPIX[0, 0], CRPIX[0, 1], CRVAL[0, 0], CRVAL[0, 1], PV[0])
                        (RA, DEC) = applytransformation(ordersol, 15. * RA, DEC, sol_astrometry_RADEC)

                        # write alert
                        writer.append({"Id": idsource,
                                       "field": field,
                                       "ccd": ccd,
                                       "ipix": ipixs[icand],
                                       "jpix": jpixs[icand],
                                       "RA": RA,
                                       "DEC": DEC,
                                       "sciMJD": MJDsci,
                                       "refMJD": MJDref,
                                       "sciepoch": sci,
                                       "refepoch": ref,
                                       "sciairmass": airmass[MJDsci],
                                       "refairmass": airmass[MJDref],
                                       "sciexptime": exptime[MJDsci],
                                       "refexptime": exptime[MJDref],
                                       "filter": filtername[MJDref],
                                       "flux": fluxes[icand],
                                       "e_flux": e_fluxes[icand],
                                       "mag": mags[icand],
                                       "e1_mag": e1_mags[icand],
                                       "e2_mag": e2_mags[icand],
                                       "prob": probs[icand],
                                       "label": labels[icand],
                                       "convref": labels[icand][-1] == "t",
                                       "npsf": 21,
                                       "psf": list(psf),
                                       "imSNR": list(imSNR),
                                       "im1": list(im1),
                                       "im2": list(im2),
                                       "imt": list(imt)})
        else:
            print("WARNING: %s does not exist" % candidatesdir)

writer.close()

