import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# Paramètres fixes
Fits1 = "71340_080302-1.fits"
Fits2 = "71340_080302-2.fits"
Fits3 = "71340_080302-3.fits"
maxdim = 4e3
PRA = "CRVAL1"
PDEC = "CRVAL2"
CPR1 = "CDELT1"
CPR2 = "CDELT2"
CPR3 = "CDELT3"
DMAX = "DATAMAX"
DMIN = "DATAMIN"
BZERO = "BZERO"
BSCALE = "BSCALE"

# Ouverture des fichiers FITS en lecture
with fits.open(Fits1) as fits1, fits.open(Fits2) as fits2, fits.open(Fits3) as fits3:
    header1 = fits1[0].header
    header2 = fits2[0].header
    header3 = fits3[0].header
    naxis = header1['NAXIS']
    bitpix = header1['BITPIX']
    naxes = [header1[f'NAXIS{i}'] for i in range(1, naxis + 1)]
    pra = header1[PRA]
    pdec = header1[PDEC]
    cpr1_1 = header1[CPR1]
    cpr2_1 = header1[CPR2]
    dmax_1 = header1[DMAX]
    dmin_1 = header1[DMIN]
    dmax_2 = header2[DMAX]
    dmin_2 = header2[DMIN]
    dmax_3 = header3[DMAX]
    dmin_3 = header3[DMIN]
    bzero = header1[BZERO]
    bscale = header1[BSCALE]

    print("Nombre d'axes:", naxis)
    print("Type de données:", bitpix)
    print("Taille des axes:", naxes)
    print(f"{PRA} = {pra}")
    print(f"{PDEC} = {pdec}")
    print(f"{CPR1}_1 = {cpr1_1}")
    print(f"{CPR2}_2 = {cpr2_1}")
    print(f"{DMAX} = {dmax_1}")
    print(f"{DMIN} = {dmin_1}")

    # Lecture des données de l'image
    imageArray1 = fits1[0].data
    imageArray2 = fits2[0].data
    imageArray3 = fits3[0].data
    
    # Opérations sur les données de l'image
    imageArray1 = imageArray1 * bscale - bzero
    imageArray1 = np.where(imageArray1 > dmax_1, dmax_1, imageArray1)
    imageArray1 = np.where(imageArray1 < dmin_1, dmin_1, imageArray1)
    
    imageArray2 = imageArray2 * bscale - bzero
    imageArray2 = np.where(imageArray2 > dmax_2, dmax_2, imageArray2)
    imageArray2 = np.where(imageArray2 < dmin_2, dmin_2, imageArray2)

    imageArray3 = imageArray3 * bscale - bzero
    imageArray3 = np.where(imageArray3 > dmax_3, dmax_3, imageArray3)
    imageArray3 = np.where(imageArray3 < dmin_3, dmin_3, imageArray3)

    # Inversion verticale de la matrice
    imageArray1 = np.flipud(imageArray1)
    imageArray2 = np.flipud(imageArray2)
    imageArray3 = np.flipud(imageArray3)

    # Coordonnées (pixels) des points sur les 3 FITS
    points1 = [(2141.8, imageArray1.shape[0]-1310.8), (2227.8, imageArray1.shape[0]-1279.8), (2258.7, imageArray1.shape[0]-1264.8), (2323.4, imageArray1.shape[0]-1230.1), (1933.6, imageArray1.shape[0]-1513.7), (2334.2, imageArray1.shape[0]-1308.5), (2411.8, imageArray1.shape[0]-1448.8), (2169.7, imageArray1.shape[0]-1503.9), (2232.3, imageArray1.shape[0]-1560.4), (2239.9, imageArray1.shape[0]-1547.6)]
    points2 = [(1877.7, imageArray2.shape[0]-1489.5), (1964.2, imageArray2.shape[0]-1458.3), (1994.1, imageArray2.shape[0]-1442.2), (2061.0, imageArray2.shape[0]-1408.8), (1669.1, imageArray2.shape[0]-1692.3), (2070.2, imageArray2.shape[0]-1487.7), (2147.7, imageArray2.shape[0]-1627.1), (1907.7, imageArray2.shape[0]-1683.2), (1969.0, imageArray2.shape[0]-1740.3), (1976.6, imageArray2.shape[0]-1726.6)]
    points3 = [(1766.9, imageArray3.shape[0]-1527.5), (1853.2, imageArray3.shape[0]-1496.7), (1884.0, imageArray3.shape[0]-1481.7), (1950.8, imageArray3.shape[0]-1447.9), (1989.8, imageArray3.shape[0]-1526.7), (1559.1, imageArray3.shape[0]-1730.8), (1796.2, imageArray3.shape[0]-1721.8), (2037.8, imageArray3.shape[0]-1666.3), (1865.9, imageArray3.shape[0]-1766.4), (1859.1, imageArray3.shape[0]-1779.8)]

    
    # Affichage de l'image
    plt.figure(figsize=(12, 10))
    plt.imshow(imageArray1, cmap='gray', extent=[0, imageArray1.shape[1], 0, imageArray1.shape[0]])
    plt.title("Image FITS 1")
    plt.show()

    plt.figure(figsize=(12, 10))
    plt.imshow(imageArray2, cmap='gray', extent=[0, imageArray2.shape[1], 0, imageArray2.shape[0]])
    plt.title("Image FITS 2")
    plt.show()

    plt.figure(figsize=(12, 10))
    plt.imshow(imageArray3, cmap='gray', extent=[0, imageArray3.shape[1], 0, imageArray3.shape[0]])
    plt.title("Image FITS 3")
    plt.show()