# generate random numbers

import numpy as np

# random field generation using Fortran GMET codes
# translation from: https://github.com/NCAR/GMET/blob/master/source/ens_generation/spcorr_grd.f90, and https://github.com/NCAR/GMET/blob/master/source/ens_generation/field_rand.f90

def LU_Decomposition(a, b):
    # replace Fortran ludcmp and lubksb
    # https://github.com/NCAR/GMET/blob/master/source/sp_regression/regression_routines.f90
    # ! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
    # ! Input:
    # !   X  = An m by n array.
    # !   TX = Precalculated transpose array of X, size n by m
    # !   Y  = An m-element vector containing the right-hand side of the linear system Ax = b.
    # ! Output:
    # !   B  = An n-element vector.

    deta = np.linalg.det(a)  # Compute the determinant of an array
    if deta == 0:
        # print('Singular matrix')
        c = 0
    else:
        ainv = np.linalg.inv(a)
        c = np.matmul(ainv, b)
    return c


def spcorr_grd(grid_lat, grid_lon, clen, outfile):

    # example input parameters
    # nspl1 = 100
    # nspl2 = 120
    # clen = 150
    # lat = np.arange(0, 10, 0.1)
    # lon = np.arange(0, 12, 0.1)
    # grid_lat = np.tile(lat[:, np.newaxis], [1, 120])
    # grid_lon = np.tile(lon[np.newaxis, :], [100, 1])
    nspl1, nspl2 = grid_lat.shape

    grid_lat = grid_lat.astype(np.float64) # important for high resolution
    grid_lon = grid_lon.astype(np.float64)
    # grid_lat = grid_lat.T # to reproduce Fortran GMET
    # grid_lon = grid_lon.T

    # ----------------------------------------------------------------------------------------
    # (0) CHECK THAT SPCORR IS NOT POPULATED ALREADY
    # ----------------------------------------------------------------------------------------


    # ----------------------------------------------------------------------------------------
    # (1) DEFINE HYPER-PARAMETERS
    # ----------------------------------------------------------------------------------------
    nnst = 10 # number of nests
    nloc = 3 # number of local points to include in the estimation

    # ----------------------------------------------------------------------------------------
    # (2) ALLOCATE SPACE FOR OUTPUT ARRAYS
    # ----------------------------------------------------------------------------------------
    # define the maximum number of previously generated points
    maxp = (nloc*2+1) ** 2
    gmsk = np.full((nspl1, nspl2), False, dtype=bool)
    spcorr_ipos = np.empty((nspl1, nspl2), dtype=object) # np.frompyfunc(list, 0, 1)(np.empty((3,2), dtype=object))
    spcorr_jpos = np.empty((nspl1, nspl2), dtype=object)
    spcorr_wght = np.empty((nspl1, nspl2), dtype=object)
    spcorr_sdev = np.empty((nspl1, nspl2), dtype=float)
    iorder = np.zeros(nspl1*nspl2, dtype=int)
    jorder = np.zeros(nspl1*nspl2, dtype=int)

    # ----------------------------------------------------------------------------------------
    # (3) LOOP THROUGH THE DIFFERENT GRID RESOLUTIONS (PROCESS COARSE RESOLUTION FIRST)
    # ----------------------------------------------------------------------------------------
    sdev = 0 # Initialize SDEV (added to account for the first point)  EÖH
    iprc = 0 # counter for the number of grid points processed
    for ires in range(nnst-1, 0-1, -1):
        incr = 2 ** ires # increment(2 ** 4 = 16, 2 ** 3 = 8, 2 ** 2 = 4, 2 ** 1 = 2, 2 ** 0 = 1)
        print('Working on Loop: ', ires)
        # ---------------------------------------------------------------------------------------
        # (4) LOOP THROUGH THE LAT-LON OF THE GRID AT A GIVEN RESOLUTION
        # ---------------------------------------------------------------------------------------
        for isp1 in range(0, nspl1, incr):
            for isp2 in range(0, nspl2, incr):
                # ! check that "current" point has not been generated yet
                if not gmsk[isp1, isp2]:
                    # ! allocate space to store the (i,j) position, and weights
                    ipos = np.zeros(maxp, dtype=int)
                    jpos = np.zeros(maxp, dtype=int)
                    wght = np.zeros(maxp, dtype=float)
                    # ! save the (i,j) position of iprc
                    iorder[iprc] = isp1
                    jorder[iprc] = isp2
                    # ! increment IPRC
                    iprc = iprc + 1
                    # ------------------------------------------------------------------------------------
                    # (5) IDENTIFY PREVIOUSLY GENERATED POINTS
                    # ------------------------------------------------------------------------------------
                    k = 0 # initialize the number of previous points generated to zero
                    # ! loop through points in the local neighbourhood
                    for jsp1 in range(max(0, isp1-(incr*nloc)), min(isp1+(incr*nloc)+1, nspl1)): #doubt
                        for jsp2 in range(max(0, isp2-(incr*nloc)), min(isp2+(incr*nloc)+1, nspl2)): #doubt
                            # ! check to see if the "local" point has been generated previously
                            if gmsk[jsp1, jsp2]:
                                ipos[k] = jsp1
                                jpos[k] = jsp2
                                k = k + 1
                    # include the (i,j) of the current point
                    ipos[k] = isp1
                    jpos[k] = isp2
                    k = k + 1
                    # ...and save the number of points
                    npts = k
                    # check that there are at least two points
                    if k >= 2:
                        # ------------------------------------------------------------------------------------
                        # (6) COMPUTE THE CORRELATION AMONG PREVIOUSLY GENERATED POINTS
                        # ------------------------------------------------------------------------------------
                        corr = np.zeros([k-1, k-1])
                        gvec = np.zeros(k-1)
                        twgt = np.zeros(k-1)
                        indx = np.zeros(k-1, dtype=int)
                        # Note that the vector of previously generated points includes the current point as its
                        # last element.  The correlation among all previously generated points are computed over
                        # elements (1...k-1) and saved in the matrix corr.  The correlation between previously
                        # generated points (1...k-1) and the current point (k) is saved in ther vector gvec.
                        for iprev in range(0, k):
                            for jprev in range(0, iprev+1):
                                if iprev == jprev:
                                    if iprev <= k-2:
                                        corr[iprev, jprev] = 1.0
                                else:
                                    lon1 = np.deg2rad(grid_lon[ipos[iprev], jpos[iprev]]) # NOTE, iprev, lon
                                    lon2 = np.deg2rad(grid_lon[ipos[jprev], jpos[jprev]]) # NOTE, jprev, lon
                                    lat1 = np.deg2rad(grid_lat[ipos[iprev], jpos[iprev]]) # NOTE, iprev, lat
                                    lat2 = np.deg2rad(grid_lat[ipos[jprev], jpos[jprev]]) # NOTE, jprev, lat
                                    # ! compute distance (km) - on the surface of a sphere
                                    dist = 6378.0 * np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2))
                                    # ! compute correlation
                                    if iprev <= k-2:
                                        # ! correlation among all previously generated points (1...k-1,1...k-1) -- corr
                                        corr[iprev, jprev] = np.exp(-(dist / clen))
                                        corr[jprev, iprev] = corr[iprev, jprev]
                                    else:
                                        # ! correlation between all previously generated points and the current point -- gvec
                                        if jprev <= k-2:
                                            gvec[jprev] = np.exp(-(dist / clen))
                        #  ------------------------------------------------------------------------------------
                        #  (7) COMPUTE THE WEIGHTS
                        #  ------------------------------------------------------------------------------------
                        #  Note that the vector of previously generated points includes the current point as its
                        #  last element.  The correlation among all previously generated points are computed over
                        #  elements (1...k-1) and saved in the matrix corr.  The correlation between previously
                        #  generated points (1...k-1) and the current point (k) is saved in ther vector gvec.
                        #  special case of the bi-variate normal
                        if k == 2:
                            wght[0] = gvec[0]
                            sdev = np.sqrt(1 - gvec[0] ** 2)
                        # all other points
                        else:
                            # temporary weight (GVEC is over-written)
                            twgt[0: k - 1] = gvec[0: k - 1]
                            # estimate weights
                            twgt[0: k - 1] = LU_Decomposition(corr, twgt)
                            # or
                            # indx, d = ludcmp(corr)
                            # twgt[0: k - 1] = lubksb(corr, indx, twgt)
                            # or
                            # from regression import linearsolver
                            # twgt = linearsolver(list(corr), np.shape(corr)[0], list(twgt))
                            # save weights and variance
                            wght[0: k - 1] = twgt[0: k - 1]
                            sdev = np.sqrt(1. - np.dot(gvec[0:k - 1], twgt[0: k - 1]))
                    gmsk[isp1, isp2] = True
                    # -------------------------------------------------------------------------------------
                    # (8) SAVE WEIGHTS IN THE SPATIAL CORRELATION STRUCTURE
                    # -------------------------------------------------------------------------------------
                    # populate the structures (-1 excludes the current (i,j) point)
                    spcorr_ipos[isp1, isp2] = ipos[0:npts-1]
                    spcorr_jpos[isp1, isp2] = jpos[0:npts-1]
                    spcorr_wght[isp1, isp2] = wght[0:npts-1]
                    spcorr_sdev[isp1, isp2] = sdev

    # save structure to output file
    np.savez_compressed(outfile, spcorr_ipos=spcorr_ipos, spcorr_jpos=spcorr_jpos, spcorr_wght=spcorr_wght, spcorr_sdev=spcorr_sdev, iorder=iorder, jorder=jorder)



def field_rand(spcorr_jpos, spcorr_ipos, spcorr_wght, spcorr_sdev, iorder, jorder, seed=np.nan):
    # ----------------------------------------------------------------------------------------
    # GET THE NUMBER OF X AND Y POINTS AND ALLOCATE SPACE FOR THE RANDOM GRID
    # ----------------------------------------------------------------------------------------

    if ~np.isnan(seed):
        np.random.seed(seed)

    nlon, nlat = spcorr_jpos.shape # nlon or nlat just represents row/col, not related to real lat/lon
    cran = np.nan * np.zeros([nlon, nlat], dtype=float)

    for igrd in range(nlon * nlat):
        ilon = iorder[igrd]
        ilat = jorder[igrd]
        # ! assign a random number to the first grid-point
        if igrd == 0:
            aran = np.random.normal(0, 1)
            cran[ilon, ilat] = aran
        # ! process gridpoints 2,...,n
        else:
            # ! get the number of "previously generated points"
            nprv = len(spcorr_wght[ilon, ilat])
            vprv = np.zeros(nprv)
            # ! build a vector of previously generated points
            for iprev in range(0, nprv):
                jlon = spcorr_ipos[ilon, ilat][iprev]
                jlat = spcorr_jpos[ilon, ilat][iprev]
                vprv[iprev] = cran[jlon, jlat] # (previously generated point)
            # ! and generate the "current" point
            aran = np.random.normal(0, 1)
            xbar = np.dot(vprv[0:nprv], spcorr_wght[ilon, ilat][0:nprv])
            cran[ilon, ilat] = xbar + spcorr_sdev[ilon, ilat] * aran

    return cran

