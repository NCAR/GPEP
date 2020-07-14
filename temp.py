rankyv = np.zeros([dnum, gnum])
for d in range(dnum):
    for g in range(gnum):
        dstn_dg = data_stn[d, g]
        if dstn_dg > -100:
            dens_dg = data_ens[d, g, :]
            rankyv[d, g] = np.sum(dens_dg < dstn_dg)