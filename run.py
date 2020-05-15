def merge_correction_stnerror(path_reastn_cv, stnlle, stndata, readata_stn, taintestindex, nearstn_locl1, nearstn_distl1,
         nearstn_locl2, nearstn_distl2, dividenum, var, hwsize, corrmode, anombound, weightmode):
    # use 2-layer cross-validation to estimate the weight and independent data of merge/correction data
    reanum, nstn, ntimes = np.shape(readata_stn)

    # initialization
    reacorr_stn = np.nan * np.zeros([reanum, nstn, ntimes], dtype=np.float32)  # corrected reanalysis data
    reamerge_weight_stn = np.nan * np.zeros([nstn, reanum])  # weight used to obtain reamerge_stn
    reamerge_stn = np.nan * np.zeros([nstn, ntimes], dtype=np.float32)  # merged reanalysis at station points

    for lay1 in range(dividenum):
        print('Correction/Merging at station points. Layer-1:', lay1)
        # extract train and test index for layer-1
        if var == 'trange':
            vari = 'tmean'  # trange and tmean have the same index
        else:
            vari = var
        trainindex1 = taintestindex[vari + '_trainindex1'][lay1, :]
        testindex1 = taintestindex[vari + '_testindex1'][lay1, :]
        stndata_trainl1 = stndata[trainindex1, :]
        stndata_testl1 = stndata[testindex1, :]
        stnlle_trainl1 = stnlle[trainindex1, :]
        stnlle_testl1 = stnlle[testindex1, :]

        # filename: save inputs for each layer-1
        file_reacorrl1 = path_reastn_cv + 'reacorr_' + var + '_layer1_' + str(lay1 + 1) + '.npz'
        if os.path.isfile(file_reacorrl1):
            datatemp = np.load(file_reacorrl1)
            reacorr_trainl1 = datatemp['reacorr']
            weight_trainl1 = datatemp['reaweight']
            del datatemp
        else:

            # layer-2: start
            reacorr_trainl1 = np.zeros([reanum, len(trainindex1), ntimes], dtype=np.float32)
            weight_trainl1 = np.zeros([len(trainindex1), reanum], dtype=np.float32)

            for lay2 in range(dividenum):
                print('Correction/Merging at station points. Layer-2:', lay2)
                # extract train and test index for layer-2 (subsets of trainindex1)
                trainindex2 = taintestindex[vari + '_trainindex2'][lay1, lay2, :]
                testindex2 = taintestindex[vari + '_testindex2'][lay1, lay2, :]
                stndata_trainl2 = stndata[trainindex2, :]
                stndata_testl2 = stndata[testindex2, :]
                stnlle_trainl2 = stnlle[trainindex2, :]
                stnlle_testl2 = stnlle[testindex2, :]

                for rr in range(reanum):
                    print('Correction/Merging at station points. Reanalysis:', rr)
                    readata_trainl2 = readata_stn[rr, trainindex2, :]
                    readata_testl2 = readata_stn[rr, testindex2, :]

                    ### calculate corrected reanalysis data
                    # calculate anomaly at the train stations
                    anom_ori = calculate_anomaly(readata_trainl2, stndata_trainl2, hwsize, corrmode,
                                                 upbound=anombound[1], lowbound=anombound[0])
                    # extrapolate the ratio to the test stations
                    anom_ext = extrapolation(anom_ori, nearstn_locl2[lay1,lay2,:],nearstn_distl2[lay1,lay2,:])
                    # correct data at the test stations
                    readata_testl2_corr = error_correction(readata_testl2, anom_ext, mode=corrmode)
                    tf, index = ismember(testindex2, trainindex1)
                    reacorr_trainl1[rr, index, :] = readata_testl2_corr

            if weightmode == 'BMA':
                for i in range(len(trainindex1)):
                    dobs = stndata_trainl1[i, :]
                    drea = reacorr_trainl1[:,i,:].T
                    w, sigma, sigma_s = bma(dobs, drea)
                    weight_trainl1[i, :] = w
            else:
                for rr in range(reanum):
                    weight_trainl1[:, rr] = calweight(stndata_trainl1, reacorr_trainl1[rr, :, :], weightmode)

            np.savez_compressed(file_reacorrl1, reacorr=reacorr_trainl1, stnlle=stnlle_trainl1,
                                reaweight=weight_trainl1)
            # layer-2: end


        # output: repeat error correction using train stations in layer-1 (as in layer-2 only 0.9*0.9=0.81 stations are used)
        # extrapolate from train stations to test stations
        for rr in range(reanum):
            readata_trainl1 = readata_stn[rr, trainindex1, :]
            readata_testl1 = readata_stn[rr, testindex1, :]
            anom_ori = calculate_anomaly(readata_trainl1, stndata_trainl1, hwsize, corrmode,
                                         upbound=anombound[1], lowbound=anombound[0])
            anom_ext = extrapolation(anom_ori, nearstn_locl1[lay1,:],nearstn_distl1[lay1,:])
            readata_testl1_corr = error_correction(readata_testl1, anom_ext, mode=corrmode)
            reacorr_stn[rr, testindex1, :] = readata_testl1_corr

        # output: extrapolate the weight from train stations to test stations (in layer-1)
        weight_testl1 = extrapolation(weight_trainl1, nearstn_locl1[lay1,:],nearstn_distl1[lay1,:])
        reamerge_weight_stn[testindex1, :] = weight_testl1

        # output: merge reanalysis products at the test stations
        nstn_testl1 = len(testindex1)
        mergedata_testl1 = np.nan * np.zeros([nstn_testl1, ntimes])
        for i in range(ntimes):
            datain = np.zeros([nstn_testl1, reanum], dtype=np.float32)
            for rr in range(reanum):
                datain[:, rr] = reacorr_stn[rr, testindex1, i]
            dataout = weightmerge(datain, weight_testl1)
            mergedata_testl1[:, i] = dataout
        reamerge_stn[testindex1, :] = mergedata_testl1

    return reamerge_stn, reamerge_weight_stn, reacorr_stn
