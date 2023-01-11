! Use this to replace the raw estimate_forcing_regression.f90 in GMET. Compiling and run. It has the same example with the Python code


! Main routine for processing station and grid predictor data and calculating the spatial regression

subroutine estimate_forcing_regression (nTotPredictors, gen_sta_weights, sta_weight_name, nwp_input_list, &
  & nDynPRedictors, nwp_vars, nwp_prcp_var, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
  & stnid, stnvar, directory, sta_limit, kfold_trials, pcp, pop, pcperr, obs_max_pcp, tmean, &
  & tmean_err, trange, trange_err, mean_autocorr, mean_tp_corr, error, &
  & pcp_2, pop_2, pcperr_2, tmean_2, tmean_err_2, trange_2, trange_err_2, use_stn_weights)

  ! ==============================================================================================
  ! This routine is called during MODE 2 usage:  creates gridded ensembles from station/point data
  ! ==============================================================================================

  use string_mod
  use utim
  use combination_routines
  use type
  implicit none

  ! ===== start interfaces =======
  interface
    subroutine read_transform_exp (ntimes, file_name, texp)
      use type
      integer (i4b), intent (in) :: ntimes
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: texp (:)
    end subroutine read_transform_exp

    subroutine read_station (stnvar, stnid, directory, st_rec, end_rec, vals, tair_vals, &
   & vals_miss, vals_miss_t, error)
      use type
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: stnid
      character (len=500), intent (in) :: directory
      integer (i4b), intent (in) :: st_rec, end_rec
      real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)
      logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
      integer, intent (out) :: error
    end subroutine read_station

    subroutine compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, &
                                   sta_limit,sta_data,tair_data, &
                                   close_meta,close_meta_t,close_loc,close_loc_t, &
                                   close_count,close_count_t,close_weights,close_weights_t,error)
      use type
      ! inputs
      character(len=500),intent(in) :: sta_weight_name     ! name of station weight binary file
      integer(I4B), intent(in)      :: ngrid               ! number of grid points
      integer(I4B), intent(in)      :: nstns               ! number of stations
      real(DP), intent(in)          :: Z(:,:)              ! grid metadata array
      real(DP), intent(in)          :: X(:,:)              ! station metadata array
      real(DP), intent(in)          :: search_distance     ! default station search distance
      integer(I4B), intent(in)      :: sta_limit           ! maximum number of stations used in regression
      real(DP), intent(in)          :: sta_data(:,:)       ! station data values for precipitation
      real(DP), intent(in)          :: tair_data(:,:,:)    ! station air temperature data
      !in/out
      real(DP), intent(inout)     :: close_meta(:,:,:)     ! 
      real(DP), intent(inout)     :: close_meta_t(:,:,:)
      integer(I4B), intent(inout) :: close_loc(:,:)        ! indices of nearest neighbors for pcp, dim (ngrid, sta_limit)
      integer(I4B), intent(inout) :: close_loc_t(:,:)      ! indices of nearest neighbors for pcp, dim (ngrid, sta_limit)
      integer(I4B), intent(inout) :: close_count(:)
      integer(I4B), intent(inout) :: close_count_t(:)
      real(DP), intent(inout)     :: close_weights(:,:)
      real(DP), intent(inout)     :: close_weights_t(:,:) 
      integer(I4B),intent(inout)  :: error
    end subroutine compute_station_weights

    subroutine write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input
      use type
      ! inputs
      character(len=500), intent(in)    :: sta_weight_name
      real(DP), intent(in)      :: close_meta(:,:,:)
      real(DP), intent(in)      :: close_meta_t(:,:,:)
      integer(I4B), intent(in)  :: close_loc(:,:)
      integer(I4B), intent(in)  :: close_loc_t(:,:)
      real(DP), intent(in)      :: close_weights(:,:)
      real(DP), intent(in)      :: close_weights_t(:,:)
      integer(I4B), intent(in)  :: close_count(:)
      integer(I4B), intent(in)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine write_station_weights

    subroutine read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output
      use type
      !input
      character(len=500), intent(in)    :: sta_weight_name
      !output
      real(DP), intent(out)      :: close_meta(:,:,:)
      real(DP), intent(out)      :: close_meta_t(:,:,:)
      integer(I4B), intent(out)  :: close_loc(:,:)
      integer(I4B), intent(out)  :: close_loc_t(:,:)
      real(DP), intent(out)      :: close_weights(:,:)
      real(DP), intent(out)      :: close_weights_t(:,:)
      integer(I4B), intent(out)  :: close_count(:)
      integer(I4B), intent(out)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine read_station_weights

    subroutine read_nwp(currentTime,nwp_timestep_file,nTotPredictors,nDynPRedictors,nwp_vars,station_grid,x,z,error)
      use string_mod
      use utim
      use type
      character (len=2000), intent (in) :: nwp_timestep_file
      character (len=100), intent (in)  :: nwp_vars(:)
      integer(i4b), intent(in)          :: nTotPredictors     !total number of predictors
      integer(i4b), intent(in)          :: nDynPRedictors        !number of NWP predictors
      integer(i4b), intent(in)          :: station_grid(:)        !nearest grid point for every station location
      real(dp), intent(in)          :: currentTime  !current timestep unix time
      real (dp), intent (inout) :: x(:,:), z(:,:) !station and grid predictor matrices
      integer, intent (inout) :: error !error integer
    end subroutine read_nwp

    subroutine station_grid_correspondence(X,Z,close_weights,close_weights_t,close_loc,close_loc_t,&
                                           & nSta,nearestGridpoint)
      use type
      real(dp), intent(in)        :: X(:,:)
      real(dp), intent(in)        :: Z(:,:)
      real(dp), intent(in)        :: close_weights(:,:)
      real(dp), intent(in)        :: close_weights_t(:,:)
      integer(I4B), intent(in)    :: close_loc(:,:)
      integer(I4B), intent(in)    :: close_loc_t(:,:)
      integer(I4B), intent(in)    :: nSta
      integer(I4B), intent(out)   :: nearestGridpoint(:)
    end subroutine station_grid_correspondence

    subroutine normalize_x (x)
      use type
      real (dp), intent (inout) :: x (:, :)
    end subroutine normalize_x

    !subroutine max_x (x, smax)
    !  use type
    !  real (dp), intent (in) :: x (:)
    !  real (dp), intent (out) :: smax
    !end subroutine max_x
    
    subroutine normalize_y (texp, y)
      use type
      real (dp), intent (in) :: texp !transform exponent
      real (dp), intent (inout) :: y (:)
    end subroutine normalize_y

    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares

    subroutine logistic_regression (x, y, tx, yp, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      integer (i4b), intent (in) :: yp (:)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine logistic_regression

    subroutine generic_corr (prcp_data, tair_data, lag, window, auto_corr, t_p_corr)
      use type
      real (dp), intent (in) :: prcp_data (:)
      real (dp), intent (in) :: tair_data (:, :)
      integer (i4b), intent (in) :: lag
      integer (i4b), intent (in) :: window
      real (dp), intent (out) :: auto_corr
      real (dp), intent (out) :: t_p_corr
    end subroutine generic_corr

    subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
      use type
      implicit none
      real (dp), intent (in) :: lat1, lon1, lat2, lon2
      real (dp), intent (out) :: dist
    end subroutine calc_distance

    subroutine heapsort (n, ra, rn)
      use type
      implicit none
      integer (i4b), intent (in) :: n
      integer (i4b), dimension (:), intent (inout) :: rn
      real (dp), dimension (:), intent (inout) :: ra
    end subroutine heapsort

    subroutine kfold_crossval(X, Y, W, kfold_trials, kfold_nsamp, n_train, max_n_test, xval_combinations, varUncert)
      use type
      implicit none
      !inputs/outputs 
      real(dp), intent(in)      :: X(:,:)                  ! full station attribute array
      real(dp), intent(in)      :: Y(:)                    ! full station variable value vector  
      real(dp), intent(in)      :: W(:,:)                  ! full station diagonal weight matrix
      integer(I4B), intent(in)  :: kfold_trials            ! number of kfold xval trials
      integer(I4B), intent(in)  :: kfold_nsamp             ! number of samples to draw from
      integer(I4B), intent(in)  :: n_train                 ! number of stations in training sample
      integer(i4b), intent(in)  :: max_n_test              ! maximum test sample size over all trials
      integer(I4B), intent(in)  :: xval_combinations(:,:)  ! array of sampling combinations (integer array indicies)
      real(sp), intent(out)     :: varUncert               ! output: uncertainty estimate from kfold trials
    end subroutine kfold_crossval
    
    !subroutine comb()? --- in a module!
  end interface
  ! =========== end interfaces, start code =============

  real (dp), intent (inout)  :: x (:, :), z (:, :)  ! station and grid point description arrays
  real (dp), intent (in)     :: maxdistance         ! max distance for weight function
  integer (i4b), intent (in) :: ngrid               ! number of grid points
  real (dp), intent (in)     :: times (:)           ! time step array

  ! AWW added next
  integer (i4b), intent (in) :: st_rec, end_rec

  character (len=100), intent (in)  :: stnid (:)           ! station id array
  character (len=100), intent (in)  :: stnvar              ! control file variables
  character (len=500), intent (in)  :: directory
  character (len=500), intent(in)   :: gen_sta_weights     ! flag for generating station weight file
  character (len=500), intent(in)   :: sta_weight_name     ! station weight file name
  character (len=500), intent(in)   :: use_stn_weights     ! flag for doing distance weighted regression
  character (len=2000),intent(in)   :: nwp_input_list      ! file containing list of NWP predictor input files
  character (len=100), intent(in)   :: nwp_vars(:)         ! list of nwp predictor variables
  character (len=100), intent(in)   :: nwp_prcp_var        ! name of nwp predictor variable for precipitation

  character (len=2000)              ::  nwp_timestep_file  ! name of a NWP predictor file

  integer(I4B), intent(inout)       :: nTotPredictors      ! number of total predictors
  integer(I4B), intent(in)          :: nDynPRedictors      ! number of NWP predictors
  integer(I4B), intent(in)          :: sta_limit           ! total number of stations in regression sample size
  integer(I4B), intent(in)          :: kfold_trials        ! number of kfold xval trials

  real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :) ! output variables for precipitation
  real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)        ! OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)      ! OLS trange estimate and error

  real (sp), allocatable, intent (out) :: tmean_2 (:, :), tmean_err_2 (:, :)    ! OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange_2 (:, :), trange_err_2 (:, :)  ! OLS trange estimate and error
  real (sp), allocatable, intent (out) :: pcp_2 (:, :), pop_2 (:, :), pcperr_2 (:, :)

  integer, intent (out)   :: error              ! integer error flag
  real (dp), intent (out) :: mean_autocorr (:)  ! mean autocorrelation from all stations over entire time period
  real (dp), intent (out) :: mean_tp_corr (:)   ! mean correlation for mean temp and precip

  ! vary at each grid point and time step
  real (dp), intent (out) :: obs_max_pcp (:, :) ! max of normalized time step precipitation

  ! Local declarations
  real (dp), allocatable :: y (:), b (:)
  real (dp), allocatable :: y_red (:)           ! reduced matrix for the predictand
  real (dp), allocatable :: x_red (:, :)        ! reduced matrix for the predictors
  real (dp), allocatable :: tmp_x_red (:, :)    ! reduced matrix for extended predictors
  real (dp), allocatable :: x_red_t (:, :)      ! reduced matrix for predictors, transposed (?)
  real (dp), allocatable :: Z_reg(:)            ! final grid point predictor vector

  ! condensed these variables into just 4 that get re-used
  real (dp), allocatable :: twx_red (:, :), tx_red (:, :)        ! reduced matricies (orig)
  real (dp), allocatable :: twx_red_2 (:, :), tx_red_2 (:, :)    ! reduced matricies (dims 2)

  real (dp), allocatable :: w_base (:, :)                        ! initial distance weight matrix
  real (dp), allocatable :: w_pcp_red (:, :), w_temp_red (:, :), w_pcp_red_vector (:)  ! reduced distance weight matricies

  real (dp), allocatable :: y_tmean (:), y_trange (:)            ! transformed station data arrays
  real (dp), allocatable :: y_tmean_red (:), y_trange_red (:)    ! transformed station data arrays
  real (dp), allocatable :: stn_prcp (:), prcp_data (:,:), tair_data (:,:,:), stn_tair (:,:)  ! orig stn data arrays
  real (dp), allocatable :: auto_corr (:)      ! lag-1 autocorrelation for stations over entire time period used
  real (dp), allocatable :: t_p_corr (:)       ! correlation between temp and precip
  integer (i4b), allocatable :: yp (:)         ! binary for logistic regression
  integer (i4b), allocatable :: yp_red (:)     ! reduced binary for logistic regression

  logical, allocatable :: stn_miss (:), stn_miss_t (:) ! missing value logical arrays

  real (dp) :: auto_corr_sum, tp_corr_sum

  real (dp) :: step_max                               ! timestep statistics  
  real (dp) :: ss_tot, ss_res                         ! r-squared and variance correction   ! not yet used

  integer (i4b) :: xsize, sta_limit2                              ! size of second dimension of input X array
  integer (i4b) :: ntimes, nstns
  integer (i4b) :: t, i, j, g, ndata, nodata, cnt
  integer (i4b) :: ndata_t, nodata_t
  integer (i4b) :: lag, window
  integer (i4b) :: auto_cnt, tp_cnt

  integer(I4B)  :: nBase                              ! number of geophysical predictors
  integer(I4B)  :: prcpPredictInd                     ! index of precipitation predictor in predictor vector
  integer(I4B),allocatable  :: noSlopePredicts(:)
  integer(I4B),allocatable  :: noPrcpPredicts(:)
  integer(I4B),allocatable  :: noPrcpNoSlopePredicts(:)
  logical                   :: drop_precip = .false.
  logical                   :: no_precip = .false.

  ! variables for tracking closest N stations for precipitation
  !integer (i4b), parameter   :: sta_limit = 30
  integer (i4b), allocatable :: close_loc (:, :)
  integer (i4b), allocatable :: close_count (:)
  real (dp), allocatable     :: close_weights (:, :)
  real (dp), allocatable     :: close_meta (:, :, :)
  real (dp)                  :: max_distance
  real (dp), parameter       :: search_distance = 1000.0

  ! variables for tracking closest N stations for temperature
  integer (i4b), allocatable :: close_loc_t (:, :)
  integer (i4b), allocatable :: close_count_t (:)
  real (dp), allocatable     :: close_weights_t (:, :)
  real (dp), allocatable     :: close_meta_t (:, :, :)
  real (dp)                  :: max_distance_t
  
  real (dp)                  :: tmp_weight, wgtsum_pcp, wgtsum_temp, errsum

  integer(I4B), allocatable :: nearestGridpoint(:)

  integer (i4b) :: slope_flag_pcp
  integer (i4b) :: slope_flag_temp

  ! variables to check for singular matrix
  real (dp), allocatable :: tmp (:, :)
  real (dp), allocatable :: vv (:),vv_temp(:)

  ! variables for timing code AJN
  integer (i4b) :: t1, t2, count_rate
  integer (i4b) :: tg1, tg2

  ! variables for keeping track of max_distance modifications
  integer (i4b), allocatable :: expand_flag (:), expand_flag_t (:)
  real (dp), allocatable     :: expand_dist (:), expand_dist_t (:)

  ! kfold cross validation variables 
  integer(i4b), allocatable :: xval_combinations(:,:)         ! array of sampling combination indices
  integer(i4b), allocatable :: sampling_order(:)              ! vector of sampling indices
  integer(i4B)              :: kfold_nsamp                    ! number of samples to draw from !AW same as nstation?
  
  integer(I4B)              :: n_total, n_test, n_train       ! renamed loop limit variables  
  integer(i4b)              :: end_test_stn, start_test_stn   ! endpoints of range of indices in test 
  integer(i4b)              :: trial_n_test, trial_n_train    ! test & train sample size in each trial (varies)
  integer(i4b)              :: max_trial_n_test               ! size of largest test sample in trials
  integer(i4b)              :: nte, ntr                       ! counters for start/end of train/test indices in array
  
  real (dp), parameter      :: transform_exp = 4.0            ! power law transform exponent; must match setting in ens_generation code

  !==============================================================!
  !                     code starts below here                   !
  !==============================================================!
	
  
  sta_limit2 = 35
  xsize = 6
  allocate (x_red(sta_limit2, xsize))
  allocate (tx_red(xsize, sta_limit2))
  allocate (w_pcp_red(sta_limit2, sta_limit2))
  allocate (twx_red(xsize, sta_limit2))
  allocate (w_pcp_red_vector(sta_limit2))
  allocate (y_red(sta_limit2))
  
  ! data
  x_red = transpose(reshape((/1.000000e+00,  3.458940e+01, -1.184547e+02,  6.367000e+02, -7.564000e+00,  2.711000e+00, &
1.000000e+00,  3.438690e+01, -1.185342e+02,  3.758000e+02, -1.046100e+01,  4.945000e+00, &
1.000000e+00,  3.449390e+01, -1.182714e+02,  8.635000e+02, -6.687000e+00,  8.084000e+00, &
1.000000e+00,  3.470500e+01, -1.184242e+02,  9.327000e+02, -4.846000e+00, -2.170000e-01, &
1.000000e+00,  3.432940e+01, -1.184006e+02,  4.724000e+02, -1.102800e+01,  7.943000e+00, &
1.000000e+00,  3.431670e+01, -1.185500e+02,  7.215000e+02, -1.107900e+01,  5.416000e+00, &
1.000000e+00,  3.440610e+01, -1.187569e+02,  2.170000e+02, -1.372900e+01, -5.280000e-01, &
1.000000e+00,  3.448830e+01, -1.181419e+02,  9.623000e+02, -4.959000e+00,  8.682000e+00, &
1.000000e+00,  3.424470e+01, -1.185250e+02,  2.890000e+02, -1.076500e+01,  4.833000e+00, &
1.000000e+00,  3.474110e+01, -1.182117e+02,  7.126000e+02, 3.260000e-01, -1.710000e-01, &
1.000000e+00,  3.474360e+01, -1.187242e+02,  1.374600e+03, -6.755000e+00, -3.813000e+00, &
1.000000e+00,  3.475000e+01, -1.187167e+02,  1.226800e+03, -3.942000e+00, -2.677000e+00, &
1.000000e+00,  3.458830e+01, -1.180939e+02,  7.943000e+02, 1.093000e+00,  5.567000e+00, &
1.000000e+00,  3.429470e+01, -1.181883e+02,  7.090000e+02, -1.306000e+01,  9.611000e+00, &
1.000000e+00,  3.462940e+01, -1.180836e+02,  7.693000e+02, 2.977000e+00,  3.458000e+00, &
1.000000e+00,  3.420060e+01, -1.183575e+02,  2.259000e+02, -1.275200e+01,  7.601000e+00, &
1.000000e+00,  3.418190e+01, -1.185744e+02,  2.408000e+02, -9.294000e+00,  3.120000e+00, &
1.000000e+00,  3.422220e+01, -1.182378e+02,  4.709000e+02, -1.413100e+01,  8.969000e+00, &
1.000000e+00,  3.418220e+01, -1.181383e+02,  3.435000e+02, -1.502000e+01,  8.760000e+00, &
1.000000e+00,  3.483280e+01, -1.188650e+02,  1.092700e+03,  2.736000e+00, -1.163000e+00, &
1.000000e+00,  3.423080e+01, -1.180711e+02,  1.740400e+03, -1.510600e+01,  1.001800e+01, &
1.000000e+00,  3.408420e+01, -1.185989e+02,  2.271000e+02, -8.123000e+00,  2.397000e+00, &
1.000000e+00,  3.406970e+01, -1.184428e+02,  1.311000e+02, -9.495000e+00,  3.607000e+00, &
1.000000e+00,  3.414810e+01, -1.181444e+02,  2.633000e+02, -1.502000e+01,  8.760000e+00, &
1.000000e+00,  3.450280e+01, -1.178944e+02,  9.296000e+02, 1.073000e+00,  3.717000e+00, &
1.000000e+00,  3.401580e+01, -1.184514e+02,  5.300000e+01, -8.343000e+00,  2.458000e+00, &
1.000000e+00,  3.405110e+01, -1.182353e+02,  7.010000e+01, -1.213700e+01,  5.152000e+00, &
1.000000e+00,  3.410610e+01, -1.181000e+02,  1.372000e+02, -1.514300e+01,  8.431000e+00, &
1.000000e+00,  3.400750e+01, -1.184997e+02,  4.300000e+00, -8.343000e+00,  2.458000e+00, &
1.000000e+00,  3.400530e+01, -1.184136e+02,  2.070000e+01, -9.192000e+00,  2.974000e+00, &
1.000000e+00,  3.502330e+01, -1.187497e+02,  4.343000e+02, 4.263000e+00,  9.616000e+00, &
1.000000e+00,  3.431750e+01, -1.178419e+02,  1.629200e+03, -9.666000e+00,  3.394000e+00, &
1.000000e+00,  3.407640e+01, -1.188803e+02,  4.877000e+02, -7.230000e+00,  3.000000e+00, &
1.000000e+00,  3.510110e+01, -1.184222e+02,  1.286300e+03, -5.493000e+00, -1.356000e+00, &
1.000000e+00,  3.504920e+01, -1.181619e+02,  8.336000e+02, -4.850000e+00, -8.311000e+00/), (/xsize, sta_limit2/) ))
                    
  w_pcp_red_vector = (/0.0388809 , 0.03832985, 0.03817355, 0.03802249, 0.03745646, &
0.03711991, 0.03588912, 0.03574807, 0.0351188 , 0.03482195, &
0.03478868, 0.03477532, 0.03415687, 0.03338446, 0.03333384, &
0.03278177, 0.03189196, 0.0317251 , 0.02656133, 0.02649087, &
0.02641929, 0.02530463, 0.02513205, 0.02471289, 0.0241944 , &
0.02079555, 0.02064448, 0.02017512, 0.02005636, 0.01975397, &
0.01841752, 0.01702782, 0.01630212, 0.01601774, 0.01559474/)

  w_pcp_red = 0.0
  do i = 1, sta_limit2
	w_pcp_red (i, i) = w_pcp_red_vector(i)
  end do
	
  y_red = (/4.416522  ,  1.9206624 ,  2.5261984 ,  3.1489668 ,  3.6673174 , &
        3.7148647 ,  2.908574  ,  1.9249749 ,  1.885819  ,  0.7568283 , &
        2.6417985 ,  1.5446663 ,  1.5836291 ,  4.1812468 ,  1.0297337 , &
        0.63115644,  2.529779  ,  2.1395035 ,  3.5515766 ,  3.2349849 , &
        4.497094  ,  0.1865406 ,  2.4116774 ,  3.3983598 ,  3.8984742 , &
        0.7568283 ,  2.0682702 ,  1.264296  ,  0.34635973,  2.6417985 , &
        5.8850327 ,  4.673559  ,  2.4834137 ,  4.097801  , -4.    /)
  

  tx_red  = transpose (x_red)
  twx_red = matmul (tx_red, w_pcp_red)
  
  print *, 'x_red:'
  do i = 1, sta_limit2
    print *, x_red(i, :)
  end do
  
  print *, 'w_pcp_red:'
  do i = 1, sta_limit2
    print *, w_pcp_red(i, :)
  end do
  
              
  call system_clock (tg1, count_rate)
  
  do i = 1, 100000
    call least_squares (x_red, y_red, twx_red, b)        ! solve regression for stations
  end do
  
  call system_clock (tg2, count_rate)
  print *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)
  print *,b
  
  ! call logistic_regression (x_red, y_red, twx_red, yp_red, b)

end subroutine estimate_forcing_regression
