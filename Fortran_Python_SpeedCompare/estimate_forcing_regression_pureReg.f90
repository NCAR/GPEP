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
	
  
  sta_limit2 = 7
  xsize = 3
  allocate (x_red(sta_limit2, xsize))
  allocate (tx_red(xsize, sta_limit2))
  allocate (w_pcp_red(sta_limit2, sta_limit2))
  allocate (twx_red(xsize, sta_limit2))
  allocate (w_pcp_red_vector(sta_limit2))
  allocate (y_red(sta_limit2))
  
  ! data
  x_red = transpose(reshape((/1.0, 1.2, 1.4, &
                    1.0, 2.2, 1.42, &
                    1.0, 4.2, 3.33, &
                    1.0, 12.0, 3.5, &
                    1.0, 0.4, 2.0, &
                    1.0, 0.523, 1.2, &
                    1.0, 5.5, 2.3/), (/xsize, sta_limit2/) ))
                    
  w_pcp_red_vector = (/0.1, 0.5, 0.23, 0.6, 0.23, 0.66, 0.9/)
  w_pcp_red = 0.0
  do i = 1, sta_limit2
	w_pcp_red (i, i) = w_pcp_red_vector(i)
  end do
	
  y_red = (/0.2, 1.2, 0.4, 0.3, 0.55, 0.9, 0.44/)
  

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
