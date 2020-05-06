
%%% JRA55
%%% Note: JRA55 can only read original global data since some parameters in the
%%% functions are fixed. Besides, the initialization time for prcp/snow and
%%% tmin/tmax/tmean is different. so the first data record in each
%%% tmin/tmax/tmean is deleted. If data are downloaded in another approach,
%%% this code must be revised because time is different. Future, more
%%% flexible settings are needed.

% Plato path
Inpath={'/datastore/GLOBALWATER/CommonData/Prcp_GMET/JRA-55-Punzip',... % precipitation path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/JRA-55-Punzip',...  % snowfall path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/JRA-55-Tunzip',...  % Tmin
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/JRA-55-Tunzip',...  % Tmax
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/JRA-55-Tunzip'};    % Tmean
FileDEMLow='../DEM_process/JRA55_DEM2.mat';
FileDEMHigh='/home/gut428/GMET/NA_basic/MERIT_DEM/NA_DEM_010deg_trim.asc';
TLRfile='../MERRA2_TLR/MERRA2_TLR.mat';

Outpath={'/home/gut428/JRA55_YearNC',...
    '/home/gut428/JRA55_YearNC',...
    '/home/gut428/JRA55_YearNC',...
    '/home/gut428/JRA55_YearNC',...
    '/home/gut428/JRA55_YearNC'};
for i=1:length(Outpath)
   if ~exist(Outpath{i},'dir')
      mkdir(Outpath{i}); 
   end
end

% basic information of target region
BasicInfo.tXll=-50;  % top right
BasicInfo.tYll=85;
BasicInfo.Xll=-180;   % bottom left
BasicInfo.Yll=5;
BasicInfo.cellsize=0.1;  % target resolution
% year=[1979,2018]; % start end year
DataInfo.prefix={'fcst_phy2m.061_tprat.reg_tl319.',...
    'fcst_phy2m.064_srweq.reg_tl319.',...
    'minmax_surf.016_tmin.reg_tl319.',...
    'minmax_surf.015_tmax.reg_tl319.',...
    'fcst_surf.011_tmp.reg_tl319.'}; % file name prefix
DataInfo.suffix={'.tang390256','.tang390256','.tang392048','.tang392048','.tang392049'}; % file name suffix
DataInfo.varname={'TPRAT_GDS4_SFC_ave3h','SRWEQ_GDS4_SFC_ave3h','TMIN_GDS4_HTGL','TMAX_GDS4_HTGL','TMP_GDS4_HTGL'}; % variable to read
DataInfo.imethod={'linear','linear','lapserate','lapserate','lapserate'}; % interpolation method
DataInfo.prefixout={'JRA55_prcp_','JRA55_snow_','JRA55_tmin_','JRA55_tmax_','JRA55_tmean_'}; % prefix for output file

parfor year=1979:2018
    f_JRA55_read(Inpath,FileDEMLow,FileDEMHigh,Outpath,year,BasicInfo,DataInfo,TLRfile);
end