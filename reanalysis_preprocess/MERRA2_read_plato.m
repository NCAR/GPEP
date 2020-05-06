clc;clear;close all

% Plato path
Inpath={'/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2',... % precipitation path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2',...  % snowfall path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2T',...  % Tmin
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2T',...  % Tmax
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2T'};    % Tmean
FileDEMLow='../DEM_process/MERRA2_DEM2.mat';
FileDEMHigh='/home/gut428/GMET/NA_basic/MERIT_DEM/NA_DEM_010deg_trim.asc';
TLRfile='../MERRA2_TLR/MERRA2_TLR.mat';

Outpath={'/home/gut428/MERRA2_YearNC',...
    '/home/gut428/MERRA2_YearNC',...
    '/home/gut428/MERRA2_YearNC',...
    '/home/gut428/MERRA2_YearNC',...
    '/home/gut428/MERRA2_YearNC'};
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
% year=[1980,2018]; % start end year
DataInfo.suffix={'.SUB','.SUB','','',''}; % file name suffix
DataInfo.varname={'PRECTOT','PRECSNO','T2MMIN','T2MMAX','T2MMEAN'}; % variable to read
DataInfo.imethod={'linear','linear','lapserate','lapserate','lapserate'}; % interpolation method
DataInfo.prefixout={'MERRA2_prcp_','MERRA2_snow_','MERRA2_tmin_','MERRA2_tmax_','MERRA2_tmean_'}; % prefix for output file

parfor year=1980:2018
    f_MERRA2_read2(Inpath,FileDEMLow,FileDEMHigh,Outpath,Outpath2,year,BasicInfo,DataInfo,TLRfile);
end
