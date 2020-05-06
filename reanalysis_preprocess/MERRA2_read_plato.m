clc;clear;close all

% Plato path
Inpath={'/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2',... % precipitation path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2T',...  % Tmin
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/MERRA2T'};    % Tmax
FileDEMLow='../DEM_process/MERRA2_DEM2.mat';


Outpath={'/home/gut428/MERRA2_day_raw',...
    '/home/gut428/MERRA2_day_raw',...
    '/home/gut428/MERRA2_day_raw'};
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
DataInfo.suffix={'.SUB','',''}; % file name suffix
DataInfo.varname={'PRECTOT','T2MMIN','T2MMAX'}; % variable to read
DataInfo.prefixout={'MERRA2_prcp_','MERRA2_tmin_','MERRA2_tmax_'}; % prefix for output file

parfor year=1980:2018
    f_MERRA2_read(Inpath,Outpath,year,DataInfo);
end
