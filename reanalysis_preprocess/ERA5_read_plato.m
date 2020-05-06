clc;clear;close

% Plato path
Inpath={'/datastore/GLOBALWATER/CommonData/Prcp_GMET/ERA5/tp',... % precipitation path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/ERA5TSF',...  % snowfall path
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/ERA5TSF',...  % Tmin
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/ERA5TSF',...  % Tmax
    '/datastore/GLOBALWATER/CommonData/Prcp_GMET/ERA5TSF'};    % Tmean

Outpath={'/home/gut428/ERA5_day_raw',...
    '/home/gut428/ERA5_day_raw',...
    '/home/gut428/ERA5_day_raw',...
    '/home/gut428/ERA5_day_raw',...
    '/home/gut428/ERA5_day_raw'};
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
% year=[1979,2018]; % start and end year
DataInfo.prefix={'ERA5_','ERA5_','ERA5_'}; % file name prefix
DataInfo.suffix={'','','','',''}; % file name suffix
DataInfo.varname={'tp','mn2t','mx2t'}; % variable to read
DataInfo.prefixout={'ERA5_prcp_','ERA5_tmin_','ERA5_tmax_'}; % prefix for output file

% parfor year=1980:2018
for year=1980:2018
    f_ERA5_read(Inpath,Outpath,year,DataInfo);
end
