clc;clear;close

% Path in Plato
Inpath='/datastore/GLOBALWATER/CommonData/Prcp_GMET/GHCN-D/ghcnd_all'; % Path of original station files
Outpath='/home/gut428/GapFill/Data/ghcn-d';  % Path to store output files
Infile_inventory='/datastore/GLOBALWATER/CommonData/Prcp_GMET/GHCN-D/documents/ghcnd-inventory.txt'; % GHCN-D station inventory
Infile_station='/datastore/GLOBALWATER/CommonData/Prcp_GMET/GHCN-D/documents/ghcnd-stations.txt'; % GHCN-D station list
Infile_mask='/home/gut428/GMET/NA_basic/MERIT_DEM/NA_DEM_010deg_trim.asc';  % only extract stations within the mask extent
Outfile_station=[Outpath,'/GaugeValid.mat'];  % File to store all information of output stations

% Basic settings
Overwrite=1; % 1: overwrite files in Outpath. Otherwise: skip existing files in Outpath.
BasicInfo.period_range=[1979,2018]; % [start year, end year]
BasicInfo.period_len=[8,100]; % the lower/upper number of years that are within period_range
BasicInfo.VarOut={'Date(yyyymmdd)','Precipitation(mm)','Tmin(C)','Tmax(C)'};
BasicInfo.VarRead={'PRCP','TMIN','TMAX'}; % variables to be read
BasicInfo.scalefactor=[0.1, 0.1, 0.1]; % scale factor of original variables
BasicInfo.missingvalue=[-9999, -9999, -9999];

% two ways for to extract rain gauges
% (1) lat/lon extent
% SR.seflag=1;
% SR.lat_range=[0,90]; % latitude range
% SR.lon_range=[-180,0]; % longitude range
% (2) spatial mask
BasicInfo.seflag=2; % indicates which type of extent definition
BasicInfo.maskfile=Infile_mask;

% start reading
f_GHCN_read(Inpath,Infile_inventory,Infile_station,BasicInfo,Outpath,Outfile_station,Overwrite);
