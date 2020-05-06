clc;clear;close
% Path in Plato
Inpath='/datastore/GLOBALWATER/CommonData/Prcp_GMET/Environment Canada/ECCC_raw';
Infile_list='/datastore/GLOBALWATER/CommonData/Prcp_GMET/Environment Canada/Station Inventory EN.csv';
Outpath='/home/gut428/GapFill/Data/eccc';
Infile_mask='/home/gut428/GMET/NA_basic/MERIT_DEM/NA_DEM_010deg_trim.asc';
Outfile_gaugeAll=[Outpath,'/GaugeValid.mat'];

% Basic settings
ECCCinfo.period_range=[1979,2018]; % [start year, end year]
ECCCinfo.period_len=[8,100]; % the least/most number of years that are within period_range
ECCCinfo.VarRead={'TotalPrecip_mm_','MinTemp__C_','MaxTemp__C_'};
ECCCinfo.scalefactor=[1, 1, 1];
ECCCinfo.VarOut={'Date(yyyymmdd)','Precipitation(mm)','Tmin(C)','Tmax(C)'};

% two ways for to extract rain gauges
% (1) lat/lon extent
% SR.seflag=1;
% SR.lat_range=[0,90]; % latitude range
% SR.lon_range=[-180,0]; % longitude range
% (2) spatial mask
ECCCinfo.seflag=2;
ECCCinfo.maskfile=Infile_mask;

% start reading
f_ECCC_read(Inpath,Infile_list,Outpath,Outfile_gaugeAll,ECCCinfo);
