function f_JRA55_read(Inpath,Outpath,year,DataInfo)
% precipitation unit: mm/day
% basic information of file characteristics and reading
% dbstop if error
prefix=DataInfo.prefix;
suffix=DataInfo.suffix;
varname=DataInfo.varname;
prefixout=DataInfo.prefixout;
varnum=length(varname);

for vv=1:varnum  % for each var, re-read the basic information
    % specific for MERRA2 since its date could change after 2014
    Indir=dir(fullfile(Inpath{vv},'*.nc'));
    if isempty(Indir)
        Indir=dir(fullfile(Inpath{vv},'*.nc4'));
    end
    yearin=cell(length(Indir),1);
    lensuff=length(suffix{vv})+9;
    for i=1:length(Indir)
        yearin{i}=Indir(i).name(end-lensuff-3:end-lensuff);
    end
    yearin=str2double(yearin);
    
    % basic information of original JRA55 data
    tempind=find(yearin==year(1)); tempind=tempind(1);
    file=[Inpath{vv},'/',Indir(tempind).name];
    latitude=ncread(file,'g4_lat_2'); % lat/lon of grid centers
    longitude=ncread(file,'g4_lon_3');
    longitude=mod(longitude+180,360)-180; % 0-360 to -180/180
    
    indlat = latitude>=0 & latitude<=90;
    indlon = longitude>=-180 & longitude<=-50;
    latitude=latitude(indlat);
    longitude=longitude(indlon);
    
    for yy=year(1):year(end)
        Outfile=[Outpath{vv},'/',prefixout{vv},num2str(yy),'.mat'];
        if ~exist(Outfile,'file')
            fprintf('JRA55 Data in process. Var %d; Year %d--%d\n',vv,yy,year(end));
            data=[];
            indyy=find(yearin==yy);
            
            for ii=1:length(indyy)
                prefixii=Indir(indyy(ii)).name(1:length(prefix{vv}));
                if strcmp(prefixii,prefix{vv})                
                    Infile=[Inpath{vv},'/',Indir(indyy(ii)).name];
                    [varint0,varmonth]=f_JRA55_VarRead(Infile,varname{vv});
                    
                    varint0=varint0(indlat,indlon,:);
                    data=cat(3,data,varint0);
                end
            end
            
            data=single(data);
            save(Outfile,'data','latitude','longitude','-v7.3');
%             f_save_nc(Outfile,data,BasicInfo,DEMHigh);
        else
            fprintf('JRA55 Data already exist. Var %d; Year %d--%d\n',vv,yy,year(end));
        end
    end
end
end

function [vardd,month]=f_JRA55_VarRead(Infile,varname)
var3h=ncread(Infile,varname); % hourly data
var3h=reshape(var3h,[size(var3h,1),size(var3h,2),size(var3h,3)*size(var3h,4)]);
var3h=permute(var3h,[2,1,3]);
initial_time0=ncread(Infile,'initial_time0');
initial_time0=initial_time0';

switch varname
    case 'TPRAT_GDS4_SFC_ave3h'
        method='mean';
    case 'SRWEQ_GDS4_SFC_ave3h'
        method='mean';
    case 'TMIN_GDS4_HTGL'
        method='min';
        var3h(:,:,1:2)=[];
        initial_time0(1,:)=[];
    case 'TMAX_GDS4_HTGL'
        method='max';
        var3h(:,:,1:2)=[];
        initial_time0(1,:)=[];
    case 'TMP_GDS4_HTGL'
        method='mean';
        var3h(:,:,1:2)=[];
        initial_time0(1,:)=[];
end
vardd=f_3h2d(var3h,method);

% month
month=nan*zeros(size(initial_time0,1)/4,1);
for i=1:length(month)
   month(i)=str2double(initial_time0(i*4,1:2)); 
end

% unit conversion
switch varname
    case {'TPRAT_GDS4_SFC_ave3h','SRWEQ_GDS4_SFC_ave3h'}
        vardd=vardd; % mm/day
    case {'TMIN_GDS4_HTGL','TMAX_GDS4_HTGL','TMP_GDS4_HTGL'}
        vardd=vardd-273.15; % K to C
end
end

function dout=f_3h2d(din,method)
% hourly data to daily data by method--mean min max sum
% unit conversion
days=size(din,3)/8;
dout=nan*zeros(size(din,1),size(din,2),days);

for dd=1:days
    vd=din(:,:,dd*8-7:dd*8); % hourly to daily
    switch method
        case 'mean'
            vd=nanmean(vd,3);
        case 'min'
            vd=min(vd,[],3);
        case 'max'
            vd=max(vd,[],3);
        case 'sum'
            vd=nansum(vd,3);
    end
    dout(:,:,dd)=vd;
end
end
