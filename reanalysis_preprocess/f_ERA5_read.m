function f_ERA5_read(Inpath,Outpath,year,DataInfo)
% read ERA5 data and conduct resolution transform
% ERA5 data is stored by month

% basic information of file characteristics and reading
prefix=DataInfo.prefix;
suffix=DataInfo.suffix;
varname=DataInfo.varname;
prefixout=DataInfo.prefixout;
varnum=length(varname);

for vv=1:varnum  % for each var, re-read the basic information
    % basic information of original ERA5 data. All data in this directory must
    % have the same information
    Infile=[Inpath{vv},'/',prefix{vv},num2str(year(1)),num2str(1,'%.2d'),suffix{vv},'.nc'];
    latitude=ncread(Infile,'latitude');
    longitude=ncread(Infile,'longitude');

    for yy=year(1):year(end)
        Outfile=[Outpath{vv},'/',prefixout{vv},num2str(yy),'.nc4'];
        if ~exist(Outfile,'file')
            fprintf('ERA5 Data in process. Var %d; Year %d--%d\n',vv,yy,year(end));
            % read, hour 2 day, unit conversion, interpolation
            daysyy=datenum(yy,12,31)-datenum(yy,1,1)+1;
            data=zeros(length(latitude),length(longitude),daysyy);
            flag=1;
            for mm=1:12
                Infile=[Inpath{vv},'/',prefix{vv},num2str(yy),num2str(mm,'%.2d'),suffix{vv},'.nc'];
                if yy==1979&&mm==1&&ismember(varname{vv},{'tp'})
                    varmd=f_ERA5_VarRead(Infile,varname{vv},1);
                else
                    varmd=f_ERA5_VarRead(Infile,varname{vv},0); 
                end

                daysmm=size(varmd,3);
                data(:,:,flag:flag+daysmm-1)=varmd;
                flag=flag+daysmm;
            end
            % save data
            save(Outfile,'data','latitude','longitude','-v7.3');
%             f_save_nc(Outfile,data,BasicInfo,DEMHigh);
        else
            fprintf('ERA5 Data already exist. Var %d; Year %d--%d\n',vv,yy,year(end));
        end
    end
end

end

function vardd=f_ERA5_VarRead(Infile,varname,flag)
varhh=ncread(Infile,varname); % hourly data
varhh=permute(varhh,[2,1,3]); % to the normal lat/lon map
if flag==1
    add=nan*zeros(size(varhh,1),size(varhh,2),7);
    varhh=cat(3,add,varhh);
end

switch varname
    case 'tp'
        method='sum';
    case 'sf'
        method='sum';
    case 'mn2t'
        method='min';
    case 'mx2t'
        method='max';
    case 't2m'
        method='mean';
end
vardd=f_h2d(varhh,method);

% unit conversion
switch varname
    case {'tp','sf'}
        vardd=vardd*1000; % m to mm
    case {'mn2t','mx2t','t2m'}
        vardd=vardd-273.15; % K to C
end
end

function dout=f_h2d(din,method)
% hourly data to daily data by method--mean min max sum
% unit conversion
days=size(din,3)/24;
dout=nan*zeros(size(din,1),size(din,2),days);

for dd=1:days
    vd=din(:,:,dd*24-23:dd*24); % hourly to daily
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