function f_JRA55_read(Inpath,FileDEMLow,FileDEMHigh,Outpath,year,BasicInfo,DataInfo,TLRfile)
% precipitation unit: mm/day
% basic information of file characteristics and reading
% dbstop if error
prefix=DataInfo.prefix;
suffix=DataInfo.suffix;
varname=DataInfo.varname;
imethod=DataInfo.imethod;
prefixout=DataInfo.prefixout;
varnum=length(varname);

% basic information of target region
tXll=BasicInfo.tXll;  % top right
if tXll<0
   tXll=tXll+360;
end
tYll=BasicInfo.tYll;
Xll=BasicInfo.Xll;   % bottom left
if Xll<0
   Xll=Xll+360;
end
Yll=BasicInfo.Yll;
cellsize=BasicInfo.cellsize;  % new resolution
Ncol=(tXll-Xll)/cellsize;
Nrow=(tYll-Yll)/cellsize;
X1=(Xll+cellsize/2):cellsize:(Xll+cellsize*Ncol-cellsize/2); % lat/lon of grid centers
Y1=(Yll+cellsize*Nrow-cellsize/2):-cellsize:(Yll+cellsize/2);
[XX1,YY1]=meshgrid(X1,Y1);

% basic variables to be saved in the output file
varsav={'tXll','tYll','Xll','Yll','Nrow','Ncol','cellsize','varname','imethod'};
for i=1:length(varsav)
    command=['info.',varsav{i},'=',varsav{i},';'];
    eval(command);
end

% read temperature lapse rate information
[TLR,TLRdate]=f_TLR_process(TLRfile,XX1,YY1);


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
    latori0=ncread(file,'g4_lat_2'); % lat/lon of grid centers
    lonori0=ncread(file,'g4_lon_3');
    latori0=sort(latori0,'descend'); lonori0=sort(lonori0,'ascend');
    [XX00,YY00]=meshgrid(lonori0,latori0);
    %%%% since JRA has irregular lat/lon, we interpolate JRA into regular grid
    latori=(90-0.5625/2):-0.5625:(-90+0.5625/2);
    lonori=0:0.5625:(360-0.5625);
    [XX0,YY0]=meshgrid(lonori,latori);
    REAinfo.Xsize=abs(lonori(2)-lonori(1));
    REAinfo.Ysize=abs(latori(2)-latori(1));
    REAinfo.xll=min(lonori)-REAinfo.Xsize/2;
    REAinfo.yll=min(latori)-REAinfo.Ysize/2;
    REAinfo.nrows=size(XX0,1);
    REAinfo.ncols=size(XX0,2);
    
    % if there exists lapserate in imethod, calcualte temperature changes in
    % each reanalysis grid pixel
    if ismember('lapserate',imethod)
        % read DEM data
        % DEMLow must have larger or equal spatial extent compared with DEMhigh,
        % otherwise it is hard to downscale
        % DEMhigh must have the same spatial extent with BasicInfo
        load(FileDEMLow,'DEM','Info');
        DEMLow=DEM;
        InfoLow=Info;
        if InfoLow.xll<0
            InfoLow.xll=InfoLow.xll+360;  % -180 180 to 0 360
        end
        
        clear DEM Info
        mm=arcgridread_tgq(FileDEMHigh);
        DEMHigh=mm.mask;
        clear mm
        % interpolate DEMLow to match reanalysis. Theoretically, DEMLow should
        % totally match reanalysis. But sometimes due to marginal differences,
        % interpolation is necessary.
        latLow=(InfoLow.yll+InfoLow.Ysize*InfoLow.nrows-InfoLow.Ysize/2):-InfoLow.Ysize:(InfoLow.yll+InfoLow.Ysize/2);
        lonLow=(InfoLow.xll+InfoLow.Xsize/2):InfoLow.Xsize:(InfoLow.xll+InfoLow.Xsize*InfoLow.ncols-InfoLow.Xsize/2);
        [XXLow,YYLow]=meshgrid(lonLow,latLow);
        DEMLow=interp2(XXLow,YYLow,DEMLow,XX0,YY0,'linear');
        
    end
    
    method=imethod{vv};
    
    for yy=year(1):year(end)
        Outfile=[Outpath{vv},'/',prefixout{vv},num2str(yy),'.nc4'];
        if ~exist(Outfile,'file')
            fprintf('JRA55 Data in process. Var %d; Year %d--%d\n',vv,yy,year(end));
            data=[];
            indyy=find(yearin==yy);
            
            for ii=1:length(indyy)
                prefixii=Indir(indyy(ii)).name(1:length(prefix{vv}));
                if strcmp(prefixii,prefix{vv})                
                    Infile=[Inpath{vv},'/',Indir(indyy(ii)).name];
                    [varint0,varmonth]=f_JRA55_VarRead(Infile,varname{vv});
                    varint1=ff_interpolate(varint0,XX00,YY00,XX0,YY0,'linear');

                    if strcmp(method,'lapserate')
                        varint2=ff_interpolate(varint1,XX0,YY0,XX1,YY1,'near');
                        % find the lapse rate for this month
                        monthU=unique(varmonth);
                        varint=nan*zeros(Nrow,Ncol,size(varint2,3));
                        for uu=1:length(monthU)
                            TLRym=TLR(:,:,TLRdate==yy*100+monthU(uu));
                            Tadd=ff_Tdownscale_lp(XX1,YY1,REAinfo,DEMHigh,DEMLow,TLRym);
                            induu=varmonth==monthU(uu);
                            varint(:,:,induu)=ff_Tdownscale_add(varint2(:,:,induu),Tadd);
                        end
                        varint(isnan(varint))=varint2(isnan(varint));  % over pixels without dem
                    else
                        varint=ff_interpolate(varint1,XX0,YY0,XX1,YY1,method);
                    end
                    data=cat(3,data,varint);
                end
            end
            
%             save(Outfile,'data','BasicInfo','REAinfo','-v7.3');
            f_save_nc(Outfile,data,BasicInfo,DEMHigh);
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
