function f_MERRA2_read(Inpath,FileDEMLow,FileDEMHigh,Outpath,year,BasicInfo,DataInfo,TLRfile)
% read MERRA2 data and conduct resolution transform

% basic information of file characteristics and reading
suffix=DataInfo.suffix;
varname=DataInfo.varname;
imethod=DataInfo.imethod;
prefixout=DataInfo.prefixout;
varnum=length(varname);

% basic information of target region
tXll=BasicInfo.tXll;  % top right
tYll=BasicInfo.tYll;
Xll=BasicInfo.Xll;   % bottom left
Yll=BasicInfo.Yll;
cellsize=BasicInfo.cellsize;  % new resolution
Ncol=(tXll-Xll)/cellsize;
Nrow=(tYll-Yll)/cellsize;
X1=(Xll+cellsize/2):cellsize:(Xll+cellsize*Ncol-cellsize/2); % lat/lon of grid centers
Y1=(Yll+cellsize*Nrow-cellsize/2):-cellsize:(Yll+cellsize/2);
[XX1,YY1]=meshgrid(X1,Y1);

% varables for saving into the structure (info)
varsav={'tXll','tYll','Xll','Yll','Nrow','Ncol','cellsize','varname','imethod'};
for i=1:length(varsav)
    command=['info.',varsav{i},'=',varsav{i},';'];
    eval(command);
end

% read temperature lapse rate information
[TLR,TLRdate]=f_TLR_process(TLRfile,XX1,YY1);


for vv=1:varnum  % for each var, re-read the basic information
    
    % specific for MERRA2 since its prefix changes over time
    Indir=dir(fullfile(Inpath{vv},'*.nc'));
    lensuff=length(suffix{vv})+3;
    if isempty(Indir)
        Indir=dir(fullfile(Inpath{vv},'*.nc4'));
        lensuff=length(suffix{vv})+4;
    end
    datein=cell(length(Indir),1);
    for i=1:length(Indir)
        datein{i}=Indir(i).name(end-lensuff-7:end-lensuff);
    end
    
    % basic information of original MERRA2 data
    Infile=[Inpath{vv},'/',Indir(1).name];
    latori=ncread(Infile,'lat');
    lonori=ncread(Infile,'lon');
    latori=sort(latori,'descend'); lonori=sort(lonori,'ascend');
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
            fprintf('MERRA2 Data in process. Var %d; Year %d--%d\n',vv,yy,year(end));
            
            Startdate=datenum(yy,1,1);  %yyyymmdd
            Enddate=datenum(yy,12,31);
            daysyy=Enddate-Startdate+1;
            data=zeros(Nrow,Ncol,daysyy);
            flagmissing=0;
            for i=Startdate:Enddate
                datei=datestr(i,'yyyymmdd');
                [~,z]=ismember(datei,datein);
                if z==0
                    flagmissing=flagmissing+1;
                    if flagmissing>=5
                       error('Missing days in a year >=5'); 
                    end
                    continue;
                end
                Infile=[Inpath{vv},'/',Indir(z).name];
                var=f_MERRA2_VarRead(Infile,varname{vv});
                if strcmp(method,'lapserate')
                    varint0=ff_interpolate(var,XX0,YY0,XX1,YY1,'near');
                    
                    % find the lapse rate for this month
                    mm=str2double(datei(5:6));
                    TLRym=TLR(:,:,TLRdate==yy*100+mm);
                    Tadd=ff_Tdownscale_lp(XX1,YY1,REAinfo,DEMHigh,DEMLow,TLRym);
                    
                    varint=ff_Tdownscale_add(varint0,Tadd);
                    varint(isnan(varint))=varint0(isnan(varint));  % over pixels without dem
                else
                    varint=ff_interpolate(var,XX0,YY0,XX1,YY1,method);
                end
                data(:,:,i-Startdate+1)=varint;
            end         
            % save data
%             save(Outfile,'data','BasicInfo','REAinfo','-v7.3');
            f_save_nc(Outfile,data,BasicInfo,DEMHigh);
        else
            fprintf('MERRA2 Data already exist. Var %d; Year %d--%d\n',vv,yy,year(end));
        end
    end
end

end

function var=f_MERRA2_VarRead(Infile,varname)
var=ncread(Infile,varname); % hourly data
var=flipud(var'); % to the normal lat/lon map
% unit conversion
switch varname
    case {'PRECTOT','PRECSNO'}
        var=var*3600*24; % mm/day
    case {'T2MMIN','T2MMAX','T2MMEAN'}
        var=var-273.15; % F to C
end
end

