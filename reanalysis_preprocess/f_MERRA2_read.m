function f_MERRA2_read(Inpath,Outpath,year,DataInfo)
% read MERRA2 data and conduct resolution transform

% basic information of file characteristics and reading
suffix=DataInfo.suffix;
varname=DataInfo.varname;
prefixout=DataInfo.prefixout;
varnum=length(varname);

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
    latitude=ncread(Infile,'lat');
    longitude=ncread(Infile,'lon');
  
    for yy=year(1):year(end)
        Outfile=[Outpath{vv},'/',prefixout{vv},num2str(yy),'.mat'];
        if ~exist(Outfile,'file')
            fprintf('MERRA2 Data in process. Var %d; Year %d--%d\n',vv,yy,year(end));
            
            Startdate=datenum(yy,1,1);  %yyyymmdd
            Enddate=datenum(yy,12,31);
            daysyy=Enddate-Startdate+1;
            data=zeros(length(latitude),length(longitude),daysyy);
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
                data(:,:,i-Startdate+1)=var;
            end         
            % save data
            data=single(data);
            save(Outfile,'data','latitude','longitude','-v7.3');
%             f_save_nc(Outfile,data,BasicInfo,DEMHigh);
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

