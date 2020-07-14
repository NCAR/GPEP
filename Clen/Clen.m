% calcualte the spatial correlation lengths for each season
function Clen(vars, monuse)
Infile='/datastore/GLOBALWATER/CommonData/GapFill_new/RawStnData/AllGauge_QC.nc4';
%Infile='/Users/localuser/Research/GapFill/AllGauge_QC.nc4';
% vars={'prcp','tmean','trange'};

if ischar(monuse)
    monuse=str2double(monuse);
end

cctype='Pearson';
leastyear=20;
outfileCC=['CC_',cctype,'_',vars,'_',num2str(monuse),'.mat'];
outfileFit=['Fit_',cctype,'.mat'];
if strcmp(vars,'prcp')
    searchradius=1000;
elseif strcmp(vars,'tmean')
    searchradius=2000;
elseif strcmp(vars,'trange')
    searchradius=1500;
end


if ~exist(outfileCC,'file')
    LLE=ncread(Infile,'LLE');
    date=ncread(Infile,'date');
    month=floor(mod(date,10000)/100);
    % read variables
    if strcmp(vars,'prcp')
        datavv=ncread(Infile,'prcp');
    end
    if strcmp(vars,'tmean')
        tmin=ncread(Infile,'tmin');
        tmax=ncread(Infile,'tmax');
        datavv=(tmin+tmax)/2;
        clear tmin tmax
    end
    if strcmp(vars,'trange')
        tmin=ncread(Infile,'tmin');
        tmax=ncread(Infile,'tmax');
        datavv=abs(tmax-tmin);
        clear tmin tmax
    end
    datavv=datavv(month==monuse,:);
    
    % exclude stations that cannot satisfy requirement
    valnum=sum(~isnan(datavv),1);
    indno=valnum<leastyear*28;
    datavv(:,indno)=[];
    LLE(indno,:)=[];
    gnum=size(LLE,1);
    
    % start calculation
    ccg=cell(gnum-1,1);
    distg=cell(gnum-1,1);
    for i=1:gnum-1
        fprintf('Variable %s--Gauge i %d--%d\n',vars,i,gnum);
        % extact qualified data to calculate cc with i
        llei=LLE(i+1:end,:);
        indexi=i+1:size(datavv,2);
        % according to distance
        disti=lldistkm(LLE(i,1),LLE(i,2),LLE(i+1:end,1),LLE(i+1:end,2));
        indno=disti>searchradius;
        indexi(indno)=[];
        llei(indno,:)=[];
        disti(indno)=[];
        if isempty(llei)
            continue;
        end
        % according to overlap data
        datavvi=datavv(:,indexi);
        cci=single(nan*zeros(size(llei,1),1));
        for j=1:length(indexi)
            dataij=[datavvi(:,j),datavv(:,i)];
            dataij(isnan(dataij(:,1))|isnan(dataij(:,2)),:)=[];
            if size(dataij,1)>=leastyear*28
                ccij=corr(dataij(:,1),dataij(:,2),'Type',cctype);
                cci(j)=ccij;
            end
        end
        
        indno=isnan(cci);
        cci(indno)=[];
        disti(indno)=[];
        disti=single(disti);
        ccg{i}=cci;
        distg{i}=disti;
    end
    
    % assemble ccv
    CC=[];
    DIST=[];
    for jj=1:length(ccg)
        ccjj=ccg{jj};
        distjj=distg{jj};
        if ~isempty(ccjj)
            CC=cat(1,CC,ccjj);
            DIST=cat(1,DIST,distjj);
        end
    end
    CC=single(CC);
    DIST=single(DIST);
    save(outfileCC,'CC','DIST','cctype','leastyear','searchradius','vars','-v7.3');
else
    load(outfileCC,'CC','DIST','cctype','leastyear','searchradius','vars');
end

% fit
if ~exist(outfileFit,'file')
    ind=isnan(CC)|isnan(DIST);
    CC(ind)=[];
    DIST(ind)=[];
    % Fit
    [xData, yData] = prepareCurveData( DIST, CC );
    
    % Set up fittype and options.
    ft1 = fittype( 'exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    ft2 = fittype( 'a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    if strcmp(vars,'prcp')
        opts.StartPoint = 150;
    else
        opts.StartPoint = 800;
    end
    
    % Fit model to data.
    [FIT1, gof1] = fit( xData, yData, ft1, opts );
    [FIT2, gof2] = fit( xData, yData, ft2, opts );
    save(outfileFit,'FIT1','FIT2');
else
    load(outfileFit,'FIT1','FIT2');
end

end
