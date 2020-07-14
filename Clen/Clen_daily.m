% calcualte the spatial correlation lengths for each season
clc;clear;close all
Infile='/datastore/GLOBALWATER/CommonData/GapFill_new/RawStnData/AllGauge_QC.nc4';
%Infile='/Users/localuser/Research/GapFill/AllGauge_QC.nc4';
vars={'prcp','tmean','trange'};
searchradius=[1000,2000,2000];
cctype='Pearson';
leastyear=20;
outfileCC=['Clen_',cctype,'_c0.mat'];
outfileFit=['Fit_',cctype,'_c0.mat'];


if ~exist(outfileCC,'file')
    
    LLE=ncread(Infile,'LLE');
    date=ncread(Infile,'date');
    month=floor(mod(date,10000)/100);
    
    gnum=size(LLE,1);
    CC=cell(length(vars),1);
    DIST=cell(length(vars),1);
    
    
    for vv=1:length(vars)
        radius=searchradius(vv);
        % read variables
        if strcmp(vars{vv},'prcp')
            datavv=ncread(Infile,'prcp');
        end
        if strcmp(vars{vv},'tmean')
            tmin=ncread(Infile,'tmin');
            tmax=ncread(Infile,'tmax');
            datavv=(tmin+tmax)/2;
        end
        if strcmp(vars{vv},'trange')
            tmin=ncread(Infile,'tmin');
            tmax=ncread(Infile,'tmax');
            datavv=abs(tmax-tmin);
        end
        
        ccv=cell(gnum-1,1);
        distv=cell(gnum-1,1);
        parfor i=1:gnum-1
            fprintf('Variable %s--Gauge i %d--%d\n',vars{vv},i,gnum);
            % extact qualified data to calculate cc with i
            llei=LLE(i+1:end,:);
            indexi=i+1:size(datavv,2);
            % according to distance
            disti=lldistkm(LLE(i,1),LLE(i,2),LLE(i+1:end,1),LLE(i+1:end,2));
            indno=disti>radius;
            indexi(indno)=[];
            llei(indno,:)=[];
            disti(indno)=[];
            if isempty(llei)
                continue;
            end
            % according to overlap data
            datavvi=datavv(:,indexi);
            temp=datavvi+repmat(datavv(:,i),1,size(datavvi,2));
            snum=sum(~isnan(temp),1);
            indno=snum<leastyear*365;
            llei(indno,:)=[];
            disti(indno)=[];
            if isempty(llei)    
               continue; 
            end

            cci=nan*zeros(size(llei,1),12);
            for j=1:size(llei,1)
                ccij=f_month_cc(datavv(:,i),datavvi(:,j),month,cctype);
                cci(j,:)=ccij;
            end
            ccv{i}=cci;
            distv{i}=disti;

        end
        
        % assemble ccv
        CCall=[];
        DISTall=[];
        for jj=1:length(ccv)
            ccjj=ccv{jj};
            distjj=distv{jj};
            nn=sum(isnan(ccjj),2);
            indnn=nn==12;
            ccjj(indnn,:)=[];
            distjj(indnn,:)=[];
            if ~isempty(ccjj)
                CCall=cat(1,CCall,ccjj);
                DISTall=cat(1,DISTall,distjj);
            end
        end
        CC{vv}=single(CCall);
        DIST{vv}=single(DISTall);
    end
    save(outfileCC,'CC','DIST','cctype','leastyear','radius','vars','-v7.3');
else
    load(outfileCC,'CC','DIST','cctype','leastyear','radius','vars');
end

% fit
if ~exist(outfileFit,'file')
    FIT=cell(length(CC),12);
    for vv=1:length(CC)
        CCv=CC{vv};
        DISTv=DIST{vv};
        for m=1:size(CCv,2)
            CCvm=CCv(:,m);
            DISTvm=DISTv;
            ind=isnan(CCvm)|isnan(DISTvm);
            CCvm(ind)=[];
            DISTvm(ind)=[];
            % Fit
            [xData, yData] = prepareCurveData( DISTvm, CCvm );

            % Set up fittype and options.
%             ft = fittype( 'a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
            ft = fittype( 'a*exp(-x/b)', 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            if strcmp(vars{vv},'prcp')
            opts.StartPoint = 150;
            else
            opts.StartPoint = 800;   
            end

            % Fit model to data.
            [fitresult, gof] = fit( xData, yData, ft, opts );
            FIT{vv,m}=fitresult;
        end
    end
    save(outfileFit,'FIT');
else
    load(outfileFit,'FIT');
end
