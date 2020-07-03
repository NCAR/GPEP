% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear
% allmode={'daily', 'climo', 'daily_anom'};
allmode={'daily'};
for am=1:length(allmode)
    TIME_MODE=allmode{am}; %daily, climo, daily_anom
    
    plato_outpath='/home/gut428/scratch/GMET/EMDNA_out';
    Path_exe='/home/gut428/scratch/GMET/EMDNA_out';
    Path_stngrid='/home/gut428/scratch/GMET/StnGridInfo';
    path_oi='/home/gut428/scratch/GMET/GMET_OIinput';
    gridinfo='gridinfo_whole.nc';
    
    if strcmp(TIME_MODE,'climo')
        name_mode='climo';
    elseif strcmp(TIME_MODE,'daily_anom')
        name_mode='anom';
    elseif strcmp(TIME_MODE,'daily')
        name_mode='daily';
    end
    Path_script=['/Users/localuser/Downloads/scripts/run_ens_',name_mode];
    
    
    year=1979:2018;
    % for downscale
    MAX_DISTANCE=100;
    % for ensemble
    if strcmp(TIME_MODE,'climo')
        hours=ones(1,12);
    else
        hours=15*ones(1,12);
    end
    % clen=[466,467,435,403,362,325,287,278,362,445,461,477];
    clen=ones(12,1)*150;
    start_ens=1;
    stop_ens=100;
    
    if ~exist(Path_script,'dir')
        mkdir(Path_script);
    end
    date=datenum(year(1),1,1):datenum(year(end),12,31);
    date=datestr(date','yyyymmdd');
    date=str2double(num2cell(date,2));
    yyyymm=floor(date/100);
    yyyy=floor(date/10000);
    mm=floor(mod(date,10000)/100);
    dd=mod(date,100);
    
    %% 3. configure files of ensemble
%     outfixd=['regress_',name_mode,'_'];
%     prefixe='config.ensemble.';
%     outfixe1=['ensemble_',name_mode,'_'];
%     outfixe2='';
    outfixd=['reg_'];
    prefixe='config.ensemble.';
    outfixe1=['ens_'];
    outfixe2='';
    for i=1:length(year)
        plato_outpathi=[plato_outpath,'/',num2str(year(i))];
        for m=1:12
            yyyymmim=year(i)*100+m;
            ind=yyyymm==yyyymmim;
            dateim=date(ind);
            
            if strcmp(TIME_MODE,'daily_anom') || strcmp(TIME_MODE,'daily')
                ntimes=length(dateim);
            elseif strcmp(TIME_MODE,'climo')
                ntimes=1;
            end
            
            outfile=[Path_script,'/',prefixe,num2str(yyyymmim),'.txt'];
            fid=fopen(outfile,'w');
            fprintf(fid,'&PARAMS\n');
            fprintf(fid,['start_time	= ',num2str(1),'\n']);
            fprintf(fid,['ntimes		= ',num2str(ntimes),'\n']);
            fprintf(fid,['start_ens	= ',num2str(start_ens),'\n']);
            fprintf(fid,['stop_ens	= ',num2str(stop_ens),'\n']);
            fprintf(fid,['clen		= ',num2str(clen(m)),'\n']);
            fprintf(fid,['grid_name	= "',Path_stngrid,'/',gridinfo,'"\n']);
            if strcmp(TIME_MODE,'daily_anom')
                fprintf(fid,['climo_path	= "',plato_outpathi,'"\n']);
            end
            fprintf(fid,['out_forc_name_base	= "',plato_outpathi,'/',outfixe1,num2str(yyyymmim),outfixe2,'"\n']);
            fprintf(fid,['in_regr_name	= "',path_oi,'/',outfixd,num2str(yyyymmim),'.nc"\n']);
            fprintf(fid,['time_mode       = "',TIME_MODE,'"\n']); % different time mode
            fprintf(fid,'/\n');
            fclose(fid);
        end
    end
    
    %% 4. Plato sbatch for ensemble
    for i=1:length(year)
        plato_outpathi=[plato_outpath,'/',num2str(year(i))];
        for m=1:12
            yyyymmim=year(i)*100+m;
            outfile=[Path_script,'/gra_ens_',num2str(yyyymmim),'.sh'];
            fidout=fopen(outfile,'w');
            fprintf(fidout,'#!/bin/bash\n');
            fprintf(fidout,['#SBATCH --job-name=Ens_',num2str(yyyymmim),'\n']);
            fprintf(fidout,'#SBATCH --account=rpp-kshook\n');
            fprintf(fidout,['#SBATCH --time=0-',num2str(hours(m)),':00:00\n']);
            fprintf(fidout,'#SBATCH --mem=30G\n');
            fprintf(fidout,['mkdir -p ',plato_outpathi,'\n']);
            
            exefile=[Path_exe,'/generate_ensemble.exe'];
            fprintf(fidout,['chmod a+x ',exefile,'\n']);
            fprintf(fidout,[exefile,' ',prefixe,num2str(yyyymmim),'.txt\n']);
%             fprintf(fidout,'rm *.out\n');
            fclose(fidout);
        end
    end
    
end