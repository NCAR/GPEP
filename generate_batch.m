% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear
% mkdir('./Sbatchscripts')

% for y=2017:2018
%     for m=1:12
%         yyyymmim=y*100+m;
%         outfile=['Plato_',num2str(yyyymmim),'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=PG_',num2str(yyyymmim),'\n']);
%         fprintf(fidout,['#SBATCH --time=0-12:00:00\n']);
%         fprintf(fidout,'#SBATCH --mem=20G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
% 
% 
%         str1=num2str(y*100+m);
% %         str2=num2str(y*10000+m*100+eomday(y,m));
%         fprintf(fidout,['srun python -u main_CAI_update.py ',str1,'\n']);
% %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%     end
% end


% for y=1979:2018
%     outfile=['Plato_',num2str(y),'.sh'];
%     fidout=fopen(outfile,'w');
%     fprintf(fidout,'#!/bin/bash\n');
%     fprintf(fidout,['#SBATCH --job-name=err_',num2str(y),'\n']);
%     fprintf(fidout,['#SBATCH --time=0-5:00:00\n']);
%     fprintf(fidout,'#SBATCH --mem=10G\n');
%     fprintf(fidout,'module load python/3.7.4\n');
% 
%     fprintf(fidout,['srun python -u main_CAI_update.py ',num2str(y),'\n']);
%     %         fprintf(fidout,'rm *.out\n');
%     fclose(fidout);
% end

% var={'prcp','tmean','trange'};
% mode={'RMSE','BMA'};
%
% flag=1;
% for i=1:3
%     for j=1:2
%         outfile=['Plato_mergecorr_',num2str(flag),'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=mercorr',num2str(flag),'\n']);
%         fprintf(fidout,['#SBATCH --time=2-0:0:0\n']);
%         fprintf(fidout,'#SBATCH --mem=20G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
%         fprintf(fidout,['srun python -u reanalysis_correction_merge.py ',var{i},' ',mode{j},'\n']);
%         %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%
%         flag=flag+1;
%     end
% end

% var={'prcp','tmean','trange'};
var={'trange'};
mode={'RMSE','BMA'};
year=1979:1:2018;
flag=1;
for i=1:1
    for j=1:2
        for y=1:length(year)
            stri=[var{i},'_',mode{j},'_',num2str(year(y))];
            outfile=['Plato_',stri,'.sh'];
            fidout=fopen(outfile,'w');
            fprintf(fidout,'#!/bin/bash\n');
            fprintf(fidout,['#SBATCH --job-name=mercorr','\n']);
            fprintf(fidout,['#SBATCH --time=0-12:0:0\n']);
            fprintf(fidout,'#SBATCH --mem=30G\n');
            fprintf(fidout,'module load python/3.7.4\n');
            
            stri=[var{i},' ',mode{j},' ',num2str(year(y)),' ',num2str(year(y)+1)];
            fprintf(fidout,['srun python -u reanalysis_correction_merge.py ',stri,'\n']);
            %         fprintf(fidout,'rm *.out\n');
            fclose(fidout);
        end
    end
end
