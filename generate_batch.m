% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear
% mkdir('./Sbatchscripts')

% for y=1979:2018
%     for m=1:12
%         yyyymmim=y*100+m;
%         outfile=['Plato_reg_',num2str(yyyymmim),'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=PG_',num2str(yyyymmim),'\n']);
%         fprintf(fidout,['#SBATCH --time=0-2:00:00\n']);
%         fprintf(fidout,'#SBATCH --mem=10G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
% 
%         str1=num2str(y*10000+m*100+1);
%         str2=num2str(y*10000+m*100+eomday(y,m));
%         fprintf(fidout,['srun python -u main_daily_run.py ',str1,' ',str2,'\n']);
% %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%     end
% end

for y=1979:2018
    outfile=['Plato_reg_',num2str(y),'.sh'];
    fidout=fopen(outfile,'w');
    fprintf(fidout,'#!/bin/bash\n');
    fprintf(fidout,['#SBATCH --job-name=PG_',num2str(y),'\n']);
    fprintf(fidout,['#SBATCH --time=1-0:00:00\n']);
    fprintf(fidout,'#SBATCH --mem=15G\n');
    fprintf(fidout,'module load python/3.7.4\n');

    fprintf(fidout,['srun python -u s3_stn_regression.py ',num2str(y),'\n']);
    %         fprintf(fidout,'rm *.out\n');
    fclose(fidout);
end


% for y=1979:2018
%     for m=1:12
%         yyyymmim=y*100+m;
%         outfile=['Plato_',num2str(yyyymmim),'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=reapop\n']);
%         fprintf(fidout,['#SBATCH --time=2-0:00:00\n']);
%         fprintf(fidout,'#SBATCH --mem=15G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
% 
% 
%         str1=num2str(y*10000+m*100+1);
%         str2=num2str(y*10000+m*100+eomday(y,m));
% %         fprintf(fidout,['srun python -u main_daily_run.py ',str1,' ',str2,'\n']);
%         fprintf(fidout,['srun python -u reanalysis_pop.py ',num2str(y),' ',num2str(m),'\n']);
% %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%     end
% end


% for y=1979:2:2018
%     outfile=['Plato_',num2str(y),'.sh'];
%     fidout=fopen(outfile,'w');
%     fprintf(fidout,'#!/bin/bash\n');
%     fprintf(fidout,['#SBATCH --job-name=y2m\n']);
%     fprintf(fidout,['#SBATCH --time=0-2:00:00\n']);
%     fprintf(fidout,'#SBATCH --mem=5G\n');
%     fprintf(fidout,'module load python/3.7.4\n');
% 
%     fprintf(fidout,['srun python -u year_2_month.py ',num2str(y),' ',num2str(y+1),'\n']);
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
% flag=1;
% for i=1:3
%     for y=1979:1:2018
%         stri=[var{i},'_',num2str(y)];
%         outfile=['Plato_',stri,'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=corrmerge','\n']);
%         fprintf(fidout,['#SBATCH --time=0-3:0:0\n']);
%         fprintf(fidout,'#SBATCH --mem=35G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
%         
%         stri=[var{i},' BMA QM ',num2str(y),' ',num2str(y)];
%         fprintf(fidout,['srun python -u reanalysis_correction_merge.py ',stri,'\n']);
%         %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%     end
% end

% vars={'prcp','tmean','trange','pop'};
% flag=1;
% for v=1:4
%     for m=1:12
%         outfile=['Plato_',num2str(flag),'.sh'];
%         fidout=fopen(outfile,'w');
%         fprintf(fidout,'#!/bin/bash\n');
%         fprintf(fidout,['#SBATCH --job-name=',vars{v},num2str(m),'\n']);
%         fprintf(fidout,['#SBATCH --time=1-00:00:00\n']);
%         fprintf(fidout,'#SBATCH --mem=35G\n');
%         fprintf(fidout,'module load python/3.7.4\n');
%         fprintf(fidout,['srun python -u observation_reanalysis_merge.py ',vars{v},' ',num2str(m),'\n']);
%         %         fprintf(fidout,'rm *.out\n');
%         fclose(fidout);
%         flag=flag+1;
%     end
% end


% time=[1:100:14610,14610];
% for i=1:length(time)-1
%     outfile=['Plato_',num2str(i),'.sh'];
%     fidout=fopen(outfile,'w');
%     fprintf(fidout,'#!/bin/bash\n');
%     fprintf(fidout,['#SBATCH --job-name=pop',num2str(i),'\n']);
%     fprintf(fidout,['#SBATCH --time=0-14:00:00\n']);
%     fprintf(fidout,'#SBATCH --mem=15G\n');
%     fprintf(fidout,'module load python/3.7.4\n');
%     fprintf(fidout,['srun python -u temprun.py ',num2str(time(i)),' ',num2str(time(i+1)),'\n']);
%     %         fprintf(fidout,'rm *.out\n');
%     fclose(fidout);
% end


