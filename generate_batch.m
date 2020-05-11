% generate sbatch scripts so that we can submit many jobs at the same time
clc;clear
% mkdir('./Sbatchscripts')

for y=1979:2018
    for m=1:12
        yyyymmim=y*100+m;
        outfile=['./Sbatchscripts/Plato_',num2str(yyyymmim),'.sh'];
        fidout=fopen(outfile,'w');
        fprintf(fidout,'#!/bin/bash\n');
        fprintf(fidout,['#SBATCH --job-name=PG_',num2str(yyyymmim),'\n']);
        fprintf(fidout,['#SBATCH --time=2-00:00:00\n']);
        fprintf(fidout,'#SBATCH --mem=20G\n');
        fprintf(fidout,'module load python/3.7.4\n');
        

        str1=num2str(y*10000+m*100+1);
        str2=num2str(y*10000+m*100+eomday(y,m));
        fprintf(fidout,['srun python main_CAI.py ',str1,' ',str2,'\n']);
%         fprintf(fidout,'rm *.out\n');
        fclose(fidout);
    end
end

