% generate scripts for batch running
% clc;clear
vars={'prcp','tmean','trange'};

for i=1:1
    for j=1:12
        outfile=['Plato_',vars{i},num2str(j),'.sh'];
        fidout=fopen(outfile,'w');
        fprintf(fidout,'#!/bin/bash\n');
        fprintf(fidout,['#SBATCH --job-name=clen\n']);
        fprintf(fidout,'#SBATCH --time=0-8:00:00\n');
        fprintf(fidout,'#SBATCH --mem=20G\n');
        fprintf(fidout,'module load matlab/R2017b\n');
        fprintf(fidout,['./Clen ',vars{i},' ',num2str(j)]);
        fclose(fidout);
    end
end
