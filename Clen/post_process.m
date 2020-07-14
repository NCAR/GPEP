% load('Fit_Pearson.mat');
load('Fit_month_Pearson.mat');
Clen=zeros(3,12); %prcp, tmean, trange
for i=1:12
    for j=1:3
        Clen(j,i)=FIT{j,i}.b;
    end
end
