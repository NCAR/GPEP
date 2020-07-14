clc;clear;close all
load('Fit_1500_Pearson.mat');
for i=1:3
    for j=1:12
        Clen(i,j)=FIT{i,j}.b;
        C0(i,j)=FIT{i,j}.a;
    end
end

figure('color','w')
subplot(2,1,1)
plot(Clen(1,:),'-*');
xlabel('Month')
ylabel('Correlation length');
title('Precipitation');

subplot(2,1,2)
plot(Clen(2,:),'-*');
xlabel('Month')
ylabel('Correlation length');
title('Temperature');