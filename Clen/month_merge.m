clc;clear;close all

for i=1:12
    file=['/Users/localuser/Downloads/Fit_Spearman_prcp_',num2str(i),'.mat'];
    load(file,'FIT1','FIT2');
    fit1_clen(i,1)=FIT1.b;
    fit2_c0(i,1)=FIT2.a;
    fit2_clen(i,1)=FIT2.b;
end

save('Clen_prcp_Spearman.mat','fit1_clen','fit2_c0','fit2_clen');