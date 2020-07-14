% calculate average auto_corr and tp_corr from all months
clc;clear
Inpath='/datastore/GLOBALWATER/CommonData/EMDNA_new/GMET_OIinput';
auto_corr=zeros(40,12);
tp_corr=zeros(40,12);
for i=1:40
    fprintf('%d\n',i);
    for j=1:12
        file=[Inpath,'/reg_',num2str((i+1978)*100+j),'.nc'];
        temp=ncread(file,'auto_corr');
        auto_corr(i,j)=temp(1);
        temp=ncread(file,'tp_corr');
        tp_corr(i,j)=temp(1);
    end
end
auto_corr_month=nanmean(auto_corr,1);
tp_corr_month=nanmean(tp_corr,1);
save('monthly_auto_tp_corr.mat','auto_corr','tp_corr','auto_corr_month','tp_corr_month');