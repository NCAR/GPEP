load('/Users/localuser/Github/GMET_py/outputs.mat')
dd2=load('/Users/localuser/Github/GMET_py/outputs-old.mat');
close all
% file='~/GMET/eCAI/GMET_test_ecai-2/outputs/regress_ts.ecai.nc';
file='~/GMET/Example_tgq/outputs_dailyGMET/regress_daily.nc';
pop2=ncread(file,'pop');
pcp2=ncread(file,'pcp');
tmean2=ncread(file,'tmean');
trange2=ncread(file,'trange');
pop2=permute(pop2,[2,1,3]);
pcp2=permute(pcp2,[2,1,3]);
tmean2=permute(tmean2,[2,1,3]);
trange2=permute(trange2,[2,1,3]);

pcp=(pcp/4+1).^4;
pcp2=(pcp2/4+1).^4;

% % stn lat/lon to row/col
% stninfo(:,2)=floor( (41 - stninfo(:,2))/0.0625 ) + 1;
% stninfo(:,3)=floor( (stninfo(:,3) + 109.5)/0.0625 ) + 1;

% tmean
var='pop';
clims=[0,1];

figure('color','w')
pic_num = 1;
for day=1:31
% day=31;
pm2=flipud( dd2.pop(:,:,day) );
pm1=flipud( pop(:,:,day) );
pm0=pcp_stn(:,day);
if strcmp(var,'pop')
    pm0(pm0>0)=1;
end

% pm2=flipud( nanmean(pcp2,3) );
% pm1=flipud( nanmean(pcp,3) );
% if strcmp(var,'pop')
%     pm0=pcp_stn;
%     pm0(pcp_stn>0)=1;
%     pm0=nanmean(pm0,2);
% else
%     pm0=nanmean(pcp_stn,2);
% end

subplot(1,3,1)
scatter(stninfo(:,3),stninfo(:,2),50,pm0,'filled');
caxis(clims);
colormap(jet)
colorbar
title(['station ',var]);

subplot(1,3,2)
imagesc(pm1)
caxis(clims);
colormap(jet)
colorbar
title(['PyGMET ',var]);

subplot(1,3,3)
imagesc(pm2)
caxis(clims);
colormap(jet)
colorbar
title(['eCAI ',var]);

pause

% drawnow;
% F=getframe(gcf);
% I=frame2im(F);
% [I,map]=rgb2ind(I,256);
% if pic_num == 1
%     imwrite(I,map,'test.gif','gif', 'Loopcount',inf,'DelayTime',0.5);
% else
%     imwrite(I,map,'test.gif','gif','WriteMode','append','DelayTime',0.5);
% end
% pic_num = pic_num + 1;
end