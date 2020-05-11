

close
shapefile='~/Google Drive/Datasets/Shapefiles/North America/From MERIT DEM/North_America.shp';
shp=shaperead(shapefile);

load testpoint
load testpoint_stnerror


subplot(3,1,1)

hold on
for q=1:length(shp)
ll1=shp(q).X; ll2=shp(q).Y; % delete Hawaii
plot(ll1,ll2,'-k');
end

colormap(jet)
scatter(x_red_use(:,3),x_red_use(:,2),50,y_tmean_red,'filled')
colorbar

plot(gridinfo_use(3),gridinfo_use(2),'r*','markersize',10)

hold off

xlim([-130,-50])
ylim([70,85])

subplot(3,1,2)

hold on
file='gridinfo_whole.nc';
dem=ncread(file,'elev');
dem=flipud(dem');
dem=dem(1:150,500:1300);
dem=flipud(dem);
imagesc(dem,'alphadata',~isnan(dem))
colorbar
colormap(jet)

subplot(3,1,3)
hold on
for q=1:length(shp)
ll1=shp(q).X; ll2=shp(q).Y; % delete Hawaii
plot(ll1,ll2,'-k');
end

scatter(stninfo(:,3),stninfo(:,2),50,tmean_err(:,1),'filled')
hold off
colorbar
caxis([-10,10])
xlim([-130,-50])
ylim([70,85])


