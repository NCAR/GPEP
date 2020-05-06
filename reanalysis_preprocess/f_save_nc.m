function f_save_nc(gridfile,data,BasicInfo,DEMHigh)
for i=1:size(data,3)
   temp=data(:,:,i);
   temp(isnan(DEMHigh))=-999;
   data(:,:,i)=temp;
end

nccreate(gridfile,'data','Datatype','single',...
'Dimensions',{'lat',size(data,1),'lon',size(data,2),'time',size(data,3)},...
'Format','netcdf4','DeflateLevel',9,'FillValue',-999);
ncwrite(gridfile,'data',data);

info=[BasicInfo.Xll, BasicInfo.Yll, ...
    BasicInfo.tXll, BasicInfo.tYll, ...
    BasicInfo.cellsize];
nccreate(gridfile,'info','Datatype','double',...
'Dimensions',{'dimi',length(info)},...
'Format','netcdf4','DeflateLevel',9,'FillValue',-999);
ncwrite(gridfile,'info',info);
ncwriteatt(gridfile,'info','description','xll yll txll tyll cellsize');
end