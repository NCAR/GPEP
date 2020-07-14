function CC=f_month_cc(data1,data2,month,cctype)
dbstop if error
data=[data1,data2];
CC=nan*zeros(1,12);
for mm=1:12
    indm=month==mm;
    dij=data(indm,:);
    dij(isnan(dij(:,1))|isnan(dij(:,2)),:)=[];
    if length(dij)>=300
        CC(mm)=corr(dij(:,1),dij(:,2),'Type',cctype);
    end
end

end