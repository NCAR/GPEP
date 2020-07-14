function CC=f_monthmon_cc(data1,data2,cctype)
dbstop if error
data=[data1,data2];
CC=nan*zeros(1,12);
for mm=1:12
    dij=data(mm:12:end,:);
    dij(isnan(dij(:,1))|isnan(dij(:,2)),:)=[];
    if length(dij)>=10
        CC(mm)=corr(dij(:,1),dij(:,2),'Type',cctype);
    end
end

end