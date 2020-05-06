function KGEgroup=f_KEG(Obs,Pre)
pre_mean = nanmean(Pre);
obs_mean = nanmean(Obs);
r = nansum((Pre - pre_mean) .* (Obs - obs_mean)) / sqrt(nansum((Pre - pre_mean).^2).*nansum((Obs - obs_mean).^2));
gamma = (std(Pre)/pre_mean) / (std(Obs) / obs_mean);
beta = nanmean(Pre)/nanmean(Obs);
KGE = 1 - sqrt((r - 1)^2 + (gamma - 1)^2 + (beta - 1)^2);
KGEgroup = [KGE,r,gamma,beta];
end
