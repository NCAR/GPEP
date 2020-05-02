pcpobs=pcpobs(pcpobs>0);
pcpest=pcpest(pcpest>0);


subplot(2,2,1)
qqplot((pcpobs.^0.25-1)*4);
title('Station observations: box-cox')

subplot(2,2,2)
qqplot(pcpobs);
title('Station observations: none')
ylim([0,200]);

subplot(2,2,3)
qqplot((pcpest.^0.25-1)*4);
title('Station estimates: box-cox')

subplot(2,2,4)
qqplot(pcpest);
title('Station estimates: none')
ylim([0,200]);