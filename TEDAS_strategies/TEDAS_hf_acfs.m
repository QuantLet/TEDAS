clc
clear
load('hedfhd_short.mat') %pre-loading data 
load('hedfhdnames_short.mat')
XRET = hedfd_short(1:end,2:end-1); 

figure, hold('on')
subplot(2,2,1)
  
autocorr(XRET(:,67))
name1 = hedfhdnames_short{67};
title(name1(13:end-17))
subplot(2,2,2)
  
autocorr(XRET(:,59))
name2 = hedfhdnames_short{59};
title(name2(13:end-17))
subplot(2,2,3)
 
autocorr(XRET(:,85))
name5 = hedfhdnames_short{85};
title(name5(13:end-20))
subplot(2,2,4)
  
autocorr(XRET(:,45))
name4 = hedfhdnames_short{45};
title(name4(13:end-17))
hold off