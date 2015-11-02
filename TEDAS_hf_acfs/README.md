
![http://quantnet.wiwi.hu-berlin.de/style/banner.png](http://quantnet.wiwi.hu-berlin.de/style/banner.png)

## ![qlogo](http://quantnet.wiwi.hu-berlin.de/graphics/quantlogo.png) **<QNAME>**


```yaml
Name of QuantLet : TEDAS_hf_acfs

Published in : Tail Event Driven Asset Allocation
  
Description : 'produces autocorrelation plots for 4 different hedge funds returns'

Keywords : 'asset, returns, log-returns, index, strategy, time-series,
autocorrelation function'

See also : 'TEDAS_dcc_garch, TEDAS_ogarch, TEDAS_qrRho, TEDAS_strategies,
TEDASstrategies2, TEDAS_gestalts'

Author : Sergey Nasekin

Submitted : 30/10/2015
  
Datafile : hedfhd_short.mat, hedfhdnames_short.mat 

Output:  autocorrelation plots for 4 different hedge funds

Example: 
```
![Picture1](TEDAS_hf_acfs.png)

```Matlab
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
```
