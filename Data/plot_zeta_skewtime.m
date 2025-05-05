clear all
close all

file_nm = "runSQG_WN_n10_Ro_2_fin"+".mat"; load(file_nm)

%%
figure(123); pplot(8,0.8,10); box on; hold on
plot(ts_ary,zeta_skew_ary); 
hold off

%%
title("Skewness of vorticity")
xlabel("$t$")
ylabel("Skew($\zeta^\mathrm{t}$)")

hold off

%%
figure(123)
savefig("figs/WN_zeta_skewtime")
