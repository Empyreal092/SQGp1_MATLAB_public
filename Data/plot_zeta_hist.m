clear all
close all

%%
case_ary = ["1" "2"];
%%
for Roi = [1 2]
    
    figure(123); pplot(8,0.8,10); box on; hold on
    
    %%
    zeta_all = [];
    for ti = 1:length(case_ary)
        file_nm = "runSQG_WN_n10_Ro_"+case_ary(Roi)+"_fin"+".mat"; load(file_nm)
        %%
%         aa_filter = ( K_ < 1/2*n/2*(1/k_scale));
        [u,v] = UV_sqgp1(b_k,Ro,ikx_,iky_,K_,IK_,aa_filter);
        u_k = fft2_n(u);
        u_k_filt = u_k; u_k_filt(K_ <= 0) = 0;
        v_k = fft2_n(v);
        v_k_filt = v_k; v_k_filt(K_ <= 0) = 0;
        zeta = ifft2_n( -u_k_filt.*iky_+v_k_filt.*ikx_,'symmetric' )*0.15;
        zeta_all = [zeta_all; zeta(:)];
    end
    
    histogram(zeta_all(:),'DisplayStyle','stairs','LineWidth',1,'Normalization','pdf','DisplayName',"Ro="+Ro); hold on
    set(gca,'YScale','log')
    
end
%%
title("PDF of $\{\varepsilon\}\zeta^\mathrm{t}/f$")
xlabel("$\{\varepsilon\}\zeta^\mathrm{t}/f$")
% ylabel("PDF")
legend show
legend('Location','northeast')
ylimm = ylim; ylim([5e-4 ylimm(2)*2])
xlim([-1.75 2.75])

hold off

%%
figure(123)
savefig("figs/WN_zeta_hist")
