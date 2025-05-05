% clear fg13
% close all

figure_default % Default Plot Values

addpath(genpath('SubPlots'))
timer_plot = tic;

if ~exist('fg13')
    fg13 = figure(13);
    fg13.WindowState = 'maximized';
end

%% Flux Calc

kek = b_k.*conj(b_k);
z_spec = real( isospectrum(kek) );

% zt_spec_av = isospectrum(zt_av);
% et_av = zt_av.*IK_;
% et_spec_av = isospectrum(et_av);
% 
% FluxZ=-cumsum(zt_spec_av);
% FluxE=-cumsum(et_spec_av);
% 
% zt_d = isospectrum(zdfac.*kek); % Dissipation density of action
% et_d = isospectrum(zdfac.*kek.*IK_); % Dissipation density of action
% 
% switch FT
%     case 'WN'
%         zt_f = isospectrum(2*f_val^2*frc_mask); % Dissipation density of action
%         et_f = isospectrum(2*f_val^2*frc_mask.*IK_); % Dissipation density of action
%     case {'instab_ring','NF','instab_fixZ'}
%         zt_f = isospectrum(zffac.*kek); % Dissipation density of action
%         et_f = isospectrum(zffac.*kek.*IK_); % Dissipation density of action
%     case 'WN_instab'
%         zt_f = isospectrum(2*f_val^2*frc_mask+zffac.*kek); 
%         et_f = isospectrum(2*f_val^2*frc_mask.*IK_*zffac.*kek.*IK_); 
% end
        

%% Plot: Real part of the Solution
subaxis(2,4,1,'ML',0.05,'MR',0.05,'MT',0.05,'MB',0.07,'SV',0.09,'SH',0.04);

b_real = ifft2_n(b_k,'symmetric');
heatmap2d(b_real,px_,py_); axis equal
b_r_max = max( [max(b_real,[],'all'),max(-b_real,[],'all'),1e-10]);
caxis([-b_r_max b_r_max]); colormap(redblue);
title("Ro="+Ro+"; $t=$"+t)

%%
subaxis(2,4,2);
[u,v] = UV_sqgp1(b_k,Ro,ikx_,iky_,K_,IK_,aa_filter);
u_k = fft2_n(u); v_k = fft2_n(v);
zeta_real = ifft2_n( -u_k.*iky_+v_k.*ikx_,'symmetric' );

zeta_real = zeta_real;
heatmap2d(zeta_real,px_,py_); axis equal
z_r_max = max( [max(zeta_real,[],'all'),max(-zeta_real,[],'all'),1e-10] );
caxis([-z_r_max z_r_max]); colormap(redblue);
title("Ro="+Ro+"; $t=$"+t)

%% Plot: spectra
subaxis(2,4,3);

z_max = floor(log(max(z_spec(2:end)))/log(10)+1.5); z_max_diff = z_max+1;

loglog(k1_,z_spec); hold on

loglog(k1_,k1_.^-(5/3)*10^(z_max_diff-2),'k')
loglog(k1_,k1_.^-(2)*10^(z_max_diff-2),'k--')
xline(n/2/k_scale); xline(n/2*(1/2)/k_scale,'--')
xline(frc_k_peak/k_scale)
xlimm = xlim; xlim([1/k_scale xlimm(2)])
% ylim([10^(-15) 10^z_max]);
hold off

%%
subaxis(2,4,4);
% %%%%
plot(ts_ary,zeta_skew_ary); hold on
hold off
% %%%
%% 
subaxis(2,4,5);
plot(ts_ary,zeta_std_ary);

%% 
subaxis(2,4,6);
% semilogx(k1_,FluxE,'Color',[0 157 115]/256); hold on
% semilogx(k1_,et_d.*k1_*k_scale,'b')
% semilogx(k1_,et_f,'r')
% xline(n/2/k_scale); xline(n/2*(2/3)/k_scale,'--')
% xlimm = xlim; xlim([1/k_scale xlimm(2)])
% yline(0)
% title("Flux $E$"); 
% hold off

%%%
[u,v] = UV_sqgp1(b_k,Ro,ikx_,iky_,K_,IK_,aa_filter);
u_k = fft2_n(u);
u_k_filt = u_k; u_k_filt(K_ <= -10) = 0;
v_k = fft2_n(v);
v_k_filt = v_k; v_k_filt(K_ <= -10) = 0;
zeta = ifft2_n( -u_k_filt.*iky_+v_k_filt.*ikx_,'symmetric' );
zeta_neg = zeta(zeta<=0); 
zeta_sym = [zeta_neg(:); -zeta_neg(:)];
% histogram(zeta_sym,'DisplayStyle','stafirs','LineWidth',0.5,'Normalization','pdf','DisplayName',"Ro="+Ro); hold on

histogram(zeta(:),'DisplayStyle','stairs','LineWidth',1,'Normalization','pdf','DisplayName',"Ro="+Ro); 
set(gca,'YScale','log')
hold off

%%
subaxis(2,4,7);
plot(cumsum(dt_ary),dt_ary)

%% Plot: Correction size
subaxis(2,4,8);
plot(ts_ary,b_mean_ary); hold on
% plot(ts_ary,Ro2_err_ary);
yline(0);
% ylim([-0.3 0.6])
hold off

%%
drawnow
t_plot = toc(timer_plot); disp("Plotting took "+t_plot+" seconds");
