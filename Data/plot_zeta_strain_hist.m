clear all
close all

%%
asp_rat = 0.55;
fig1 = figure(111); pplot(8,asp_rat,10); box on; hold on
fig2 = figure(222); pplot(8,asp_rat,10); box on; hold on

%%
case_ary = ["1" "2"];

%% 
zeta_all = [];
strain_all = [];
div_all = [];

file_nm = "runSQG_WN_n10_Ro_2_fin"+".mat"; load(file_nm)

[u,v] = UV_sqgp1(y1,Ro,ikx_,iky_,K_,IK_,aa_filter);
u_k = fft2_n(u); v_k = fft2_n(v);
zeta = ifft2_n( -u_k.*iky_+v_k.*ikx_,'symmetric' )*0.15;
str1 = ifft2_n( u_k.*ikx_-v_k.*iky_,'symmetric' );
str2 = ifft2_n( v_k.*ikx_+u_k.*iky_,'symmetric' );
strain = sqrt(str1.^2+str2.^2)*0.15;

zeta_all = [zeta_all; zeta(:)];
strain_all = [strain_all; strain(:)];

%%
div = ifft2_n( u_k.*ikx_+v_k.*iky_,'symmetric' )*0.15;
div_all = [div_all; div(:)];

% Create buns and obtain count
x_bin_num = 100;
Xedges=linspace( -2 , 3 ,x_bin_num+1);
Yedges=linspace( 0 , 3 ,x_bin_num/2+1);
histN = histcounts2(zeta_all,strain_all,Xedges,Yedges);
% use bin edges to find the center of each bin in X and Y
[X,Y] = meshgrid(mean([Xedges(1:end-1);Xedges(2:end)],1),...
    mean([Yedges(1:end-1);Yedges(2:end,1)],1));

xmax = max(X(:)); xmin = min(X(:)); 
ymax = max(Y(:)); ymin = min(Y(:)); 

%%
% Create countor, transponsing N so rows correspond to Y and columns correspond to X
figure(111)
% contourf(X,Y,N',50)

histogram2(zeta_all,strain_all,[x_bin_num x_bin_num/2],'DisplayStyle','tile','Normalization','pdf','EdgeColor','none')
axis equal

plot(0:3,0:3,'k--','LineWidth',1)
plot([xmin:0 0],-[xmin:0 0],'k--','LineWidth',1)

ylabel("$\{\varepsilon\}\sigma^t/f$"); xlabel("$\{\varepsilon\}\zeta^t/f$")
axis([xmin xmax ymin ymax])
set(gca,'ColorScale','log'); box on

cb = colorbar();
cmocean('amp')

title("Surface Vorticity-Strain JPDF")

x1=get(gca,'position');
x=get(cb,'Position');
x(3)=0.03;
set(cb,'Position',x)
set(gca,'position',x1)

%%
tic

meanbin_div = nan(size(histN));
for xi = 1:length(Xedges)-1
    x_low_filt = zeta_all >= Xedges(xi);
    x_high_filt = zeta_all < Xedges(xi+1);
    x_filt = and(x_low_filt,x_high_filt);
    for yi = 1:length(Yedges)-1
        if histN(xi,yi) < 2
            meanbin_div(xi,yi) = 0;
        else
            y_low_filt = strain_all >= Yedges(yi);
            y_high_filt = strain_all < Yedges(yi+1);
            y_filt = and(y_low_filt,y_high_filt);

            tot_filt = and(x_filt,y_filt);

            meanbin_div(xi,yi) = mean(div_all(tot_filt))/Ro;
%             pdf_bin(xi,yi) = sum(tot_filt)/length(div_all);
%             if sum(tot_filt) == 1
%                 pdf_bin(xi,yi) = 0;
%             end
        end
    end
end

toc

%%
div_pad = -floor( (log(min(abs(meanbin_div(meanbin_div~=0))))/log(10))/1.8 ); 

data_log = nan(size(meanbin_div));
data_log(meanbin_div>0) = log(meanbin_div(meanbin_div>0))/log(10)+div_pad;
data_log(meanbin_div<0) = -log(-meanbin_div(meanbin_div<0))/log(10)-div_pad;
data_log(isnan(data_log)) = 0;

figure(222)

colorplot = pcolor(X,Y,data_log'); hold on
set(colorplot, 'EdgeColor', 'none')
axis equal

plot(0:10,0:10,'k--','LineWidth',1)
plot([xmin:0 0],-[xmin:0 0],'k--','LineWidth',1)

axis([xmin xmax ymin ymax])
cmocean('balance')
colormap(); box on
ylabel("$\{\varepsilon\}\sigma^t/f$"); xlabel("$\{\varepsilon\}\zeta^t/f$")

b_r_max = max(abs(data_log(:))); caxis([-b_r_max b_r_max]); 
%%
cb = colorbar('Ticks',[-3,-2,-1,0,1,2,3]); 
cb.Label.Interpreter = 'latex';
cblab = cb.TickLabels; 
cblab_ary = str2double(cblab);
cblablog_ary = ["","","","0","","",""];
cblablog_ary(cblab_ary<0) = "$-10^{"+( -cblab_ary(cblab_ary<0)-div_pad )+"}$";
cblablog_ary(cblab_ary>0) = "$10^{"+( cblab_ary(cblab_ary>0)-div_pad )+"}$";
cb.TickLabels = cblablog_ary;
cb.TickLabelInterpreter = 'latex';

title("Conditional Mean of $\{\varepsilon\}\delta^t/f$")
% cb.Label.String = "$$"; cb.Label.Interpreter = 'latex'; cb.FontSize = 10;

x1=get(gca,'position');
x=get(cb,'Position');
x(3)=0.03;
set(cb,'Position',x)
set(gca,'position',x1)
%%
figure(111)
savefig("figs/WN_zetastrain_hist")
figure(222)
savefig("figs/WN_zetastrain_div_hist")


