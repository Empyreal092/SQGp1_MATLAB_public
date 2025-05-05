clear all
close all

%%
case_ary = ["1" "2"]; 
time_str = 500;
%%

fig=figure(123); hold on; box on; pplot(15,0.5,10)
subaxis(1,2,1,'ML',0.08,'MR',0.13,'MT',0.1,'MB',0.12,'SV',0.03,'SH',0.02);

subaxis(1,2,1)
title("$$\mathrm{SQG}, t="+time_str+"$")
ylabel("$y$")
xlabel("$x$")
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)
subaxis(1,2,2)
title("$\mathrm{SQG}^{+1}, \mathrm{Ro}="+"0.15"+"$")
xlabel("$x$")
set(gca,'yticklabel',[])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10)

%%

for Roi = 1:length(case_ary)
    file_nm = "runSQG_WN_n10_Ro_"+case_ary(Roi)+"_fin"+".mat"; load(file_nm)
    %%
%     aa_filter = ( K_ < 1/2*n/2*(1/k_scale));

    subaxis(1,2,Roi); hold on
    [u,v] = UV_sqgp1(y1,Ro,ikx_,iky_,K_,IK_,aa_filter);
    u_k = fft2_n(u); v_k = fft2_n(v);
    zeta = ifft2_n( -u_k.*iky_+v_k.*ikx_,'symmetric' )*0.15;

    colorplot = pcolor(px_,py_,zeta);
    set(colorplot, 'EdgeColor', 'none')
    axis equal

    z_r_max = 1.5;
    caxis([-z_r_max z_r_max]); 
    cmapb = cmocean('balance'); 
    colormap(cmapb)
    
end

%%
h1 = axes(fig,'visible','off'); 
h1.Title.Visible = 'on'; h1.XLabel.Visible = 'on'; h1.YLabel.Visible = 'on';
title(h1,"");
c = colorbar(h1,'Position',[0.89 0.12 0.03 0.78]);  % attach colorbar to h
c.LineWidth = 1;
c.Label.String = "$\{\varepsilon\}\zeta^\mathrm{t}/f$"; c.Label.Interpreter = 'latex'; c.FontSize = 10;
caxis(h1,[-z_r_max,z_r_max]); 

%%
savefig("figs/WN_zeta_2d")

