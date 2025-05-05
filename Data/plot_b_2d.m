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
    b_real = ifft2_n(b_k,'symmetric');
    colorplot = pcolor(px_,py_,b_real-mean(b_real(:)));
    set(colorplot, 'EdgeColor', 'none')
    axis equal
    b_r_max = 4.0;
    caxis([-b_r_max b_r_max]); 
    cmt = cmocean('thermal');  
    colormap(cmt)
    
end

%%
h1 = axes(fig,'visible','off'); 
h1.Title.Visible = 'on'; h1.XLabel.Visible = 'on'; h1.YLabel.Visible = 'on';
title(h1,"");
c = colorbar(h1,'Position',[0.89 0.12 0.03 0.78]);  % attach colorbar to h
c.LineWidth = 1;
c.Label.String = "$(b^\mathrm{t}-\overline{b^\mathrm{t}})/(UfL/H)$"; c.Label.Interpreter = 'latex'; c.FontSize = 10;
caxis(h1,[-b_r_max,b_r_max]); 

%%
savefig("figs/WN_b_2d")

