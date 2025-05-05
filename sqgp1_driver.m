% clear all
% close all

%% Numerical solution of the SQG+ equation
% For detail, see note

%% Initialization
timer_wholescipt = tic; % Start the timer for the whole scipt
addpath(genpath('SubRoutine'))  

if ~exist('script_name','var') script_name = mfilename; end
mkdir("Run_Data/"+script_name)

% rng(111);

%% Determine if this is an new run
if ~exist('end_time','var') end_time = 50; end

if ~exist('RHS_1','var')
    new_run=true; tnew = [0 end_time];
    t=tnew(1);
    t_end = tnew(2);
    new_run_msg = "This is a NEW run till " + t_end; disp(new_run_msg)
else
    new_run=false; tconti = [t end_time];
    t_end = tconti(2);
    new_run_msg = "This is a CONT run, from t="+t+" to "+t_end; disp(new_run_msg)
end

%% SQG+1 Parameters.
% Rossby number
if ~exist('Ro','var') Ro = 0.1; end
disp("Ro = "+Ro);

%% Numerical Settings
time_step_method = "IF-RK4";
if ~exist('aa_on','var') aa_on = true; end

if ~exist('log_n','var') log_n = 8; end
n=2^log_n; % Number of gridpoints

%%% Timestep specifications
if ~exist('dt_tune','var') dt_tune = 1; end

disp("log(n) = "+log_n+", dt_tune = "+dt_tune);

%% Run Time Plot Settings
if ~exist('doruntimeplots','var') doruntimeplots = true; end
disp("do runtimeplots = "+doruntimeplots);

%% IC, Frc, Diss types
% Initial condition
if ~exist('frc_k_peak','var') frc_k_peak = 5; end
if ~exist('IC','var') 
    IC = 'zero';
end

% Choice of forcing
if ~exist('FT','var')
    % FT = 'NF';
    FT='WN'; frc_k = [-1 0 1]+frc_k_peak; ZfluxTheory = 2e-3;
end

% Choice of dissipation
if ~exist('diss_switch_IR','var') diss_switch_IR = true; end
diss_switch_UV = true; UV_diss_pos = log_n;

%% Numerical setup details
% physical grid
k_scale = frc_k_peak/1.6;
dx = 2*pi/n*k_scale;
L = 2*pi*k_scale;
[px_,py_] = meshgrid([0:dx:(L-dx)]-L/2,[0:dx:(L-dx)]-L/2);

% Fourier grid
[kx_,ky_] = meshgrid([0:n/2 1-n/2:-1]/k_scale,[0:n/2 1-n/2:-1]/k_scale);
%%% Some variance of k for future calculations
ikx_ = 1i*kx_;  % only ever used for computing horizontal derivatives
iky_ = 1i*ky_;

k1_ = [0:n/2-1]/k_scale;
ik1_ = 1./k1_; ik1_(1) = 0;
K_ = sqrt(kx_.^2+ky_.^2); Kint_ = K_*k_scale;
IK_ = 1./K_; IK_(K_==0)=0; % for division by K_ w/o nans

% Anti-aliasing filter
if aa_on
    aa_filter = ( K_ < 1/2*n/2*(1/k_scale));
    disp("anti-aliasing: ON");
else
    aa_filter = ( K_ < 2*n/k_scale);
    disp("anti-aliasing: OFF");
end

%% Samping Settings
dts = 0.1*k_scale;
dtp = 1*k_scale;   % Plotting time interval for run time diags
dtsave = 10*k_scale;   % Plotting time interval for run time diags

tsrate = 1/200; % Decay rate for autoregressive sampling.  Fast rate is 10x that.
tsfac = 1-exp(-tsrate*dts); tsfacf= 1-exp(-10*tsrate*dts);

%% IC, Frc, Diss types details
switch IC  
    case 'zero'
        bk_0 = zeros(n,n);
        IC_msg = "IC: "+IC; disp(IC_msg);
end

%%%%%%%%%%%% Dissipation %%%%%%%%%%%%
if ~diss_switch_IR
    disp("IR Dissipation: OFF")
    nu_IR = 0;
else
    disp("IR Dissipation: ON")
    nu_IR = 2e-3;
end
if ~diss_switch_UV
    disp("UV Dissipation: OFF")
    nu_UV=0;
else
    disp("UV Dissipation: ON")
    if ~exist('nu_UV_const','var')
        nu_UV_const = 1e-12;
    end
	nu_UV=nu_UV_const*(2^8/2^UV_diss_pos)^8*(k_scale^8);
end

%%%%%%%%%%%% Forcing %%%%%%%%%%%%
switch FT
    case 'NF'
        aff = zeros(n,n); zffac = aff;
        forcing_msg = "Forcing: "+FT; disp(forcing_msg);
    case 'WN'
        aff = zeros(n,n); zffac = aff; 
        frc_mask = zeros(n,n);
        for K = [frc_k frc_k(2)]
            frc_mask = frc_mask + 1/2*((floor(Kint_+.5)-K)==0);
        end
        totnum_frck = sum(frc_mask.^2,'all');
        f_val=sqrt(ZfluxTheory/(totnum_frck)); % not multiply 2 is because real field is required
        forcing_msg = "Forcing: "+FT+" @ "+mat2str(frc_k)+" w/ XFlux: " + ZfluxTheory; disp(forcing_msg);
end

%% New Run Setup
if(new_run)
    b_k=bk_0;
    
    %% Samples Setup
    b_mean_ary = nan(1,round(tnew(2)/dts));
    zeta_std_ary = nan(1,round(tnew(2)/dts));
    zeta_skew_ary = nan(1,round(tnew(2)/dts));

    ts_ary = nan(1,round(tnew(2)/dts));
    dts_i = 1;
    
    dt_ary = [];
    
    %% time in loop
    % t has been set by now, so:
    tp   = t + dtp;  % next plotting time
    ts   = t + dts;  % next sampling time
%     tint = t + dtint;    % next integer sampling time
    tsave = t + dtsave;    % next save time
    
    timer_step = tic; % the first step time counter starts here, otherwise it starts in the plotting loop
else
    timer_step = tic; % the first step time counter starts here, otherwise it starts in the plotting loop
end

%%
dth = 1e-16;
while t <= t_end-dth
    if sum(isnan(b_k),'all') > 0 disp("NaN detected"); break; end
    %% RK4, with forcing and dissipation
    y1 = b_k;
    [RHS_1,u,v] = RHS_sqgp1(y1,Ro,ikx_,iky_,K_,IK_,aa_filter);    
    
    %% CFL
    maxU = max( [max(abs(u),[],'all') max(abs(v),[],'all')] );
    
    dt = dt_tune*(dx/maxU); dt = min( max(dt,1e-4*k_scale), 1e-1*k_scale);
    if (t+dt >= ts-dth) dt = ts-t; end
    dth = dt/2;
    
    %% Forcing, and dissipation
    diss1=exp(-(nu_IR+nu_UV*K_.^8)*dt); dissh=exp(-(nu_IR+nu_UV*K_.^8)*dth);
    forc1 = exp(zffac*dt); forch = exp(zffac*dth);
    itfm1 = diss1.*forc1; itfmh = dissh.*forch; itfm0 = 1;
    
    dt_ary = [dt_ary dt];
    %%
    y2 = itfmh.*b_k+dth*itfmh.*RHS_1;
    RHS_2 = RHS_sqgp1(y2,Ro,ikx_,iky_,K_,IK_,aa_filter);
    y3 = itfmh.*b_k+dth*itfm0.*RHS_2;
    RHS_3 = RHS_sqgp1(y3,Ro,ikx_,iky_,K_,IK_,aa_filter);
    y4 = itfm1.*b_k+dt *itfmh.*RHS_3;
    RHS_4 = RHS_sqgp1(y4,Ro,ikx_,iky_,K_,IK_,aa_filter);
    
    RHSm = (itfm1.*RHS_1 + 2*itfmh.*RHS_2 + 2*itfmh.*RHS_3 + itfm0.*RHS_4)/6;
    
    b_k = itfm1.*b_k + dt*RHSm;
    switch FT
        case {'NF'}
            % do nothing
            b_k = b_k;
        case {'WN'}
            f_dcorr = f_val*sqrt(dt)*diss1.*frc_mask;
            b_k = b_k+f_dcorr.*(rand_conjsym_for_ifft2(n));
    end
    t=t+dt;
    
    %% Time to sample
    if (t >= ts-dth)
        b_mean_ary(dts_i) = b_k(1,1);
        %%%
        [u,v] = UV_sqgp1(b_k,Ro,ikx_,iky_,K_,IK_,aa_filter);
        u_k = fft2_n(u); v_k = fft2_n(v);
        zeta_real = ifft2_n( -u_k.*iky_+v_k.*ikx_,'symmetric' );
        zeta_std_ary(dts_i) = std(zeta_real(:));
        zeta_skew_ary(dts_i) = mean(zeta_real.^3,'all')/std(zeta_real(:))^3;
        %%%
        ts_ary(dts_i) = t;
        ts = ts + dts; dts_i = dts_i + 1;
    end
    
    %% Runtime Plotting
    if doruntimeplots && t >= tp-dth
        t_onestep = toc(timer_step); disp("Between update took "+t_onestep+" seconds, now @ t = "+t);
        timer_step = tic;
        tp = tp + dtp;
        
        RuntimePlotting_sqgp1
    elseif not(doruntimeplots) && t >= tp-dth
        t_onestep = toc(timer_step);
        disp("Between update took "+t_onestep+" seconds, now @ t = "+t);
        timer_step = tic;
        tp = tp + dtp;
    end
    
    if t>= tsave-dth
        clear fg13
        save("Run_Data/"+script_name+"/"+script_name+"_t"+round(t)+".mat")
        tsave = tsave+dtsave;
    end
    
end

% disp("Autostop would have stopped the run at "+auto_stop_time)
tend = toc(timer_wholescipt); disp("Whole scipt took "+tend/60+" minutes or "+tend/(60*60)+" hours")

clear fg13
% script_name = mfilename;
save("Run_Data/"+script_name+"/"+script_name+"_fin.mat")
disp("END OF SCRIPT")


