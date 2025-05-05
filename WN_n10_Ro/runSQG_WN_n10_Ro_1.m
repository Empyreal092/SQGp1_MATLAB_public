full_path_name = mfilename('fullpath');
start_pos = strfind(full_path_name,'runSQG');

script_name = full_path_name(start_pos:end);

doruntimeplots = false;

Ro = 0;
log_n = 9;

frc_k_peak = 5;
end_time = frc_k_peak*100;

nu_UV_const = 1e-12;

%%
run("../sqgp1_driver")