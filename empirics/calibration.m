%%%%%%%CALIBRATION GREEN AMBIGUITY%%%%%%%%%

%Standard Parameters
rho_y = 2; %elasticity of substitution between green and dirty outputs (Papageorgiu et al. 2017)
%rho_y = 0.7; %cobb-douglas case
nu = 0.7; %share of dirty output
alpha = 0.3; %capital share for dirty production
sigma_h = 1; %inverse of Frisch elasticity/disutility of labor
delta = 0.025;  %depreciation rate (for now the same for green and dirty production)
%delta = 1;
beta = 0.99; %discount factor
betad = 0.98;
betag = betad;
%betag = 0.943;
betab = 0.985;
gamma_i = 3;
gammap=2; %risk aversion
gammad=1;
gammag=1;
gammab=1;
kappa = 0.1;
md = 0.9;
mg = 0.9;
phi=0.25;

%Environmental Parameters
epsilon = 0.7; %sensitivity of emissions to dirty production (Heutel 2012)
delta_m = 1-0.5^(1/(83*4)); %quarterly decay rate
%delta_m = 0.9;
m_bar = 581; %pre-industrial concentration of carbon
e_star=2.013261219344034; %damage 2100 5% output loss
%e_star = 0.506484709506948; %damage 2020 0.5% output loss
xi=4.223836122353320e-05; %damage 2100 5% output loss
%xi = 1.640734330436812e-05; %damage 2020 0.5% output loss

%Parameters for exogenous processes
rho_ad = 0.9; %autoregressive parameter for dirty technological progress
rho_ag = 0.9; %autoregressive parameter for green technological progress
rho_amb = 0.9; %autoregressive parameter for ambiguity progress
rho_fd = 0.9;
rho_fg = 0.9;
% sig_ad = 0.005; %standard deviation for dirty technological progress
% sig_ag = 0.005; %standard deviation for green technological progress
% sig_amb = 0.1; %standard deviation for ambiguity progress

params = [rho_y, nu, alpha, sigma_h, delta, beta, phi, betad, betag, betab, kappa, ...
         mg, md, gammap, gammad, gammag, gammab, epsilon, gamma_i, delta_m, m_bar, e_star, xi, ...
         rho_ad, rho_ag, rho_amb, rho_fd, rho_fg];
     
save ('par_calibration.mat', 'params')

save ('params_calibration.mat', 'rho_y', 'nu', 'alpha', 'sigma_h','delta', 'beta', 'phi', ...
    'betad', 'betag', 'betab', 'kappa', 'mg', 'md', 'gammap', 'gammad', 'gammag' , 'gammab', ...
    'epsilon', 'gamma_i', 'delta_m', 'm_bar', 'e_star', 'xi', 'rho_ad', 'rho_ag', 'rho_amb', 'rho_fd', 'rho_fg')

