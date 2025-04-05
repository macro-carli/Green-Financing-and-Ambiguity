%%%%%%%STEADY_STATE GREEN AMBIGUITY%%%%%%%%%
run calibration

par_calibration = load('par_calibration.mat');
params = par_calibration.params;

%y_ss = 1;
%h_ss = 0.17;
%m_ss = 891.3378;
%m_ss = 1650.537;
%e_ss = delta_m*(m_ss-m_bar);

amb_ss = 0.008; %to match highest loan rate spread 2007-2017
%amb = 0;
ag_ss = exp(-amb_ss/(1-rho_ag));
ad_ss = 1;

r_ss = beta^(-1);

% %Spread loan rate/interest rate target = 2;
% spread_ss = 2;
% targ_z = r_ss+spread_ss;
% 0.6 is ss target for charge-off rate 2007-2017
% 
% scale_par_fd = log(0.6)*(1-rho_fd)/(-ad_ss); 
% %scale_par_fd =1;
% zfd =(-ad_ss*scale_par_fd/(1-rho_fd));
% fd_ss = 1*exp(zfd);
% 
% zfg=(-ag_ss*scale_par_fd/(1-rho_fg));
% fg_ss = 1*exp(zfg);

fd_ss = 0.35;
fg_ss = 0.35;

zd_ss = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fd_ss);
zg_ss = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fg_ss);

ezd_ss = zd_ss*(1-fd_ss);
ezg_ss = zg_ss*(1-fg_ss);

qd_ss = 1;
qg_ss = 1;

chi = 1.8049; %damage 2100 5% output loss

taud = 0;

params = [params amb_ss ad_ss ag_ss r_ss fd_ss fg_ss zd_ss zg_ss qd_ss qg_ss chi taud];

%Initial Values (computed with fall method starting from yd, yg and
%implying ad, ag)

yd0 = 1.3; 
yg0 = 0.8;
pd0 = 1.06;
pg0 = 0.88;
kg0 = 4.67;
wg0 = 1.11;
hg0 = 0.45;
kd0 = 9.09;
wd0 = 1.54;
hd0 = 0.63;
e0 = 0.5195;
y0 = 0.5;
lambda0 = 0.41;
c0 = 1.57;
invd0 = 0.23;
invg0 = 0.12;
m0 = 1650.537;
l0 = 12.2;
lg0 = 4.14;
ld0 = 8.06;
b0 = 11;
taug0 = 0;

X0 = [yd0 yg0 pd0 pg0 kg0 wg0 hg0 kd0 wd0 hd0 e0 y0 lambda0 c0 invd0 invg0 m0 l0 lg0 ld0 b0 taug0];


options = optimset('Display','off', 'TolFun', 1e-12, 'TolX', 1e-12);
[sol,Fval,exitflag] = fsolve('ss_fin_sens_solver',X0,options,params);

% options = optimset('Display','off', 'TolFun', 1e-12, 'TolX', 1e-12);
% sol=fsolve(@(x) ss_solver_sens(x, rho_y, nu, alpha1, alpha2, alpha, sigma_h, delta, epsilon, gammap,...
%                  m_bar, xi, chi, ad_ss, ag_ss, rd_ss, rg_ss, delta_m, e_star, tau), [0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;1;0.5;0.5;0.5;0.5;0.5], options);
%                 %m_bar, xi, m_ss, e_ss, ad_ss, ag_ss, rd_ss, rg_ss, px), [0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;0.5;1;0.5;0.5;0.5;0.5;0.5], options);
%                 %m_bar, xi, m_ss, e_ss, ad_ss, ag_ss, rd_ss, rg_ss, px), [0.7;0.3;1;1;2.8;0.4;0.5;6.6;0.6;0.8;0.17;1;1.3;0.8;0.16;0.07;0.5], options);
%                 %m_bar, xi, m_ss, e_ss, ad_ss, ag_ss, rd_ss, rg_ss, px), [0.5;0.5;1.05;0.9;4;0.15;0.2;4;0.14;0.18;0.2;1;1.25;0.8;0.1;0.1;0.5], options);

yd_ss = sol(1); 
yg_ss = sol(2);
pd_ss = sol(3);
pg_ss = sol(4);
kg_ss = sol(5);
wg_ss = sol(6);
hg_ss = sol(7);
kd_ss = sol(8);
wd_ss = sol(9);
hd_ss = sol(10);
e_ss = sol(11);
y_ss = sol(12);
lambda_ss = sol(13);
c_ss = sol(14);
invd_ss = sol(15); 
invg_ss = sol(16);
m_ss = sol(17);
l_ss = sol(18);
lg_ss = sol(19);
ld_ss = sol(20);
b_ss = sol(21);
taug_ss = sol(22);

h_ss = hd_ss + hg_ss;

dam_ss = exp(-xi*(m_ss-m_bar));

welf_ss = (c_ss^(1-gammap)/(1-gammap)- hg_ss^(1+sigma_h)/(1+sigma_h) - hd_ss^(1+sigma_h)/(1+sigma_h))/(1-beta);
%welf_ss = (log(c_ss)- hg_ss^(1+sigma_h)/(1+sigma_h) - hd_ss^(1+sigma_h)/(1+sigma_h))/(1-beta);


% save ('ss_calibration.mat', 'm_ss', 'e_ss', 'amb_ss', 'ag_ss','ad_ss', 'r_ss', 'rg_ss', 'rd_ss',...
%     'yd_ss', 'yg_ss', 'pd_ss', 'pg_ss', 'kg_ss', 'wg_ss', 'hg_ss', 'kd_ss', 'wd_ss', 'hd_ss',...
%     'chi', 'y_ss', 'lambda_ss', 'c_ss', 'invd_ss', 'invg_ss', 'x_ss', 'h_ss')
