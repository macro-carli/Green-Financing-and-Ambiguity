%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ys,params,check]=ambig_envir_fin_steadystate(ys,exe,M_,options_)
% compute s.s.
%global M_ options_
check = 0;

%the following statements ensure that the variables are not interpreted by
%MATLAB as functions...
beta=NaN;
alpha=NaN;
% Here we load the values of the deep parameters in a loop.
Np = M_.param_nbr;                                            
for i = 1:Np
    paramname = M_.param_names{i};
    eval([ paramname ' = M_.params(' int2str(i) ');']);
end

% read in instrument values
% for ii = 1:size(options_.instruments,1)
%   eval([options_.instruments{ii} ' = ys(strmatch(options_.instruments{ii},M_.endo_names,''exact'')) ;']);
% end
% initialize indicator
check = 0;

% m = 1795.377;
% e = delta_m*(m-m_bar)-e_star;
% %chi = 0.342945253112513;
% 
% amb = 0.001;
% fd = 0.002;
% fg = 0.002 + amb;
% %fd = (1/(1-rho_fd));
% %fg = (1/(1-rho_fg));
% r = 1/beta;
% % ezd = kappa*(1/betab) + (1-kappa)*(1/beta); %effective gross dirty rate
% % zd = ezd/(1-fd); %contractual gross dirty rate
% % ezg = ezd; %effective gross green rate (no arbitrage)
% % zg = ezg/(1-fg); %contractual gross green rate
% zd = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fd);
% zg = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fg);
% 
% disc=1;
% qd=1;
% qg=1;
% yd=2.1;
% yg=0.9;
% y = (nu^(1/rho_y)*yd^((rho_y-1)/rho_y) + (1-nu)^(1/rho_y)*yg^((rho_y-1)/rho_y))^(rho_y/(rho_y-1));
% pd = ((nu)^(1/rho_y))*((yd/y)^(-1/rho_y));
% pg = ((1-nu)^(1/rho_y))*((yg/y)^(-1/rho_y));
% kd = (alpha*pd*yd)/(1/betad - (1-delta) - md/(zd*(1-fd)*betad) +md);
% kg = (alpha*pg*yg)/(1/betag - (1-delta) - mg/(zg*(1-fg)*betag) +mg);
% invd = delta*kd;
% invg= delta*kg;
% lg = mg*(qg/(zg*(1-fg)))*kg;
% ld = md*(qd/(zd*(1-fd)))*kd;
% l = ld+lg;
% b = (1-kappa)*l;
% % pib = zd*(1-fd)*ld + zg*(1-fg)*lg + b - r*b - l;
% % pid = pd*yd + ld - invd - (1-alpha)*pd*yd - zd*(1-fd)*ld;
% % pig = pg*yg + lg - invg - (1-alpha)*pg*yg - zg*(1-fg)*lg;
% c = y-invg-invd;
% lambda = c^(-gammap);
% hd = ((1-alpha)*(pd*yd*lambda))^(1/(sigma_h+1));
% hg = ((1-alpha)*(pg*yg*lambda))^(1/(sigma_h+1));
% h = hg+hd;
% wd = ((1-alpha)*pd*yd)/hd;
% wg = ((1-alpha)*pg*yg)/hg;
% ad = yd/(exp(-xi*(m-m_bar))*(kd^(alpha))*(hd^(1-alpha)));
% ag = yg/(exp(-xi*(m-m_bar))*(kg^(alpha))*(hg^(1-alpha)));
% chi = e/(yd^(epsilon));
% 
% 
% welf = (c^(1-gammap)/(1-gammap) - hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))/(1-beta);
% %welf = (log(c) - hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))/(1-beta);
% tau =1;


params = [rho_y, nu, alpha, sigma_h, delta, beta, phi, betad, betag, betab, kappa, ...
         mg, md, gammap, gammad, gammag, gammab, epsilon, gamma_i, delta_m, m_bar, xi, ...
         rho_ad, rho_ag, rho_amb, rho_fd, rho_fg];

zamb = 0;
amb = 0.008; %to match highest loan rate spread 2007-2017
%amb = 0;
ag = exp(-amb/(1-rho_ag));
ad = 1;

fin = 1;

r = beta^(-1);

% %Spread loan rate/interest rate target = 2;
% spread = 2;
% targ_z = r+spread;
% 0.6 is ss target for charge-off rate 2007-2017

%scale_par_fd = log(0.6)*(1-rho_fd)/(-ad); 
scale_par_fd =1;
%zfd =(-ad*scale_par_fd/(1-rho_fd));
%fd = 1*exp(zfd);
fd = 0.35;

%zfg=(-ag*scale_par_fd/(1-rho_fg));
%fg = 1*exp(zfg);
fg = 0.35;

zd = (kappa*weightd*(1/betab) + (1-kappa*weightd)*(1/beta))/(1-fd);
zg = (kappa*weightg*(1/betab) + (1-kappa*weightg)*(1/beta))/(1-fg);

ezd=(1-fd)*zd;
ezg=(1-fg)*zg;

spreadd=zd/r;
spreadg=zg/r;
diffz=zg/zd;
kappa_flu = kappa;

qd = 1;
qg = 1;

disc = 1;

taud = 0;

chi = 1.8049;

params = [params amb ad ag r fd fg zd zg qd qg chi taud weightd weightg];

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
m0 = 1795.377;
l0 = 12.2;
lg0 = 4.14;
ld0 = 8.06;
b0 = 11;
taug0 = 0;

X0 = [yd0 yg0 pd0 pg0 kg0 wg0 hg0 kd0 wd0 hd0 e0 y0 lambda0 c0 invd0 invg0 m0 l0 lg0 ld0 b0 taug0];


options = optimset('Display','off', 'TolFun', 1e-12, 'TolX', 1e-12);
[sol,Fval,exitflag] = fsolve('ss_fin_sens_solver_mod',X0,options,params);

yd = sol(1); 
yg = sol(2);
pd = sol(3);
pg = sol(4);
kg = sol(5);
wg = sol(6);
hg = sol(7);
kd = sol(8);
wd = sol(9);
hd = sol(10);
e = sol(11);
y = sol(12);
lambda = sol(13);
c = sol(14);
invd = sol(15); 
invg = sol(16);
m = sol(17);
l = sol(18);
lg = sol(19);
ld = sol(20);
b = sol(21);
taug = sol(22);

h = hd + hg;

welf = (c^(1-gammap)/(1-gammap)- hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))/(1-beta);
%welf = (log(c)- hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))/(1-beta);

scale_par_fd_idx = find(strcmp(M_.param_names, 'scale_par_fd'));
M_.params(scale_par_fd_idx)=scale_par_fd;
eval(['scale_par_fd = M_.params(' int2str(scale_par_fd_idx) ');']);

chi_idx = find(strcmp(M_.param_names, 'chi'));
M_.params(chi_idx)=chi;
eval(['chi = M_.params(' int2str(chi_idx) ');']);

ad_idx = find(strcmp(M_.param_names, 'ad_ss'));
M_.params(ad_idx)=ad;
eval(['ad_ss = M_.params(' int2str(ad_idx) ');']);

ag_idx = find(strcmp(M_.param_names, 'ag_ss'));
M_.params(ag_idx)=ag;
eval(['ag_ss = M_.params(' int2str(ag_idx) ');']);

fd_idx = find(strcmp(M_.param_names, 'fd_ss'));
M_.params(fd_idx)=fd;
eval(['fd_ss = M_.params(' int2str(fd_idx) ');']);

fg_idx = find(strcmp(M_.param_names, 'fg_ss'));
M_.params(fg_idx)=fg;
eval(['fg_ss = M_.params(' int2str(fg_idx) ');']);

amb_idx = find(strcmp(M_.param_names, 'amb_ss'));
M_.params(amb_idx)=amb;
eval(['amb_ss = M_.params(' int2str(amb_idx) ');']);

% ld_idx = find(strcmp(M_.param_names, 'ld_ss'));
% M_.params(ld_idx)=ld;
% eval(['ld_ss = M_.params(' int2str(ld_idx) ');']);
% 
% lg_idx = find(strcmp(M_.param_names, 'lg_ss'));
% M_.params(lg_idx)=lg;
% eval(['lg_ss = M_.params(' int2str(lg_idx) ');']);


params=NaN(Np,1);
for iter = 1:length(M_.params) %update parameters set in the file
  eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
end

eval([ 'params(' num2str(chi_idx) ') = ' M_.param_names{chi_idx} ';' ])
eval([ 'params(' num2str(ad_idx) ') = ' M_.param_names{ad_idx} ';' ])
eval([ 'params(' num2str(ag_idx) ') = ' M_.param_names{ag_idx} ';' ])
eval([ 'params(' num2str(fd_idx) ') = ' M_.param_names{fd_idx} ';' ])
eval([ 'params(' num2str(fg_idx) ') = ' M_.param_names{fg_idx} ';' ])
eval([ 'params(' num2str(amb_idx) ') = ' M_.param_names{amb_idx} ';' ])
% eval([ 'params(' num2str(ld_idx) ') = ' M_.param_names{ld_idx} ';' ])
% eval([ 'params(' num2str(lg_idx) ') = ' M_.param_names{lg_idx} ';' ])

NumberOfEndogenousVariables = M_.orig_endo_nbr; %auxiliary variables are set automatically
for ii = 1:NumberOfEndogenousVariables
  varname = M_.endo_names{ii};
  eval(['ys(' int2str(ii) ') = ' varname ';']);
end
