%%%GREEN FINANCE AND AMBIGUITY%%%%%
%%%MARCO CARLI%%%%%

var yd yg pd pg kg wg hg kd wd hd y lambda c invd invg inv r qd qg m e h ag amb ad zamb
    welf tau fd fg zd zg ezd ezg lg ld l b disc disp;
var ycycle ccycle icycle loancycle;

varexo sh_ag sh_amb sh_ad sh_fg sh_fd;

parameters rho_y nu alpha sigma_h delta beta phi betad betag betab kappa mg md
           gammap gammad gammag gammab epsilon gamma_i delta_m m_bar e_star xi
           rho_ad rho_ag rho_amb rho_fd rho_fg amb_ss ad_ss ag_ss fd_ss fg_ss chi scale_par_fd;

load params_calibration

%Macro and financial parameters
set_param_value('rho_y',rho_y);
set_param_value('nu',nu);
set_param_value('alpha',alpha);
set_param_value('sigma_h',sigma_h);
set_param_value('delta',delta);
set_param_value('beta',beta);
set_param_value('betad',betad);
set_param_value('betag',betag);
set_param_value('betab',betab);
set_param_value('kappa',kappa);
set_param_value('mg',mg);
set_param_value('md',md);
set_param_value('gammap',gammap);
set_param_value('gammad',gammad);
set_param_value('gammag',gammag);
set_param_value('gammab',gammab);
set_param_value('gamma_i',gamma_i);
set_param_value('epsilon',epsilon);
set_param_value('rho_ag',rho_ag);
set_param_value('rho_ad',rho_ad);
set_param_value('rho_amb',rho_amb);
set_param_value('rho_fg',rho_fg);
set_param_value('rho_fd',rho_fd);
set_param_value('phi',phi);



%Climate module
set_param_value('delta_m', delta_m);
set_param_value('m_bar', m_bar);
set_param_value('e_star', e_star); 
set_param_value('xi', xi);
%set_param_value('chi', chi);



model;
c^(-gammap) = lambda;
disc = lambda(+1)/lambda;
hg^sigma_h = lambda*wg;
hd^sigma_h = lambda*wd;
beta*disc*r = 1;
kg = (1 - delta) * kg(-1) + (1-(gamma_i/2)*((invg/invg(-1) - 1)^2))*invg;
kd = (1 - delta) * kd(-1) + (1-(gamma_i/2)*((invd/invd(-1) - 1)^2))*invd;
lg = mg*(qg(+1)/(zg*(1-fg(+1))))*kg;
ld = md*(qd(+1)/(zd*(1-fd(+1))))*kd;
betag*disc*(zg*(1-fg(+1)))-1 = betag*disc*((zg*(1-fg(+1)))/(mg*qg(+1)))*((alpha*pg(+1)*yg(+1))/kg + qg(+1)*(1-delta)) - (qg*zg*(1-fg(+1)))/(mg*qg(+1));
betad*disc*(zd*(1-fd(+1)))-1 = betad*disc*((zd*(1-fd(+1)))/(md*qd(+1)))*((alpha*pd(+1)*yd(+1))/kd + qd(+1)*(1-delta)) - (qd*zd*(1-fd(+1)))/(md*qd(+1));
qg = 1 + qg*(gamma_i/2)*(invg/invg(-1)-1)^2 + qg*gamma_i*(invg/invg(-1)-1)*(invg/invg(-1)) - gamma_i*betag*(disc)*qg(+1)*((invg(+1)/invg -1))*((invg(+1)^2)/(invg^2)); 
qd = 1 + qd*(gamma_i/2)*(invd/invd(-1)-1)^2 + qd*gamma_i*(invd/invd(-1)-1)*(invd/invd(-1)) - gamma_i*betad*(disc)*qd(+1)*((invd(+1)/invd -1))*((invd(+1)^2)/(invd^2)); 
yg=(exp(-xi*(m-m_bar)))*ag*(kg(-1)^alpha)*(hg^(1-alpha));
yd=(exp(-xi*(m-m_bar)))*ad*(kd(-1)^alpha)*(hd^(1-alpha));
wg=(1-alpha)*pg*yg/hg;
wd=(1-alpha)*pd*yd/hd;
e=chi*(yd)^(epsilon);
m-m_bar=(1-delta_m)*(m(-1)-m_bar)+e+e_star;
yd=nu*(tau*pd)^(-rho_y)*y;
yg=(1-nu)*(pg)^(-rho_y)*y;
y=(nu^(1/rho_y)*(yd)^((rho_y-1)/rho_y)+(1-nu)^(1/rho_y)*(yg)^((rho_y-1)/rho_y))^(rho_y/(rho_y-1));
l = ld+lg;
b = (1-kappa)*l;
betab*disc*(zg*(1-fg(+1))) = kappa + (1-kappa)*betab*r*disc;
betab*disc*(zd*(1-fd(+1)))= kappa + (1-kappa)*betab*r*disc;
h=hg+hd;
y=c+invg+invd;
log(ad) = rho_ad*log(ad(-1)) + sh_ad;
log(ag) = rho_ag*log(ag(-1)) - amb(-1) + sh_ag;
%log(amb/amb_ss) = rho_amb*log(amb(-1)/amb_ss) + sh_amb;
amb =amb_ss*exp(zamb);
zamb = rho_amb*zamb(-1) + sh_amb;
fg = 0.6 + sh_fg;
fd = 0.6 + sh_fd;
welf = (c^(1-gammap)/(1-gammap) - hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))+ beta*welf(+1);
tau = (yd/STEADY_STATE(yd))^0;
ezd=(1-fd)*zd;
ezg=(1-fg)*zg;

[name='dispersion observed']
disp = 0.24*2*amb;

[name='GDP observed']
%ycycle = log(y)-log(y(-1));
%ycycle = log(y)-log(STEADY_STATE(y));
ycycle = log(y);

[name='c observed']
%ccycle = log(c)-log(c(-1));
ccycle = log(c)-log(steady_state(c));

inv = invd + invg;
[name='i observed']
icycle = log(inv)-log(steady_state(inv));
%icycle = log(inv)-log(inv(-1));

%[name='fd observed']
%fdcycle = log(fd)-log(fd_ss);
%fdcycle = log(fd)-log(fd(-1));

%[name='fg observed']
%fgcycle = log(fg)-log(fg_ss);
%fgcycle = log(fg)-log(fg(-1));

[name='loan observed']
loancycle = log(l)-log(steady_state(l));
%loancycle = log(l)-log(l(-1));

end;

resid(1);
steady;

check;
model_diagnostics;

shocks;
var sh_amb = 0.16^2;
var sh_ag = 0.013^2;
var sh_ad = 0.013^2;
end;

%stoch_simul(order=2,irf=41, pruning, nograph);
stoch_simul(order=2,pruning,nograph,periods=1000,hp_filter=1600);

mean_welfare=oo_.mean(strmatch('welf',M_.endo_names,'exact'))
steady_state_welfare=oo_.dr.ys(strmatch('welf',M_.endo_names,'exact'))

mean_disp_data= 0.004
mean_disp=oo_.mean(strmatch('disp',M_.endo_names,'exact'))
std_disp_data = 0.0015
std_disp=oo_.var((strmatch('disp',oo_.var_list,'exact')),(strmatch('disp',oo_.var_list,'exact')))^0.5
%autoc_disp=autocorr(disp);

c_ss=oo_.dr.ys(strmatch('c',M_.endo_names,'exact'))
hd_ss=oo_.dr.ys(strmatch('hd',M_.endo_names,'exact'))
hg_ss=oo_.dr.ys(strmatch('hg',M_.endo_names,'exact'))


%welf_ss = (c_ss^(1-gammap)/(1-gammap) - hg_ss^(1+sigma_h)/(1+sigma_h) - hd_ss^(1+sigma_h)/(1+sigma_h))/(1-beta)

fac_ss = (1/(1-beta))*(1/(1-gammap))*(c_ss)^(1-gammap)
fah_ss = (1/(1-beta))*((1/(1+sigma_h))*(hd_ss^(1+sigma_h) + hg_ss^(1+sigma_h)))
welf_cost = (1-((mean_welfare + fah_ss)/fac_ss)^(1/(1-gammap)))*100


if options_.order==2
cond_welfare=oo_.dr.ys(strmatch('welf',M_.endo_names,'exact'))+oo_.dr.ghs2(oo_.dr.inv_order_var(strmatch('welf',M_.endo_names,'exact')))/2
    sigma_h_idx = find(strcmp(M_.param_names, 'sigma_h'));
    sigma_h=M_.params(sigma_h_idx);
    beta_idx = find(strcmp(M_.param_names, 'beta'));
    beta=M_.params(beta_idx);
    %welfare_cost=(1-exp((1-beta)*(oo_.dr.ghs2(oo_.dr.inv_order_var(strmatch('welf',M_.endo_names,'exact')))/2)))*100
    %welfare_cost=(1-(1+((1-beta)*oo_.dr.ghs2(oo_.dr.inv_order_var(strmatch('welf',M_.endo_names,'exact')))/2)/oo_.steady_state(strmatch('welf',M_.endo_names,'exact')))^(1/(sigma_l*(1-gam))))*100
mean_welfare=oo_.mean(strmatch('welf',M_.endo_names,'exact'))
 
%mean_c=(oo_.mean(strmatch('c',M_.endo_names,'exact'))-oo_.dr.ys(strmatch('c',M_.endo_names,'exact')))*100/oo_.steady_state(strmatch('c',M_.endo_names,'exact'))
%std_c=oo_.var((strmatch('c',oo_.var_list,'exact')),(strmatch('c',oo_.var_list,'exact')))^0.5

steady_state_welfare=oo_.dr.ys(strmatch('welf',M_.endo_names,'exact'))

fac_ss = (1/(1-beta))*(1/(1-gammap))*(c_ss)^(1-gammap)
fah_ss = (1/(1-beta))*((1/(1+sigma_h))*(hd_ss^(1+sigma_h) + hg_ss^(1+sigma_h)))
welf_cost = (1-((mean_welfare + fah_ss)/fac_ss)^(1/(1-gammap)))*100

%welfare_cost_bis=(1-exp(mean_welfare - steady_state_welfare))*100
end


%Matching first-order moments
%see file data_calibration
i_gdp_ratio_data=0.2646
i_gdp_ratio=oo_.dr.ys(strmatch('inv',M_.endo_names,'exact'))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))

%Matching second-order moments
std_y_data=0.0101  %0.0121 nominal 
std_y=oo_.var((strmatch('ycycle',M_.endo_names,'exact')),(strmatch('ycycle',M_.endo_names,'exact')))^0.5
std_c_y_ratio_data=0.8497
std_c_y_ratio=oo_.var((strmatch('ccycle',M_.endo_names,'exact')),(strmatch('ccycle',M_.endo_names,'exact')))^0.5/oo_.var((strmatch('ycycle',M_.endo_names,'exact')),(strmatch('ycycle',M_.endo_names,'exact')))^0.5
std_i_y_ratio_data=4.6719
std_i_y_ratio=oo_.var((strmatch('icycle',M_.endo_names,'exact')),(strmatch('icycle',M_.endo_names,'exact')))^0.5/oo_.var((strmatch('ycycle',M_.endo_names,'exact')),(strmatch('ycycle',M_.endo_names,'exact')))^0.5
corr_c_y_data=0.9224
corr_c_y=corr(ccycle, ycycle)
corr_i_y_data=0.8903
corr_i_y=corr(icycle, ycycle)
autoc_y_data=0.9004
autoc_y=autocorr(ycycle)
autoc_c_data=0.8644
autoc_c=autocorr(ccycle)
autoc_i_data=0.8842
autoc_i=autocorr(icycle)
autoc_chargeoff_data=0.8721
%autoc_fd=autocorr(fdcycle)
%autoc_fg=autocorr(fgcycle)

mean_c=(oo_.mean(strmatch('c',M_.endo_names,'exact'))-oo_.dr.ys(strmatch('c',M_.endo_names,'exact')))*100/oo_.steady_state(strmatch('c',M_.endo_names,'exact'));
std_c=oo_.var((strmatch('c',oo_.var_list,'exact')),(strmatch('c',oo_.var_list,'exact')))^0.5;
std_y=oo_.var((strmatch('y',oo_.var_list,'exact')),(strmatch('y',oo_.var_list,'exact')))^0.5;
std_y_net=oo_.var((strmatch('y_net',oo_.var_list,'exact')),(strmatch('y_net',oo_.var_list,'exact')))^0.5;
std_invd=oo_.var((strmatch('invd',oo_.var_list,'exact')),(strmatch('invd',oo_.var_list,'exact')))^0.5;
std_invg=oo_.var((strmatch('invg',oo_.var_list,'exact')),(strmatch('invg',oo_.var_list,'exact')))^0.5;


mean_disp_data= 0.004
mean_disp=oo_.mean(strmatch('disp',M_.endo_names,'exact'))
std_disp_data = 0.0015
std_disp=oo_.var((strmatch('disp',oo_.var_list,'exact')),(strmatch('disp',oo_.var_list,'exact')))^0.5
%autoc_disp=autocorr(disp);

