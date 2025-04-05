%%%GREEN FINANCE AND AMBIGUITY%%%%%
%%%MARCO CARLI%%%%%

format long

var yd yg pd pg kg wg hg kd wd hd y lambda c invd invg r qd qg m e h ag amb zamb ad
    welf taud taug fd fg zd zg ezd ezg lg ld l b disc spreadd spreadg diffz kappa_flu weightd_flu weightg_flu md_flu mg_flu;

varexo sh_ag sh_amb sh_ad sh_fg sh_fd;

parameters rho_y nu alpha sigma_h delta beta phi betad betag betab kappa weightd weightg mg md
           gammap gammad gammag gammab epsilon gamma_i delta_m m_bar xi psi
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
%set_param_value('weightg',weightg);
%set_param_value('weightd',weightd);
weightd=1;
weightg=1;
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
%set_param_value('e_star', e_star); 
set_param_value('xi', xi);
%set_param_value('chi', chi);


welf_mat=[];
%psi_vec=[-10:0.01:-9]';
%psi_vec=[-10:1:10]';
psi_vec=0;

for jjj=1:size(psi_vec,1);
psi=psi_vec(jjj,1);


model;
c^(-gammap) = lambda;
disc = lambda(+1)/lambda;
hg^sigma_h = lambda*wg;
hd^sigma_h = lambda*wd;
beta*disc*r = 1;
kg = (1 - delta) * kg(-1) + (1-(gamma_i/2)*((invg/invg(-1) - 1)^2))*invg;
kd = (1 - delta) * kd(-1) + (1-(gamma_i/2)*((invd/invd(-1) - 1)^2))*invd;
lg = mg_flu*(qg(+1)/(zg*(1-fg(+1))))*kg;
ld = md_flu*(qd(+1)/(zd*(1-fd(+1))))*kd;
betag*disc*(zg*(1-fg(+1)))-1 = betag*disc*((zg*(1-fg(+1)))/(mg_flu*qg(+1)))*((alpha*pg(+1)*yg(+1))/kg + qg(+1)*(1-delta)) - (qg*zg*(1-fg(+1)))/(mg_flu*qg(+1));
betad*disc*(zd*(1-fd(+1)))-1 = betad*disc*((zd*(1-fd(+1)))/(md_flu*qd(+1)))*((alpha*pd(+1)*yd(+1))/kd + qd(+1)*(1-delta)) - (qd*zd*(1-fd(+1)))/(md_flu*qd(+1));
qg = 1 + qg*(gamma_i/2)*(invg/invg(-1)-1)^2 + qg*gamma_i*(invg/invg(-1)-1)*(invg/invg(-1)) - gamma_i*betag*(disc)*qg(+1)*((invg(+1)/invg -1))*((invg(+1)^2)/(invg^2)); 
qd = 1 + qd*(gamma_i/2)*(invd/invd(-1)-1)^2 + qd*gamma_i*(invd/invd(-1)-1)*(invd/invd(-1)) - gamma_i*betad*(disc)*qd(+1)*((invd(+1)/invd -1))*((invd(+1)^2)/(invd^2)); 
yg=(exp(-xi*(m-m_bar)))*ag*(kg(-1)^alpha)*(hg^(1-alpha));
yd=(exp(-xi*(m-m_bar)))*ad*(kd(-1)^alpha)*(hd^(1-alpha));
wg=(1-alpha)*pg*yg/hg;
wd=(1-alpha)*pd*yd/hd;
e=chi*(yd)^(epsilon);
m-m_bar=(1-delta_m)*(m(-1)-m_bar)+e;
yd=nu*((1+taud)*pd)^(-rho_y)*y;
yg=(1-nu)*((1-taug)*pg)^(-rho_y)*y;
y=(nu^(1/rho_y)*(yd)^((rho_y-1)/rho_y)+(1-nu)^(1/rho_y)*(yg)^((rho_y-1)/rho_y))^(rho_y/(rho_y-1));
l = ld+lg;
b = (1-kappa_flu*weightd_flu)*ld + (1-kappa_flu*weightg_flu)*lg;
betab*disc*(zg*(1-fg(+1))) = kappa_flu*weightg_flu + (1-kappa_flu*weightg_flu)*betab*r*disc;
betab*disc*(zd*(1-fd(+1)))= kappa_flu*weightd_flu + (1-kappa_flu*weightd_flu)*betab*r*disc;
h=hg+hd;
y=c+invg+invd;
log(ad) = rho_ad*log(ad(-1)) + sh_ad;
log(ag) = rho_ag*log(ag(-1)) - amb(-1) + sh_ag;
%log(amb/amb_ss) = rho_amb*log(amb(-1)/amb_ss) + sh_amb;
amb =amb_ss*exp(zamb);
zamb = rho_amb*zamb(-1) + sh_amb;
fg = 0.35 + sh_fg;
fd = 0.35 + sh_fd;
welf = (c^(1-gammap)/(1-gammap) - hg^(1+sigma_h)/(1+sigma_h) - hd^(1+sigma_h)/(1+sigma_h))+ beta*welf(+1);
ezd=(1-fd)*zd;
ezg=(1-fg)*zg;
spreadd=zd/r;
spreadg=zg/r;
diffz=zg/zd;

%taud=0;
taud = (e/STEADY_STATE(e))^(0)-1; %psi_e=5.37 
taug*pg*yg = taud*pd*yd;

%kappa_flu = (kappa)*(l/l(-1))^(psi); 
kappa_flu = kappa*(l/steady_state(l))^(0); %psi_ld=4 %psi_lg=11 %psi_l=7.11 %psi_qd=1 %psi_qg=1 %psi_y=26

%weightd_flu = (weightd)*(l/l(-1))^(psi); 
weightd_flu = (weightd)*(ld/steady_state(ld))^(0); %psi_ld=4.25 %psi_lg=7 %psi_l=7 %psi_qd=1

%weightg_flu = (weightg)*(qg/qg(-1))^(0); 
weightg_flu = (weightg)*(lg/steady_state(lg))^(0); %psi_ld=5 %psi_lg=20.46 %psi_l=9 %psi_qg=2

%md_flu = md*(yd/yd(-1))^(psi);
md_flu = md*(yd/steady_state(yd))^(0); %psi_yd=-9.65

%mg_flu = mg*(yg/yg(-1))^(psi);
mg_flu = mg*(yg/steady_state(yg))^(psi); %psi_yg=-3.93 %psi_y=-18

end;

resid;
steady;

check;
model_diagnostics;

shocks;
var sh_amb = 0.16^2;
%var sh_ag = 0.013^2;
%var sh_ad = 0.013^2;
end;


%stoch_simul(order=2,irf=41,pruning);
stoch_simul(order=1,irf=41,pruning,nograph);
%stoch_simul(order=2,pruning, nograph);
%stoch_simul(order=1,irf=41,pruning, nograph);
%stoch_simul(order=2,pruning,periods=2000,drop=1000,simul_replic=1,nograph);


%ss_noamb = oo_.dr.ys;
%save ss_noamb ss_noamb;
load ss_noamb;

%ss_amb_nopol = oo_.dr.ys;
%save ss_amb_nopol ss_amb_nopol;
load ss_amb_nopol;

ss_amb = oo_.dr.ys;
save ss_amb ss_amb;
%load ss_amb;


if options_.order==2

load mean_c_base.mat
load mean_hd_base.mat
load mean_hg_base.mat
load mean_e_base.mat
load mean_welfare_base.mat
load mean_welfare2_base.mat
load steady_state_welfare_base.mat
load c_ss_baseline.mat
load hd_ss_baseline.mat
load hg_ss_baseline.mat

mean_welfare=oo_.mean(strmatch('welf',M_.endo_names,'exact'))
mean_welfare2=oo_.dr.ys(strmatch('welf',M_.endo_names,'exact'))+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(strmatch('welf',M_.endo_names,'exact')))
steady_state_welfare=oo_.dr.ys(strmatch('welf',M_.endo_names,'exact'))

c_ss=oo_.dr.ys(strmatch('c',M_.endo_names,'exact'));
hd_ss=oo_.dr.ys(strmatch('hd',M_.endo_names,'exact'));
hg_ss=oo_.dr.ys(strmatch('hg',M_.endo_names,'exact'));
e_ss=oo_.dr.ys(strmatch('e',M_.endo_names,'exact'))
mean_c=oo_.mean(strmatch('c',M_.endo_names,'exact'))
mean_yd=oo_.mean(strmatch('yd',M_.endo_names,'exact'))
mean_yg=oo_.mean(strmatch('yg',M_.endo_names,'exact'))
mean_kd=oo_.mean(strmatch('kd',M_.endo_names,'exact'))
mean_kg=oo_.mean(strmatch('kg',M_.endo_names,'exact'))
mean_hd=oo_.mean(strmatch('hd',M_.endo_names,'exact'))
mean_hg=oo_.mean(strmatch('hg',M_.endo_names,'exact'))
mean_e=oo_.mean(strmatch('e',M_.endo_names,'exact'))

std_c=oo_.var((strmatch('c',M_.endo_names,'exact')),(strmatch('c',M_.endo_names,'exact')))^0.5;
std_hd=oo_.var((strmatch('hd',M_.endo_names,'exact')),(strmatch('hd',M_.endo_names,'exact')))^0.5;
std_hg=oo_.var((strmatch('hg',M_.endo_names,'exact')),(strmatch('hg',M_.endo_names,'exact')))^0.5;
std_e=oo_.var((strmatch('e',M_.endo_names,'exact')),(strmatch('e',M_.endo_names,'exact')))^0.5;
std_kd=oo_.var((strmatch('kd',M_.endo_names,'exact')),(strmatch('kd',M_.endo_names,'exact')))^0.5;
std_kg=oo_.var((strmatch('kg',M_.endo_names,'exact')),(strmatch('kg',M_.endo_names,'exact')))^0.5;
std_zd=oo_.var((strmatch('zd',M_.endo_names,'exact')),(strmatch('zd',M_.endo_names,'exact')))^0.5;
std_zg=oo_.var((strmatch('zg',M_.endo_names,'exact')),(strmatch('zg',M_.endo_names,'exact')))^0.5;
std_fd=oo_.var((strmatch('fd',M_.endo_names,'exact')),(strmatch('fd',M_.endo_names,'exact')))^0.5;
std_fg=oo_.var((strmatch('fg',M_.endo_names,'exact')),(strmatch('fg',M_.endo_names,'exact')))^0.5;
std_welf=oo_.var((strmatch('welf',M_.endo_names,'exact')),(strmatch('welf',M_.endo_names,'exact')))^0.5;


fac_ss = (1/(1-beta))*(1/(1-gammap))*(c_ss)^(1-gammap);
fah_ss = (1/(1-beta))*((1/(1+sigma_h))*(hd_ss^(1+sigma_h) + hg_ss^(1+sigma_h)));
welf_cost_flu = (1-((mean_welfare + fah_ss)/fac_ss)^(1/(1-gammap)))*100


end


welf_mat(jjj,1) = welf_cost_flu;
%welf_mat(jjj,1) = mean_welfare;
end

welf_opt = min(welf_mat)
[idx_psi] = find(welf_mat==welf_opt);
psi_opt = psi_vec(idx_psi, 1)



%/*
yd_sh_amb2 = (yd_sh_amb./ss_amb(strmatch('yd',M_.endo_names,'exact'))).*ss_noamb(strmatch('yd',M_.endo_names,'exact'));
yg_sh_amb2 = (yg_sh_amb./ss_amb(strmatch('yg',M_.endo_names,'exact'))).*ss_noamb(strmatch('yg',M_.endo_names,'exact'));
pd_sh_amb2 = (pd_sh_amb./ss_amb(strmatch('pd',M_.endo_names,'exact'))).*ss_noamb(strmatch('pd',M_.endo_names,'exact'));
pg_sh_amb2 = (pg_sh_amb./ss_amb(strmatch('pg',M_.endo_names,'exact'))).*ss_noamb(strmatch('pg',M_.endo_names,'exact'));
kg_sh_amb2 = (kg_sh_amb./ss_amb(strmatch('kg',M_.endo_names,'exact'))).*ss_noamb(strmatch('kg',M_.endo_names,'exact'));
wg_sh_amb2 = (wg_sh_amb./ss_amb(strmatch('wg',M_.endo_names,'exact'))).*ss_noamb(strmatch('wg',M_.endo_names,'exact'));
hg_sh_amb2 = (hg_sh_amb./ss_amb(strmatch('hg',M_.endo_names,'exact'))).*ss_noamb(strmatch('hg',M_.endo_names,'exact'));
kd_sh_amb2 = (kd_sh_amb./ss_amb(strmatch('kd',M_.endo_names,'exact'))).*ss_noamb(strmatch('kd',M_.endo_names,'exact'));
wd_sh_amb2 = (wd_sh_amb./ss_amb(strmatch('wd',M_.endo_names,'exact'))).*ss_noamb(strmatch('wd',M_.endo_names,'exact'));
hd_sh_amb2 = (hd_sh_amb./ss_amb(strmatch('hd',M_.endo_names,'exact'))).*ss_noamb(strmatch('hd',M_.endo_names,'exact'));
y_sh_amb2 = (y_sh_amb./ss_amb(strmatch('y',M_.endo_names,'exact'))).*ss_noamb(strmatch('y',M_.endo_names,'exact'));
lambda_sh_amb2 = (lambda_sh_amb./ss_amb(strmatch('lambda',M_.endo_names,'exact'))).*ss_noamb(strmatch('lambda',M_.endo_names,'exact'));
c_sh_amb2 = (c_sh_amb./ss_amb(strmatch('c',M_.endo_names,'exact'))).*ss_noamb(strmatch('c',M_.endo_names,'exact'));
invd_sh_amb2 = (invd_sh_amb./ss_amb(strmatch('invd',M_.endo_names,'exact'))).*ss_noamb(strmatch('invd',M_.endo_names,'exact'));
invg_sh_amb2 = (invg_sh_amb./ss_amb(strmatch('invg',M_.endo_names,'exact'))).*ss_noamb(strmatch('invg',M_.endo_names,'exact'));
r_sh_amb2 = (r_sh_amb./ss_amb(strmatch('r',M_.endo_names,'exact'))).*ss_noamb(strmatch('r',M_.endo_names,'exact'));
ad_sh_amb2 = (ad_sh_amb./ss_amb(strmatch('ad',M_.endo_names,'exact'))).*ss_noamb(strmatch('ad',M_.endo_names,'exact'));
ag_sh_amb2 = (ag_sh_amb./ss_amb(strmatch('ag',M_.endo_names,'exact'))).*ss_noamb(strmatch('ag',M_.endo_names,'exact'));
m_sh_amb2 = (m_sh_amb./ss_amb(strmatch('m',M_.endo_names,'exact'))).*ss_noamb(strmatch('m',M_.endo_names,'exact'));
e_sh_amb2 = (e_sh_amb./ss_amb(strmatch('e',M_.endo_names,'exact'))).*ss_noamb(strmatch('e',M_.endo_names,'exact'));
h_sh_amb2 = (h_sh_amb./ss_amb(strmatch('h',M_.endo_names,'exact'))).*ss_noamb(strmatch('h',M_.endo_names,'exact'));
ag_sh_amb2 = (ag_sh_amb./ss_amb(strmatch('ag',M_.endo_names,'exact'))).*ss_noamb(strmatch('ag',M_.endo_names,'exact'));
%zamb_sh_amb2 = (zamb_sh_amb./ss_amb(strmatch('zamb',M_.endo_names,'exact'))).*ss_noamb(strmatch('zamb',M_.endo_names,'exact'));
ad_sh_amb2 = (ad_sh_amb./ss_amb(strmatch('ad',M_.endo_names,'exact'))).*ss_noamb(strmatch('ad',M_.endo_names,'exact'));
qd_sh_amb2 = (qd_sh_amb./ss_amb(strmatch('qd',M_.endo_names,'exact'))).*ss_noamb(strmatch('qd',M_.endo_names,'exact'));
qg_sh_amb2 = (qg_sh_amb./ss_amb(strmatch('qg',M_.endo_names,'exact'))).*ss_noamb(strmatch('qg',M_.endo_names,'exact'));
fd_sh_amb2 = (fd_sh_amb./ss_amb(strmatch('fd',M_.endo_names,'exact'))).*ss_noamb(strmatch('fd',M_.endo_names,'exact'));
fg_sh_amb2 = (fg_sh_amb./ss_amb(strmatch('fg',M_.endo_names,'exact'))).*ss_noamb(strmatch('fg',M_.endo_names,'exact'));
fd_sh_amb2 = (fd_sh_amb./ss_amb(strmatch('fd',M_.endo_names,'exact'))).*ss_noamb(strmatch('fd',M_.endo_names,'exact'));
fg_sh_amb2 = (fg_sh_amb./ss_amb(strmatch('fg',M_.endo_names,'exact'))).*ss_noamb(strmatch('fg',M_.endo_names,'exact'));
zd_sh_amb2 = (zd_sh_amb./ss_amb(strmatch('zd',M_.endo_names,'exact'))).*ss_noamb(strmatch('zd',M_.endo_names,'exact'));
zg_sh_amb2 = (zg_sh_amb./ss_amb(strmatch('zg',M_.endo_names,'exact'))).*ss_noamb(strmatch('zg',M_.endo_names,'exact'));
ezd_sh_amb2 = (ezd_sh_amb./ss_amb(strmatch('ezd',M_.endo_names,'exact'))).*ss_noamb(strmatch('ezd',M_.endo_names,'exact'));
ezg_sh_amb2 = (ezg_sh_amb./ss_amb(strmatch('ezg',M_.endo_names,'exact'))).*ss_noamb(strmatch('ezg',M_.endo_names,'exact'));
ld_sh_amb2 = (ld_sh_amb./ss_amb(strmatch('ld',M_.endo_names,'exact'))).*ss_noamb(strmatch('ld',M_.endo_names,'exact'));
lg_sh_amb2 = (lg_sh_amb./ss_amb(strmatch('lg',M_.endo_names,'exact'))).*ss_noamb(strmatch('lg',M_.endo_names,'exact'));
l_sh_amb2 = (l_sh_amb./ss_amb(strmatch('l',M_.endo_names,'exact'))).*ss_noamb(strmatch('l',M_.endo_names,'exact'));
b_sh_amb2 = (b_sh_amb./ss_amb(strmatch('b',M_.endo_names,'exact'))).*ss_noamb(strmatch('b',M_.endo_names,'exact'));
spreadd_sh_amb2 = (spreadd_sh_amb./ss_amb(strmatch('spreadd',M_.endo_names,'exact'))).*ss_noamb(strmatch('spreadd',M_.endo_names,'exact'));
spreadg_sh_amb2 = (spreadg_sh_amb./ss_amb(strmatch('spreadg',M_.endo_names,'exact'))).*ss_noamb(strmatch('spreadg',M_.endo_names,'exact'));
diffz_sh_amb2 = (diffz_sh_amb./ss_amb(strmatch('diffz',M_.endo_names,'exact'))).*ss_noamb(strmatch('diffz',M_.endo_names,'exact'));
%*/



//save results
vec_ss = [ss_noamb(strmatch('yd',M_.endo_names,'exact')) ss_noamb(strmatch('yg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('pd',M_.endo_names,'exact')) ss_noamb(strmatch('pg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('kg',M_.endo_names,'exact')) ss_noamb(strmatch('wg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('hg',M_.endo_names,'exact')) ss_noamb(strmatch('kd',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('wd',M_.endo_names,'exact')) ss_noamb(strmatch('hd',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('y',M_.endo_names,'exact'))  ss_noamb(strmatch('lambda',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('c',M_.endo_names,'exact'))  ss_noamb(strmatch('invd',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('invg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('r',M_.endo_names,'exact')) ss_noamb(strmatch('ad',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('ag',M_.endo_names,'exact')) ss_noamb(strmatch('m',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('e',M_.endo_names,'exact')) ss_noamb(strmatch('h',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('ag',M_.endo_names,'exact')) ss_noamb(strmatch('amb',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('ad',M_.endo_names,'exact')) ss_noamb(strmatch('zamb',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('qd',M_.endo_names,'exact')) ss_noamb(strmatch('qg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('fd',M_.endo_names,'exact')) ss_noamb(strmatch('fg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('fd',M_.endo_names,'exact')) ss_noamb(strmatch('fg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('zd',M_.endo_names,'exact')) ss_noamb(strmatch('zg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('ezd',M_.endo_names,'exact')) ss_noamb(strmatch('ezg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('ld',M_.endo_names,'exact')) ss_noamb(strmatch('lg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('l',M_.endo_names,'exact')) ss_noamb(strmatch('b',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('spreadd',M_.endo_names,'exact')) ss_noamb(strmatch('spreadg',M_.endo_names,'exact')) ...
          ss_noamb(strmatch('diffz',M_.endo_names,'exact')) oo_.dr.ys(strmatch('md_flu',M_.endo_names,'exact')) ...
          oo_.dr.ys(strmatch('mg_flu',M_.endo_names,'exact')) oo_.dr.ys(strmatch('kappa_flu',M_.endo_names,'exact')) ...
          oo_.dr.ys(strmatch('weightd_flu',M_.endo_names,'exact')) oo_.dr.ys(strmatch('weightg_flu',M_.endo_names,'exact'))];


mat_sh_amb2 = [yd_sh_amb2 yg_sh_amb2 pd_sh_amb2 pg_sh_amb2 kg_sh_amb2 wg_sh_amb2 hg_sh_amb2 kd_sh_amb2 wd_sh_amb2 ... 
              hd_sh_amb2 y_sh_amb2 lambda_sh_amb2 c_sh_amb2 invd_sh_amb2 invg_sh_amb2 r_sh_amb2 ad_sh_amb2 ag_sh_amb2 ...
              m_sh_amb2 e_sh_amb2 h_sh_amb2 ag_sh_amb2 amb_sh_amb ad_sh_amb2 zamb_sh_amb qd_sh_amb2 qg_sh_amb2 ...
              fd_sh_amb2 fg_sh_amb2 fd_sh_amb2 fg_sh_amb2 zd_sh_amb2 zg_sh_amb2 ezd_sh_amb2 ezg_sh_amb2 ....
              ld_sh_amb2 lg_sh_amb2 l_sh_amb2 b_sh_amb2 spreadd_sh_amb2 spreadg_sh_amb2 diffz_sh_amb2 ....
              md_flu_sh_amb mg_flu_sh_amb kappa_flu_sh_amb weightd_flu_sh_amb weightg_flu_sh_amb];


/*
IRFs_sh_amb_w=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_w(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_w(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_w_dyn_fin_results', 'IRFs_sh_amb_w')
*/


/*
IRFs_sh_amb=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_dyn_fin_results', 'IRFs_sh_amb')
*/






/*
IRFs_sh_amb_kappa=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_kappa(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_kappa(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_kappa_dyn_fin_results', 'IRFs_sh_amb_kappa')
*/



/*
IRFs_sh_amb_kappa_lg=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_kappa_lg(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_kappa_lg(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_kappa_lg_dyn_fin_results', 'IRFs_sh_amb_kappa_lg')
*/



/*
IRFs_sh_amb_weightd=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_weightd(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_weightd(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_weightd_dyn_fin_results', 'IRFs_sh_amb_weightd')
*/



/*
IRFs_sh_amb_weightg=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_weightg(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_weightg(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_weightg_dyn_fin_results', 'IRFs_sh_amb_weightg')
*/


/*
IRFs_sh_amb_md=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_md(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_md(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_md_dyn_fin_results', 'IRFs_sh_amb_md')
*/



/*
IRFs_sh_amb_mg=[];

for i = 1:size(mat_sh_amb2,2)
IRFs_sh_amb_mg(:,i)=(mat_sh_amb2(:,i)/vec_ss(i))*100;
if i == 16 | i == 17 | i == 18 | i == 23 | i == 25 | i == 28 | i == 29 | i == 30 | i == 31 | i == 32 | i == 33 | i == 34 | i == 35 | i == 40 | i == 41 | i == 42 | i == 43 | i == 44 | i == 45 | i == 46 | i == 47
IRFs_sh_amb_mg(:,i)=(mat_sh_amb2(:,i))*100;
end
end

save('amb_mg_dyn_fin_results', 'IRFs_sh_amb_mg')
*/









/*
IRFs_sh_amb_w=[];
IRFs_sh_amb_w(:,1)=yd_sh_amb2*100/ss_noamb(strmatch('yd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,2)=yg_sh_amb2*100/ss_noamb(strmatch('yg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,3)=pd_sh_amb2*100/ss_noamb(strmatch('pd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,4)=pg_sh_amb2*100/ss_noamb(strmatch('pg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,5)=kg_sh_amb2*100/ss_noamb(strmatch('kg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,6)=wg_sh_amb2*100/ss_noamb(strmatch('wg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,7)=hg_sh_amb2*100/ss_noamb(strmatch('hg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,8)=kd_sh_amb2*100/ss_noamb(strmatch('kd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,9)=wd_sh_amb2*100/ss_noamb(strmatch('wd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,10)=hd_sh_amb2*100/ss_noamb(strmatch('hd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,11)=y_sh_amb2*100/ss_noamb(strmatch('y',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,12)=lambda_sh_amb2*100/ss_noamb(strmatch('lambda',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,13)=c_sh_amb2*100/ss_noamb(strmatch('c',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,14)=invd_sh_amb2*100/ss_noamb(strmatch('invd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,15)=invg_sh_amb2*100/ss_noamb(strmatch('invg',M_.endo_names,'exact'));;
IRFs_sh_amb_w(:,16)=r_sh_amb2*100;
IRFs_sh_amb_w(:,17)=zad_sh_amb2*100;
IRFs_sh_amb_w(:,18)=zag_sh_amb2*100;
IRFs_sh_amb_w(:,19)=m_sh_amb2*100/ss_noamb(strmatch('m',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,20)=e_sh_amb2*100/ss_noamb(strmatch('e',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,21)=h_sh_amb2*100/ss_noamb(strmatch('h',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,22)=ag_sh_amb2*100/ss_noamb(strmatch('ag',M_.endo_names,'exact'));
%IRFs_sh_amb_w(:,23)=amb_sh_amb*100/ss_noamb(strmatch('amb',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,23)=amb_sh_amb*100;
IRFs_sh_amb_w(:,24)=ad_sh_amb2*100/ss_noamb(strmatch('ad',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,25)=zamb_sh_amb*100;
IRFs_sh_amb_w(:,26)=qd_sh_amb2*100/ss_noamb(strmatch('qd',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,27)=qg_sh_amb2*100/ss_noamb(strmatch('qg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,28)=fd_sh_amb2*100;
IRFs_sh_amb_w(:,29)=fg_sh_amb2*100;
IRFs_sh_amb_w(:,30)=zfd_sh_amb2*100;
IRFs_sh_amb_w(:,31)=zfg_sh_amb2*100;
IRFs_sh_amb_w(:,32)=zd_sh_amb2*100;
IRFs_sh_amb_w(:,33)=zg_sh_amb2*100;
IRFs_sh_amb_w(:,34)=ezd_sh_amb2*100;
IRFs_sh_amb_w(:,35)=ezg_sh_amb2*100;
IRFs_sh_amb_w(:,36)=ld_sh_amb2*100/ss_noamb(strmatch('ld',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,37)=lg_sh_amb2*100/ss_noamb(strmatch('lg',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,38)=l_sh_amb2*100/ss_noamb(strmatch('l',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,39)=b_sh_amb2*100/ss_noamb(strmatch('b',M_.endo_names,'exact'));
IRFs_sh_amb_w(:,40)=spreadd_sh_amb2*100;
IRFs_sh_amb_w(:,41)=spreadg_sh_amb2*100;
IRFs_sh_amb_w(:,42)=diffz_sh_amb2*100;


save('amb_w_dyn_fin_results','IRFs_sh_amb_w')
*/


/*
IRFs_sh_amb=[];
IRFs_sh_amb(:,1)=yd_sh_amb2*100/ss_noamb(strmatch('yd',M_.endo_names,'exact'));
IRFs_sh_amb(:,2)=yg_sh_amb2*100/ss_noamb(strmatch('yg',M_.endo_names,'exact'));
IRFs_sh_amb(:,3)=pd_sh_amb2*100/ss_noamb(strmatch('pd',M_.endo_names,'exact'));
IRFs_sh_amb(:,4)=pg_sh_amb2*100/ss_noamb(strmatch('pg',M_.endo_names,'exact'));
IRFs_sh_amb(:,5)=kg_sh_amb2*100/ss_noamb(strmatch('kg',M_.endo_names,'exact'));
IRFs_sh_amb(:,6)=wg_sh_amb2*100/ss_noamb(strmatch('wg',M_.endo_names,'exact'));
IRFs_sh_amb(:,7)=hg_sh_amb2*100/ss_noamb(strmatch('hg',M_.endo_names,'exact'));
IRFs_sh_amb(:,8)=kd_sh_amb2*100/ss_noamb(strmatch('kd',M_.endo_names,'exact'));
IRFs_sh_amb(:,9)=wd_sh_amb2*100/ss_noamb(strmatch('wd',M_.endo_names,'exact'));
IRFs_sh_amb(:,10)=hd_sh_amb2*100/ss_noamb(strmatch('hd',M_.endo_names,'exact'));
IRFs_sh_amb(:,11)=y_sh_amb2*100/ss_noamb(strmatch('y',M_.endo_names,'exact'));
IRFs_sh_amb(:,12)=lambda_sh_amb2*100/ss_noamb(strmatch('lambda',M_.endo_names,'exact'));
IRFs_sh_amb(:,13)=c_sh_amb2*100/ss_noamb(strmatch('c',M_.endo_names,'exact'));
IRFs_sh_amb(:,14)=invd_sh_amb2*100/ss_noamb(strmatch('invd',M_.endo_names,'exact'));
IRFs_sh_amb(:,15)=invg_sh_amb2*100/ss_noamb(strmatch('invg',M_.endo_names,'exact'));;
IRFs_sh_amb(:,16)=r_sh_amb2*100;
IRFs_sh_amb(:,17)=zad_sh_amb2*100;
IRFs_sh_amb(:,18)=zag_sh_amb2*100;
IRFs_sh_amb(:,19)=m_sh_amb2*100/ss_noamb(strmatch('m',M_.endo_names,'exact'));
IRFs_sh_amb(:,20)=e_sh_amb2*100/ss_noamb(strmatch('e',M_.endo_names,'exact'));
IRFs_sh_amb(:,21)=h_sh_amb2*100/ss_noamb(strmatch('h',M_.endo_names,'exact'));
IRFs_sh_amb(:,22)=ag_sh_amb2*100/ss_noamb(strmatch('ag',M_.endo_names,'exact'));
%IRFs_sh_amb(:,23)=amb_sh_amb*100/0.05;
IRFs_sh_amb(:,23)=amb_sh_amb*100;
IRFs_sh_amb(:,24)=ad_sh_amb2*100/ss_noamb(strmatch('ad',M_.endo_names,'exact'));
IRFs_sh_amb(:,25)=zamb_sh_amb*100;
IRFs_sh_amb(:,26)=qd_sh_amb2*100/ss_noamb(strmatch('qd',M_.endo_names,'exact'));
IRFs_sh_amb(:,27)=qg_sh_amb2*100/ss_noamb(strmatch('qg',M_.endo_names,'exact'));
IRFs_sh_amb(:,28)=fd_sh_amb2*100;
IRFs_sh_amb(:,29)=fg_sh_amb2*100;
IRFs_sh_amb(:,30)=zfd_sh_amb2*100;
IRFs_sh_amb(:,31)=zfg_sh_amb*100;
IRFs_sh_amb(:,32)=zd_sh_amb2*100;
IRFs_sh_amb(:,33)=zg_sh_amb2*100;
IRFs_sh_amb(:,34)=ezd_sh_amb2*100;
IRFs_sh_amb(:,35)=ezg_sh_amb2*100;
IRFs_sh_amb(:,36)=ld_sh_amb2*100/ss_noamb(strmatch('ld',M_.endo_names,'exact'));
IRFs_sh_amb(:,37)=lg_sh_amb2*100/ss_noamb(strmatch('lg',M_.endo_names,'exact'));
IRFs_sh_amb(:,38)=l_sh_amb2*100/ss_noamb(strmatch('l',M_.endo_names,'exact'));
IRFs_sh_amb(:,39)=b_sh_amb2*100/ss_noamb(strmatch('b',M_.endo_names,'exact'));
IRFs_sh_amb(:,40)=spreadd_sh_amb2*100;
IRFs_sh_amb(:,41)=spreadg_sh_amb2*100;
IRFs_sh_amb(:,42)=diffz_sh_amb2*100;


save('amb_dyn_fin_results','IRFs_sh_amb')
*/
