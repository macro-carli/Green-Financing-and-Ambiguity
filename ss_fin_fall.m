run calibration

m = 1650.537;
%chi = 0.342945253112513;
e = delta_m*(m-m_bar)-e_star;

zamb=0;
amb = 0.005; %0.005
%amb = 0.00001;
zag = (-amb/(1-rho_ag));
ag = 1*exp(zag);
zad = 0;
ad = 1*exp(zad);

r = beta^(-1);

% %Spread loan rate/interest rate target = 2;
% spread_ss = 2;
% targ_z = r_ss+spread_ss;

%scale_par_fd = log(0.77)*(1-rho_fd)/(-ad); 
scale_par_fd =1;
zfd =(-ad*scale_par_fd/(1-rho_fd));
fd = 1*exp(zfd);

zfg=(-ag*scale_par_fd/(1-rho_fg));
fg = 1*exp(zfg);

zd = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fd);
zg = (kappa*(1/betab) + (1-kappa)*(1/beta))/(1-fg);


disc=1;
discd=1;
discg=1;
discb=1;
qd=1;
qg=1;
yd=1.236399093224243;
yg=0.763509963014383;
y = (nu^(1/rho_y)*yd^((rho_y-1)/rho_y) + (1-nu)^(1/rho_y)*yg^((rho_y-1)/rho_y))^(rho_y/(rho_y-1));
pd = ((nu)^(1/rho_y))*((yd/y)^(-1/rho_y));
pg = ((1-nu)^(1/rho_y))*((yg/y)^(-1/rho_y));
kd = (alpha*pd*yd)/(1/betad - (1-delta) - md/(zd*(1-fd)*betad) +md);
kg = (alpha*pg*yg)/(1/betag - (1-delta) - mg/(zg*(1-fg)*betag) +mg);
invd = delta*kd;
invg= delta*kg;
lg = mg*(qg/(zg*(1-fg)))*kg;
ld = md*(qd/(zd*(1-fd)))*kd;
l = ld+lg;
b = (1-kappa)*l;
%pib = zd*(1-fd)*ld + zg*(1-fg)*lg + b - r*b - l;
%pid = pd*yd + ld - invd - (1-alpha)*pd*yd - zd*(1-fd)*ld;
%pig = pg*yg + lg - invg - (1-alpha)*pg*yg - zg*(1-fg)*lg;
%c = y-pib-pig-pid-invg-invd;
c = y-invg-invd;
%lambdab = pib^(-gammab);
%lambdad = pid^(-gammad);
%lambdag = pig^(-gammag);
lambda = c^(-gammap);
hd = ((1-alpha)*(pd*yd*lambda))^(1/(sigma_h+1));
hg = ((1-alpha)*(pg*yg*lambda))^(1/(sigma_h+1));
h = hg+hd;
wd = ((1-alpha)*pd*yd)/hd;
wg = ((1-alpha)*pg*yg)/hg;
ad = yd/(exp(-xi*(m-m_bar))*(kd^(alpha))*(hd^(1-alpha)));
ag = yg/(exp(-xi*(m-m_bar))*(kg^(alpha))*(hg^(1-alpha)));
amb=log(ag)*(1-rho_ag);
chi = e/(yd^(epsilon));


