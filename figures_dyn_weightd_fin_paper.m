%%%%%%%FIGURES GREEN AMBIGUITY%%%%%%%%%
load amb_w_dyn_fin_results;
load amb_dyn_fin_results;
%load amb_w_kappa2_dyn_fin_results;
load amb_weightd_dyn_fin_results;

TT=size(IRFs_sh_amb_w,1)-1;
time=0:1:TT;

iyd=1;
iyg=2;
ipd=3;
ipg=4;
ikg=5;
iwg=6;
ihg=7;
ikd=8;
iwd=9;
ihd=10;
iy=11;
ilambda=12;
ic=13;
iinvd=14;
iinvg=15;
ir=16;
izad=17;
izag=18;
im=19;
ie=20;
ih=21;
iag=22;
iamb=23;
iad=24;
izamb=25;
iqd=26;
iqg=27;
ifd=28;
ifg=29;
izfd=30;
izfg=31;
izd=32;
izg=33;
inzd=34;
inzg=35;
ild=36;
ilg=37;
il=38;
ib=39;
ispreadd=40;
ispreadg=41;
idiffz=42;
imd=43;
img=44;
ikappa=45;
iweightd=46;
iweightg=47;

z_vec = 1:41;
IRFs_sh_amb(:,iag)=0*z_vec';

%%%%%%%%AMBIGUITY COMPARISON%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%AMBIGUITY COMPARISON%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(5,3,1)
plot(time, IRFs_sh_amb(:,iy),'b', time, IRFs_sh_amb_weightd(:,iy),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Total Output, $Y$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,2)
plot(time, IRFs_sh_amb(:,iyg),'b', time, IRFs_sh_amb_weightd(:,iyg),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Green Output, $Y^g$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,3)
plot(time, IRFs_sh_amb(:,iyd),'b', time, IRFs_sh_amb_weightd(:,iyd),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Dirty Output, $Y^d$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,4)
plot(time, IRFs_sh_amb(:,ikg),'b', time, IRFs_sh_amb_weightd(:,ikg),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Green Capital, $K^g$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,5)
plot(time, IRFs_sh_amb(:,ikd),'b', time, IRFs_sh_amb_weightd(:,ikd),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Dirty Capital, $K^d$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,6)
plot(time, IRFs_sh_amb(:,ic),'b', time, IRFs_sh_amb_weightd(:,ic),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Consumption, $C$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,7)
plot(time, IRFs_sh_amb(:,ihd),'b', time, IRFs_sh_amb_weightd(:,ihd),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Dirty Hours, $H^d$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,8)
plot(time, IRFs_sh_amb(:,ihg),'b', time, IRFs_sh_amb_weightd(:,ihg),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Green Hours, $H^g$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,9)
plot(time, IRFs_sh_amb(:,ib),'b', time, IRFs_sh_amb_weightd(:,ib),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Bonds, $B$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,10)
plot(time, IRFs_sh_amb(:,il),'b', time, IRFs_sh_amb_weightd(:,il),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Loans, $L$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,11)
plot(time, IRFs_sh_amb(:,ild),'b', time, IRFs_sh_amb_weightd(:,ild),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Dirty Loans, $L^d$', 'interpreter', 'latex','FontSize',10);
% plot(time, IRFs_sh_amb(:,iwd),'b', time, IRFs_sh_amb_weightd(:,iwd),'m--','LineWidth', 1);
% ylabel('\%','interpreter', 'latex','FontSize',10)
% title('Dirty Wages, $W^d$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,12)
plot(time, IRFs_sh_amb(:,ilg),'b', time, IRFs_sh_amb_weightd(:,ilg),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Green Loans, $L^g$', 'interpreter', 'latex','FontSize',10);
% plot(time, IRFs_sh_amb(:,iwg),'b', time, IRFs_sh_amb_weightd(:,iwg),'m--','LineWidth', 2);
% ylabel('\%','interpreter', 'latex','FontSize',10)
% title('Green Wages, $W^g$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,13)
plot(time, IRFs_sh_amb(:,ie),'b', time, IRFs_sh_amb_weightd(:,ie),'m--','LineWidth', 1);
ylabel('\%','interpreter', 'latex','FontSize',10)
title('Emissions, $E$', 'interpreter', 'latex','FontSize',10);
subplot(5,3,14)
plot(time, IRFs_sh_amb(:,iamb),'b', time, IRFs_sh_amb_weightd(:,iamb),'m--','LineWidth', 1);
ylabel('p.p.','interpreter', 'latex','FontSize',10)
title('Ambiguity, $a$', 'interpreter', 'latex','FontSize',10);
%legend('Undistorted','Worst case', 'Orientation', 'Vertical', 'Location', 'North','interpreter', 'latex','FontSize',10);
subplot(5,3,15)
plot(time, IRFs_sh_amb(:,iweightd),'b', time, IRFs_sh_amb_weightd(:,iweightd),'m--','LineWidth', 1);
ylabel('p.p.','interpreter', 'latex','FontSize',10)
title('Weight on Dirty Loans, $w^d$', 'interpreter', 'latex','FontSize',10);
hl=legend('Baseline','Rule on $w^d_t$', 'Orientation', 'Horizontal', 'Location', 'Best','interpreter', 'latex','FontSize',10);
newPosition = [0.52 0.05 0 0];
newUnits = 'normalized';
set(hl,'Position', newPosition,'Units', newUnits);
% print('Graphs\IRF_fin_weightd_paper', '-depsc')
