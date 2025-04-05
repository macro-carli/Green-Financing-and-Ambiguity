clear

load('LoanRateSpread.mat')
load('ChargeoffRatesonBusinessLoans.mat')
load('IndustrialProduction.mat')
load('CapacityUtilization.mat')
load('CILoans.mat')


LOAN_RATE = LoanRateSpread.Loan_Rate;
FEDFUNDS = LoanRateSpread.Fed_Funds_Rate;
SPREAD = LoanRateSpread.Rate_Spread;
CHARGE_OFF = ChargeoffRatesonBusinessLoans.Charge_Off_Rate;
DATES = ChargeoffRatesonBusinessLoans.Date;
INDPRO = IndustrialProduction.INDPRO;
CAP_UTI = CapacityUtilization.Cap_Uti;
LOANS = CILoans.Loans;

mean_loan_rate = mean(LOAN_RATE); %4.2

lCHARGE_OFF = log(CHARGE_OFF);
lINDPRO = log(INDPRO);
lCAP_UTI = log(CAP_UTI);
lLOANS = log(LOANS);

lINDPRO = lINDPRO(265:end-1);
lCAP_UTI = lCAP_UTI(73:end-1);
lLOANS = lLOANS(49:end-1);

lCHARGE_OFF = lCHARGE_OFF - mean(lCHARGE_OFF);
lINDPRO = lINDPRO - mean(lINDPRO);
lCAP_UTI = lCAP_UTI - mean(lCAP_UTI);
lLOANS = lLOANS - mean(lLOANS);

%plot(DATES, lINDPRO, DATES, lCHARGE_OFF)


[lCHARGE_OFF_trend, lCHARGE_OFF_cycle] = hp(lCHARGE_OFF, 1600);
[lINDPRO_trend, lINDPRO_cycle] = hp(lINDPRO, 1600);
[lCAP_UTI_trend, lCAP_UTI_cycle] = hp(CAP_UTI, 1600);
[lLOANS_trend, lLOANS_cycle] = hp(CAP_UTI, 1600);
lCHARGE_OFF_cycle_bp = bpass(lCHARGE_OFF, 6, 32);
lINDPRO_cycle_bp = bpass(lINDPRO, 6, 32);
lCAP_UTI_cycle_bp = bpass(lCAP_UTI, 6, 32);
lLOANS_cycle_bp = bpass(lLOANS, 6, 32);

std_indpro = std(lINDPRO_cycle);
std_chargeoff = std(lCHARGE_OFF_cycle);
std_loans = std(lLOANS_cycle);

%plot(DATES, lINDPRO_cycle, DATES, lCHARGE_OFF_cycle)

maxlags=6;
[autocorr_chargeoff,lags]=xcorr(lCHARGE_OFF_cycle, maxlags, 'coeff');

maxlags=6;
[corr_chargeoff,lags]=xcorr(lINDPRO_cycle_bp, lCHARGE_OFF_cycle_bp, maxlags, 'coeff');
stem(lags,corr_chargeoff)

% maxlags=6;
% [corr_caputi,lags]=xcorr(lCAP_UTI_cycle_bp, lCHARGE_OFF_cycle_bp, maxlags, 'coeff');
% stem(lags,corr_caputi)


