%%%% DISPERSION DATA %%%%

%Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "RGDP";
opts.DataRange = "A11:G233";

% Specify column names and types
opts.VariableNames = ["Survey_DateT", "Var2", "Var3", "Var4", "RGDP_P25_D3T1", "RGDP_P75_D3T1", "RGDP_D3T1"];
opts.SelectedVariableNames = ["Survey_DateT", "RGDP_P25_D3T1", "RGDP_P75_D3T1", "RGDP_D3T1"];
opts.VariableTypes = ["string", "char", "char", "char", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, ["Survey_DateT", "Var2", "Var3", "Var4"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Survey_DateT", "Var2", "Var3", "Var4"], "EmptyFieldRule", "auto");

% Import the data
DispersionRGDP = readtable("C:\Users\marco\Desktop\Economics and Finance (PhD in)\Research2\Ambig_enivr_codes\main\empirics_disp\Dispersion_3.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts


%%
Date_times = DispersionRGDP.(1);
dates = datetime(Date_times,'InputFormat','yyyyQQQ');

lRGDP_P25 = DispersionRGDP.(2);
lRGDP_P75 = DispersionRGDP.(3);
ldisp = DispersionRGDP.(4);

%Convert to loglevels
l1RGDP_P25 = lRGDP_P25/100;
l1RGDP_P75 = lRGDP_P75/100;
l1disp = l1RGDP_P75 - l1RGDP_P25;

%Convert to levels (billions of dollars)
RGDP_P25 = exp(lRGDP_P25/100);
RGDP_P75 = exp(lRGDP_P75/100);
disp = RGDP_P75 - RGDP_P25;

%Restrict sample
dates_start = find(dates==datetime('01-Jan-1985'));
dates_end = find(dates==datetime('01-Oct-2019'));
dates_85 = dates(dates_start:dates_end);

lRGDP_P25_85 = lRGDP_P25(dates_start:dates_end);
lRGDP_P75_85 = lRGDP_P75(dates_start:dates_end);
ldisp_85 = ldisp(dates_start:dates_end);

l1RGDP_P25_85 = l1RGDP_P25(dates_start:dates_end);
l1RGDP_P75_85 = l1RGDP_P75(dates_start:dates_end);
l1disp_85 = l1disp(dates_start:dates_end);

RGDP_P25_85 = RGDP_P25(dates_start:dates_end);
RGDP_P75_85 = RGDP_P75(dates_start:dates_end);
disp_85 = disp(dates_start:dates_end);

[~,l1disp_85_cycle] = hpfilter(l1disp_85,1600);


mean_ldisp = mean(ldisp)
mean_l1disp = mean(l1disp)
mean_disp = mean(disp)

std_ldisp = std(ldisp)
std_l1disp = std(l1disp)
std_disp = std(disp)

mean_ldisp_85 = mean(ldisp_85)
mean_l1disp_85 = mean(l1disp_85)
mean_disp_85 = mean(disp_85)
mean_l1disp_85_cycle = mean(l1disp_85_cycle)

std_ldisp_85 = std(ldisp_85)
std_l1disp_85 = std(l1disp_85)
std_disp_85 = std(disp_85)
std_l1disp_85_cycle = std(l1disp_85_cycle)

%autocorr(l1disp_85)
figure
plot(dates_85, l1disp_85)
print('disp', '-depsc')
