clc; clear all; close all;
%% stats for unitcounts
tbl = readtable('C:\Users\mmelin\Downloads\acute_chronic_unit_count.csv');

% Convert predictors to categorical where appropriate
tbl.is_chronic = categorical(tbl.is_chronic);
tbl.insertion_site = categorical(tbl.insertion_site);
tbl.session_num = categorical(tbl.session_num);

% fixed effect: is_chronic, random effects: insertion_site, timepoint, with nested session_num
sua_lme = fitlme(tbl, 'single_units ~ is_chronic + (1|timepoint) + (1|insertion_site) + (1|insertion_site:session_num)');
mua_lme = fitlme(tbl, 'multi_units ~ is_chronic + (1|timepoint) + (1|insertion_site) + (1|insertion_site:session_num)');

% Show results
disp(sua_lme);
disp(mua_lme);

%compare(lme1, lme2)

%% now run the dredge stats
tbl = readtable('C:\Users\mmelin\Downloads\acute_chronic_dredge.csv');

% Convert predictors to categorical where appropriate
tbl.is_chronic = categorical(tbl.is_chronic);
tbl.insertion_site = categorical(tbl.insertion_site);
tbl.session_num = categorical(tbl.session_num);

% fixed effect: is_chronic, random effects: insertion_site, timepoint, with nested session_num
lme = fitlme(tbl, 'session_drift ~ is_chronic + (1|timepoint) + (1|insertion_site) + (1|insertion_site:session_num)');

% Show results
disp(lme);
