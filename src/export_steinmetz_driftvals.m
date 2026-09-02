%%
% Modified code originally from https://github.com/nsteinme/steinmetz-et-al-2021
% to extract the slope of the linear fits to unit count
%%

clear all
close all
clc
addpath(genpath('/home/mmelin/steinmetz-et-al-2021'))
%data directory, MUST extract zip first
processed_folder = '/home/mmelin/steinmetz-et-al-2021/fig2/data';

%% linear fit


tp = list_files(processed_folder,'*.mat');
mm = 1;
ax = [];
lab_list = {'haesler','hantman','lee','moser','o_keefe','carandini'};
clabs = lines(6);
time_start = 3;
ev_corr = [];
yield_corr = [];
ev_slope = [];
yield_slope = [];
c_sig = [];
for ii = 1:length(tp)
    load(tp{ii});
    [ndays,nprobes] = size(dat);
    
    for ip = 1:nprobes
        dd = vertcat(dat(:,ip).date);
        
        tmpx = dd-min(dd)+time_start ;
        tmpy = vertcat(dat(:,ip).mev);
        idx = ~isnan(tmpx)&~isnan(tmpy)&tmpx<=60;
        
        
        log_tmpy = log10(tmpy);
        lab_id = find(strcmp(dat(1,ip).lab,lab_list));

        p1 = polyfit(tmpx(idx),log_tmpy(idx),1);
        ev_slope = [ev_slope,p1(1)];
        c_sig = [c_sig;clabs(lab_id,:)];
        

        [R,P]=corrcoef(tmpx(idx),log_tmpy(idx));
        ev_corr = [ev_corr,P(2,1)];
       
        
        dd = vertcat(dat(:,ip).date);
        tmpx = dd-min(dd)+time_start;
        tmpy = vertcat(dat(:,ip).mgood);
        idx = ~isnan(tmpx)&~isnan(tmpy);
               
        log_tmpy = log10(tmpy);
        dat(:,ip).name;
        lab_id = find(strcmp(dat(1,ip).lab,lab_list));

        
        p1 = polyfit(tmpx(idx),log_tmpy(idx),1);
        yield_slope = [yield_slope,p1(1)];
        
        [R,P]=corrcoef(tmpx(idx),log_tmpy(idx));
        yield_corr = [yield_corr,P(2,1)];
        
        
    end
    
end

save("steinmetz_slope.mat",'yield_slope')

