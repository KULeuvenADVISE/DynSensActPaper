clc; clear all;

db = 'SINS';
scale = 1; %Lifetime/Energy [0/1]
chains{1} = {'framelogMel__none','CNN_64_128_FC64_DO0_2'};
chains{2} = {'MFCCddd__mustd','DNN_64_DO0_2'};
% dirs
dirs.bbpath = 'WASN-EM';
dirs.configs = 'configs';
dirs.chains = 'chains';
addpath(genpath(dirs.bbpath));
dirs.gen = fullfile(dirs.configs);
dirs.sens = fullfile(dirs.configs);
dirs.proc = fullfile(dirs.chains);
dirs.comm = fullfile(dirs.configs);

% overall setting
gen.method = ['general_' db]; 
sens.method = 'sensing';
proc{1}.method = 'FE/ADC';
comm.method = 'wireless'; 

%% Communication unit
figure;
comm.N_R_list = [8 0];
for c2=1:length(comm.N_R_list)
    comm.N_T_list = round((2.^(1:25)));
    ow.ignore_energy = [1 1 0];
    for c=1:length(comm.N_T_list)
        ow.N_T = comm.N_T_list(c);
        ow.N_R = comm.N_R_list(c2);
        fc = wasn_wrapper(dirs,gen,sens,proc,comm,ow);
        LT1(c,1) = fc.LT;
        E(c,1) = fc.E_sep(1,3);
    end
    comm.N_T_list = comm.N_T_list; %model based on data of 10s, normalize to 10s
    if scale==0
        plot(comm.N_T_list,LT1,'Color',[0 0 0],'LineWidth',1); hold on; % LT line
    else
        plot(comm.N_T_list,E,'Color',[0 0 0],'LineWidth',1); hold on; % LT line
    end
end

%% Sensing
gen.method = ['general_' db]; 
sens.method = 'sensing';
proc{1}.method = 'FE/ADC';
comm.method = 'wireless'; 
ow.ignore_energy = [0 1 1];
fc = wasn_wrapper(dirs,gen,sens,proc,comm,ow);
if scale==0
    plot([comm.N_T_list(1) comm.N_T_list(end)],[fc.LT fc.LT],'-.','Color',[0 0 0],'LineWidth',1);
else
    plot([comm.N_T_list(1) comm.N_T_list(end)],[fc.E fc.E],'-.','Color',[0 0 0],'LineWidth',1);
end
    
%% Processing unit
ls = {'--','-.'};
for k=1:length(chains)
    % overall setting
    proc{1}.method = ['FE/' chains{k}{1}];
    proc{2}.method = ['NN/' chains{k}{2}];

    %% Bits versus life time
    % primary energy consumption line
    ow = struct;
    ow.ignore_energy = [1 0 1];
    fc = wasn_wrapper(dirs,gen,sens,proc,comm,ow);
    hold on;  
    if scale==0
        plot([comm.N_T_list(1) comm.N_T_list(end)],[fc.LT fc.LT],'--','Color',[0 0 0],'LineWidth',1);
    else
        plot([comm.N_T_list(1) comm.N_T_list(end)],[fc.E fc.E],'--','Color',[0 0 0],'LineWidth',1);
    end
end

%% Draw text/line
if scale==0
    set(gca, 'XScale', 'log'); hold on;
    max_b = 10^4;
    audio_bits = 12*16000*15; % do audio bits
    classf_bits = 10 * 16; % do classf bits
    plot([audio_bits audio_bits],[max_b min(LT1)],':','Color',[0 0 0]+0.5,'LineWidth',1);
    plot([classf_bits classf_bits],[max_b min(LT1)],':','Color',[0 0 0]+0.5,'LineWidth',1);
    xlim([0 max(comm.N_T_list)]);
    %ylim([0 max_b]);
    xlabel('Amount of transmitted information N_T (byte/estimate)');
    ylabel('Lifetime (days)');
else
    set(gca, 'YScale', 'log'); hold on;
    set(gca, 'XScale', 'log'); hold on;
    max_b = 10^4;
    audio_bits = 12*16000*15; % do audio bits
    classf_bits = 10 * 16; % do classf bits
    plot([audio_bits audio_bits],[max_b min(E)],':','Color',[0 0 0]+0.5,'LineWidth',1);
    plot([classf_bits classf_bits],[max_b min(E)],':','Color',[0 0 0]+0.5,'LineWidth',1);
    xlim([min(comm.N_T_list) max(comm.N_T_list)]);
    ylim([10^-1 max_b]);
    xlabel('Amount of transmitted information N_T (byte/estimate)');
    ylabel('Energy consumption (mJ)');
    txt = 'Classification';
    text(5,max_b/2,txt,'BackgroundColor',[0.9 0.9 0.9],'FontSize',12)
    txt = 'Feature extraction';
    text(4000,max_b/2,txt,'BackgroundColor',[0.9 0.9 0.9],'FontSize',12)
    txt = 'Raw audio';
    text(3300000,max_b/2,txt,'BackgroundColor',[0.9 0.9 0.9],'FontSize',12)
    % Draw line spec
    txt = '$\mathcal{E}_T, \mathcal{E}_R=1B$';
    text(3,0.35,txt,'interpreter','latex', 'FontSize', 12);
    txt = '$\mathcal{E}_T, \mathcal{E}_R=0B$';
    text(3,0.18,txt,'interpreter','latex', 'FontSize', 12);
    txt = '$\mathcal{E}_{P,high}$';
    text(3,140,txt,'interpreter','latex', 'FontSize', 12);
    txt = '$\mathcal{E}_{P,low}$';
    text(3,50,txt,'interpreter','latex', 'FontSize', 12);
    txt = '$\mathcal{E}_S$';
    text(3,1,txt,'interpreter','latex', 'FontSize', 12);
end
print -depsc2 -r864 -painters test.eps