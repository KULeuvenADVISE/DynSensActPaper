function fc = wasn_wrapper(dirs,gen,sens,proc,comm,ow)
%% Code for acquiring energy consumption for a single estimate (network node)
verbose = 0; %0: show final output, 1: show all

% -- Init system -- %
gen = general_loadparam(gen.method,dirs.gen); % get general parameters (fs, µC, ..)

% -- Sensing layer -- %
if verbose, disp('%%%% Sensing %%%%'); end;
sens.conf = sensing_loadparam(sens.method, dirs.sens, gen); % get sensing parameters
[sens.cons, sens.os] = sensing_consumption(sens.conf, gen); % [mJ]

% -- Processing layer -- %
input_shape = sens.os; prev_chain = []; E_proc = []; name_proc = cell(0); %init as empty at first run
for p=1:length(proc) % for each processing chain
    if verbose, disp(['%%%% Processing layer: ' proc{p}.method ' %%%%']); end;
    % get values
    proc{p}.conf = proc_loadparam(proc{p}.method, dirs.proc, input_shape, gen); % get params
    [proc{p}.os, proc{p}.ops, proc{p}.par] = proc_info(proc{p}.conf, gen, input_shape,verbose); % get output shape, ops and nr. parameters
    [proc{p}.ma, proc{p}.ms_o, proc{p}.ms_p] = memo_acc_stor(proc{p}.os, proc{p}.par, proc{p}.ops, proc{p}.conf, prev_chain, gen,verbose); % get memory used in storing ops/params
    [proc{p}.cons.all,proc{p}.cons.op,proc{p}.cons.ms_o,proc{p}.cons.ms_p,proc{p}.cons.ma] = bits_to_energy(proc{p}.ma, proc{p}.ms_o, proc{p}.ms_p, proc{p}.ops, gen); % energy consumption
    % keep from previous chain
    input_shape = proc{p}.os(end,:);
    prev_chain = [prev_chain; proc{p}.conf];
end

% -- Communication layer -- %
if verbose, disp('%%%% Communication %%%%'); end;
comm.conf = comm_loadparam(comm.method, dirs.comm,gen); % get comm params
if ~isfield(ow,'N_T')
    comm.conf.N_T = prod(input_shape)*gen.S; % informative bits to communicate
else
    comm.conf.N_T = ow.N_T;
end
if ~isfield(ow,'N_R')
    comm.conf.N_R = 0; % informative bits to communicate
else
    comm.conf.N_R = ow.N_R;
end
comm.cons = comm_consumption(comm.conf,sens.conf,gen);
comm.cons = sum(comm.cons);

% -- Final output -- %
if isfield(ow,'ignore_energy')
    if ow.ignore_energy(1)
    	sens.cons = 0;
    end
    if ow.ignore_energy(2)
        for p=1:length(proc)
            proc{p}.cons.all = 0;
        end
    end
    if ow.ignore_energy(3)
        comm.cons = 0;
    end
end
fc = prep_chain_info(sens,proc,comm,gen); % get everything in proper format to plot/print
fc.NT = comm.conf.N_T;
end