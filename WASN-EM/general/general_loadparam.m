function general_param = general_loadparam(general_str,param_dir,gen_param)
% ----------------------------------------------------------------------
%  Function to read the general config from the yaml file
%  Author: Gert Dekkers, KU Leuven
% ----------------------------------------------------------------------
% Syntaxis: general_param = general_loadparam(feat_str,gen_param)
% Inputs:
% (1) general_str       general config str str, refers to yaml file [str]
% (2) param_dir      	directory of string if where the params are located - optional [str]
% (3) gen_param      	gen_param, that could be used as constant value overload in the yaml file
% Outputs:
% (1) general_param     general config parameters [struct]      

if nargin~=3, gen_param = struct; end; % no parameter overloading

if isempty(param_dir) % if empty, use default
    param_dir = fullfile('general','params');
end

general_param = process_yaml(general_str,param_dir,gen_param);
