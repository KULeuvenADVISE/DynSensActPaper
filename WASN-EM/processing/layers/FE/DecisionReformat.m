% ----------------------------------------------------------------------
%  Decision format behaviour model
%
%  Not described in document. Refers to methods to reduce data
%  representation of the decision output. E.g. label, ranked, posterior,
%  softVQ, ...
%  Used in:
%  G. Dekkers et al, “Sensor fusion and selection in acoustic sensor network for
%  classification of daily activities"
%
%  Author: Gert Dekkers, KU Leuven
% ----------------------------------------------------------------------
% Syntaxis: [output_shape, complexity, nr_parameters] = DecisionReformat(pp,gp,input_shape)
% Inputs:
% (1) pp                output shape of the layer given an input shape and parameters
% (2) gp                complexity of the layer given an input shape and parameters
% (3) input_shape       number of parameters of the layer given an input shape and parameters
% Outputs:
% (1) output_shape      output shape of the layer given an input shape and parameters
% (2) complexity        complexity of the layer given an input shape and parameters
% (3) nr_parameters     number of parameters of the layer given an input shape and parameters
%
% Usage example (chain):
%   - class_name: Aggregation
%     config:
%       type: [mean,std] % use mean and std as output / other options could be [mean], [std]

function [output_shape, complexity, nr_parameters] = DecisionReformat(pp,gp,input_shape)
    % var inits
    output_shape = zeros(1,gp.nr_dimensions);
    complexity = zeros(1,gp.nr_arop);
    nr_parameters = zeros(1,1);
    
    % Choose DecisionReformat
    % Current assumption: no storage/computational cost
    switch pp.type
        case 'posterior'
            output_shape = input_shape*pp.precision/gp.S;
        case 'posteriorVQ'
            output_shape([gp.chid gp.featid gp.frameid]) = [pp.VQ_bits/gp.S 1 1];
        case 'ranked'
            rankbits = nextpow2(factorial(gp.CLASS_COUNT)/(factorial(pp.rank_size)*factorial(gp.CLASS_COUNT-pp.rank_size))*factorial(pp.rank_size));
            output_shape([gp.chid gp.featid gp.frameid]) = [rankbits/gp.S 1 1];
        case 'rankpost'
            rankbits = nextpow2(factorial(gp.CLASS_COUNT)/(factorial(pp.rank_size)*factorial(gp.CLASS_COUNT-pp.rank_size))*factorial(pp.rank_size));
            posteriorbits = pp.precision*pp.rank_size;
            output_shape([gp.chid gp.featid gp.frameid]) = [(posteriorbits+rankbits)/gp.S 1 1];
        case 'rankpostVQ'   
            rankbits = nextpow2(factorial(gp.CLASS_COUNT)/(factorial(pp.rank_size)*factorial(gp.CLASS_COUNT-pp.rank_size))*factorial(pp.rank_size));
            output_shape([gp.chid gp.featid gp.frameid]) = [(pp.VQ_bits+rankbits)/gp.S 1 1];
    end
end

