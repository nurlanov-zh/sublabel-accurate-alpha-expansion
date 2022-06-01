% Reproducing the results from 
% Table 1: Comparison of various methods in terms of the accuracy and 
% efficiency for robust image denoising; 
% Table 2: Energies of discretized sublabel-accurate solutions;
% Figure 1; Figure 2.

run_previous_methods = false;

discrete_method = 'GCO';            % 'GCO' or 'exact'
sublabel_method = 'QL';             % 'QL' or 'LM' or 'QM' 
discretize_sublabel = true;

coef_imresize = 0.5;
lmb = 0.6;
%labels = [5	10	15 20	30	40	50	70	100 120 150 200 256];

labels = [5 10];
num_labels = length(labels);

trunc_smooth = -1;      % default: no truncation in smoothness term 
k_smooth = 1;           

%%% For non-convex priors experiment, see eq. (10):
%%% Truncated Linear smoothness cost
%{
trunc_smooth = 0.6;
k_smooth = 1;
lmb = 0.6;
%}

%%% Truncated Quadratic smoothness cost
%{
trunc_smooth = 0.7;
k_smooth = 2;
lmb = 3;
%}

timings_discrete = [];
energies_discrete = [];

timings_ours = [];
energies_ours = [];
energies_raw_ours = [];
non_convex_terms = [];

energies_ours_discretized = [];

for L = labels
    [time_elapsed_GCO, energy_to_compare_GCO, ...
    time_elapsed_sublabel, energy_to_compare_sublabel, ...
    energy_raw_sublabel, count_non_convex, energy_to_compare_dissublabel ...
    ] = discrete_plus_refine_denoising( L, coef_imresize, lmb, trunc_smooth, ...
    k_smooth, discrete_method, sublabel_method, discretize_sublabel);

    timings_discrete(end+1) = time_elapsed_GCO;
    energies_discrete(end+1) = energy_to_compare_GCO;
    
    timings_ours(end+1) = time_elapsed_sublabel;
    energies_ours(end+1) = energy_to_compare_sublabel;
    energies_raw_ours(end+1) = energy_raw_sublabel;
    non_convex_terms(end+1) = count_non_convex;
    
    energies_ours_discretized(end+1) = energy_to_compare_dissublabel;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Previous sublabel-accurate methods by Pock and Mollenhoff:
if run_previous_methods
    fprintf('\n___________________________________________________________\n');
    fprintf('Previous sublabel-accurate methods by Pock and Mollenhoff.\n');
    fprintf('To run them make sure you have valid CUDA device!\n')
    pause(3);

    timings_Pock = [];
    energies_Pock = [];
    energies_Pock_discretized = [];

    timings_Mollenhoff = [];
    energies_Mollenhoff = [];
    energies_Mollenhoff_discretized = [];

    for L = labels
        [time_Pock, energy_Pock, energy_Pock_discretized ...
            ] = Pock_denoising(L, coef_imresize, lmb, discretize_sublabel);

        [time_Moll, energy_Moll, energy_Moll_discretized ...
            ] = Mollenhoff_denoising(L, coef_imresize, lmb, discretize_sublabel);

        timings_Pock(end+1) = time_Pock;
        energies_Pock(end+1) = energy_Pock;
        energies_Pock_discretized(end+1) = energy_Pock_discretized;

        timings_Mollenhoff(end+1) = time_Moll;
        energies_Mollenhoff(end+1) = energy_Moll;
        energies_Mollenhoff_discretized(end+1) = energy_Moll_discretized;

    end
end




