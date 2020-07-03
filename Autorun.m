%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This file currently allows creation of meta-models from datasets
%% using BioGP and EvoNN, and the optimization of the meta-models
%% using EvoNN, BioGP, or RVEA.
%%
%% To be added : Multiple configuration support.
%%             : Creation of NN models using standard methods
%%             : Creation of NN and GP models using RVEA
%%             : Inclusion of other standard modeling methods
%%             : GUI
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% The parameters of various algorithms involved can be changed
%% from the Configuration.m file. Default values for those
%% parameters are saved in Default.mat
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs are saved in the 'Output' folder. Temp folder stores
%% temporary outputs to be used by the program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Autorun()
clc
do_training = false;
do_optimization = true;

Training_Algorithms = {'EvoDN2'};  %{'EvoNN' 'BioGP'}
Optimization_Algorithms = {'cRVEA'}; %{'RVEA' 'cRVEA'}

Problems = {'test1'};
in_index = [1:12];  %in_index = [a:b ; c:d; e:f] a:b for Problem 1, c:d for problem 2, and so on;
out_index = [13:14];

use_defaults = false; % if true, Default.mat will be used, otherwise Configuration.m will be used.
multi_config = false;  %future work

if use_defaults
    parameters = importdata('Default.mat');
elseif ~multi_config
    parameters = Configuration();
%else
    %multiconfig support
end

parameters.in_index = in_index;
parameters.out_index = out_index;
if do_training
    for Prob = 1:length(Problems)
        for param = 1:length(parameters)
            for Algo = 1:length(Training_Algorithms)
                oldpath = path;
                addpath(genpath([pwd '\' Training_Algorithms{Algo}]));
                disp(parameters(param));
                Train(Problems{Prob},parameters(param));
                path(oldpath);
                pause(5);
                close all;
            end
        end
    end
end

if do_optimization
    for Prob = 1:length(Problems)
        for param = 1:length(parameters)
            for Algo = 1:length(Training_Algorithms)
                for opt = 1:length(Optimization_Algorithms)
                    oldpath = path;
                    addpath(genpath([pwd '\' Training_Algorithms{Algo}]));
                    addpath(genpath([pwd '\' Optimization_Algorithms{opt}]));
                    savedir = fullfile(pwd,'Output',Problems{Prob},Training_Algorithms{Algo},parameters(param).name);
                    addpath(savedir);
                    Opt(Problems{Prob},Training_Algorithms{Algo},parameters(param),savedir);
                    path(oldpath);
                    pause(5);
                    close all;
                end
            end
        end
    end
end
end
    