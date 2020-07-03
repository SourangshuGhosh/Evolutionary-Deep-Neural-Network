function output = Configuration()
%Configuration - Return or create a configuration
%
% Syntax: output = Configuration()
%
% Long description
output.name = 'tataprob'; %configuration name
save_configuration = true; 
%%EvoNN Training Configs=============================
Evotrain.subsets = 1; Evotrain.overlap = 1;      %number of partitions of datafile and overlap b/w them
Evotrain.generations = 100;       % 10 max generations for evolution
Evotrain.nonodes = 5;           %maximum number of nodes
Evotrain.Prey_popsize = 500;      %500 Initial popsize
Evotrain.no_Prey_preferred = 300; %500 Desired popsize
Evotrain.no_new_Prey = 200;       %500 new prey introduced every KillInterval
Evotrain.Predator_popsize = 50;  %Number of Predators
Evotrain.no_x = 50;              %lattice size (no of rows)
Evotrain.no_y = 50;              %lattice size (no of cols)
Evotrain.ploton = 50;            %set 0 for no plots or 1 for plots at every generation
%======================================================
output.EvoTrain = Evotrain;
%% EVoDN2 ============================================
%EvDtrain.in_index = [1:30];
%EvDtrain.out_index = [31:32];
EvDtrain.generations = 50; % 10 max generations for evolution
EvDtrain.Pop_str{1}{1} = [1:6];           %maximum number of nodes
EvDtrain.Pop_str{1}{2} = [6 9];
EvDtrain.Pop_str{2}{1} = [4:9];           %maximum number of nodes
EvDtrain.Pop_str{2}{2} = [6 8];
EvDtrain.Pop_str{3}{1} = [7:12];           %maximum number of nodes
EvDtrain.Pop_str{3}{2} = [6 8];
%EvDtrain.Pop_str{4}{1}= [2:5];
%EvDtrain.Pop_str{4}{2}= [4 9];
% EvDtrain.Pop_str{5}{1}= [9:12];
% EvDtrain.Pop_str{5}{2}= [4 3]
EvDtrain.Prey_popsize = 80;      %500 Initial popsize
EvDtrain.no_Prey_preferred = 100; %500 Desired popsize
EvDtrain.no_new_Prey = 50;       %500 new prey introduced every KillInterval
EvDtrain.Predator_popsize = 40;  %Number of Predators
EvDtrain.no_x = 20;              %lattice size (no of rows)
EvDtrain.no_y = 20;              %lattice size (no of cols)
EvDtrain.ploton = 50;%set 0 for no plots or 1 for plots at every generation
%%==========================================================================================
output.EvDtrain=EvDtrain;
%%BioGP Training Config===============================
Biotrain.evo_type = 2; %set 1 for only Biobj evolution and 2 for first single obj followed by Biobj evolution
Biotrain.max_depth = 6;%max depth to which a tree grows
Biotrain.max_roots = 6;%max subtrees that a tree grows
Biotrain.subsets = 1;
Biotrain.overlap = 1;
Biotrain.generation1 = 10;           %generations for single obj evolution set Biotrain.evo_type = 2
Biotrain.generations = 15;          %max generations for evolution (initial 20)
Biotrain.maxrank = 30;               %maxrank retained at KillInterval
Biotrain.KillInterval = 10;          %Interval at which bad preys are eliminated
Biotrain.Prey_popsize = 500;         %Initial popsize
Biotrain.no_Prey_preferred = 300;    %Desired popsize
Biotrain.no_new_Prey = 200;          %new prey introduced every KillInterval
Biotrain.Predator_popsize = 60;     %Number of Predators
Biotrain.no_x = 60;                  %lattice size (no of rows)
Biotrain.no_y = 60;                  %lattice size (no of cols)
Biotrain.tour_size = 5;              %tournament size for single objective GP
Biotrain.ploton = 1;                 %set 0 for no plots or 1 for plots at every generation
Biotrain.lambda = 0.5;               %parameter to evaluate complexity of tree = lambda*depth + (1-lambda)*nodes
Biotrain.err_red_lim = 1e-3;         %any subtree contibuting less than this value for error reduction is eliminated
%====================================================
output.Biotrain = Biotrain;


%%EvoNN Optimization Configs=========================
EvoOpt.obj(1) = 1 ;  %set 1 for min and -1 for max
EvoOpt.obj(2) = -1 ;  %set 1 for min and -1 for max
EvoOpt.Prey_popsize = 500;         %Initial popsize
EvoOpt.no_Prey_preferred = 500;    %Desired popsize
EvoOpt.no_new_Prey = 200;          %new prey introduced every KillInterval
EvoOpt.Predator_popsize = 50;     %Number of Predators 100
EvoOpt.no_generations = 10;       %max generations
EvoOpt.P_move_prey = 0.5;          %Prob with which a Prey moves
EvoOpt.P_mut = 0.5;                %prob of choosing a prey for mutation
%Prob of Xover is 1 for every Prey
EvoOpt.F_bad = 1e5;                %fitness assigned to preys performing badly
%2D-lattice
EvoOpt.no_x = 50;                  %lattice size (no of rows) 50
EvoOpt.no_y = 50;                  %lattice size (no of cols) 50
EvoOpt.KillInterval = 4;          %Interval at which bad preys are eliminated
EvoOpt.maxrank = 20;               %maxrank retained at KillInterval
EvoOpt.ploton = 1;                 %set 0 for no plots or 1 for plots at every generation
%constraints
EvoOpt.useConstraints = true; 
EvoOpt.LB_F = [18.27 1.89];              %set upper bound for F1 & F2 respectively
EvoOpt.UB_F = [29.05 12.88];                   %set lower bound for F1 & F2 respectively
%=================================================================
output.EvoOpt = EvoOpt;

%%BioGP Optimization Configs=========================
BioOpt.obj(1) = 1 ;  %set 1 for min and -1 for max
BioOpt.obj(2) = -1 ;  %set 1 for min and -1 for max
BioOpt.Prey_popsize = 500;         %Initial popsize
BioOpt.no_Prey_preferred = 500;    %Desired popsize
BioOpt.no_new_Prey = 200;          %new prey introduced every KillInterval
BioOpt.Predator_popsize = 50;     %Number of Predators 100
BioOpt.no_generations = 10;       %max generations
BioOpt.P_move_prey = 0.5;          %Prob with which a Prey moves
BioOpt.P_mut = 0.5;                %prob of choosing a prey for mutation
%Prob of Xover is 1 for every Prey
BioOpt.F_bad = 1e5;                %fitness assigned to preys performing badly
%2D-lattice
BioOpt.no_x = 50;                  %lattice size (no of rows) 50
BioOpt.no_y = 50;                  %lattice size (no of cols) 50
BioOpt.KillInterval = 4;          %Interval at which bad preys are eliminated
BioOpt.maxrank = 20;               %maxrank retained at KillInterval
BioOpt.ploton = 1;                 %set 0 for no plots or 1 for plots at every generation
BioOpt.useConstraints = true; 
BioOpt.LB_F = [18.27 1.89];              %set upper bound for F1 & F2 respectively
BioOpt.UB_F = [29.05 12.88];                   %set lower bound for F1 & F2 respectively
%=================================================================
output.BioOpt = BioOpt;

%%RVEA Optimization Configs===========================
RVEAopt.Generations = 120;
RVEAopt.p1p2 = num2cell([16 0]); %%[p1 p2] define the number of reference vectors. p1 is the number of divisions along an axis
RVEAopt.N = 10;  %%defines the population size.
RVEAopt.alpha = 2; % the parameter in APD, the bigger, the faster RVEA converges
RVEAopt.fr = 0.05; % frequency to call reference vector
%=====================================================
output.RVEAopt = RVEAopt;

%%cRVEA Optimization Configs===========================
cRVEAopt.obj = [1 1] ;  %set 1 for min and -1 for max
cRVEAopt.Generations = 100;
cRVEAopt.p1p2 = num2cell([20 0]); %%[p1 p2] define the number of reference vectors. p1 is the number of divisions along an axis
cRVEAopt.N = 20;  %%defines the population size.
cRVEAopt.alpha = 0.05; % the parameter in APD, the bigger, the faster cRVEA converges
cRVEAopt.fr = 0.1; % frequency to call reference vector

cRVEAopt.eqCon{1} = ''; %equality constraints(f(Var,Obj)=0)
cRVEAopt.ieqCon{1} = '';
%cRVEAopt.ieqCon{1} = '1.57-obj1';
cRVEAopt.ieqCon{1} = 'obj1'; %inequality contraints(f(Var,Obj)>0)
%cRVEAopt.ieqCon{3} = '68.85-obj2';
cRVEAopt.ieqCon{2} = 'obj2';
%cRVEAopt.ieqCon{5} = '2.49-obj3';
%cRVEAopt.ieqCon{6} = 'obj3-.70';
%cRVEAopt.ieqCon{7} = '2270-obj4';
%cRVEAopt.ieqCon{8} = 'obj4-2210';
%cRVEAopt.ieqCon{8} = 'obj8';

%cRVEAopt.ieqCon{6} = 'obj3-130';
% %=====================================================
output.cRVEAopt = cRVEAopt;


if save_configuration
save([output.name '.mat'], 'output');
end
end