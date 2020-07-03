%function autorun
warning off
%myFun - Description
%
% Syntax autorun
%
% This will be replaced when code added to mastercode
	Problem_name = 'ZDT1_100';
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	EvDtrain.in_index = [1:30];
	EvDtrain.out_index = [31:32];
	EvDtrain.generations = 20; % 10 max generations for evolution
	EvDtrain.Pop_str{1}{1} = [1:15];           %maximum number of nodes
	EvDtrain.Pop_str{1}{2} = [15 15];
    EvDtrain.Pop_str{2}{1} = [15:30];           %maximum number of nodes
	EvDtrain.Pop_str{2}{2} = [16 15];
	EvDtrain.Prey_popsize = 80;      %500 Initial popsize
	EvDtrain.no_Prey_preferred = 100; %500 Desired popsize
	EvDtrain.no_new_Prey = 50;       %500 new prey introduced every KillInterval
	EvDtrain.Predator_popsize = 40;  %Number of Predators
	EvDtrain.no_x = 20;              %lattice size (no of rows)
	EvDtrain.no_y = 20;              %lattice size (no of cols)
	EvDtrain.ploton = 1;            %set 0 for no plots or 1 for plots at every generation
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	parameters.EvDtrain = EvDtrain;
    parameters.EvDtrain.name = Problem_name;
    
	Train(Problem_name,parameters);
%end
