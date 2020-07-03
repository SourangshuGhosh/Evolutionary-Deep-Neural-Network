function Train (Problem_name, parameters)
%PP_GP_SUBSETS for training GP models on dataset
%Created by Brijesh Kumar Giri
%   2011-2012
%In the plot F1 (x-axis) shows the complexity and
%F2 (y-axis) shows the training error
%====Instruction on creating datafile for the code=========================
% if a .xls file for dataset is available then go to the respective directory
% and double click it. Matlab will load data from the .xls file in a matrix of
% type double in a variable named data. the headers in .xls file will loaded
% matrix of type cell in a variable named textdata. Make sure that the xls
% contains two rows of headers, 1st row gives a detailed name of variable and
% 2nd row gives names to be used as terminal nodes (keep it as simple as
% possible, best would be to give an alphabet followed by a number)
% next execute the following line in commamnd window:
% data = [textdata; num2cell(data)];
% now save variable data in a .mat file by executing the line below:
% save Nmae_of _file data
%try this with example demo.xls
%also check out data_demo.amt and other demo data files

clc
%=====================DND==================================================
global Setslog figure_handle

figure_handle = [];
RandStream('mt19937ar','seed', sum(100*clock));
Setslog = [];
%=====================Define parameters=======================================
filename = [Problem_name '.xls'];         %Data file
in_index = parameters.in_index;           %independent variables column no.
out_index = parameters.out_index;         %dependent variable column no.
savedir = fullfile(pwd,'Output',Problem_name,'BioGP',parameters.name);
mkdir(savedir);

Setslog.evo_type = parameters.Biotrain.evo_type;  %set 1 for only Biobj evolution and 2 for first single obj followed by Biobj evolution

Setslog.in_index = in_index;            %input variable cols in datafile
%Setslog.out_index = out_index;          %output variable col in datafile

[DataSet,paraname,DATA] = xlsread(filename);
Setslog.DATA = DATA;
Setslog.paraname= paraname;
Setslog.DataSet = DataSet;

F1_set = {};                            %single arity functions
F2_set = {'add' 'sub' 'mul' 'rdiv'};    %double arity functions
F3_set = {};                            %triple arity functions
%for higher arity functions define F4_set, F5_set and so on till the arity
%desired and include them in F.set sequentially. for ex if 5 arity function
%is desired define F4_set ={}; F5_set = {'your_func'}; and then modify
%F.set = {F1_set F2_set F3_set F4_set F5_set};
F.set = {F1_set F2_set F3_set};
j = length(F.set); temp = {}; F.sets = [];
for i = 1:j
    F.sets = [F.sets i*ones(1,length(F.set{i}))];
    temp = [temp F.set{i}];
end
F.set = temp;
Setslog.F = F;

T1_set = paraname(2,in_index); %keep terminals name in datafile as simple as possible see run_BioGP and demo datafiles
T2_set = {};                  %ERC empheral random floating point constant and ADF's with 0 arity
T.set = [T1_set T2_set];
Setslog.T = T;

Setslog.max_depth = parameters.Biotrain.max_depth;              %max depth to which a tree grows
Setslog.max_roots = parameters.Biotrain.max_roots;              %max subtrees that a tree grows

subsets = parameters.Biotrain.subsets; overlap = parameters.Biotrain.overlap;           % previous 2 and 20 number of partitions of datafile and overlap b/w them
set_size = length(DataSet(:,1)); Setslog.set_size = set_size;
Setslog.subsets = subsets;
Setslog.overlap = overlap;
subset_size = round((set_size + (subsets-1)*overlap) / subsets);
Setslog.subset_size = subset_size;

Setslog.generation1 = parameters.Biotrain.generation1;           %generations for single obj evolution set Setslog.evo_type = 2
Setslog.generations = parameters.Biotrain.generations;          %max generations for evolution (initial 20)
Setslog.maxrank = parameters.Biotrain.maxrank;               %maxrank retained at KillInterval
Setslog.KillInterval = parameters.Biotrain.KillInterval;          %Interval at which bad preys are eliminated
Setslog.Prey_popsize = parameters.Biotrain.Prey_popsize;         %Initial popsize
Setslog.no_Prey_preferred = parameters.Biotrain.no_Prey_preferred;    %Desired popsize
Setslog.no_new_Prey = parameters.Biotrain.no_new_Prey;          %new prey introduced every KillInterval
Setslog.Predator_popsize = parameters.Biotrain.Predator_popsize;     %Number of Predators
Setslog.no_x = parameters.Biotrain.no_x;                  %lattice size (no of rows)
Setslog.no_y = parameters.Biotrain.no_y;                  %lattice size (no of cols)
Setslog.tour_size = parameters.Biotrain.tour_size;              %tournament size for single objective GP
Setslog.ploton = parameters.Biotrain.ploton;                 %set 0 for no plots or 1 for plots at every generation
Setslog.lambda = parameters.Biotrain.lambda;               %parameter to evaluate complexity of tree = lambda*depth + (1-lambda)*nodes
Setslog.err_red_lim = parameters.Biotrain.err_red_lim;         %any subtree contibuting less than this value for error reduction is eliminated
%================================DND==========================================
for out = out_index
Setslog.out_index = out;
if Setslog.ploton
    figure_handle = [figure(1) figure(2) figure(3) figure(4)];
    scrsz = get(0,'ScreenSize');
    set(figure_handle(1), 'OuterPosition', [0*scrsz(3) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(1))
    set(figure_handle(2), 'OuterPosition', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(2))
    set(figure_handle(3), 'OuterPosition', [0*scrsz(3) 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(3))
    set(figure_handle(4), 'OuterPosition', [scrsz(3)/2 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(4))
    pause(1)
end
warning('off', 'all')
eval(['delete ' savedir '\Y' num2str(out-out_index(1)+1) '.mat']);
warning('on', 'all')
%============================================

%choose datasets for training
if subsets == 1
    training_subsets = 1;
    no_run = 1; Setslog.no_run = no_run;
else
    training_subsets = [eye(subsets,subsets); ones(1,subsets)];
    no_run = 1+subsets; Setslog.no_run = no_run;
end
Setslog.chosen_training_sets = training_subsets;
%============================================

%creating training data
data_index = (1:1:length(DataSet))';
for i = 1:no_run
    chosen_data_index = [];
    from = 1; to = from + subset_size - 1;
    for j = 1:subsets-1
        if training_subsets(i,j) == 1
            chosen_data_index = [chosen_data_index; data_index(from:to,:)];
        end
        from = to - overlap + 1; to = from + subset_size - 1;
    end
    
    if training_subsets(i,subsets) == 1
        chosen_data_index = [chosen_data_index; data_index(from:length(data_index(:,1)),:)];
    end
    chosen_data_index = unique(chosen_data_index);
    Setslog.dataset(i).data_index = chosen_data_index;
    Setslog.dataset(i).in = Setslog.DataSet(chosen_data_index,in_index);
    Setslog.dataset(i).out = Setslog.DataSet(chosen_data_index,Setslog.out_index);
end
%====================================

%train GP tree
for i = 1:no_run
    if i < no_run
        fprintf('\nTraining Dataset %d\n\n', i);
    else
        fprintf('\nTraining Whole Dataset\n\n');
    end
    PP_GP(i);
    pause(5)
end
fprintf('\n\nTrainining Subsets\n');
disp(training_subsets);
%====================================
eval(['save ' savedir '\Y' num2str(out-out_index(1)+1) '.mat Setslog'])
copyfile([pwd '\BioGP\evaluate_obj.m'], savedir);
close all;
autooutput(Setslog,savedir);
close all;
svr(Setslog,[savedir '\svr_Y' num2str(out-out_index(1)+1)],'trend.mat');
close all;
end
end