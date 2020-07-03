function Train(Problem_name, parameters)
%PP_NNGA_SUBSETS for training Neural Net models on dataset
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
savedir = fullfile(pwd,'Output',Problem_name,'EvoNN',parameters.name);
mkdir(savedir);

Xmin = eps; Setslog.Xmin = Xmin;          %normalization range for variables
Xmax = 1; Setslog.Xmax = Xmax;

Setslog.in_index = in_index;

[DataSet,paraname,DATA] = xlsread(filename);
Setslog.DATA = DATA;
Setslog.paraname= paraname;
Setslog.DataSet = DataSet;

subsets = parameters.EvoTrain.subsets; overlap = parameters.EvoTrain.overlap;           %number of partitions of datafile and overlap b/w them
set_size = length(DataSet(:,1)); Setslog.set_size = set_size;
Setslog.subsets = subsets;
Setslog.overlap = overlap;
subset_size = round((set_size + (subsets-1)*overlap) / subsets);
Setslog.subset_size = subset_size;

Setslog.generations = parameters.EvoTrain.generations;   % 10 max generations for evolution
Setslog.nonodes = parameters.EvoTrain.nonodes;           %maximum number of nodes
Setslog.noinnodes = length(Setslog.in_index);
Setslog.nooutnodes = 1;
Setslog.maxrank = 20;                           %maxrank retained at KillInterval
Setslog.KillInterval = 5;                       %Interval at which bad preys are eliminated
Setslog.Prey_popsize = parameters.EvoTrain.Prey_popsize; %500 Initial popsize
Setslog.no_Prey_preferred = parameters.EvoTrain.no_Prey_preferred;   %500 Desired popsize
Setslog.no_new_Prey = parameters.EvoTrain.no_new_Prey;               % 500 new prey introduced every KillInterval
Setslog.Predator_popsize = parameters.EvoTrain.Predator_popsize;     %Number of Predators
Setslog.no_x = parameters.EvoTrain.no_x;                 %lattice size (no of rows)
Setslog.no_y = parameters.EvoTrain.no_y;                 %lattice size (no of cols)
Setslog.ploton = parameters.EvoTrain.ploton;             %set 0 for no plots or 1 for plots at every generation
%======================================================================
for out = out_index
Setslog.out_index = out;
if Setslog.ploton
    figure_handle = [figure(1) figure(2) figure(3) figure(4)];
    scrsz = get(0,'ScreenSize');
    set(figure_handle(1), 'OuterPosition', [0*scrsz(3) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(1))
    set(figure_handle(2), 'OuterPosition', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(2))
    set(figure_handle(3), 'OuterPosition', [0*scrsz(3) 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(3))
    set(figure_handle(4), 'OuterPosition', [scrsz(3)/2 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(4))
    pause(20)
    
end
    


%scale the data
DataSet_sc = [];
for i=1:length(DataSet(1,:))
    Data_min(1,i) = min(DataSet(:,i));
    Data_max(1,i) = max(DataSet(:,i));
    DataSet_sc = [DataSet_sc Xmin+ (DataSet(:,i)-Data_min(1,i))/(Data_max(1,i)-Data_min(1,i))*(Xmax-Xmin)];
end
Setslog.DataSet_sc = DataSet_sc;
Setslog.Data_min = Data_min;
Setslog.Data_max = Data_max;

eval(['delete ' savedir '\Y' num2str(out-out_index(1)+1) '.mat'])
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
    Setslog.dataset(i).in = Setslog.DataSet_sc(chosen_data_index,in_index);
    Setslog.dataset(i).out = Setslog.DataSet_sc(chosen_data_index,Setslog.out_index);
end
%====================================

%train neural nets
for i = 1:no_run
    if i < no_run
        fprintf('\nTraining Dataset %d\n\n', i);
    else
        fprintf('\nTraining Whole Dataset\n\n');
    end
    PP_NNGA(i);
    pause(5)
end
fprintf('\n\nTrainining Subsets\n');
disp(training_subsets);
%====================================
disp(['save ' savedir '\Y' num2str(out-out_index(1)+1) '.mat Setslog'])
eval(['save ' savedir '\Y' num2str(out-out_index(1)+1) '.mat Setslog'])
save([savedir '\parameters.m'],'parameters');
copyfile([pwd '\EvoNN\evaluate_obj.m'], savedir);

close all;
autooutput(Setslog,savedir);
close all;
svr(Setslog,[savedir '\svr_Y' num2str(out-out_index(1)+1)],'trend.mat');
close all;
end
end