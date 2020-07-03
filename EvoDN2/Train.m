function Train(Problem_name, param)
	%myFun - Description
	%
	% Syntax: Train(input)
	%
	% Trains Deep Neural Networks using PPGA
	clc
	%=====================DND==================================================
	global figure_handle parameters
    
	figure_handle = [];
	RandStream('mt19937ar','seed', sum(100*clock));
	parameters = [];
    parameters = param;
	%=====================Define parameters=======================================
	filename = [Problem_name '.xlsx']; %Data file
    
	in_index = parameters.in_index;           %independent variables column no.
    out_index = parameters.out_index;         %dependent variable column no.
    parameters = parameters.EvDtrain;
	parameters.in_index = in_index;           %independent variables column no.
    parameters.out_index = out_index;
	savedir = fullfile(pwd,'Output',Problem_name,'EvoDN2',param.name); % remove parameter_name to match mastercode% SwagataRoy
	mkdir(savedir);

	Xmin = eps; parameters.Xmin = Xmin;          %normalization range for variables
	Xmax = 1; parameters.Xmax = Xmax;

	parameters.in_index = in_index;

    [DataSet,paraname,DATA] = xlsread(filename);
	parameters.DATA = DATA;
	parameters.paraname= paraname;
	parameters.DataSet = DataSet;
    parameters.noinnodes = length(in_index);
	parameters.nooutnodes = 1;
    parameters.maxrank = 20;
    parameters.KillInterval =10;
    subsets = 1; overlap = 0;           %number of partitions of datafile and overlap b/w them
	set_size = length(DataSet(:,1)); parameters.set_size = set_size;
	parameters.subsets = subsets;
	parameters.overlap = overlap;
	subset_size = round((set_size + (subsets-1)*overlap) / subsets);
	parameters.subset_size = subset_size;

	%======================================================================
	for out = out_index
		parameters.out_index = out;
		if parameters.ploton
			figure_handle = [figure(1) figure(2) figure(3) figure(4)];
			scrsz = get(0,'ScreenSize');
			set(figure_handle(1), 'OuterPosition', [0*scrsz(3) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(1))
			set(figure_handle(2), 'OuterPosition', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(2))
			set(figure_handle(3), 'OuterPosition', [0*scrsz(3) 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(3))
			set(figure_handle(4), 'OuterPosition', [scrsz(3)/2 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf(figure_handle(4))
			pause(40)
		end

		%scale the data
		DataSet_sc = [];
		for i=1:length(DataSet(1,:))
			Data_min(1,i) = min(DataSet(:,i));
			Data_max(1,i) = max(DataSet(:,i));
			DataSet_sc = [DataSet_sc Xmin+ (DataSet(:,i)-Data_min(1,i))/(Data_max(1,i)-Data_min(1,i))*(Xmax-Xmin)];
		end
		parameters.DataSet_sc = DataSet_sc;
		parameters.Data_min = Data_min;
		parameters.Data_max = Data_max;

		eval(['delete ' savedir '\Y' num2str(out-out_index(1)+1) '.mat'])
		%============================================

		%choose datasets for training

		if subsets == 1
			training_subsets = 1;
			no_run = 1; parameters.no_run = no_run;
		else
			training_subsets = [eye(subsets,subsets); ones(1,subsets)];
			no_run = 1+subsets; parameters.no_run = no_run;
		end
		parameters.chosen_training_sets = training_subsets;
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
			parameters.dataset(i).data_index = chosen_data_index;
			parameters.dataset(i).in = parameters.DataSet_sc(chosen_data_index,in_index);
			parameters.dataset(i).out = parameters.DataSet(chosen_data_index,parameters.out_index);
		end
		%====================================

		%train neural nets
		for i = 1:no_run
			if i < no_run
				fprintf('\nTraining Dataset %d\n\n', i);
			else
				fprintf('\nTraining Whole Dataset\n\n');
			end
			PPGA();
			pause(5)
		end
		fprintf('\n\nTrainining Subsets\n');
		disp(training_subsets);
		%====================================
		eval(['save ' savedir '\Y' num2str(out-out_index(1)+1) '.mat parameters'])
        disp(['save ' savedir '\Y' num2str(out-out_index(1)+1) '.mat parameters'])
        save([savedir '\parameters.m'],'parameters');
		copyfile([pwd '\EvoDN2\evaluate_obj.m'], savedir);
		close all;
		autooutput(parameters,savedir);
        close all;
        svr(parameters,[savedir '\svr_Y' num2str(out-out_index(1)+1)],'trend.mat');
        close all;
	end
end