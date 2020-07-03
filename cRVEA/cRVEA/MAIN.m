%The main file of constrained RVEA (cRVEA)
function MAIN(Problem,M,parameters,savedir)
%clc;
format compact;tic;
disp(Problem);
%basic settings
[Generations,N,p1,p2] = P_settings('cRVEA',Problem,M,parameters);
Evaluations = Generations*N; % max number of fitness evaluations
alpha = parameters.cRVEAopt.alpha; % the parameter in APD, the bigger, the faster RVEA converges
fr = parameters.cRVEAopt.fr; % frequency to call reference vector
obj = parameters.cRVEAopt.obj;
eqCon = parameters.cRVEAopt.eqCon;
ieqCon = parameters.cRVEAopt.ieqCon;
FE = 0; % fitness evaluation counter
FinalSelection = 0; % flat to indicate if it is the late generation

for i = 1:M
    obj_val(i) = loadwlog(importdata([savedir '\Y' num2str(i) '.mat']),0,0, Problem); %see PPopt for reference.
end


%reference vector initialization
[N,Vs] = F_weight(p1,p2,M);
for i = 1:N
    Vs(i,:) = Vs(i,:)./norm(Vs(i,:));
end;
V = Vs;
Generations = floor(Evaluations/N);

%calculat neighboring angle for angle normalization
cosineVV = V*V';
[scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
acosVV = acos(scosineVV(:,2));
refV = (acosVV);

%population initialization
[Population,Boundary,Coding] = P_objective('init',obj_val,M,N);
FunctionValue = P_objective('value',obj_val,M,Population,obj,eqCon,ieqCon);

for Gene = 0 : Generations - 1
    if(Gene == Generations - 1)
        FinalSelection = 1;
    end;

    %random mating and reproduction
    [MatingPool] = F_mating(Population);
    Offspring = P_generator(MatingPool,Boundary,Coding,N);  FE = FE + size(Offspring, 1);
    Population = [Population; Offspring];
    FunctionValue = [FunctionValue; P_objective('value',obj_val,M,Offspring,obj,eqCon,ieqCon);];
    
    %APD based selection
    theta0 =  (Gene/(Generations))^alpha*(M);

    [Selection] = F_select(FunctionValue,V, theta0, refV, FinalSelection);
    Population = Population(Selection,:);
    FunctionValue = FunctionValue(Selection,:);

    %reference vector adaption
    if(mod(Gene, ceil(Generations*fr)) == 0)
        %update the reference vectors
        Zmin = min( FunctionValue(:, 1:end - 1),[],1);
        Zmax = max( FunctionValue(:, 1:end - 1),[],1);
        V = Vs;
        V = V.*repmat((Zmax - Zmin)*1.0,N,1);
        for i = 1:N
            V(i,:) = V(i,:)./norm(V(i,:));
        end;
        %update the neighborning angle value for angle normalization
        cosineVV = V*V';
        [scosineVV, neighbor] = sort(cosineVV, 2, 'descend');
        acosVV = acos(scosineVV(:,2));
        refV = (acosVV); 
    end;

%    clc; fprintf('Progress %4s%%\n',num2str(round(Gene/Generations*100,-1)));
    clc; fprintf('Progress %4s%%\n',num2str(round(Gene/Generations*100)));
end;
P_output(Population,toc,'cRVEA',Problem,M, obj_val,savedir,obj,eqCon,ieqCon);
end

function wlog = loadwlog(Setslog, setno, index,Problem)
switch Problem
    case 'EvoNN'
        wlog.nonodes = Setslog.nonodes;
        wlog.noinnodes = Setslog.noinnodes;
        wlog.nooutnodes = Setslog.nooutnodes;
        wlog.Xmin = Setslog.Xmin;
        wlog.Xmax = Setslog.Xmax;
        f_index = Setslog.out_index;
        wlog.ymin = Setslog.Data_min(1,f_index);
        wlog.ymax = Setslog.Data_max(1,f_index);
        wlog.xmin = Setslog.Data_min(1,Setslog.in_index);
        wlog.xmax = Setslog.Data_max(1,Setslog.in_index);
        if index == 0, index = Setslog.dataset(1).pareto.select; end
        if setno == 0, setno = Setslog.no_run; end
        Net = Setslog.dataset(setno).pareto.P(index);
        wlog.w = Net.w;
        wlog.W = Net.W;
        wlog.in_index = Setslog.in_index;
        wlog.out_index = f_index;
    case 'BioGP'
        wlog.T = Setslog.T;
        f_index = Setslog.out_index; in_index = Setslog.in_index;
        wlog.DataSet = Setslog.DataSet;
        wlog.ymin = min(Setslog.DataSet(:,f_index));
        wlog.ymax = max(Setslog.DataSet(:,f_index));
        wlog.xmin = min(Setslog.DataSet(:,in_index));
        wlog.xmax = max(Setslog.DataSet(:,in_index));
        if index == 0, index = 1; end
        if setno == 0, setno = Setslog.no_run; end
        wlog.func_tree = Setslog.dataset(setno).pareto.P(index);
        wlog.in_index = in_index;
        wlog.out_index = f_index;
    case 'EvoDN2'
        setno = 1;
        select = 1;
       % wlog.nonodes = Setslog.nonodes;
        wlog.noinnodes = Setslog.noinnodes;
        wlog.nooutnodes = Setslog.nooutnodes;
        wlog.Xmin = Setslog.Xmin;
        wlog.Xmax = Setslog.Xmax;
        f_index = Setslog.out_index;
        wlog.ymin = Setslog.Data_min(1,f_index);
        wlog.ymax = Setslog.Data_max(1,f_index);
        wlog.xmin = Setslog.Data_min(1,Setslog.in_index);
        wlog.xmax = Setslog.Data_max(1,Setslog.in_index);
        wlog.in = Setslog.dataset(setno).in;
        wlog.out = Setslog.dataset(setno).out;
        wlog.Pop_str = Setslog.Pop_str;
        wlog.net = (Setslog.dataset.pareto.P(select));
        wlog.endnet = (Setslog.dataset.pareto.P(select).endnet);
        num_subnets = length(wlog.Pop_str);
        wlog.in_end = [];
        %select = Setslog.dataset.pareto.select;
        
        outmax = max(wlog.out);
        outmin = min(wlog.out);
        no_points = length(wlog.out);
%         for subnet = 1:num_subnets
%             wlog.in = Setslog.dataset(setno).in;
%             wlog.in = wlog.in(:,wlog.Pop_str{subnet}{1});
%             no_layer = length(wlog.Pop_str{subnet}{2});
%             wlog.snet = wlog.net.subnet{subnet};
%             wlog.in_end = [wlog.in_end wlog.in];
%         end
%          wlog.Beta = Setslog.Beta{select};
         in_index = Setslog.in_index;
        wlog.in_index = Setslog.in_index;
        wlog.out_index = f_index;
         
end
end
