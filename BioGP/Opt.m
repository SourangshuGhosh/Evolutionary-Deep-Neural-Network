function Opt(Problem,Training_Algorithm, parameters, savedir)
if ~strcmp(Training_Algorithm,'BioGP')
return;
end
%Predator-Prey GA for Multi-Objective Optimization
clc
global LB UB obj LB_F UB_F
global Setslog figure_handle no_x no_y lattice
global generation F_bad setno

RandStream('mt19937ar','seed', sum(100*clock));
figure_handle = [];

% plst = {'\adf' '\PP_util' '\util'};
% for i = 1:length(plst)
%     path([pwd plst{i}], path);
% end
%===============User Input=======================================
obj(1) = parameters.BioOpt.obj(1);  %set 1 for min and -1 for max
obj(2) = parameters.BioOpt.obj(2);  %set 1 for min and -1 for max
LB_F = parameters.BioOpt.LB_F;
UB_F = parameters.BioOpt.UB_F;
Setslog = {};
for count = 1:length(parameters.out_index)
Setslog(count) = {[savedir '\Y' num2str(count) '.mat']};
end
% obj(1) = -1;             %set 1 for min and -1 for max
% obj(2) = -1;             %set 1 for min and -1 for max
% Setslog = {};
% Setslog(1) = {'Setslog2-3-6_data_num.mat'};
% Setslog(2) = {'Setslog2-3-7_data_num.mat'};
% if constraint is applicable
%for constraint function introduce Setslog(3) and so on
% Setslog(3) = {'Setslog2-10-8_data_mod_smb1.mat'};
% Setslog(4) = {'Setslog2-10-9_data_mod_smb1.mat'};
%loadwlog loads the tree from Setslog file, the prototype:
%loadwlog(Setslog, setno, index);
%setno is the training subset used and index is tree index obtained from
%the pareto front of training, 1 being the tree with lowest error seeting
%setno =0 takes the last training subset and index = 0 sets the default
%tree chosen by BioGP training
wlog(1) = loadwlog(importdata(Setslog{1}), 0, 0);
wlog(2) = loadwlog(importdata(Setslog{2}), 0, 0);
% wlog(3) = loadwlog(importdata(Setslog{3}), 0, 0);
% wlog(4) = loadwlog(importdata(Setslog{4}), 0, 0);
%define extra wlog for checking the contraint and edit the contra function
%as desired
novar = length(wlog(1).in_index);
LB = []; UB = []; svstr =[savedir '\BioGP_Pareto.mat'];
for i = 1:2
    LB = [LB; wlog(i).xmin];
    UB = [UB; wlog(i).xmax];
end
LB = max(LB); UB = min(UB);

Prey_popsize = parameters.BioOpt.Prey_popsize;         %Initial popsize
no_Prey_preferred = parameters.BioOpt.no_Prey_preferred;    %Desired popsize
no_new_Prey = parameters.BioOpt.no_new_Prey;          %new prey introduced every KillInterval
Predator_popsize = parameters.BioOpt.Predator_popsize;     %Number of Predators 100

no_generations = parameters.BioOpt.no_generations;       %max generations
P_move_prey = parameters.BioOpt.P_move_prey;          %Prob with which a Prey moves
P_mut = parameters.BioOpt.P_mut;                %prob of choosing a prey for mutation
%Prob of Xover is 1 for every Prey
F_bad = parameters.BioOpt.F_bad;                %fitness assigned to preys performing badly
%2D-lattice
no_x = parameters.BioOpt.no_x;                  %lattice size (no of rows) 50
no_y = parameters.BioOpt.no_y;                  %lattice size (no of cols) 50
KillInterval = parameters.BioOpt.KillInterval;          %Interval at which bad preys are eliminated
maxrank = parameters.BioOpt.maxrank;               %maxrank retained at KillInterval
ploton = parameters.BioOpt.ploton;                 %set 0 for no plots or 1 for plots at every generation

%=============================DND=========================================
Setslog = [];
Setslog.generations = no_generations; Setslog.ploton = ploton;
Setslog.Predator_popsize = Predator_popsize; Setslog.Prey_popsize = Prey_popsize;
if Setslog.ploton
    figure_handle = [figure(1) figure(2) figure(3) figure(4)];
    scrsz = get(0,'ScreenSize'); setno = [];
    set(figure_handle(1), 'OuterPosition', [0*scrsz(3) scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf
    set(figure_handle(2), 'OuterPosition', [scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]); clf
    set(figure_handle(3), 'OuterPosition', [0*scrsz(3) 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf
    set(figure_handle(4), 'OuterPosition', [scrsz(3)/2 0*scrsz(4) scrsz(3)/2 scrsz(4)/2]); clf
    pause(1)
end
%=============================DND=========================================
%**************************************************
%   INITIALIZATION
lattice = zeros(no_x+2,no_y+2);
%initialization of Prey_pop
Prey = rand(Prey_popsize,novar);
for i = 1:novar
    Prey(:,i) = Prey(:,i)*(UB(i)-LB(i))+LB(i);
end
%location of prey
for i = 1:length(Prey)
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = i;
    end
end
%initialization of Predator_pop
Predators = (0:1/(Predator_popsize-1):1)'; %value to weight F1
%location of predators
for i = 1:length(Predators(:,1))
    [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
    if ~isempty(emptyx > 0)
        j = ceil(rand*(length(emptyx)));
        lattice(emptyx(j)+1,emptyy(j)+1) = -i; %Predator indentified with neg values
    end
end
MovePrey(0, 0, 0);
%**************************************************
%   GENERATIONS
for generation = 1:no_generations  
    Prey_new = [];
    len = length(Prey(:,1));
    % Move of prey /10 trials /1step
    for i = 1:len
        if rand < P_move_prey
            %Identify location
            [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
            xpos = xpos+1; ypos = ypos+1;
            for j = 1:10 %10 trial for free spot
                dx = round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                dy = round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                if lattice(xpos+dx,ypos+dy) == 0
                    lattice(xpos,ypos) = 0; %removing prey i
                    MovePrey(xpos+dx, ypos+dy, i); break
                end
            end
        end
    end
    
    % Breeding of prey
    for i = 1:len
        [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == i);
        xpos = xpos+1; ypos = ypos+1;
        moore = lattice(xpos-1:xpos+1,ypos-1:ypos+1);
        %Remove the i-prey
        moore(2,2) = 0;
        [matex,matey] = find(moore >= 1 & moore <= len);
        if ~isempty(matex)
            parent2 = ceil(rand*length(matex));
            parent2 = lattice(xpos-2+matex(parent2),ypos-2+matey(parent2));
            %**********Crossover*********************
            Offspring = Prey([i;parent2],:);
            Offspring = XOVER(Offspring);
            %**********Mutations*********************
            for j = 1:2
                if rand < P_mut
                    Offspring(j,:) = MUTE(Offspring(j,:), generation);
                end
            end
            %Random location of offspring /10 trials
            %Placing offsprings
            for l = 1:2
                for j = 1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos) == 0
                        Prey = [Prey; Offspring(l,:)];
                        lattice(xpos,ypos) = length(Prey);
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end
    
    F1 = obj(1)*evalF(Prey, wlog, 1);
    F2 = obj(2)*evalF(Prey, wlog, 2);
    if parameters.BioOpt.useConstraints
        Fout = Contraints('BioGP',Prey, obj(1)*F1, obj(2)*F2, wlog);
        F1 = Fout(:,1); F2 = Fout(:,2);
    end
    F1 = obj(1)*F1;
    F2 = obj(2)*F2;
    [fonrank front] = NONDOM_SORT([F1 F2]);
    
    % Removing all except rank <= maxrank
    if generation/KillInterval == round(generation/KillInterval)
        if generation == no_generations
            maxrank = 0;
        else
            Prey_new = rand(no_new_Prey,novar);
            for i = 1:novar
                Prey_new(:,i) = Prey_new(:,i)*(UB(i)-LB(i))+LB(i);
            end
        end
        indfr = find(fonrank > maxrank);
        F1(indfr) = F_bad+eps; F2(indfr) = F_bad+eps;
        [Prey F1 F2] = KillBadPrey(Prey, F1, F2);
        [fonrank front] = NONDOM_SORT([F1 F2]);
    end
    
    % Move of predators /killing
    crodit = CROW_SORT([F1 F2], front);
    f1 = (F1 - min(F1))./(max(F1) - min(F1));
    f2 = (F2 - min(F2))./(max(F2) - min(F2));
    PredMoves = floor((length(Prey) - no_Prey_preferred)/Predator_popsize);
    fprintf('\nGeneration %i: Predatorpop %i PredMoves %i\n',generation,length(Predators(:,1)),PredMoves);
    fprintf('Preypop before: %i; Preypop after: ', length(Prey(:,1)))
    for i = 1:length(Predators)
        for k = 1:PredMoves
            [xpos, ypos] = find(lattice(2:no_x+1,2:no_y+1) == -i);
            xpos = xpos+1; ypos = ypos+1;
            [matex,matey] = find(lattice(xpos-1:xpos+1,ypos-1:ypos+1) > 0);
            if length(matex) > 1 %prey available
                rows=[];
                for t = 1:length(matex)
                    rows = [rows; lattice(matex(t)+xpos-2,matey(t)+ypos-2)];
                end
                f = (Predators(i)*f1(rows) + (1-Predators(i))*f2(rows)).*(front(rows)-1); % elitism for front = 1
                if length(unique(front(rows))) == 1
                    f = (f+1)./(1+crodit(rows));
                end
                [~,pos] = max(f);
                j = rows(pos(1)); lattice(xpos,ypos) = 0; %removing predator i
                [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == j);
                xpos = xpos+1; ypos = ypos+1;
                Prey = MovePredator(Prey, xpos, ypos, i, j);
                f1(j) = []; f2(j) = []; front(j) = []; crodit(j) = [];
                F1(j) = []; F2(j) = []; fonrank(j) = [];
            else %Only move
                for j=1:10 %10 trial for free spot
                    dx=round(rand(1,1)*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                    dy=round(rand(1,1)*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                    if lattice(xpos+dx,ypos+dy)==0
                        lattice(xpos,ypos)=0; %removing predator i
                        Prey = MovePredator(Prey, xpos+dx, ypos+dy, i, inf);
                        break
                    end
                end
            end
        end
    end
    fprintf('%i\n' , length(Prey(:,1)))

    PlotLattice
    PlotPareto(obj(2)*F2, obj(1)*F1, fonrank)
    
    % Placing new Prey
    for i = 1:length(Prey_new)
        for k = 1:10 %10 trial for free spot
            [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
            if ~isempty(emptyx > 0)
                j = ceil(rand*length(emptyx));
                lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey(:,1)) + 1;
                Prey = [Prey; Prey_new(i,:)]; break
            end
        end
    end
end

F1 = obj(1)*evalF(Prey, wlog, 1);
F2 = obj(2)*evalF(Prey, wlog, 2);
if parameters.BioOpt.useConstraints
    Fout = Contraints('BioGP',Prey, obj(1)*F1, obj(2)*F2, wlog);
    F1 = Fout(:,1); F2 = Fout(:,2);
end
F1 = obj(1)*F1;
F2 = obj(2)*F2;
[fonrank, front] = NONDOM_SORT([F1 F2]);
PlotPareto(obj(2)*F2, obj(1)*F1, fonrank)

Prey = Prey(front == 1,:); F1 = F1(front == 1); F2 = F2(front == 1);
f = [(1:length(F1))' F1 F2]; f = sortrows(f, 3); f = sortrows(f, 2); i = 2;
while i ~= length(f(:,1))
    if f(i,2) == f(i-1,2), f(i,:) = []; else i = i+1; end
end
pareto.P = Prey(f(:,1),:);
pareto.F1 = obj(1)*F1(f(:,1));
pareto.F2 = obj(2)*F2(f(:,1));
pareto.cons = [];
if length(wlog) > 2
    for i = 3:length(wlog)
        pareto.cons = [pareto.cons evalF(pareto.P, wlog, i)];
    end
end
eval(['save ' svstr ' pareto']);
saveas(figure_handle(4), [savedir '\BioGP_Pareto'], 'jpg')
saveas(figure_handle(4), [savedir '\BioGP_Pareto'], 'fig')
% for i = 1:length(plst)
%     rmpath([pwd plst{i}]);
% end
end

function wlog = loadwlog(Setslog, setno, index)
wlog.T = Setslog.T;
f_index = Setslog.out_index; 
in_index = Setslog.in_index;
wlog.ymin = min(Setslog.DataSet(:,f_index));
wlog.ymax = max(Setslog.DataSet(:,f_index));
wlog.xmin = min(Setslog.DataSet(:,in_index));
wlog.xmax = max(Setslog.DataSet(:,in_index));
if index == 0, index = 1; end
if setno == 0, setno = Setslog.no_run; end
wlog.func_tree = Setslog.dataset(setno).pareto.P(index);
wlog.in_index = in_index;
wlog.out_index = f_index;
end

function out = evalF(Prey, wlog, oi)

func_tree = wlog(oi).func_tree; 
for i = 1:length(wlog(oi).in_index)
    eval([wlog(oi).T.set{i} '= Prey(:,i);']);
end
coeff = func_tree.bias; fval = [];
for i = 1:func_tree.no_roots
    tree_val = eval(cat(2,func_tree.root(i).tree{:}));
    fval = [fval tree_val.*ones(length(Prey(:,1)),1)];
    coeff = [coeff; func_tree.root(i).w];
end
fval = [ones(length(Prey(:,1)),1) fval]*coeff;
out = fval;

end