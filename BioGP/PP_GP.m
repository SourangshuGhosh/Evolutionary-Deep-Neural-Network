function Prey = PP_GP(set_no)
%Created by Brijesh Kumar Giri
%M.Tech Project

%General PP parameters
global generation
global lattice no_x no_y
global Setslog setno

%General GP parameters
global F T in_index out_index
global max_depth max_arity min_arity
global F_bad

%initialize parameters
global choose_term

%selection parameters
global tour_size

%xover parameters
global choose_xfunc Px_stand Px_fair P_xover
%mutation parameters
global Pm_stand Pm_small Pm_mono P_mut

%General PP parameters
setno = set_no;                                 %training dataset number
no_generations = Setslog.generations;           %max generations
Prey_popsize = Setslog.Prey_popsize;            %Initial popsize
no_Prey_preferred = Setslog.no_Prey_preferred;  %Desired popsize
no_new_Prey = Setslog.no_new_Prey;              %new prey introduced every KillInterval
Predator_popsize = Setslog.Predator_popsize;    %Number of Predators
P_move_prey = 0.5;                              %Prob with which a Prey moves

%2D-lattice
no_x = Setslog.no_x;                            %lattice size (no of rows)
no_y = Setslog.no_y;                            %lattice size (no of cols)
KillInterval = Setslog.KillInterval;            %Interval at which bad preys are eliminated
maxrank = Setslog.maxrank;                      %maxrank retained at KillInterval

%General GP parameters
F = Setslog.F;                                  %Function set
T = Setslog.T;                                  %Terminal set
in_index = Setslog.in_index;                    %input variable cols in datafile
out_index = Setslog.out_index;                  %output variable col in datafile

max_depth = Setslog.max_depth;                  %max depth to which a tree grows
max_roots = Setslog.max_roots;                  %max subtrees that a tree grows
max_arity = max(F.sets);                        %max arity of function set
min_arity = min(F.sets);                        %min arity of function set
F_bad = 1e6;                                    %fitness assigned to preys performing badly

%initialize parameters
choose_term = 0.2;                              %prob of choosing a node to be a terminal

%selection parameters
tour_size = Setslog.tour_size;                  %tournament size for single objective GP

%xover parameters
choose_xfunc = 0.9;                             %prob of choosing a function node for Xover
P_xover = 1;                                    %Xover prob of Preys
Px_stand = 0.5; Px_stand = P_xover*Px_stand;    %fraction of standard xover
Px_fair = P_xover - Px_stand;                   %fraction of Hieght fair xover

%mutation parameters
P_mut = 0.3;                                    %prob of choosing a prey for mutation
Pm_stand = 1/3; Pm_stand = Pm_stand*P_mut;      %fraction of standard mutation
Pm_small = 1/3; Pm_small = Pm_small*P_mut;      %fraction of small mutation
Pm_mono = P_mut - Pm_stand - Pm_small;          %fraction of monoparental mutation

%**************************************************************************
%   INITIALIZATION
lattice = zeros(no_x+2,no_y+2);
%initialization of Prey_pop
Prey = initialize(Prey_popsize, max_roots, max_depth);
[Prey F1 F2] = GPevaltree(Prey, 0);
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

%**************************************************************************

%**************************************************************************
% GENERATIONS
for generation = 1:no_generations
    if generation <= Setslog.generation1 && Setslog.evo_type == 2
        Prey = tournament_select(Prey);
        Prey = xover(Prey, 0);
        Prey = mutate(Prey);
        [Prey F1 F2]  = GPevaltree(Prey, 0); a = find(F1 == min(F1));
        disp(['Gen ' num2str(generation) '; Best Ind fit: ' num2str(F1(a(1)))])
        PlotPareto(F1, F2, [])
        continue
    end
    
    Prey_new = [];
    len = length(Prey);
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
            Offspring = Prey([i;parent2]);
            Offspring = xover(Offspring, 0);
            %**********Mutations*********************
            for j = 1:2
                Offspring(j) = mutate(Offspring(j));
            end
            %Random location of offspring /10 trials
            %Placing offsprings
            for l = 1:2
                for j = 1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos) == 0
                        Prey = [Prey; Offspring(l)];
                        lattice(xpos,ypos) = length(Prey);
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end
    
    [Prey F1 F2] = GPevaltree(Prey, 0);
    [fonrank front] = NONDOM_SORT([F1 F2]);
    
    % Removing all except rank <= maxrank
    if generation/KillInterval == round(generation/KillInterval)
        if generation == no_generations
            maxrank = 0;
        else
            Prey_new = initialize(no_new_Prey, max_roots, max_depth);
        end
        indfr = find(fonrank > maxrank);
        F1(indfr) = F_bad+eps; F2(indfr) = F_bad+eps;
        [Prey F1 F2] = KillBadPrey(Prey, F1, F2);
        [fonrank front] = NONDOM_SORT([F1 F2]);
    end
    
    % Move of predators /killing
    crodit = CROW_SORT([F1 F2], front);
    F1 = (F1 - min(F1))./(max(F1) - min(F1));
    F2 = (F2 - min(F2))./(max(F2) - min(F2));
    PredMoves = floor((length(Prey) - no_Prey_preferred)/Predator_popsize);
    fprintf('\nGeneration %i: Predatorpop %i PredMoves %i\n',generation,length(Predators(:,1)),PredMoves);
    fprintf('Preypop before: %i; Preypop after: ', length(Prey))
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
                f = (Predators(i)*F1(rows) + (1-Predators(i))*F2(rows)).*(front(rows)-1); % elitism for front = 1
                if length(unique(front(rows))) == 1
                    f = (f+1)./(1+crodit(rows));
                end
                [~,pos] = max(f);
                j = rows(pos(1)); lattice(xpos,ypos) = 0; %removing predator i
                [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == j);
                xpos = xpos+1; ypos = ypos+1;
                Prey = MovePredator(Prey, xpos, ypos, i, j);
                F1(j) = []; F2(j) = []; front(j) = []; crodit(j) = []; fonrank(j) = [];
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
    fprintf('%i\n' , length(Prey))
    
    [Prey F1 F2] = GPevaltree(Prey, 0);
    PlotLattice
    PlotPareto(F1, F2, fonrank)
    
    % Placing new Prey
    for i = 1:length(Prey_new)
        for k = 1:10 %10 trial for free spot
            [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
            if ~isempty(emptyx > 0)
                j = ceil(rand*length(emptyx));
                lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey) + 1;
                Prey = [Prey; Prey_new(i)]; break
            end
        end
    end
    
end

% [Prey F1 F2] = GPevaltree(Prey, 0);
[Prey F1 F2] = GPevaltree(Prey, 1);
[fonrank, front] = NONDOM_SORT([F1 F2]);
PlotPareto(F1, F2, fonrank)

Prey = Prey(front == 1); F1 = F1(front == 1); F2 = F2(front == 1);
f = [(1:length(F1))' F1 F2]; f = sortrows(f, 3); f = sortrows(f, 2); i = 2;
while i ~= length(f(:,1))
    if f(i,2) == f(i-1,2), f(i,:) = []; else i = i+1; end
end
Setslog.dataset(setno).pareto.P = Prey(f(:,1));
Setslog.dataset(setno).pareto.F1 = F1(f(:,1));
Setslog.dataset(setno).pareto.F2 = F2(f(:,1));

end