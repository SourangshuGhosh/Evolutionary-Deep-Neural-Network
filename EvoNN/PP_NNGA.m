function PP_NNGA(train_set_no)

global nonodes noinnodes nooutnodes
global lattice no_x no_y
global generation
global Setslog setno F_bad

setno = train_set_no;

%NNGA definitions
nonodes = Setslog.nonodes; % %maximum number of nodes
noinnodes = Setslog.noinnodes;
nooutnodes = Setslog.nooutnodes;
P_omit_max = 0.99; %After ranking max value for random generated P_omit
P_omit_min = 0.3; %-''- min value 
P_node_xover = 0.8; %Prob with which a node is exchanged
P_mutation = 0.3;   %Prob with which a connection mutates
Mut_alfa = 0.7;     %Mutation parameter
wlow = -5; %Lower bound for randomly generated weights w
whigh = 5; %Upper bound for randomly generated weights w

%PP definitions
Prey_popsize = Setslog.Prey_popsize;                %Initial popsize
no_Prey_preferred = Setslog.no_Prey_preferred;      %Desired popsize
Predator_popsize = Setslog.Predator_popsize;
no_new_Prey = Setslog.no_new_Prey;
no_generations = Setslog.generations;
P_move_prey = 0.3;
%2D-lattice
no_x = Setslog.no_x;
no_y = Setslog.no_y;
KillInterval = Setslog.KillInterval;
maxrank = Setslog.maxrank;

P_omit_match = P_omit_min;
F_bad = 1e6;        %fitness assigned to preys performing badly

%**************************************************
%   INITIALIZATION
lattice = zeros(no_x+2,no_y+2);
%initialization of Prey_pop
Prey = rand(Prey_popsize,nonodes,noinnodes+1)*(whigh-wlow)+wlow;
%Eliminating some matches
for i = 1:Prey_popsize
    for j = 1:nonodes
        for k = 1:noinnodes
            if rand(1,1) < P_omit_match
                Prey(i,j,k+1) = 0;                
            end
        end
    end
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
for generation=1:no_generations
    
    Prey_new = [];
    len = length(Prey(:,1,1));
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
            %Exchanging "node-packages"
            Offsprng = Prey([i;parent2],:,:);
            for j = 1:nonodes
                if rand < P_node_xover
                    Offsprng_tmp = Offsprng(1,j,:);
                    Offsprng(1,j,:) = Offsprng(2,j,:);
                    Offsprng(2,j,:) = Offsprng_tmp;
                end
            end
            %**********Mutations*********************
            for l = 1:2
                for j = 1:nonodes
                    for k = 1:noinnodes
                        if Offsprng(l,j,k) ~= 0 && rand < P_mutation
                            %Find randomly two other individuals with current match active
                            alternatives = find(Prey(:,j,k) ~= 0);
                            alternatives(alternatives == i) = [];
                            if length(alternatives) >= 2
                                select1 = ceil(rand*length(alternatives));
                                tmp = select1;
                                select1 = alternatives(select1);
                                alternatives(tmp) = [];
                                select2 = ceil(rand*length(alternatives));
                                select2 = alternatives(select2);
                                %Self adapting mutation
                                Offsprng(l,j,k) = Offsprng(l,j,k)+Mut_alfa*(1-generation/no_generations)*(Prey(select1,j,k)-Prey(select2,j,k));
                            end
                        end
                    end
                end
            end
            %Random location of offspring /10 trials
            %Placing offspring
            for l = 1:2
                for j=1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos)==0
                        Prey = [Prey; Offsprng(l,:,:)];
                        lattice(xpos,ypos) = length(Prey(:,1,1));
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end    

    [F1 F2 UW IC] = EvalF1F2(Prey); 
    [fonrank front] = NONDOM_SORT([F1 F2]);
    
    % Removing all except rank n
    if generation/KillInterval == round(generation/KillInterval)
        if generation < no_generations
            Prey_new = rand(no_new_Prey,nonodes,noinnodes+1)*(whigh-wlow)+wlow;
            P_omit_match = (P_omit_max-P_omit_min)*(generation/(no_generations-1))+P_omit_min;
            for i = 1:no_new_Prey
                for j = 1:nonodes
                    for k = 1:noinnodes
                        if rand < P_omit_match
                            Prey_new(i,j,k+1) = 0;
                        end
                    end
                end
            end
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
                    dx=round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                    dy=round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
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
    
    [F1 F2 UW IC] = EvalF1F2(Prey);
    PlotLattice
    PlotPareto(F1, F2, fonrank)
    
    % Placing new Prey
    for i = 1:length(Prey_new)
        for k = 1:10 %10 trial for free spot
            [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
            if ~isempty(emptyx > 0)
                j = ceil(rand*length(emptyx));
                lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey(:,1,1)) + 1;
                Prey = [Prey; Prey_new(i,:,:)]; break
            end
        end
    end

end

[F1 F2 UW IC] = EvalF1F2(Prey);
[fonrank, front] = NONDOM_SORT([F1 F2]);
PlotPareto(F1, F2, fonrank); figure(4)

Prey = Prey(front == 1,:,:); F1 = F1(front == 1); F2 = F2(front == 1);
UW = UW(front == 1); IC = IC(front == 1);
f = [(1:length(F1))' F1 F2]; f = sortrows(f, 3); f = sortrows(f, 2); i = 2;
while i ~= length(f(:,1))
    if f(i,2) == f(i-1,2), f(i,:) = []; else i = i+1; end
end
F1 = F1(f(:,1)); F2 = F2(f(:,1)); UW = UW(f(:,1));
IC = IC(f(:,1)); [~,i] = min(IC); hold on
Setslog.dataset(setno).pareto.select = i;
plot(F2(i), F1(i), 'dr', 'LineWidth', 2)

P = [];
for i = 1:length(f(:,1))
    for j = 1:nonodes
        for k = 1:noinnodes+1
            P(i).w(j,k) = Prey(f(i,1),j,k);
        end
    end
    P(i).W = UW(i).W; P(i).IC = IC(i);
end

Setslog.dataset(setno).pareto.P = P;
Setslog.dataset(setno).pareto.F1 = F1;
Setslog.dataset(setno).pareto.F2 = F2;

end
%-------------------------------------------------
% Evaluate F1 and F2
%-------------------------------------------------
function[F1 F2 UW IC] = EvalF1F2(Prey)
global nonodes noinnodes F_bad

for i = 1:length(Prey(:,1,1))
    for j = 1:nonodes
        for k = 1:noinnodes+1
            w(j,k) = Prey(i,j,k);
        end
    end
    [fval,W,InfoC] = NNevalnet(w);
    F1(i) = fval;
    if isnan(F1(i)) || isinf(F1(i))
        F1(i) = F_bad+eps; F2(i) = F_bad+eps;
    end
    F2(i) = length(find(w(:,2:noinnodes+1))); %Does not include bias parameters!
    UW(i).W = W;
    IC(i) = InfoC;
    %Should somehow be able to retrieve the different networks...
end
F1 = F1'; F2 = F2'; IC = IC'; UW = UW';
end