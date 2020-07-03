function PPGA()
global no_x no_y lattice F_bad
global parameters
%%Descrition: Use PPGA to create DNN models
%NNGA definitions
Pop_str = parameters.Pop_str; %Cell defining the input varables and 
								%structure of various subnets
P_omit_max = 0.99;  %After ranking max value for random generated P_omit
P_omit_min = 0.3;   %-''- min value 

%PP definitions
Prey_popsize = parameters.Prey_popsize;                %Initial popsize
no_Prey_preferred = parameters.no_Prey_preferred;      %Desired popsize
Predator_popsize = parameters.Predator_popsize;
newPrey_popsize = parameters.no_new_Prey;
no_generations = parameters.generations;
P_move_Prey = 0.3;

%2D-lattice
no_x = parameters.no_x;
no_y = parameters.no_y;

%Killing of bad Prey
KillInterval = parameters.KillInterval;
maxrank = parameters.maxrank;
P_omit_match = P_omit_min;
F_bad = 1e6;%fitness assigned to Preys performing badly

lattice = zeros(no_x+2,no_y+2);[Prey, FVal] = create_Population(Prey_popsize, Pop_str);

%Placement of Prey
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

for num_Gen = 1:no_generations
    Pop_size = length(Prey);
    Prey_new = [];
    % Calculate Standard deviation in population as mutval here%

    % Movement of all population members
    for Pop_index = 1 : Pop_size
        if rand < P_move_Prey
            %Identify location
            [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == Pop_index);
            xpos = xpos+1; ypos = ypos+1;
            for j = 1:10 %10 trial for free spot
                dx = round(rand*(3-eps)-1.5); %-1 to the left, 0 no move, 1 to the right
                dy = round(rand*(3-eps)-1.5); %-1 down, 0 no move, 1 up
                if lattice(xpos+dx,ypos+dy) == 0
                    lattice(xpos,ypos) = 0; %removing prey i
                    MovePrey(xpos+dx, ypos+dy, Pop_index); 
                    break
                end
            end
        end
    end

    % Breeding of all population members
    for Pop_index = 1 : Pop_size
        [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == Pop_index);
        xpos = xpos+1; ypos = ypos+1;
        moore = lattice(xpos-1:xpos+1,ypos-1:ypos+1);
        %Remove the i-prey
        moore(2,2) = 0;
        [matex,matey] = find(moore >= 1 & moore <= Pop_size);
        if ~isempty(matex)
            parent2 = ceil(rand*length(matex));
            parent2 = lattice(xpos-2+matex(parent2),ypos-2+matey(parent2));
            % Write the following function
            [Offsprng, FVal_offsp] = create_offspring(Pop_index, parent2, ...
                                                     Prey, Pop_str, num_Gen, ...
                                                     no_generations); 
            %Random placement of offsprings /10 trials
            for l = 1:2
                for j=1:10 %10 trial for free spot
                    xpos = round(rand(1,1)*(no_x-1)+1)+1;
                    ypos = round(rand(1,1)*(no_y-1)+1)+1;
                    if lattice(xpos,ypos)==0
                        Prey = [Prey Offsprng(l)];
                        FVal = [FVal ; FVal_offsp(l,:)];
                        lattice(xpos,ypos) = length(Prey);
                        break
                    end
                end
            end
        end
        MovePrey(0, 0, 0);
    end

    [fonrank, front] = NONDOM_SORT([FVal]);  % Check this function

    % Removing all except rank n
    if num_Gen/KillInterval == round(num_Gen/KillInterval)
        if num_Gen < no_generations
            % Generate new prey to balance the population
            [Prey_new, FVal_new] = create_Population(newPrey_popsize, Pop_str);
        end
        indfr = find(fonrank > maxrank);
        FVal(indfr,:) = F_bad+eps;
        [Prey FVal] = KillBadPrey(Prey, FVal);  % Write this function
        [fonrank front] = NONDOM_SORT([FVal]);  % Write this function
    end

    crodit = CROW_SORT([FVal], front); %   Write this function
    FValKill = [];
    FValKill(:,1) = (FVal(:,1) - min(FVal(:,1)))./(max(FVal(:,1)) - min(FVal(:,1)));
    FValKill(:,2) = (FVal(:,2) - min(FVal(:,2)))./(max(FVal(:,2)) - min(FVal(:,2)));

    % Move of predators /killing
    PredMoves = floor((length(Prey) - no_Prey_preferred)/Predator_popsize);
    fprintf('\nGeneration %i: Predatorpop %i PredMoves %i\n',num_Gen,length(Predators(:,1)),PredMoves);
    fprintf('Preypop before: %i; Preypop after: ', length(Prey))
    for i = 1:Predator_popsize
        for k = 1:PredMoves
            [xpos, ypos] = find(lattice(2:no_x+1,2:no_y+1) == -i);
            xpos = xpos+1; ypos = ypos+1;
            [matex,matey] = find(lattice(xpos-1:xpos+1,ypos-1:ypos+1) > 0);
            if length(matex) > 1 %prey available
                rows=[];
                for t = 1:length(matex)
                    rows = [rows; lattice(matex(t)+xpos-2,matey(t)+ypos-2)];
                end

                % Calculating Fitness Value for all Prey near Predator i, with Elitism
                f = (FValKill(rows,:)*[Predators(i);1-Predators(i)]).*(front(rows)-1);

                % Using Crowding in case all Prey have same fitness
                if length(unique(front(rows))) == 1
                    f = (f+1)./(1+crodit(rows));
                end

                % Killing Prey
                [~,pos] = max(f);
                j = rows(pos(1)); lattice(xpos,ypos) = 0; %removing predator i
                [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == j);
                xpos = xpos+1; ypos = ypos+1;
                Prey = MovePredator(Prey, xpos, ypos, i, j);
                FValKill(j,:) = []; FVal(j,:) = [];
                front(j) = []; crodit(j) = []; 
                fonrank(j) = []; %CHECK!!!!! => checked, works => see MovePredator
            else    %Only move
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
    PlotLattice % Write this function
    PlotPareto(FVal, fonrank, num_Gen)   % Write this function

    %Random placement of New Prey /10 trials
    for i = 1:length(Prey_new)
        for k = 1:10 %10 trial for free spot
            [emptyx,emptyy] = find(lattice(2:no_x+1,2:no_y+1) == 0);
            if ~isempty(emptyx > 0)
                j = ceil(rand*length(emptyx));
                lattice(emptyx(j)+1,emptyy(j)+1) = length(Prey) + 1;
                Prey = [Prey Prey_new(i)];
                FVal = [FVal; FVal_new(i,:)];
                break
            end
        end
    end
end

% Choosing only Prey on the Pareto Front
[fonrank, front] = NONDOM_SORT([FVal]);
PlotPareto(FVal, fonrank, num_Gen); figure(4)
Prey = Prey(front==1);
FVal = FVal(front==1,:);

% Removing Weakly dominated Prey
f = [(1:length(FVal(:,1)))' FVal]; 
f = sortrows(f, 3); f = sortrows(f, 2); 
[~,index] = sortrows([Prey.err].'); Prey = Prey(index); clear index;
i = 2;
while i ~= length(f(:,1))
    if f(i,2) == f(i-1,2)
		f(i,:) = []; 
        Prey(i) = [];
    else 
        i = i+1; 
    end
end

FVal = f(:,2:3);

% Selecting Prey based on Corrected Akaike Information criteria
% for i = 1:length(Prey)
%     IC = Prey(i).IC;
% end
% [~,i] = min(IC); 
hold on;
% parameters.dataset(setno).pareto.select = i;
 plot(FVal(1,2), FVal(1,1), 'dr', 'LineWidth', 2)

% Saving Output
setno = 1;
parameters.dataset(setno).pareto.P = Prey;
% parameters.dataset(setno).pareto.info = IC;
parameters.dataset(setno).pareto.FVal = FVal;
