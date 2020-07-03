function [Offsprng, F_val] = create_offspring(i,parent2,Prey,Pop_str,generation,no_generations)
%myFun - Description
%
% Syntax: children = create_offspring(Prey, i,parent2,P_node_xover,P_mutation,Mut_alfa)
%
% Long description
P_node_xover = 0.8; %Prob with which a node is exchanged
P_mutation = 0.3;   %Prob with which a connection mutates
Mut_alfa = 0.7;     %Mutation parameter
P_kill_connection = 0.1;
Offsprng(1).subnet = Prey(i).subnet;
Offsprng(2).subnet = Prey(parent2).subnet;
cross_technique = 'short';
mut_technique = 'short';
F_val = [];
for subnet = 1:length(Offsprng(1).subnet)
	sub1 = Offsprng(1).subnet{subnet};
	sub2 = Offsprng(2).subnet{subnet};
	randsub = ceil(rand(2,1)*length(Prey));
	sub4 = Prey(randsub(1)).subnet{subnet};
	sub3 = Prey(randsub(2)).subnet{subnet};
	for layer = 1:length(sub1)
	    %**********Crossover*********************
		%Exchanging "node-packages"
		switch cross_technique
		case 'short'
			connections = numel(sub1{layer});
			exchange = binornd(connections,P_node_xover);
			exchange = randperm(connections,exchange);
			Offsprng_tmp = sub1{layer};
			sub1{layer}(exchange) = sub2{layer}(exchange);
			sub2{layer}(exchange) = Offsprng_tmp(exchange);
		end

		%**********Mutations*********************
		switch mut_technique
		case 'short'
			connections = numel(sub1{layer});
			mutate = binornd(connections,P_mutation);
			mutate = randperm(connections,mutate);
			sub1{layer}(mutate) = sub1{layer}(mutate) + Mut_alfa * ...
                                    (1-generation/no_generations) * ...
                                    (sub3{layer}(mutate) - sub4{layer}(mutate));
			mutate = binornd(connections,P_mutation);
			mutate = randperm(connections,mutate);
            sub2{layer}(mutate) = sub2{layer}(mutate) + Mut_alfa * ...
                                    (1-generation/no_generations) * ...
                                    (sub3{layer}(mutate) - sub4{layer}(mutate));
		end
	end
	Offsprng(1).subnet{subnet} = sub1;
	[err, complexity, endnet] = DNevalnet(Offsprng(1),Pop_str);
    Offsprng(1).err = err;
    Offsprng(1).complexity = complexity;
    Offsprng(1).endnet = endnet;
    F_val = [err, complexity];
	Offsprng(2).subnet{subnet} = sub2;
    [err, complexity, endnet] = DNevalnet(Offsprng(2),Pop_str);
    Offsprng(2).err = err;
    Offsprng(2).complexity = complexity;
    Offsprng(2).endnet = endnet;
    F_val = [F_val; [err, complexity]];
end