function [Population, FVal] = create_Population(Pop_size, Pop_str)
%create_Poulation - Creates cell array of Population of subnet NNs
%
% Syntax: Population = create_Poulation(Pop_size, Pop_str)
%
% Pop_size descries the population size
% Pop_str is a cell array with n cells. Each cell describes a subnet.
% Pop_str{i}{1} is a vector of the column indices of the input
% Pop_str{i}{2} is a vector describing the number of nodes in each layer

% Population is an struct array of size Pop_size
% Struct Population contains: Cell array- subnet (DNN connections and bias values)
%                             2D matrix- endnet (Final layer that aggregates all subnets)
%                             double- complexity (complexity of the DNN)
%                             double- err (RMSE of the DNN on training data)
%
% Cell arrray subnet has size n (see up). Each element is itself a cell array.
% Population(i).subnet{j}{k} gives the 2D matrix of connections in the k-th layer
%                            of the j-th subnet of the i-th population member.

    P_omit_match = 0.3; % Probability with which a connection is randomly killed.
    whigh = 5;  %Max value for a connection
    wlow = -5;  %min value for a connection

    Population = [];
    FVal = []; %contains Error and complexity
    Num_subnets = length(Pop_str);
    for Pop_index = 1:Pop_size
        %Population(Pop_index).structure = Pop_str;
        for s_index = 1:Num_subnets
            for layer = 1:length(Pop_str{s_index}{2})-1
                in = Pop_str{s_index}{2}(layer);    % number of input nodes
                out = Pop_str{s_index}{2}(layer+1); % number of output nodes
                % random genration of layer
                net = rand(in+1,out)*(whigh-wlow)+wlow;   % +1 for bias node
                % Killing some connections
                for i = 1:in+1
                    for j = 1:out
                        if rand(1,1) < P_omit_match
                            net(i,j) = 0;
                        end
                    end
                end
                Population(Pop_index).subnet{s_index}{layer} = net;
            end
        end
        [err, complexity, endnet] = DNevalnet(Population(Pop_index), Pop_str);
        FVal = [FVal; [err, complexity]];
        Population(Pop_index).err = err;
        Population(Pop_index).complexity = complexity;
        Population(Pop_index).endnet = endnet;
    end
end