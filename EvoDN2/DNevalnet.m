function [errorVal,F2,Beta,InfoC] = DNevalnet(net, Pop_str)
%DNevalnet - Description
%
% Syntax: [fval,W,InfoC] = DNevalnet(net)
%
% Long description
    global parameters
    %global setno
    setno=1;
	out = parameters.dataset(setno).out;
	outmax = max(out);
	outmin = min(out);
	no_points = length(out);
    num_subnets = length(Pop_str);
    in_end = [];
    for subnet = 1:num_subnets
        in = parameters.dataset(setno).in;
        in = in(:,Pop_str{subnet}{1});
        no_layer = length(Pop_str{subnet}{2});
        complexity{subnet} = ones(1,size(in,2));
        %k = 0;
        snet = net.subnet{subnet};
        cnet = snet; %for complexity calculations
        for layer = 1:no_layer - 1  %-1?
            %cnet{layer} = logical(cnet{layer});
            in = in*snet{layer}(2:end,:);
            bias = snet{layer}(1,:);
            d=size(in);
            for i=1:d(1)
                for j=1:d(2)
                    in(i,j) = in(i,j)+bias(j);
                end
            end
            in = 1./(1+exp(-in));
            complexity{subnet} = complexity{subnet}*abs(cnet{layer}(2:end,:));
            %k = k+ length(find(net{layer}));
        end
        complexity{subnet} = sum(complexity{subnet});
        in_end = [in_end in];
    end
	%% LLSQ to solve in*Beta = out
	Beta = in_end\out;
	modelout = in_end*Beta;

	errorVal = sqrt(sum(((modelout-out)/(outmax-outmin)).^2)/no_points);
	F2 = sum(cell2mat(complexity));

	% n = no_points;
	% k = k+length(find(Beta));
	% rss = sum(sum((modelout-out).^2));
	% InfoC = 2*k+n*log(rss/n);
	% InfoC = InfoC+(2*k*(k+1)/(n-k-1));
end