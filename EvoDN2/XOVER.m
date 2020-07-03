function pop_new = XOVER(pop)

global LB UB

sbxn = 2;  %SBX parameter
P_point = 0.2;  %fraction of population to undergo point Xover
novar = length(pop(1,:));
dox = zeros(floor(length(pop(:,1))/2), 2);
k = 1;
pop_new = pop;

index = 1:length(pop(:,1));
while ~isempty(index)
    if length(index) >= 2
        for j = 1:2
            select_index = ceil(rand*length(index));
            dox(k,j) = index(select_index);
            index(select_index) = [];
        end
    else
        pop_new(2*k-1,:) = pop(index,:);
        index = [];
    end
    k = k + 1;
end

for i = 1:length(dox(:,1))
    p1 = dox(i,1); p2 = dox(i,2);
    if rand < P_point
        xovrpnt = round(rand*(novar-2))+1;
        pop_new(2*i-1,:) = [pop(p1,1:xovrpnt) pop(p2,xovrpnt+1:novar)];
        pop_new(2*i,:) = [pop(p2,1:xovrpnt) pop(p1,xovrpnt+1:novar)];
    else
        for t=1:novar
            u = rand;
            if u <= 0.5
                beta = power(2*u, 1/(sbxn+1));
            else
                beta = power(1/(2*(1-u)), 1/(sbxn+1));
            end
            pop_new(2*i-1,t) = 0.5*((1+beta)*pop(p1,t) + (1-beta)*pop(p2,t));
            if pop_new(2*i-1,t) < LB(t)
                pop_new(2*i-1,t) = LB(t) + eps;
            elseif pop_new(2*i-1,t) > UB(t)
                pop_new(2*i-1,t) = UB(t) - eps;
            end
            pop_new(2*i,t) = 0.5*((1 - beta)*pop(p1,t) + (1+beta)*pop(p2,t));
            if pop_new(2*i,t) < LB(t)
                pop_new(2*i,t) = LB(t) + eps;
            elseif pop_new(2*i,t) > UB(t)
                pop_new(2*i,t) = UB(t) - eps;
            end
        end
    end
end

end