function crowd = CROW_SORT(pop, front)
%CROW crowding distance calculations
%Created by Brijesh Kumar Giri
crowd = zeros(size(pop(:,1)));
maxf = zeros(size(pop(1,:))); minf = maxf;
for i = 1:length(pop(1,:))
    maxf(i) = max(pop(:,i));
    minf(i) = min(pop(:,i));
end

for i = 1:length(pop(1,:))
    pop(:,i) = pop(:,i) / (maxf(i) - minf(i));
end

for i = min(front):max(front)
    l = find(front == i);
    pop_temp = [l pop(l,:)];
    pop_temp = sortrows(pop_temp, 2);
    l = pop_temp(:,1); pop_temp(:,1) = [];
    
    crowd(l(1)) = inf;
    for j = 2:length(l)-1
        crowd(l(j)) = abs(pop(l(j+1),1) - pop(l(j-1),1)) + abs(pop(l(j+1),2) - pop(l(j-1),2));
    end
    crowd(l(length(l))) = inf;
end
    
end

