function [rank front] = NONDOM_SORT(pop)

fonrank = zeros(size(pop(:,1)));
S = cell(size(pop(:,1)));
front = zeros(size(pop(:,1)));

for i = 1:length(pop(:,1))
    inda = find(pop(:,1) > pop(i,1) & pop(:,2) > pop(i,2));
    indb = find(pop(:,1) == pop(i,1) & pop(:,2) > pop(i,2));
    indc = find(pop(:,1) > pop(i,1) & pop(:,2) == pop(i,2));
    S(i) = {[inda; indb; indc]};
    
    inda = find(pop(:,1) < pop(i,1) & pop(:,2) < pop(i,2));
    indb = find(pop(:,1) == pop(i,1) & pop(:,2) < pop(i,2));
    indc = find(pop(:,1) < pop(i,1) & pop(:,2) == pop(i,2));
    fonrank(i) = length(inda) + length(indb) + length(indc);
end

rank = fonrank;
frontc = 1;
P = find(fonrank == 0);
fonrank(P) = -1;
front(P) = frontc;

while ~isempty(P)
    for i = 1:length(P)
        pop_dom = S{P(i)};
        for j = 1:length(pop_dom)
            fonrank(pop_dom(j)) = fonrank(pop_dom(j))- 1;
        end
    end
    frontc = frontc + 1;
    P = find(fonrank == 0);
    fonrank(P) = -1;
    front(P) = frontc;
end

end