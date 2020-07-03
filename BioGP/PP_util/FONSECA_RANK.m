function fonrank = FONSECA_RANK(pop, maxrank)

numobs = length(pop(:,1));
fonrank = ones(size(pop(:,1)));
for i = 1:numobs
    inda = find(pop(:,1) < pop(i,1) & pop(:,2) < pop(i,2));
    numa = length(inda);
    if numa < maxrank
        indb = find(pop(:,1) == pop(i,1) & pop(:,2) < pop(i,2));
        numb = length(indb);
        if numa + numb < maxrank
            indc = find(pop(:,1) < pop(i,1) & pop(:,2) == pop(i,2));
            numc = length(indc);
            if numa + numb + numc < maxrank
                fonrank(i) = numa + numb + numc;
            else
                fonrank(i) = maxrank;
            end
        else
            fonrank(i) = maxrank;
        end
    else
        fonrank(i) = maxrank;
    end
end

end