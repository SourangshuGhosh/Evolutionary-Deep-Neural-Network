function ind_new = MUTE(ind, gen)

global LB UB

P_point = ind; %Probability that a gene is mutated
ind_new = P_point;
for k=1:length(ind(1,:))
    if rand < P_point
        ind_new(k) = ind_new(k) + 0.5*(UB(k) - LB(k))*randn;% / (0.1*gen);
        if ind_new(k) < LB(k)
            ind_new(k) = LB(k) + eps;
        elseif ind_new(k) > UB(k)
            ind_new(k) = UB(k) - eps;
        end
    end
end

end