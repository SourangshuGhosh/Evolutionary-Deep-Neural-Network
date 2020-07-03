function pop = tournament_select(pop)

global tour_size

trees = struct2cell(pop); name = fieldnames(pop);
fit = cell2mat(trees(strcmp(name, 'fval'),:));

play_tourl = [];
for i = 1:tour_size+1
    index = 1:length(pop);
    if i == tour_size+1
        index = play_tourl;
    end
    k = 1;
    while ~isempty(index)
        play_tour = [];
        if length(index) >= tour_size
            for j = 1:tour_size
                select_index = ceil(rand*length(index));
                play_tour = [play_tour index(select_index)];
                index(select_index) = [];
            end
            win = play_tour(fit(play_tour) == min(fit(play_tour)));
            tour(i).pop(k) = pop(win(1));
            k = k + 1;
        else
            play_tourl = [play_tourl index];
            index = [];
        end
    end
end

tour = struct2cell(tour(:));
if isempty(tour{1})
    pop = cell2mat(tour(length(tour)));
else
    pop = cell2mat(tour);
end
pop = pop';

end