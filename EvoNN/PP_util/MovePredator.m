function Prey = MovePredator(Prey, xpos, ypos, i, j)
%-------------------------------------------------
% Move Predator
%-------------------------------------------------
global lattice no_x no_y
%xpos   position of prey
%ypos   -"-
%i      predator to move
%j      prey to kill
if j < inf
    %removing prey j from population
    Prey(j,:,:) = [];    
    [renamex,renamey] = find(lattice > j);
    for t = 1:length(renamex)
        lattice(renamex(t),renamey(t)) = lattice(renamex(t),renamey(t))-1;
    end
end

if xpos > 1 && xpos < no_x+2 && ypos > 1 && ypos < no_y+2
    lattice(xpos,ypos) = -i;
elseif xpos == 1
    if ypos == 1
        lattice(no_x+1,no_y+1) = -i;
    elseif ypos == no_y+2
        lattice(no_x+1,2) = -i;
    else
        lattice(no_x+1,ypos) = -i;
    end
elseif xpos == no_x+2
    if ypos == 1
        lattice(2,no_y+1) = -i;
    elseif ypos == no_y+2
        lattice(2,2) = -i;
    else
        lattice(2,ypos) = -i;
    end
elseif ypos == 1
    lattice(xpos,no_y+1) = -i;
else
    lattice(xpos,2) = -i;
end

lattice(2:no_y+1,1) = lattice(2:no_y+1,no_x+1);
lattice(2:no_y+1,no_x+2) = lattice(2:no_y+1,2);
lattice(1,2:no_x+1) = lattice(no_y+1,2:no_x+1);
lattice(no_y+2,2:no_x+1) = lattice(2,2:no_x+1);
lattice(1,1) = lattice(no_x+1,no_y+1);
lattice(1,no_x+2) = lattice(no_x+1,2);
lattice(no_y+2,1) = lattice(2,no_x+1);
lattice(no_y+2,no_x+2) = lattice(2,2);

end