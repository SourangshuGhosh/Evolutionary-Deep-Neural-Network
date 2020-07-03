function MovePrey(xpos, ypos, i)
%-------------------------------------------------
% Move Prey 
%-------------------------------------------------
%i      Prey to move
%xpos   new position for Prey
%ypos   -"-

global lattice no_x no_y

if i > 0
    if xpos > 1 && xpos < no_x+2 && ypos > 1 && ypos < no_y+2
        lattice(xpos,ypos) = i;
    elseif xpos == 1
        if ypos == 1
            lattice(no_x+1,no_y+1) = i;
        elseif ypos == no_y+2
            lattice(no_x+1,2) = i;
        else
            lattice(no_x+1,ypos) = i;
        end
    elseif xpos == no_x+2
        if ypos == 1
            lattice(2,no_y+1) = i;
        elseif ypos == no_y+2
            lattice(2,2) = i;
        else
            lattice(2,ypos) = i;
        end
    elseif ypos == 1
        lattice(xpos,no_y+1) = i;
    else
        lattice(xpos,2) = i;
    end
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