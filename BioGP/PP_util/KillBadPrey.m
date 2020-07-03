function [Prey F1 F2] = KillBadPrey(Prey, F1, F2)
%-------------------------------------------------
% Kill bad Prey 
%-------------------------------------------------

global lattice F_bad no_x no_y

PreyIndex = (1:length(Prey(:,1)))';

a = find(F1 >= F_bad);
b = find(F2 >= F_bad);
rem = unique([a;b]);

for i = 1:length(rem)            
    [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == rem(i));
    lattice(xpos+1,ypos+1) = 0;
end

Prey(rem,:)=[];
PreyIndex(rem)=[];
F1(rem)=[];
F2(rem)=[];

if isempty(PreyIndex)
    disp('Eliminated all Prey as useless')
    pause
end
for i = 1:length(PreyIndex)            
    [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) == PreyIndex(i));
    lattice(xpos+1,ypos+1) = i;
end

MovePrey(0, 0, 0);

end