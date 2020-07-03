function PlotLattice
%-------------------------------------------------
% Plot Lattice
%-------------------------------------------------
global lattice no_x no_y figure_handle Setslog

if Setslog.ploton
    h = figure_handle(2);
    set(0,'CurrentFigure',h)
    [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) > 0);
    plot(xpos,ypos,'b*');
    hold on
    [xpos,ypos] = find(lattice(2:no_x+1,2:no_y+1) < 0);
    plot(xpos,ypos,'r*');
    xlim([0 no_x+1]); ylim([0 no_y+1]);
    axis off
    set(h, 'Color', 'w')
    pause(0.1)
    hold off
end

end