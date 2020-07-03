function err_table = autooutput(Setslog,savedir)

global figure_han saveon ploton

saveon = 1;
ploton = 1;

figure_han = [figure(1) figure(2)];
scrsz = get(0, 'Screensize');
set(figure_han(1), 'OuterPosition', scrsz, 'Color', 'w'); clf
set(figure_han(2), 'OuterPosition', scrsz, 'Color', 'w'); clf

%Setslog = importdata(Setslog);
yl = Setslog.paraname{1,Setslog.out_index};
err_table = [];
for i = 1:Setslog.no_run
    pause(2)
    fprintf('Training data %d\n', i)
    err_table = [err_table; get_trained_data(Setslog, i, yl,savedir)];
end
disp(['ERROR TABLE for training of ' yl ' by EvoNN:'])
disp(err_table)

end