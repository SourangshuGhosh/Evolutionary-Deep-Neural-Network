function err_table = autooutput(Setslog,savedir)

global ploton saveon fig

ploton = 1;
saveon = 1;

fig = [figure(1) figure(2)];
scrsz = get(0, 'Screensize');
set(fig(1), 'OuterPosition', scrsz, 'Color', 'w'); clf
set(fig(2), 'OuterPosition', scrsz, 'Color', 'w'); clf

%Setslog = importdata(Setslog);
yl = Setslog.paraname{1, Setslog.out_index};
err_table = [];
for i = 1:Setslog.no_run
    pause(2)
    fprintf('Training data %d\n', i)
    err_table = [err_table; get_trained_data(Setslog, i, yl,savedir)];
end
disp(['ERROR TABLE for training of ' yl ' by BioGP:'])
disp(err_table)

end