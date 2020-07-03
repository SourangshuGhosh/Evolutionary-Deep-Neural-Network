%function plottree('Setslog2-3-7_GA_Data.mat', 2, 4)
function plottree(Setslog, setno, tree_index)
%PLOTTREE
% Summary of this function goes here
%   Detailed explanation goes here
clc
global max_arity

Setslog = importdata(Setslog);
ind = Setslog.dataset(setno).pareto.P(tree_index);

F = Setslog.F;
max_arity = max(F.sets);

disp(['Bias: ' num2str(ind.bias)])
for i = 1:ind.no_roots
    plotroot(ind.root(i), i, ['root ' num2str(i)]);
    disp(['Root' num2str(i) ' Weight: ' num2str(ind.root(i).w) '; Error Reduction Ratio: ' num2str(ind.root(i).err_red)])
end

end

function plotroot(ind, index, fig_name)
%PLOTROOT
% Summary of this function goes here
%   Detailed explanation goes here

global max_arity

figure(index); clf
set(gcf, 'Name', fig_name)
axes
set(gca, 'visible', 'off', 'box', 'off')
set(gcf, 'color', 'w')

r = ind.tree_index;
d_max = 7;

d = ind.depth;

if d <= d_max
    d_max = d; if d_max <= 1, d_max = 2; end
end
xlim([0 max_arity^(d_max-1)]); ylim([1 d_max]); hold on

dist = cell(d_max,1); diff = 0.025;
for i = d_max:-1:1
    if i == d_max
        dist{i} = (0:1:max_arity^(i-1)-1);
    else
        dist{i} = [];
        k = 1;
        chd = dist{i+1};
        m = 0;
        while k <= length(chd)
            dist{i} = [dist{i} (chd(k) + chd(k+1))/2]; m = m+1;
            par = dist{i}; par = par(m);
            ypar = ((d_max-i+1)-diff*d_max);
            ychd = ((d_max-(i+1)+1)+diff*d_max);
            if ~isempty(find(r == max_arity^(i)+k-1,1))
                plot([par chd(k)],[ypar ychd],'k','Linewidth',1.5)
            end
            if ~isempty(find(r == max_arity^(i)+k,1))
                plot([par chd(k+1)],[ypar ychd],'k','Linewidth',1.5)
            end
            k = k+2;
        end
    end
end

halign = 'center'; lstyle = '-'; ecolor = 'red';
k = 1;
for i = 1:d_max
    par = dist{i};
    for j = max_arity^(i-1):max_arity^i-1
        txt = find(r == j);
        x = par(j-max_arity^(i-1)+1);
        y = (d_max-i+1);
        if i == d_max && length(txt) > 1
            sname = [num2str(index) num2str(k)];
            text(x, y, ['s' sname], 'HorizontalAlignment',halign, 'LineStyle',lstyle,  'EdgeColor',ecolor);
            subtree = snatch_tree(ind, txt);
            subtree = rename_ind(subtree, d_max, '1');
            plotroot(subtree, str2num(sname), ['root ' sname])
            k = k+1;
            figure(index)
        elseif ~isempty(txt)
            txt = txt(1);
            if strcmp(ind.tree{txt}, 'add(')
                txt = '+';
            elseif strcmp(ind.tree{txt}, 'mul(')
                txt = '*';
            elseif strcmp(ind.tree{txt}, 'sub(')
                txt = '-';
            elseif strcmp(ind.tree{txt}, 'rdiv(')
                txt = '/';
            elseif strcmp(ind.tree{txt}, 'exp(')
                txt = 'e^';
            elseif strcmp(ind.tree{txt}, 'sin(')
                txt = 'sin';
            elseif strcmp(ind.tree{txt}, 'rlog(')
                txt = 'log';
            elseif strcmp(ind.tree{txt}, 'rsqrt(')
                txt = 'sqrt';
            elseif strcmp(ind.tree{txt}, 'powsq(')
                txt = '^2';
            else
                txt = ind.tree{txt};
            end
            text(x, y, txt, 'HorizontalAlignment',halign, 'LineStyle',lstyle,  'EdgeColor',ecolor);
        end
    end
end

end