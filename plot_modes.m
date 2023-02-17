load('/Users/carlos/ownCloud/Code/IRE_WSA/Results/Modes_RIS_MIMO_Nty8Ntz8Nry8Nrz8D6lt2lr2Nrisx50K100000_final.mat')
max_abs = 0;
for i =1:64
    mode = ObAtRIS_OptPhase_K100000{i};
    a(i) = sum(abs(mode(:)).^2);
end

[~,ind] = maxk(a,6);

for i =1:6
    mode = ObAtRIS_OptPhase_K100000{ind(i)};

    % Abs
    fig1 = figure(1);
    subplot(3,2,i)
    max_abs = max(max_abs, max(abs(mode)));
    surf(reshape(abs(mode),[50,50]))
    axis([1,50,1,50])
    view(2)
    xlabel('x-index RIS','Interpreter', 'latex', 'FontSize', 10)
    ylabel('y-index RIS','Interpreter', 'latex', 'FontSize', 10)
    title(sprintf('\\textbf{Mode %i}', i),'Interpreter', 'latex', 'FontSize', 10)


    % Phase
    fig2 = figure(2);
    subplot(3,2,i)
    surf(reshape(angle(mode)/pi,[50,50]))
    axis([1,50,1,50])
    view(2)
    xlabel('x-index RIS','Interpreter', 'latex', 'FontSize', 10)
    ylabel('y-index RIS','Interpreter', 'latex', 'FontSize', 10)
    title(sprintf('\\textbf{Mode %i}', i),'Interpreter', 'latex', 'FontSize', 10)
end

for i = 1:6
    figure(1)
    subplot(3,2,i)
    clim manual
    clim([0 max_abs]);

    figure(2)
    subplot(3,2,i)
    clim manual
    clim([-1 1]);
end

% Abs color bar
h = axes(fig1,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
caxis(h,[0,max_abs]);             % set colorbar limits

% Phase color bar
h = axes(fig2,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
caxis(h,[-1,1]);             % set colorbar limits
