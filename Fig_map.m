clear, close all
addpath('Functions/')

load('Results/Map_RIS_MIMO_Nty8Ntz8Nry8Nrz8Nrisx32_perpendicular_fine_v2.mat')
load('Results/Map_RIS_MIMO_Nty8Ntz8Nry8Nrz8Nrisx32_perpendicular_not_opt_fine.mat', 'output_fix')


rate_ND_PSWF = output_ND_PSWF_WA./output_ND_num_WA;
rate_FOC_num = output_FOC_num_WA./output_ND_num_WA;
rate_FOC_PSWF = output_FOC_PSWF_WA./output_ND_num_WA;

figure(1),surf(X,Y, reshape(rate_ND_PSWF, size(X))),view(2),caxis([0.9 1]);
shading interp; colorbar
% exportgraphics(figure(1),'Map_ND_PSWF.pdf','ContentType','vector')
figure(2),surf(X,Y, reshape(rate_FOC_num, size(X))),view(2),caxis([0.9 1]);
shading interp; colorbar
% exportgraphics(figure(2),'Map_FOC.pdf','ContentType','vector')
figure(3),surf(X,Y, reshape(rate_FOC_PSWF, size(X))),view(2),caxis([0.9 1]);
shading interp; colorbar
% exportgraphics(figure(3),'Map_FOC_PSWF.pdf','ContentType','vector')
%% Map with subfig opt
fig1 = figure(10);
delta = lambda/2;

% Abs
subplot(3,1,1)
h=pcolor(X,Y, reshape(rate_ND_PSWF, size(X))); 
hold on
set(h,'edgecolor','none')
plot([r0_TX(1)+delta*1-delta/2*(Nt_y+1), r0_TX(1)+delta*Nt_y-delta/2*(Nt_y+1)],...
    r0_TX(2)*ones(1,2),'k','LineWidth',2);clim([0.9 1])
plot(r0_RIS(1)*ones(1,2)+0.01,...
    [r0_RIS(2)+delta*1-delta/2*(Nris_x+1), r0_RIS(2)+delta*Nris_x-delta/2*(Nris_x+1)],...
    'r','LineWidth',3);
% shading interp;
xlabel('x-axis (m)','Interpreter', 'latex', 'FontSize', 14)
ylabel('y-axis (m)','Interpreter', 'latex', 'FontSize', 14)
% title(sprintf('\\textbf{Mode %i}', 1),'Interpreter', 'latex', 'FontSize', 10)

subplot(3,1,2)
h=pcolor(X,Y, reshape(rate_FOC_num, size(X)));
set(h,'edgecolor','none')
hold on
xlabel('x-axis (m)','Interpreter', 'latex', 'FontSize', 14)
ylabel('y-axis (m)','Interpreter', 'latex', 'FontSize', 14)
plot([r0_TX(1)+delta*1-delta/2*(Nt_y+1), r0_TX(1)+delta*Nt_y-delta/2*(Nt_y+1)],...
    r0_TX(2)*ones(1,2),'k','LineWidth',2);
plot(r0_RIS(1)*ones(1,2)+0.01,...
    [r0_RIS(2)+delta*1-delta/2*(Nris_x+1), r0_RIS(2)+delta*Nris_x-delta/2*(Nris_x+1)],...
    'r','LineWidth',3);
clim([0.9 1])
% shading interp;

subplot(3,1,3)
h=pcolor(X,Y, reshape(rate_FOC_PSWF, size(X)));
set(h,'edgecolor','none')
hold on
clim([0.9 1])
xlabel('x-axis (m)','Interpreter', 'latex', 'FontSize', 14)
ylabel('y-axis (m)','Interpreter', 'latex', 'FontSize', 14)
% plot(r0_TX(1),r0_TX(2),'b^','MarkerFaceColor','k');
plot([r0_TX(1)+delta*1-delta/2*(Nt_y+1), r0_TX(1)+delta*Nt_y-delta/2*(Nt_y+1)],...
    r0_TX(2)*ones(1,2),'k','LineWidth',2);
plot(r0_RIS(1)*ones(1,2)+0.01,...
    [r0_RIS(2)+delta*1-delta/2*(Nris_x+1), r0_RIS(2)+delta*Nris_x-delta/2*(Nris_x+1)],...
    'r','LineWidth',3);
% shading interp;
% Abs color bar
h = axes(fig1,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
caxis(h,[0.9,1]);             % set colorbar limits

%% Map with subfig non opt
fig1 = figure(20);
delta = lambda/2;

rate_rnd = output_rnd./output_ND_num_WA;
rate_fix = output_fix./output_ND_num_WA;

% Abs
subplot(2,1,1)
h=pcolor(X,Y, reshape(rate_rnd, size(X))); 
hold on
set(h,'edgecolor','none')
plot([r0_TX(1)+delta*1-delta/2*(Nt_y+1), r0_TX(1)+delta*Nt_y-delta/2*(Nt_y+1)],...
    r0_TX(2)*ones(1,2),'k','LineWidth',2);clim([0 1])
plot(r0_RIS(1)*ones(1,2)+0.01,...
    [r0_RIS(2)+delta*1-delta/2*(Nris_x+1), r0_RIS(2)+delta*Nris_x-delta/2*(Nris_x+1)],...
    'r','LineWidth',3);
% shading interp;
xlabel('x-axis (m)','Interpreter', 'latex', 'FontSize', 14)
ylabel('y-axis (m)','Interpreter', 'latex', 'FontSize', 14)
% title(sprintf('\\textbf{Mode %i}', 1),'Interpreter', 'latex', 'FontSize', 10)

subplot(2,1,2)
h=pcolor(X,Y, reshape(rate_fix, size(X)));
set(h,'edgecolor','none')
hold on
xlabel('x-axis (m)','Interpreter', 'latex', 'FontSize', 14)
ylabel('y-axis (m)','Interpreter', 'latex', 'FontSize', 14)
plot([r0_TX(1)+delta*1-delta/2*(Nt_y+1), r0_TX(1)+delta*Nt_y-delta/2*(Nt_y+1)],...
    r0_TX(2)*ones(1,2),'k','LineWidth',2);
plot(0.01+r0_RIS(1)*ones(1,2)+0.01,...
    [r0_RIS(2)+delta*1-delta/2*(Nris_x+1), r0_RIS(2)+delta*Nris_x-delta/2*(Nris_x+1)],...
    'r','LineWidth',3);
clim([0 1])
% shading interp;

% Abs color bar
h = axes(fig1,'visible','off');
c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);  % attach colorbar to h
caxis(h,[0,1]);             % set colorbar limits


%% CCDF
[f_ND_PSWF,x_ND_PSWF] = ecdf(output_ND_PSWF_WA);
[f_ND_num,x_ND_num] = ecdf(output_ND_num_WA);
[f_FOC_num,x_FOC_num] = ecdf(output_FOC_num_WA);
[f_FOC_PSWF,x_FOC_PSWF] = ecdf(output_FOC_PSWF_WA);
[f_rnd_RIS,x_rnd_RIS] = ecdf(output_rnd);
[f_fix_RIS,x_fix_RIS] = ecdf(output_fix);


figure(4)
plot(x_ND_num,1-f_ND_num, 'k-x','MarkerIndices',1:300:length(f_ND_num)),hold on
plot(x_ND_PSWF,1-f_ND_PSWF,'b-o','MarkerIndices',1:300:length(f_ND_num))
plot(x_FOC_num,1-f_FOC_num,'r-^','MarkerIndices',1:300:length(f_ND_num))
plot(x_FOC_PSWF,1-f_FOC_PSWF, 'g-*','MarkerIndices',1:300:length(f_ND_num))
plot(x_rnd_RIS,1-f_rnd_RIS,'m-diamond','MarkerIndices',1:300:length(f_ND_num))
plot(x_fix_RIS,1-f_fix_RIS,'c-square','MarkerIndices',1:300:length(f_ND_num))
legend('ND-num','ND-PWSF','FOC-num','FOC-PSWF',...
    'Random RIS configuration','Uniform RIS configuration','Interpreter', 'latex', 'FontSize', 12)
xlabel('Achievable rate (bps/Hz)','Interpreter', 'latex', 'FontSize', 12)
ylabel('C-CDF in log-scale','Interpreter', 'latex', 'FontSize', 12)
% exportgraphics(figure(4),'Complementary_distribution.pdf','ContentType','vector')


figure(5)
semilogx(x_ND_num,1-f_ND_num, 'k-x','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50]),hold on
semilogx(x_ND_PSWF,1-f_ND_PSWF,'b-o','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50])
semilogx(x_FOC_num,1-f_FOC_num,'r-^','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50])
semilogx(x_FOC_PSWF,1-f_FOC_PSWF, 'g-*','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50])
semilogx(x_rnd_RIS,1-f_rnd_RIS,'m-diamond','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50])
semilogx(x_fix_RIS,1-f_fix_RIS,'c-square','MarkerIndices',[1:300:length(f_ND_num), length(f_ND_num)-100,length(f_ND_num)-50])
legend('ND-num','ND-PWSF','FOC-num','FOC-PSWF',...
    'Random RIS configuration','Uniform RIS configuration','Interpreter', 'latex', 'FontSize', 12)
xlabel('Achievable rate (bps/Hz)','Interpreter', 'latex', 'FontSize', 12)
ylabel('C-CDF','Interpreter', 'latex', 'FontSize', 12)
axis([0,290,0.5e-3,1])
grid on
% exportgraphics(figure(5),'Complementary_distribution_log_scale.pdf','ContentType','vector')
