
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                                Plot Settings                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Line_width      = 2.5;
clr_red         = [0.90,0.15,0.35]; % Red
clr_blue        = [0.20,0.35,0.75]; % Blue
clr_green       = [0.00,0.75,0.25]; % Green
clr_orange      = [0.95,0.55,0.25]; % Orange
clr_purple      = [0.40,0.20,0.95]; % Purple
clr_gray        = [0.50,0.50,0.55]; % Gray
clr_black       = [0.00,0.00,0.00]; % Gray
FontWeight      = 'normal'; 
Font_Titles     = 10.5; 
Font_axis       = 11;
x_fig = 500; y_fig = 200; w_fig = 700; h_fig = 800;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%                          FE Channel Gain = 6X                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       FE Channel Output Noise [V^2/Hz] at Different Peaking Times       %
%                    Tp = 0.6us, 1.2us, 2.4us, 3.6us                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data      = csvread('outnoise_tps.csv',2);
freq      = data(:,1);
von_0u6s  = data(:,2);
von_1u2s  = data(:,3);
von_2u4s  = data(:,4);
von_3u6s  = data(:,5);

figure(1)
set(gcf,'Position',[x_fig y_fig w_fig h_fig])
subplot(4,1,1),semilogx(freq, von_0u6s, '-' ,'color',clr_red    , 'LineWidth',Line_width); hold on;
    title ('FE Output Noise'  ,'FontSize',Font_Titles,'FontWeight','bold');
    legend('tp = 0.6us'   ,'Location','northeast'); grid on;
    ylabel('Vnoise [V^2/Hz]'  ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 60e-12]);
subplot(4,1,2),semilogx(freq, von_1u2s, '-' ,'color',clr_green  , 'LineWidth',Line_width); hold on;
    legend('tp = 1.2us'   ,'Location','northeast'); grid on;
    ylabel('Vnoise [V^2/Hz]'  ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 60e-12]);
subplot(4,1,3),semilogx(freq, von_2u4s, '-' ,'color',clr_orange , 'LineWidth',Line_width); hold on;
    legend('tp = 2.4us'   ,'Location','northeast'); grid on;
    ylabel('Vnoise [V^2/Hz]'  ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 60e-12]);
subplot(4,1,4),semilogx(freq, von_3u6s, '-' ,'color',clr_blue   , 'LineWidth',Line_width); hold on;
    legend('tp = 3.6us'   ,'Location','northeast'); grid on;
    ylabel('Vnoise [V^2/Hz]'  ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 60e-12]);
    xlabel('Freq [Hz]'    ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%            FE Channel ENC [e-] at Different Peaking Times               %
%                      Tp = 0.6us, 1.2us, 2.4us, 3.6us                    %                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data      = csvread('enc_tps.csv',2);
freq      = data(:,1);
enc_0u6s  = data(:,2);
enc_1u2s  = data(:,3);
enc_2u4s  = data(:,4);
enc_3u6s  = data(:,5);

figure(2)
set(gcf,'Position',[x_fig y_fig w_fig h_fig])
subplot(4,1,1),semilogx(freq, enc_0u6s, '-' ,'color',clr_red    , 'LineWidth',Line_width); hold on;
    title ('FE TOTAL ENC' ,'FontSize',Font_Titles,'FontWeight','bold');
    legend('tp = 0.6us'   ,'Location','northwest'); grid on;
    ylabel('ENC [e-]'     ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 200]);
subplot(4,1,2),semilogx(freq, enc_1u2s, '-' ,'color',clr_green  , 'LineWidth',Line_width); hold on;
    legend('tp = 1.2us'   ,'Location','northwest'); grid on;
    ylabel('ENC [e-]'     ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 200]);
subplot(4,1,3),semilogx(freq, enc_2u4s, '-' ,'color',clr_orange , 'LineWidth',Line_width); hold on;
    legend('tp = 2.4us'   ,'Location','northwest'); grid on;
    ylabel('ENC [e-]'     ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 200]);
subplot(4,1,4),semilogx(freq, enc_3u6s, '-' ,'color',clr_blue   , 'LineWidth',Line_width); hold on;
    legend('tp = 3.6us'   ,'Location','northwest'); grid on;
    ylabel('ENC [e-]'     ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
    xlim([0 1e9]); ylim([0 200]);
    xlabel('Freq [Hz]'    ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)

% ENC Vs Tp
% -------------------------------------------------------------------------
figure(3)
enc_tot_0u6s = max(enc_0u6s);
enc_tot_1u2s = max(enc_1u2s);
enc_tot_2u4s = max(enc_2u4s);
enc_tot_3u6s = max(enc_3u6s);
enc_tot_all  = [enc_tot_0u6s enc_tot_1u2s enc_tot_2u4s enc_tot_3u6s];
peak_time    = [0.6 1.2 2.4 3.6];

plot(peak_time, enc_tot_all, '-' ,'color',clr_black, 'LineWidth',Line_width) , hold on;
legend('Cdet = 40pF' ,'Location','northeast'); grid on;
title ('Total ENC' ,'FontSize' ,Font_Titles*1.25,'FontWeight','bold');
ylabel('ENC [e-]'  ,'FontSize' ,Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
xlabel('Peaking Time [us]'   ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
xmin = min(peak_time); xmax = max(peak_time);
xlim([xmin xmax]); ylim([100 150]);
grid on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%            FE Channel ENC [e-] Vs Detector Capacitance (Cdet)           %
%                      Tp = 0.6us, 1.2us, 2.4us, 3.6us                    %                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_0u6s      = csvread('enc_cdet_0u6s.csv',1);
data_1u2s      = csvread('enc_cdet_1u2s.csv',1);
data_2u4s      = csvread('enc_cdet_2u4s.csv',1);
data_3u6s      = csvread('enc_cdet_3u6s.csv',1);

cdet           = data_0u6s(:,1);
cdet_scale     = cdet/1e-12;

enc_cdet_0u6s  = data_0u6s(:,2);
enc_cdet_1u2s  = data_1u2s(:,2);
enc_cdet_2u4s  = data_2u4s(:,2);
enc_cdet_3u6s  = data_3u6s(:,2);

figure(4)
plot(cdet_scale, enc_cdet_0u6s, '-' ,'color',clr_red  , 'LineWidth',Line_width) , hold on;
plot(cdet_scale, enc_cdet_1u2s, '-' ,'color',clr_green, 'LineWidth',Line_width) , hold on;
plot(cdet_scale, enc_cdet_2u4s, '-' ,'color',clr_orange, 'LineWidth',Line_width) , hold on;
plot(cdet_scale, enc_cdet_3u6s, '-' ,'color',clr_blue, 'LineWidth',Line_width) , hold on;
legend('tp = 0.6us','tp = 1.2us','tp = 2.4us','tp = 3.6us'   ,'Location','northwest'); grid on;
title ('Total ENC' ,'FontSize' ,Font_Titles*1.25,'FontWeight','bold');
ylabel('ENC [e-]'  ,'FontSize' ,Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
xlabel('Detector Capacitance [pF]'   ,'FontSize',Font_Titles,'FontAngle','normal','FontWeight',FontWeight)
xlim([0 250]); ylim([50 750]);
grid on;

