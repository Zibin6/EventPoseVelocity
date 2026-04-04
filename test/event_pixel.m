clc;
clear all;
close all;
addpath(genpath('..\'));

maxTime = 0.1;
samplNum = 20;
eventNumLv = [3,6,10,15,20,30];
lineNum = 20;
noiseLv = 1;
timeNoise= 1 *10^-3;
trialNum = 10000;

errR = zeros(length(eventNumLv), length(noiseLv), trialNum);
errT = zeros(length(eventNumLv), length(noiseLv), trialNum);
errw = zeros(length(eventNumLv), length(noiseLv), trialNum);
errv = zeros(length(eventNumLv), length(noiseLv), trialNum);


for ii = 1:length(eventNumLv)
    eventNum = eventNumLv(ii);
    for jj = 1:length(noiseLv)
        noise = noiseLv(jj);
        for kk = 1:trialNum
            [Line, evt, K, Rgt, Tgt, wgt, vgt] = generate_eventline(lineNum, samplNum, eventNum, maxTime, noise, timeNoise);

            [R_linear, T_linear] = AbsLin(Line,evt(:,1));
            [R_poly, T_poly] = AbsPol(Line, evt(:,1), Rgt);

            errR(ii,jj,kk) = cal_rotation_err(R_linear, Rgt);
            errT(ii,jj,kk) = norm(Tgt - T_linear) / norm(Tgt) * 100;

            errR2(ii,jj,kk) = cal_rotation_err(R_poly, Rgt);
            errT2(ii,jj,kk) = norm(Tgt - T_poly) / norm(Tgt) * 100;

            [w_linear, v_linear] = VelLin(Line, evt, Rgt, Tgt);
            [w_opti, v_opti] = VelOpt(Line, evt, Rgt, Tgt);

            errw(ii,jj,kk) = norm(wgt - w_linear) / (norm(wgt) + norm(w_linear)) * 100;
            errv(ii,jj,kk) = norm(vgt - v_linear) / (norm(vgt) + norm(v_linear)) * 100;

            errw2(ii,jj,kk) = norm(wgt - w_opti) / (norm(wgt) + norm(w_opti)) * 100;
            errv2(ii,jj,kk) = norm(vgt - v_opti) / (norm(vgt) + norm(v_opti)) * 100;
        end
    end
end


fillcolor1(1,:) = [88, 159, 243]./255;
fillcolor1(2,:) = [249, 65, 65]./255;  
fillcolor2(1,:) = [55, 171, 120]./255; 
fillcolor2(2,:) = [243, 177, 105]./255;


n = length(eventNumLv);
interval_size = 8;      
box_width = 1.2;         
position_1 = [0:interval_size:interval_size*(n-1)] - 0.6;  
position_2 = [0:interval_size:interval_size*(n-1)] + 0.6;  


FontSize_Label = 14;    
FontSize_Title = 16;   
FontSize_Tick = 14;    
LineWidth_Box = 1.0;   

figure('Position', [100, 100, 1600, 400]);

subplot_position = [0.05, 0.15, 0.18, 0.75;  
                   0.3, 0.15, 0.18, 0.75;  
                   0.55, 0.15, 0.18, 0.75; 
                   0.8, 0.15, 0.18, 0.75];  

subplot('Position', subplot_position(1,:));
boxplot(squeeze(errR(:,1,:))', 'positions', position_1, 'width', box_width, 'symbol', '');
hold on;
boxplot(squeeze(errR2(:,1,:))', 'positions', position_2, 'width', box_width, 'symbol', '');

for i = 1:2
    groupObj = findobj('Tag', 'boxplot');
    allLineObj = groupObj(3-i).Children;
    outLineObj = findobj(allLineObj, 'Tag', 'Outliers');
    set(outLineObj, 'MarkerEdgeColor', fillcolor1(i,:))
    boxLineObj = findobj(allLineObj, 'Tag', 'Box');
    medLineObj = findobj(allLineObj, 'Tag', 'Median');

    for j = 1:n
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(1),boxLineObj(j).YData(1),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor1(i,:), 'LineWidth', LineWidth_Box);
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(2),boxLineObj(j).YData(2),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor1(i,:), 'LineWidth', LineWidth_Box);
    end
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize_Tick, 'YScale', 'log', ...
    'XTick', [0:interval_size:interval_size*(n-1)], 'LineWidth', 1, 'FontWeight', 'bold', 'YLim', [1e-2 1e1]);
xticklabels(eventNumLv);
xlabel('Number of Events', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');
ylabel('Rotation Error (deg)', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');

subplot('Position', subplot_position(2,:));
boxplot(squeeze(errT(:,1,:))', 'positions', position_1, 'width', box_width, 'symbol', '');
hold on;
boxplot(squeeze(errT2(:,1,:))', 'positions', position_2, 'width', box_width, 'symbol', '');

for i = 1:2
    groupObj = findobj('Tag', 'boxplot');
    allLineObj = groupObj(3-i).Children;
    outLineObj = findobj(allLineObj, 'Tag', 'Outliers');
    set(outLineObj, 'MarkerEdgeColor', fillcolor1(i,:))
    boxLineObj = findobj(allLineObj, 'Tag', 'Box');
    medLineObj = findobj(allLineObj, 'Tag', 'Median');

    for j = 1:n
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(1),boxLineObj(j).YData(1),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor1(i,:), 'LineWidth', LineWidth_Box);
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(2),boxLineObj(j).YData(2),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor1(i,:), 'LineWidth', LineWidth_Box);
    end
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize_Tick, 'YScale', 'log', ...
    'XTick', [0:interval_size:interval_size*(n-1)], 'LineWidth', 1, 'FontWeight', 'bold', 'YLim', [1e-2 1e2]);
xticklabels(eventNumLv);
xlabel('Number of Events', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');
ylabel('Translation Error (%)', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');

subplot('Position', subplot_position(3,:));
boxplot(squeeze(errw(:,1,:))', 'positions', position_1, 'width', box_width, 'symbol', '');
hold on;
boxplot(squeeze(errw2(:,1,:))', 'positions', position_2, 'width', box_width, 'symbol', '');

for i = 1:2
    groupObj = findobj('Tag', 'boxplot');
    allLineObj = groupObj(3-i).Children;
    outLineObj = findobj(allLineObj, 'Tag', 'Outliers');
    set(outLineObj, 'MarkerEdgeColor', fillcolor2(i,:))
    boxLineObj = findobj(allLineObj, 'Tag', 'Box');
    medLineObj = findobj(allLineObj, 'Tag', 'Median');

    for j = 1:n
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(1),boxLineObj(j).YData(1),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor2(i,:), 'LineWidth', LineWidth_Box);
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(2),boxLineObj(j).YData(2),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor2(i,:), 'LineWidth', LineWidth_Box);
    end
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize_Tick, 'YScale', 'log', ...
    'XTick', [0:interval_size:interval_size*(n-1)], 'LineWidth', 1, 'FontWeight', 'bold', 'YLim', [1e-1 1e2]);
xticklabels(eventNumLv);
xlabel('Number of Events', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');
ylabel('Angular Velocity Error (%)', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');

subplot('Position', subplot_position(4,:));
boxplot(squeeze(errv(:,1,:))', 'positions', position_1, 'width', box_width, 'symbol', '');
hold on;
boxplot(squeeze(errv2(:,1,:))', 'positions', position_2, 'width', box_width, 'symbol', '');

for i = 1:2
    groupObj = findobj('Tag', 'boxplot');
    allLineObj = groupObj(3-i).Children;
    outLineObj = findobj(allLineObj, 'Tag', 'Outliers');
    set(outLineObj, 'MarkerEdgeColor', fillcolor2(i,:))
    boxLineObj = findobj(allLineObj, 'Tag', 'Box');
    medLineObj = findobj(allLineObj, 'Tag', 'Median');

    for j = 1:n
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(1),boxLineObj(j).YData(1),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor2(i,:), 'LineWidth', LineWidth_Box);
        patch(boxLineObj(j).XData, [medLineObj(j).YData(1),boxLineObj(j).YData(2),boxLineObj(j).YData(2),medLineObj(j).YData(1),medLineObj(j).YData(1)], ...
            fillcolor2(i,:), 'LineWidth', LineWidth_Box);
    end
end

set(gca, 'FontName', 'Times New Roman', 'FontSize', FontSize_Tick, 'YScale', 'log', ...
    'XTick', [0:interval_size:interval_size*(n-1)], 'LineWidth', 1, 'FontWeight', 'bold', 'YLim', [1e-1 1e2]);
xticklabels(eventNumLv);
xlabel('Number of Events', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');
ylabel('Linear Velocity Error (%)', 'FontName', 'Times New Roman', 'FontSize', FontSize_Label, 'FontWeight', 'bold');
