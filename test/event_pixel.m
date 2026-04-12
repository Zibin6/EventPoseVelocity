%% Monte Carlo Simulation: Analysis of the number of events
%  Author: Zibin Liu and Shunkun Liang
%  Date: 2026-04-11

clc; clear; close all;
addpath(genpath('..\'));

%% --- 1. Parameter Configurations ---
maxTime    = 0.1;           
samplNum   = 20;            
lineNum    = 20;           
trialNum   = 10000;         

% Experimental variables
eventNumLv = [3, 6, 10, 15, 20, 30];
noiseLv    = 1;                     
timeNoise  = 1e-3;                  

% Pre-allocate memory for speed (Crucial for performance!)
numLv = length(eventNumLv);
numNz = length(noiseLv);

[errR, errT, errR2, errT2] = deal(zeros(numLv, numNz, trialNum));
[errw, errv, errw2, errv2] = deal(zeros(numLv, numNz, trialNum));

%% --- 2. Main Simulation Loop ---
fprintf('Simulation started. Total levels: %d\n', numLv);

for ii = 1:numLv
    eventNum = eventNumLv(ii);
    
    for jj = 1:numNz
        noise = noiseLv(jj);
        
        % Progress display
        fprintf('Processing Event Number: %d, Noise: %.2f...\n', eventNum, noise);
        
        parfor kk = 1:trialNum  % Optional: Change 'for' to 'parfor' if you have Parallel Computing Toolbox
            %% 2.1 Data Generation
            [Line, evt, K, Rgt, Tgt, wgt, vgt] = generate_eventline(...
                lineNum, samplNum, eventNum, maxTime, noise, timeNoise);

            %% 2.2 Absolute Pose Estimation (Linear vs. Polynomial)
            % Linear Method
            [R_linear, T_linear] = AbsLin(Line, evt(:,1));
            errR(ii,jj,kk) = cal_rotation_err(R_linear, Rgt);
            errT(ii,jj,kk) = norm(Tgt - T_linear) / norm(Tgt) * 100;

            % Polynomial Method
            [R_poly, T_poly] = AbsPol(Line, evt(:,1), Rgt);
            errR2(ii,jj,kk) = cal_rotation_err(R_poly, Rgt);
            errT2(ii,jj,kk) = norm(Tgt - T_poly) / norm(Tgt) * 100;

            %% 2.3 Velocity Estimation (Linear vs. Optimization)
            % Linear Velocity
            [w_linear, v_linear] = VelLin(Line, evt, Rgt, Tgt);
            errw(ii,jj,kk) = norm(wgt - w_linear) / (norm(wgt) + norm(w_linear)) * 100;
            errv(ii,jj,kk) = norm(vgt - v_linear) / (norm(vgt) + norm(v_linear)) * 100;

            % Optimized Velocity
            [w_opti, v_opti] = VelOpt(Line, evt, Rgt, Tgt);
            errw2(ii,jj,kk) = norm(wgt - w_opti) / (norm(wgt) + norm(w_opti)) * 100;
            errv2(ii,jj,kk) = norm(vgt - v_opti) / (norm(vgt) + norm(v_opti)) * 100;
        end
    end
end

fprintf('Simulation completed successfully.\n');

%% --- 3. Result Visualization (Template) ---
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
