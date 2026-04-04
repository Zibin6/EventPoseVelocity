function [Line, evt, K, R0, t0, w, v] = generate_eventline(lineNum, samplNum, eventNum, maxTime, pixelNoise, timeNoise)
    
    if nargin < 5
        pixelNoise = 0;
    end
    if nargin < 6
        timeNoise = 0;
    end

    w = 640;
    h = 480;
    cx = 320;
    cy = 240;
    fx = 800;
    fy = 800;
    K = [fx, 0, cx; 0, fy, cy; 0, 0, 1];

    lineNum2 = 2 * lineNum;

    u_d = round(w * rand(1, lineNum2) - w/2);
    v_d = round(h * rand(1, lineNum2) - h/2);
    uv = [u_d + cx; v_d + cy; ones(1, lineNum2)];

    minDepth = 5;
    maxDepth = 10;
    depth = minDepth + (maxDepth - minDepth) * rand(1, lineNum2);

    R0 = randR();
    t0 = minDepth/2 * rand(3, 1);

    linePts3d = R0' * (depth .* (K\uv) - t0);
    linePts3d = [linePts3d; ones(1, lineNum2)];

    v = (rand(3, 1) - 0.5) * 10;
    w = (rand(3, 1) - 0.5) * 0.25;

    for ii = 1:lineNum
        Line{ii}.start = linePts3d(1:3, 2*ii-1);
        Line{ii}.end = linePts3d(1:3, 2*ii);
        Line{ii}.direct = Line{ii}.start - Line{ii}.end;
        Line{ii}.normal = cross(Line{ii}.start, Line{ii}.direct);
    end

    Ts = [0 ; maxTime * (rand(samplNum - 1,1))];
    Ts = sort(Ts);

    for ii = 1:samplNum
        dt = Ts(ii);
        deltR = expmap(w * dt);
        deltt = v * dt;
        Rii = deltR * R0;
        tii = deltR * t0 + deltt;

        dt = dt + timeNoise*randn;
        for jj = 1:lineNum
            eventNumjj = eventNum;
            linePoints = pointsOnLineSegment(Line{jj}.start, Line{jj}.end, eventNumjj);
            eventPoints = K * (Rii * linePoints + tii);
            eventPoints = eventPoints ./ eventPoints(3,:);
            eventPoints(1:2,:) = eventPoints(1:2,:) + pixelNoise * randn(2,eventNumjj);
            eventPoints = inv(K) * eventPoints;

            ex = eventPoints(1,:);
            ey = eventPoints(2,:);
            eline = polyfit(ex, ey, 1);
            eDir = [eline(2)/eline(1);eline(2);0];
            eDir = eDir./norm(eDir);
            lineDirs = repmat(eDir, 1, eventNumjj);

            eventNormal = cross(eventPoints, lineDirs);
            eventNormal = eventNormal ./ vecnorm(eventNormal,2);
            eventNormalMean1 = mean(eventNormal,2);

            evt{jj,ii} = struct('time', dt, 'events', eventPoints, 'eventNormal', eventNormalMean1);
        end
    end
end

function R = randR()
    theta_x = 2 * rand(1) - 1;
    theta_y = 2 * rand(1) - 1;
    theta_z = 2 * rand(1) - 1;
    r_x = [1, 0, 0; 0, cos(theta_x), -sin(theta_x); 0, sin(theta_x), cos(theta_x)];
    r_y = [cos(theta_y), 0, sin(theta_y); 0, 1, 0; -sin(theta_y), 0, cos(theta_y)];
    r_z = [cos(theta_z), -sin(theta_z), 0; sin(theta_z), cos(theta_z), 0; 0, 0, 1];
    R = r_z * r_y * r_x;
end

function R = expmap(omega)
    theta = norm(omega);
    if theta < 1e-10
        R = eye(3);
        return;
    end
    K = skew(omega);
    R = eye(3) + (sin(theta) / theta) * K + ((1 - cos(theta)) / (theta^2)) * K^2;
end

function K = skew(v)
    K = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function P = pointsOnLineSegment(P1, P2, num)
    t = rand(num,1);
    P = (1 - t').*P1 + t'.*P2;
end