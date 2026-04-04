function [w1, v1] = VelOpt(Line, evt, R, t, w0, v0)

if nargin < 5 || isempty(w0) || isempty(v0)
    x0 = zeros(6, 1);
else
    x0 = [w0; v0];
end


options = optimoptions('lsqnonlin', 'Display', 'none');
options.Algorithm = 'levenberg-marquardt';
options.FunctionTolerance = 1e-12;
options.StepTolerance = 1e-12;

costFun = @(x) loss_function(x, Line, evt, R, t);
[x_optim] = lsqnonlin(costFun, x0, [], [], options);
w1 = x_optim(1:3);
v1 = x_optim(4:6);
end

function cost = loss_function(x, Line, evt, R0, t0)
w = x(1:3); v = x(4:6);
samplNum = size(evt, 2);
lineNum = size(evt, 1);
cost = [];
for ii = 1:samplNum
    time = evt{1,ii}.time;
    dR = expmap(w * time);
    dt = v * time;
    Rii = dR * R0;
    tii = dR * t0 + dt;
    for jj = 1:lineNum
        eventPoints = evt{jj,ii}.events;
        err = eventPoints'*(Rii * Line{jj}.normal + skew(tii) * Rii * Line{jj}.direct);
        cost = [cost;err];
    end
end
end

