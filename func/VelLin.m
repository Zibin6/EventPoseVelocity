function [w,v] = VelLin(Line, evt, R, T)
A = zeros(numel(Line) * size(evt, 2), 3);
B = zeros(numel(Line) * size(evt, 2), 1);
idx = 1;

for ii = 1:numel(Line)
    d = Line{ii}.direct;
    for jj = 1:size(evt, 2)
        t = evt{ii,jj}.time;
        nj = evt{ii,jj}.eventNormal;
        dc = R * d;
        K = kron(dc, nj)';
        A(idx, :) = [K(6)*t-K(8)*t, K(7)*t-K(3)*t, K(2)*t-K(4)*t];
        B(idx) = -K(1)-K(5)-K(9);

        idx = idx + 1;
    end
end

w = A \ B;

A = zeros(numel(Line) * size(evt, 2), 3);
B = zeros(numel(Line) * size(evt, 2), 1);
idx = 1;

for ii = 1:numel(Line)
    m = Line{ii}.normal;
    d = Line{ii}.direct;

    for jj = 1:size(evt, 2)
        evtjj = evt{ii,jj}.events;
        t = evt{ii,jj}.time;
        for kk = 1:size(evtjj, 2)
            h = evtjj(:,kk);
            Rj = expmap(w*t);
            mm = Rj * R * m;
            a = Rj * T;
            b = Rj * R * d;

            A(idx, :) = [b(2)*t*h(3) - b(3)*t*h(2), b(3)*t*h(1) - b(1)*t*h(3), b(1)*t*h(2) - b(2)*t*h(1)];
            B(idx) = -(h' * (mm + cross(a, b)));

            idx = idx + 1;
        end
    end
end

v = A \ B;

end

function err_w = Loss_w(x, Line, evt, R)
w = x;
idx = 1;
for ii = 1:numel(Line)
    d = Line{ii}.direct;
    for jj = 1:size(evt, 2)
        t = evt{ii,jj}.time;
        nj = evt{ii,jj}.eventNormal;
        dc = R * d;
        Rj = expmap(w*t);
        err_w(idx) = nj'*Rj*dc;
        idx = idx + 1;
    end
end
end
