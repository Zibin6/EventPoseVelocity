function [R, T] = AbsLin(Line, evt)

line_num = length(Line);
event_num = size(evt{1,1}.events,2);
row_idx = 1;
for ii = 1:line_num
    m1 = Line{ii}.normal(1);
    m2 = Line{ii}.normal(2);
    m3 = Line{ii}.normal(3);
    v1 = Line{ii}.direct(1);
    v2 = Line{ii}.direct(2);
    v3 = Line{ii}.direct(3);

    events = evt{ii}.events;
    for j = 1:event_num

        h1 = events(1,j);
        h2 = events(2,j);

        dii = [v1;v2;v3];
        nii = evt{ii}.eventNormal;

        KK = [dii' .* nii(1), dii' .* nii(2), dii' .* nii(3)];
        AA = [h1*v1, h1*v2, h1*v3, h2*v1, h2*v2, h2*v3, v1, v2, v3, h1*m1, h1*m2, h1*m3, h2*m1, h2*m2, h2*m3, m1, m2, m3];
        AK(2*row_idx-1:2*row_idx,:) = [AA;zeros(1,9),KK];
        row_idx = row_idx + 1;
    end
end

AE = AK(:, 1:9); AR = AK(:, 10:18);
ARP = pinv(AR);
AA = (AR*ARP - eye(2*event_num*line_num))*AE;
[~, ~, V_E] = svd(AA);
e = V_E(:, end);
E_temp = [e(1:3),e(4:6),e(7:9)];
E = E_temp';

[U_R, ~, V_R] = svd(E);

W = [0, -1, 0; 1, 0, 0; 0, 0, 1];

Ra = U_R * W * V_R';
Rb = U_R * W' * V_R';

if det(Ra) < 0
    Ra = -Ra;
end
if det(Rb) < 0
    Rb = -Rb;
end

Rs(:,:,1) = Ra;
Rs(:,:,2) = Rb;

ts = zeros(3,2);
errs = zeros(1,2);
for ii = 1:2
    Rii = Rs(:,:,ii);
    C = [];b = [];
    for jj = 1:size(Line,2)
        Ld = Line{jj}.direct;
        Lm = Line{jj}.normal;
        Rd = Rii*Ld;
        Rm = Rii*Lm;
        evtjj = evt{jj}.events;
        for ll = 1:size(evtjj,2)
            ell = evtjj(:,ll);
            ckk = [Rd(2) - Rd(3)*ell(2), Rd(3)*ell(1) - Rd(1), Rd(1)*ell(2) - Rd(2)*ell(1)];
            C = [C;ckk];
            b = [b;-1*[Rm(3) + Rm(1)*ell(1) + Rm(2)*ell(2)]];
        end
    end
    tii = C\b;
    ts(:,ii) = tii;
    errs(ii) = norm(C*tii-b);
end
[~, idx] = min(errs);
R = Rs(:,:,idx);
T = ts(:,idx);
end