function [hL] = lin_eu(hlin,d,F1,FD)
% mod eul

% dh =d.A*[hlin(1)-d.h1_0,hlin(2)-d.h2_0]'+d.B*(F1-d.F1_0)+d.E*(FD-d.FD);
% dh1Lk1=dh(1,1);
% dh1Lk2=dh(2,1);
% hL= hlin + 0.5*d.Tp*(d.F0+d.A*[(hlin(1)+0.5*dh1Lk1)-d.h1_0,(hlin(2)+0.5*dh1Lk2)-d.h2_0]'+d.B*(F1-d.F1_0)+d.E*(FD-d.FD));
dh = d.Tp*(d.A*[hlin(1)-d.h1_0,hlin(2)-d.h2_0]'+d.B*(F1-d.F1_0)+d.E*(FD-d.FD));
dh05 = 0.5*d.Tp*(d.A*[(hlin(1)+0.5*dh(1))-d.h1_0,(hlin(2)+0.5*dh(2))-d.h2_0]'+d.B*(F1-d.F1_0)+d.E*(FD-d.FD));
hL = hlin + dh05;
end