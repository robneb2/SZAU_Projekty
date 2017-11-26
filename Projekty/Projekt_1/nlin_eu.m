function [h] = nlin_eu(h1,h2,d,F1,FD_in)
% mod Euler

dh1_eu=d.Tp*(F1+FD_in-d.alfa1*sqrt(h1))/(3*d.C1*h1^2);
h(1)= h1 + (d.Tp/2)*((F1+FD_in-d.alfa1*sqrt(h1+0.5*dh1_eu))/(3*d.C1*(h1+0.5*dh1_eu)^2));

dh2_eu=d.Tp*(d.alfa1*sqrt(h1)-d.alfa2*sqrt(h2))/d.A2;
h(2)= h2 + (d.Tp/2)*(((d.alfa1*sqrt(h1+0.5*dh1_eu)-d.alfa2*sqrt(h2+0.5*dh2_eu))/d.A2));
    

% h(1)= h1 + d.Tp*(F1+d.FD-d.alfa1*sqrt(h1))/(3*d.C1*h1^2);
% h(2)= h2 + d.Tp*(d.alfa1*sqrt(h1)-d.alfa2*sqrt(h2))/d.A2;
%  



end