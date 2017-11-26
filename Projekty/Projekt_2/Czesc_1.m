clear;clc;close all;
clf;hold all;grid on
% Stale
A2 = 310;
C1 = 0.75;
alfa1 = 16;
alfa2 = 19;
tau=150;

% punkty poczatkowe
F1_0 = 51;
FD = 7;
h2_0 = ((F1_0+FD)/alfa2)^2;     
h1_0 = ((F1_0+FD)/alfa1)^2;     

%Okres probkowania.
Tp = 1;

%% struktura danych w celu ułatwienia przekazywania do funkcji
dane.A2 = A2;
dane.C1 = C1;
dane.alfa1 = alfa1;
dane.alfa2 = alfa2;
dane.FD = FD;
dane.F1_0 = F1_0;
dane.Tp = Tp;
dane.h1_0=h1_0;
dane.h2_0=h2_0;

%% wektory danych
czas_symulacji = 15000;
l_iter = floor(czas_symulacji/Tp);
% iter_opoznienia = opoznienie/Tp;

% wektory dla modelu nieliniowego
h2_eu = h2_0 * ones(l_iter, 1);
h1_eu = h1_0*ones(l_iter, 1);
h1=h1_0*ones(l_iter, 1);
h2=h2_0 * ones(l_iter, 1);

% wektory dla modelu liniowego
h2lin_Eu = h2_0 * ones(l_iter, 1);
h1lin_Eu =h1_0*ones(l_iter, 1);
h1lin=h1_0*ones(l_iter, 1)';
h2lin=h2_0*ones(l_iter, 1)';
hlin = [h1lin;h2lin];
% hlin = zeros(2,l_iter+1);

% wektory dla sterowania (wielkości wejściowej, w celu realizacji
% opóźnienia)
F1 = zeros(l_iter+1,1);
F1_in = zeros(l_iter+1,1);
FDin=zeros(l_iter+1,1);
FD_in=zeros(l_iter+1,1);

A11 = -(4*F1_0+4*FD-3*alfa1*sqrt(h1_0))/(6*C1*h1_0^3);
A12 = 0;
A21 = alfa1/(2*A2*sqrt(h1_0));
A22 = -alfa2/(2*A2*sqrt(h2_0));

A = [A11,A12;A21,A22];
B = [1/(3*C1*h1_0^2),0]';
E = [1/(3*C1*h1_0^2),0]';
C = [1 0];
dane.A = A;
dane.B = B;
dane.C = C;
dane.E = E;

%% Symulacja obiektów

F1wy=zeros(1,7);
FDwy=zeros(1,7);
tabdF1={-30, -20, -10, 0, 10, 20, 30};
tabdFD={-6, -4, -2, 0, 2, 4, 6};
plotcolor={[1 0.9882 0.5333],[0.3608 0.9961 0.3804],[0.9961 0.9804 0.2784],[0.0941 0.9961 0.1176],[0.9804 0.9529 0.0431],[0.0039 0.7725 0.0235],[0.9882 0.0118 0.0118],[0.2235 0.0588 1],[0.9882 0.2941 0.2941],[0.3882 0.2588 1],[0.9922 0.5569 0.5569],[0.5373 0.4431 1],[0.92 0.4 0.6],[0.7 0.443 1],[0.8 0.3 0.7],[0.4 0.2 0.6]};
czas=0:Tp:czas_symulacji;
for m = 1 : 7
dF1 =0;%tabdF1{m};
dFD =tabdFD{m};


for k = 1:l_iter
    
    F1_in(k) = F1_0;
    FD_in(k)=FD;
    
    % skok sterowania w 2000 kroku
    if k>2000
        F1_in(k) = F1_0+dF1;
        FD_in(k)=FD+dFD;
    end
    
    %--Realizacja opznienia
    if(k - tau/Tp > 0)
        F1(k) = F1_in(k - floor(tau/Tp));
%         FD_in(k) = FD_in(k - floor(tau/Tp));
    else
        FD_in(k)=FD;
        F1(k) = F1_0;
    end
    
    %zmodyfikowany Euler nieliniowy
    [h] = nlin_eu(h1(k),h2(k),dane,F1(k),FD_in(k));
    h1(k+1) = h(1);
    h2(k+1) = h(2);
        
    %zlinearyzowany
%   hlin(:,k+1) = hlin(:,k) + Tp*(A*[hlin(1,k)-h1_0,hlin(2,k)-h2_0]'+B*(F1(k)-F1_0)+E*(FD-FD));
    [h] = lin_eu([h1lin(k),h2lin(k)]',dane,F1(k),FD_in(k));
    h1lin(k+1) = h(1);
    h2lin(k+1) = h(2);

    F1wy(m)=F1_in(k);
    FDwy(m)=FD_in(k);
end
 stairs(czas,h2,'Color',plotcolor{2*m},'LineWidth',1.5);
    grid on; hold on
    stairs(czas,h2lin,'Color',plotcolor{2*m+1},'LineWidth',1.5);
    ylabel('h_2[cm]'); xlabel('czas[s]');
    
end

%% rysowanie wykresow
legend(sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(1),FDwy(1)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(1),FDwy(1)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(2),FDwy(2)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(2),FDwy(2)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(3),FDwy(3)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(3),FDwy(3)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(4),FDwy(4)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(4),FDwy(4)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(5),FDwy(5)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(5),FDwy(5)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(6),FDwy(6)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(6),FDwy(6)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(7),FDwy(7)),sprintf('Model liniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(7),FDwy(7)),'Location','eastoutside');

figure
stairs(czas,F1_in);
grid on
ylabel('F1[cm^3/s]'); xlabel('czas[s]');






