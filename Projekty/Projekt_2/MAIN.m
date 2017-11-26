clear;clc;close all;

%% Generuje skok jednostkowy i oblicza parametry DMC
% Znienić procedurę liczenia po zaimplementowaniu (nlin na lin)

%% Parametry regulatora
Nu = 170;
N = 240;
D = N;
lambda = 700;

param.Nu = Nu;
param.N = N;
param.D = D;
param.lambda = lambda;

% Stale
A2 = 310;
C1 = 0.75;
alfa1 = 16;
alfa2 = 19;
tau = 150;

% punkty poczatkowe
F1_0 = 51;
FD = 7;
h2_0 = ((F1_0+FD)/alfa2)^2;     
h1_0 = ((F1_0+FD)/alfa1)^2;     

%Okres probkowania.
Tp = 1;

%Model
A11 = -(4*F1_0+4*FD-3*alfa1*sqrt(h1_0))/(6*C1*h1_0^3);
A12 = 0;
A21 = alfa1/(2*A2*sqrt(h1_0));
A22 = -alfa2/(2*A2*sqrt(h2_0));

A = [A11,A12;A21,A22];
B = [1/(3*C1*h1_0^2),0]';
E = [1/(3*C1*h1_0^2),0]';
C = [1 0];

% parametry symulacji
czas_symulacji = 20000;
l_iter = floor(czas_symulacji/Tp);

%% struktura danych w celu ułatwienia przekazywania do funkcji
dane.A2 = A2;
dane.C1 = C1;
dane.alfa1 = alfa1;
dane.alfa2 = alfa2;
dane.tau = tau;
dane.FD = FD;
dane.F1_0 = F1_0;
dane.Tp = Tp;
dane.h1_0 = h1_0;
dane.h2_0 = h2_0;
dane.A = A;
dane.B = B;
dane.C = C;
dane.E = E;
dane.l_iter = l_iter;

% punkty poczatkowe
h2_0 = ((F1_0+FD)/alfa2)^2;     %Wyznaczono z definicji punktu stałego x*
h1_0 = ((F1_0+FD)/alfa1)^2;     %Wyznaczono z definicji punktu stałego x*

[ke,Ku] = DMC_param(param,dane);

%%	Odpowiedz ukladu symulowanego z regulacji DMC
Y = zeros(1,l_iter);                            %Wektor odpowiedzi uk³adu, czyszczony dla nowych wyznaczen
% YzadDMC = ((h2_0+2)*ones(1,N))';                         %Wektor wartoci zadanych                  -- zadaæ
Yhat = zeros(1,N)';                              %Wektor predykowany odpowiedzi uk³adu      -- odpowied wyznaczona

% Yzad = [h2_0*ones(1,500), YzadDMC(1)*ones(1,l_iter-499)];
kroczek = 2000;
Yzad = [(h2_0)*ones(1,kroczek),(h2_0+1)*ones(1,l_iter+kroczek-1)];
% Yzad = (h2_0)*ones(1,l_iter+1);
%% wektory danych
l_iter = floor(czas_symulacji/Tp);

% wektory dla modelu liniowego
h2lin_eu = h2_0 * ones(l_iter, 1);
h1lin_eu = h1_0*ones(l_iter, 1);
h1lin = h1_0 * ones(l_iter, 1);
h2lin = h2_0 * ones(l_iter+1, 1);

% wektory dla sterowania (wielkości wejściowej, w celu realizacji
% opóźnienia)
F1 = F1_0*ones(l_iter+1,1);
F1_in = F1_0*ones(l_iter+1,1);
dU =  (zeros(1, Nu))';
dUp = (zeros(1,(D-1)))';
h2 = h2_0 * ones(l_iter, 1);
h1 = h1_0*ones(l_iter, 1);

krok = 8000;
% FD = [(FD)*ones(1,krok),(FD+2)*ones(1,l_iter+krok-1)];
FD = [(FD)*ones(1,l_iter+1)];
for i = 1:(l_iter)  
    % wyznaczenie sterowania
    ek = Yzad(i) - h2(i);
    
    dukk = ke*ek - Ku*dUp;
    for j = (D-1):-1:2
        dUp(j) = dUp(j-1);
    end
    dUp(1) = dukk;

    F1_in(i+1) = F1_in(i) + dukk;

    % symulacja
    
    %--Realizacja opznienia
    if(i - tau/Tp > 0)
        F1(i) = F1_in(i - floor(tau/Tp));
    else
        
        F1(i) = F1_0;
    end

    %symulacja modelu nieliniowego
    h = nlin_eu(h1(i),h2(i),dane,F1(i),FD(i));
    h1(i+1) = h(1);
    h2(i+1) = h(2);
    
end

% wykresy
czas=0:Tp:czas_symulacji;
Y = Y(1:end);          %Odpowied bez "0"
figure
stairs(h2)%czas,h2)
hold on
stairs(Yzad')%czas,Yzad')
hold off
% stairs(czas,F1)
% hold off
