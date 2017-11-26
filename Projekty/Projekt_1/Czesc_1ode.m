clear;clc;

% Sta³e
A2 = 310;
C1 = 0.75;
alfa1 = 16;
alfa2 = 19;

%% punkty poczatkowe

F1_0 = 51;
FD = 7;

h2_0 = 9.3186;
tau=150;

%% wektory danych
Tp = 5;
czas_symulacji = 1000;
l_iter = czas_symulacji/Tp;
iter_opoznienia = tau/Tp;
tspan=0:Tp:czas_symulacji;
% h2 = h2_0.* ones(l_iter, 1);
% h1 = ones(l_iter, 1);
h1_0=0;
F1=F1_0;
% F1=(F1_0.*ones(l_iter,1));

%% ode45
 
[t,h] = ode45(@(t,h)liczeniestanuciag(t,h,F1_0,FD,C1,A2,alfa1,alfa2),tspan,[h1_0;h2_0]);
h1=h(:,1);
h2=h(:,2);

%% rysowanie wykresów

figure
stairs(t,h2);
ylabel('Wysokosc'); xlabel('czas');





