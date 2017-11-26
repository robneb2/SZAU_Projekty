clear;clc;close all;
%% Takagi-Sugemo dla wejścia 40-70[cm^3/s]
char_statyczna = 1;     % czy wyznaczyc char. stat.

% liczba modeli lokalnych
l_modeli = 4;           % 2-5 -> zmienia

% dodatkowe parametry
min_F1 = 35;
max_F1 = 65;
zakres = max_F1-min_F1;             % zakres zmian wielkos�ci wyjsciowej
delta_Fin = zakres/(l_modeli+1);    % do wyznaczenia punktow pracy
delta_Fin_przed = zakres/l_modeli;  % do wyznaczanie przedzialoww 
%
% Stale
A2 = 310;
C1 = 0.75;
alfa1 = 16;
alfa2 = 19;
tau = 150;
FD = 7; 

%% Okres probkowania.
Tp = 1;
czas_symulacji = 20000;
l_iter = floor(czas_symulacji/Tp);

%% struktura danych w celu ułatwienia przekazywania do funkcji
dane.A2 = A2;
dane.C1 = C1;
dane.alfa1 = alfa1;
dane.alfa2 = alfa2;
dane.tau = tau;
dane.FD = FD;
dane.Tp = Tp;
dane.l_iter = l_iter;

%% punkty pracy
F1_0 = (min_F1+delta_Fin)*ones(l_modeli,1);
h2_0 = zeros(l_modeli,1);
h1_0 = zeros(l_modeli,1);
for i = 1:l_modeli
    if i<l_modeli
        F1_0(i+1) = F1_0(i)+delta_Fin;
    end
    h2_0(i) = ((F1_0(i)+FD)/alfa2)^2;     
    h1_0(i) = ((F1_0(i)+FD)/alfa1)^2; 
end

% Modele (parametry)
A = zeros(2,2,l_modeli);
B = zeros(2,l_modeli);
E = zeros(2,l_modeli);
C = zeros(2,l_modeli); 
for i = 1:l_modeli
    A11 = -(4*F1_0(i)+4*FD-3*alfa1*sqrt(h1_0(i)))/(6*C1*sqrt(h1_0(i)));
    A12 = 0;
    A21 = alfa1/(2*A2*sqrt(h1_0(i)));
    A22 = -alfa2/(2*A2*sqrt(h2_0(i)));

    A(:,:,i) = [A11,A12;A21,A22];
    B(:,i) = [1/(3*C1*h1_0(i)^2),0]';
    E(:,i) = [1/(3*C1*h1_0(i)^2),0]';
    C(:,i) = [1 0];
end

if char_statyczna == 1
    Fin = 30:70;
    figure;plot(Fin,((Fin+FD)/alfa2).^2);
    title('Charakterystyka statyczna');
    xlabel('F1');ylabel('h2');grid on;axis equal;clear Fin;
end

% wektory dla matlaba
h1linfuz = zeros(l_iter,l_modeli);
h2linfuz = zeros(l_iter,l_modeli);
h1lin = zeros(l_iter,1);
h2lin = zeros(l_iter,1);
mi = zeros(l_modeli,1);
F1_in = 51*ones(l_iter,1);

F1_in = 35:1:65;
mi = zeros(length(F1_in),l_modeli);

deltaF1_0 = 1;
for k = 1:length(F1_in)%l_iter
    % funkcje przynależności
    for j = 1:l_modeli
        if j == 1
            mi(1) = 1-1/(1+exp(-1*(F1_in(k)-(min_F1+delta_Fin_przed))));
        elseif j == l_modeli
            mi(l_modeli) = 1/(1+exp(-1*(F1_in(k)-(max_F1-delta_Fin_przed))));
        else
            mi(j) = 1/(1+((F1_in(k)-(min_F1+((j-1)*(zakres/(l_modeli-1)))))/4)^(2*3));
        end
    end    
    
%     % regularyzacja wag
%     if sum(mi) ~= 1
%         mi = mi/sum(mi);
%     end
% 
%     % symulacja modeli
%     F1_in(k) = 51;
%     
    % skok sterowania w 1 kroku
    if k>0
        F1_in(k) = 51+deltaF1_0;
    end

    %--Realizacja opznienia
    if(k - tau/Tp > 0)
        F1(k) = F1_in(k - floor(tau/Tp));
    else
        F1(k) = 51;
    end
%     %symulacja modelu zlinearyzowanego
%     for j = 1:l_modeli    
%         dane.A = A(:,:,j);
%         dane.B = B(:,j);
%         dane.C = C(:,j);
%         dane.E = E(:,j);
%         dane.F1_0 = F1_0(j);
%         dane.h1_0 = h1_0(j);
%         dane.h2_0 = h2_0(j);
%         [h] = lin_eu([h1linfuz(k,j),h2linfuz(k,j)]',dane,F1(k),dane.FD);
%         h1linfuz(k+1,j) = h(1);
%         h2linfuz(k+1,j) = h(2);
%     end
%     h1lin(k+1) = mi'*h1linfuz(k+1,:)';
%     h2lin(k+1) = mi'*h2linfuz(k+1,:)';
end
% czas=0:Tp:czas_symulacji;
for i = 1:l_modeli
    
    plot(F1_in,mi(i,1))
    hold on
end
% figure
% plot(czas,h2lin)
    
    