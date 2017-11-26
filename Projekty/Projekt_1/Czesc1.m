clear;clc;close all;

%% Generuje skok jednostkowy i oblicza parametry DMC
% Znienić procedurę liczenia po zaimplementowaniu (nlin na lin)

% Parametry regulatora
Nu = 170;
N = 240;
D = N;
lambda = 700;

% wartosc skoku
deltaF1_0 = 1;      %skok o wartość 1 od punktu pracy

% Stale
A2 = 310;
C1 = 0.75;
alfa1 = 16;
alfa2 = 19;
tau=150;

% punkty poczatkowe
F1_0 = 51;
FD = 7;
h2_0 = ((F1_0+FD)/alfa2)^2;     %Wyznaczono z definicji punktu stałego x*
h1_0 = ((F1_0+FD)/alfa1)^2;     %Wyznaczono z definicji punktu stałego x*

%Okres probkowania.
Tp = 1;

A11 = -(4*F1_0+4*FD-3*alfa1*sqrt(h1_0))/(6*C1*h1_0^3);
A12 = 0;
A21 = alfa1/(2*A2*sqrt(h1_0));
A22 = -alfa2/(2*A2*sqrt(h2_0));

A = [A11,A12;A21,A22];
B = [1/(3*C1*h1_0^2),0]';
E = [1/(3*C1*h1_0^2),0]';
C = [1 0];

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
dane.A = A;
dane.B = B;
dane.C = C;
dane.E = E;

%% wektory danych
czas_symulacji = 40000; 
l_iter = floor(czas_symulacji/Tp);

% wektory dla modelu liniowego
h2lin_eu = h2_0 * ones(l_iter, 1);
h1lin_eu =h1_0*ones(l_iter, 1);
h1lin=h1_0*ones(l_iter, 1)';
h2lin=h2_0*ones(l_iter, 1)';
hlin = [h1lin;h2lin];

% wektory dla sterowania (wielkości wejściowej, w celu realizacji
% opóźnienia)
F1 = zeros(l_iter+1,1);
F1_in = zeros(l_iter+1,1);

%% Symulacja obiektów
for k = 1:l_iter
    
    F1_in(k) = F1_0;
    
    % skok sterowania w 1 kroku
    if k>0
        F1_in(k) = F1_0+deltaF1_0;
    end
    
    %--Realizacja opznienia
    if(k - tau/Tp > 0)
        F1(k) = F1_in(k - floor(tau/Tp));
    else
        
        F1(k) = F1_0;
    end

    %symulacja modelu zlinearyzowanego
    %zlinearyzowany
%     hlin(:,k+1) = hlin(:,k) + Tp*(A*[hlin(1,k)-h1_0,hlin(2,k)-h2_0]'+B*(F1(k)-F1_0)+E*(FD-FD));

    [h] = lin_eu([h1lin(k),h2lin(k)]',dane,F1(k),FD);
    h1lin(k+1) = h(1);
    h2lin(k+1) = h(2);
end
% h1lin = hlin(1,:);
% h2lin = hlin(2,:);

czas=0:Tp:czas_symulacji;
figure
stairs(h2lin)

%% Wyznaczanie skoku jednostkowego
Y = h2lin-min(h2lin);   % przesuniecie do 0 układu współrzędnych
Y = Y/deltaF1_0;        % przeskalowanie odpowiedzi (ma być odpowiedz na skok
          
%% Wyznaczanie macierzy Mp
Mp = zeros(N,D-1);          %macierz przeszla
for i = 1:N
    for j=1:D-1
        if i+j<=D
            Mp(i,j)=Y(i+j)-Y(j);        
        else
            Mp(i,j)=Y(D)-Y(j);
        end
    end
end

%% Wyznaczanie macierzy M
M = zeros(N,Nu);
for i=1:N
    for j=1:Nu
        if (i>=j)
            M(i,j) = Y(i-j+1);        % kolejne +1 bo przesuniecie, docelowo usun�c
        end
    end
end

%% Wyznaczanie macierzy K
K = ((M'*M+eye(Nu)*lambda)^(-1))*M';
K1 = K(1,:);
ke = sum(K1);
Ku = K1*Mp;

%%	Odpowiedz ukladu symulowanego z regulacji DMC
Y = zeros(1,l_iter);                            %Wektor odpowiedzi uk³adu, czyszczony dla nowych wyznaczen
% YzadDMC = ((h2_0+2)*ones(1,N))';                         %Wektor wartoci zadanych                  -- zadaæ
Yhat = zeros(1,N)';                              %Wektor predykowany odpowiedzi uk³adu      -- odpowied wyznaczona

% Yzad = [h2_0*ones(1,500), YzadDMC(1)*ones(1,l_iter-499)];
% Yzad = [(h2_0)*ones(1,1000),(h2_0+1)*ones(1,l_iter+999)];
Yzad = (h2_0)*ones(1,l_iter+1);
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

kroczek = 2000;


krok = 4000;
FD0=FD;
Yzadwy=zeros(1,7);
FDwy=zeros(1,7);
F1wy=zeros(1,7);
tabdYzad={-6, -4, -2, 0, 2, 4, 6};
tabdFD={-6, -4, -2, 0, 2, 4, 6};
plotcolor={[1 0.9882 0.5333],[0.3608 0.9961 0.3804],[0.9961 0.9804 0.2784],[0.0941 0.9961 0.1176],[0.9804 0.9529 0.0431],[0.0039 0.7725 0.0235],[0.9882 0.0118 0.0118],[0.2235 0.0588 1],[0.9882 0.2941 0.2941],[0.3882 0.2588 1],[0.9922 0.5569 0.5569],[0.5373 0.4431 1],[0.92 0.4 0.6],[0.7 0.443 1],[0.8 0.3 0.7],[0.4 0.2 0.6]};
figure
for m = 1 : 7
czas=0:Tp:czas_symulacji;
dYzad =tabdYzad{m};
dFD =0;%tabdFD{m};
Yzad = [(h2_0)*ones(1,kroczek),(h2_0+dYzad)*ones(1,l_iter+kroczek-1)];
FD = [(FD0)*ones(1,krok),(FD0+dFD)*ones(1,l_iter+krok-1)];
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

    
    h = nlin_eu(h1(i),h2(i),dane,F1(i),FD(i));
    h1(i+1) = h(1);
    h2(i+1) = h(2);
    
    Yzadwy(m)=Yzad(k);
    FDwy(m)=FD(krok+1);
    F1wy(m)=F1(k);
end
    %zmiana yzad
%     stairs(czas,h2,'Color',plotcolor{2*m},'LineWidth',1.5);
%     grid on; hold on
%     czas=1:Tp:length(Yzad);
%     stairs(czas,Yzad','Color',plotcolor{2*m+1},'LineWidth',1.5);
%     ylabel('h_2[cm]'); xlabel('czas[s]');
    %zmiana fd
%     stairs(czas,h2,'Color',plotcolor{2*m},'LineWidth',1.5);
%     grid on; hold on
%     if m==1
%     czas=1:Tp:length(Yzad);
%     stairs(czas,Yzad','Color',plotcolor{2},'LineWidth',1.5);
%     end
%     ylabel('h_2[cm]'); xlabel('czas[s]');
    %sterowanie
    stairs(czas,F1,'Color',plotcolor{2*m},'LineWidth',1.5)
    hold on
%     if m==1
%     czas=1:Tp:length(Yzad);
%     stairs(czas,Yzad','Color',plotcolor{2},'LineWidth',1.5);
%     end
    ylabel('F1[cm^3/s]'); xlabel('czas[s]');
    hold on
    
end

% wykresy
Y = Y(1:end);          %Odpowied bez "0"
%zmianayzad
% legend(sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(1),FDwy(1)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(1)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(2),FDwy(2)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(2)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(3),FDwy(3)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(3)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(4),FDwy(4)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(4)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(5),FDwy(5)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(5)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(6),FDwy(6)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(6)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(7),FDwy(7)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(7)),'Location','eastoutside');
%    
%zmianafd
% legend(sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(1),FDwy(1)),sprintf('Wartosc zadana Y^z^a^d=%dcm',Yzadwy(1)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(2),FDwy(2)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(3),FDwy(3)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(4),FDwy(4)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(5),FDwy(5)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(6),FDwy(6)),...
%        sprintf('Model nieliniowy Y^z^a^d=%dcm,FD=%dcm^3/s',Yzadwy(7),FDwy(7)),'Location','eastoutside');   
%sterowanie dla zmian fd
% legend(sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(1),FDwy(1)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(2),FDwy(2)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(3),FDwy(3)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(4),FDwy(4)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(5),FDwy(5)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(6),FDwy(6)),...
%        sprintf('Model nieliniowy F1=%dcm^3/s,FD=%dcm^3/s',F1wy(7),FDwy(7)),'Location','eastoutside');   
%    
%sterowanie dla zmian yzad
legend(sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(1),Yzadwy(1)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(2),Yzadwy(2)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(3),Yzadwy(3)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(4),Yzadwy(4)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(5),Yzadwy(5)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(6),Yzadwy(6)),...
       sprintf('Model nieliniowy F1=%dcm^3/s,Y^z^a^d=%dcm',F1wy(7),Yzadwy(7)),'Location','eastoutside');   
   


