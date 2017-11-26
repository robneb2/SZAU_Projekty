function [ke,Ku] = DMC_param(p,d)
dane = d;

% wektory dla modelu liniowego
h2lin_eu = d.h2_0 * ones(d.l_iter, 1);
h1lin_eu = d.h1_0*ones(d.l_iter, 1);
h1lin = d.h1_0*ones(d.l_iter, 1)';
h2lin = d.h2_0*ones(d.l_iter, 1)';

% wektory dla sterowania (wielkości wejściowej, w celu realizacji
% opóźnienia)
F1 = zeros(d.l_iter+1,1);
F1_in = zeros(d.l_iter+1,1);

% wartosc skoku
deltaF1_0 = 1;      %skok o wartość 1 od punktu pracy

%% Symulacja obiektów
for k = 1:d.l_iter
    
    F1_in(k) = d.F1_0;
    
    % skok sterowania w 1 kroku
    if k>0
        F1_in(k) = d.F1_0+deltaF1_0;
    end
    
    %--Realizacja opznienia
    if(k - d.tau/d.Tp > 0)
        F1(k) = F1_in(k - floor(d.tau/d.Tp));
    else
        
        F1(k) = d.F1_0;
    end

    %symulacja modelu zlinearyzowanego
    [h] = lin_eu([h1lin(k),h2lin(k)]',dane,F1(k),d.FD);
    h1lin(k+1) = h(1);
    h2lin(k+1) = h(2);
end

%% Wyznaczanie skoku jednostkowego
Y = h2lin-min(h2lin);   % przesuniecie do 0 układu współrzędnych
Y = Y/deltaF1_0;        % przeskalowanie odpowiedzi (ma być odpowiedz na skok
          
%% Wyznaczanie macierzy Mp
Mp = zeros(p.N,p.D-1);          %macierz przeszla
for i = 1:p.N
    for j=1:p.D-1
        if i+j<=p.D
            Mp(i,j)=Y(i+j)-Y(j);        
        else
            Mp(i,j)=Y(p.D)-Y(j);
        end
    end
end


%% Wyznaczanie macierzy M
M = zeros(p.N,p.Nu);
for i=1:p.N
    for j=1:p.Nu
        if (i>=j)
            M(i,j) = Y(i-j+1);        % kolejne +1 bo przesuniecie, docelowo usun�c
        end
    end
end

%% Wyznaczanie macierzy K
K = ((M'*M+eye(p.Nu)*p.lambda)^(-1))*M';
K1 = K(1,:);
ke = sum(K1);
Ku = K1*Mp;
end
