function [ S, dot_S, ddot_S, T ] = kubischer_spline( W_stuetz, T_ges, delta_T )
% Erzeugt aus N_I Stuetzvektoren in W_stuetz eine Trajektorie in Form kubischer Splines
% S         := Trajektorie auf Positionsebene
% dot_S     := Trajektorie auf Geschwindigkeitsebene
% ddot_S    := Trajektorie auf Beschleunigungsebene
% T         := Zeitvektor der Trajektorie

% W_stuetz  := Stuetzpunkte
% T_ges     := Dauer der Bewegung/Interpolation
% delta_T   := Taktzeit

% Erzeugt aus n Stuetzvektoren p_i und einer Zeit T eine Trajektorie aus
%stueckweise zusammengesetzten Polynomen dritten Grades (Spline)


%% --- ARBEITSBEREICH: ------------------------------------------------
%% Dimensionen pruefen/festlegen
     N_Q = size( W_stuetz,1 );
     N_I = size( W_stuetz,2 );
     N_T = floor(T_ges/( N_I - 1 )/delta_T)*( N_I - 1 );

%T_ges neu setzen, damit T_ges exakt bei N_T erreicht ist
     T_ges = N_T*delta_T;
%% --- ENDE ARBEITSBEREICH --------------------------------------------

% Aequidistanter h-Vektor
h       = T_ges / ( N_I - 1 ) * ones( 1, N_I - 1 );

%% --- ARBEITSBEREICH: ------------------------------------------------
     S = [];
     dot_S = [];
     ddot_S = [];
     T= 0:delta_T:T_ges ; % gesamte Zeit wird erneut berechnet
                          % hier:=5.94s insgesamt 

    for i = 1:N_I-1
         [s,dot_s,ddot_s] = kubischer_spline_skalar(W_stuetz, i, N_Q, N_I, h);
         S = [S s];
         dot_S = [dot_S dot_s];
         ddot_S = [ddot_S ddot_s];
    end
    
%% --- ENDE ARBEITSBEREICH --------------------------------------------

%% Erzeuge Spline fuer jede Komponente

end


%% Hilfsfunktion

%% --- ARBEITSBEREICH: ------------------------------------------------
function [s, dot_s, ddot_s] = kubischer_spline_skalar(p, teil, N_Q,N_I,h)
%% Variablen fuer Gleichungssystem A* ddot_p = r anlegen
    A = zeros(N_I,N_I);
    r = zeros(N_I,N_Q);

% Eigentliche Spline Koeffizienten
    a = zeros(N_Q,1);
    b = zeros(N_Q,1);
    c = zeros(N_Q,1);
    d = zeros(N_Q,1);

%% Erstelle A-Matrix und r-Vektor
% Erstelle A-Matrix
    for i = 1:N_I
        if i == 1
            A(i,i) = 2*h(i);
            A(i,i+1) = h(i);
        elseif i == N_I
            A(i,i) = 2*h(i-1);
            A(i,i-1) = h(i-1);
        else
            A(i,i-1) = h(i-1);
            A(i,i) = 2*(h(i)+h(i-1));
            A(i,i+1) = h(i);
        end
    end
    
% Erstelle r-Vektor 
    for i = 1: N_I   %1:12
        for k = 1:N_Q %1:3 x y z drei Freiheitsgrad
            if i == 1
                r(i,k) = 6*(p(k,i+1)-p(k,i))/h(i);
            elseif i == 12
                r(i,k) = -6*(p(k,i)-p(k,i-1))/h(i-1);
            else
                r(i,k) = -6*(p(k,i)-p(k,i-1))/h(i-1) + 6*(p(k,i+1)-p(k,i))/h(i);
            end
        end
    end
%% Gleichungssystem loesen
    ddot_p = A\r;

%% Koeffizientenberechnung
    for k = 1:N_Q     %1:3 x y z drei Freiheitsgrad
            a(k) = (ddot_p(teil+1,k)-ddot_p(teil,k))/(6*h(teil));
            b(k) = ddot_p(teil,k)/2;
            c(k) = (p(k,teil+1)-p(k,teil))/h(teil) - h(teil)*(ddot_p(teil+1,k)+2*ddot_p(teil,k))/6;
            d(k) = p(k,teil);
    end

%% Spline generieren
    s      = zeros(N_Q,54);
    dot_s  = zeros(N_Q,54);
    ddot_s = zeros(N_Q,54);
% Schleife ueber alle Intervalle
   
      
if teil==N_I-1
    N=55;
else
    N=54;
end
        
t = 0: 0.01: 0.01*N;
for i = 1:N
       for k = 1:N_Q
           s(k,i)      = a(k)*t(i)^3 + b(k)*t(i)^2 + c(k)*t(i) + d(k);
           dot_s(k,i)  = 3*a(k)*t(i)^2 + 2*b(k)*t(i) + c(k);
           ddot_s(k,i) = 6*a(k)*t(i) + 2*b(k);
       end
end

%% --- ENDE ARBEITSBEREICH --------------------------------------------


end


