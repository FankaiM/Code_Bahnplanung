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
    T_ges = N_T*delta_T; % hier T_ges = 5.94
%% --- ENDE ARBEITSBEREICH --------------------------------------------

% Aequidistanter h-Vektor
h       = T_ges / ( N_I - 1 ) * ones( 1, N_I - 1 );

%% --- ARBEITSBEREICH: ------------------------------------------------
    S = [];
    dot_S = [];
    ddot_S = [];
    T = 0:delta_T:T_ges;

    for i = 1:N_I-1
         [s,dot_s,ddot_s] = kubischer_spline_skalar(W_stuetz,i,h);
         S = [S;s];
         dot_S = [dot_S;dot_s];
         ddot_S = [ddot_S;ddot_s];
    end
    
    S = S';
    dot_S = dot_S';
    ddot_S = ddot_S';
%% --- ENDE ARBEITSBEREICH --------------------------------------------

%% Erzeuge Spline fuer jede Komponente

end


%% Hilfsfunktion

%% --- ARBEITSBEREICH: ------------------------------------------------
function [ s, dot_s,ddot_s] = kubischer_spline_skalar(p,teil,h)
%% Variablen fuer Gleichungssystem A* ddot_p = r anlegen
    A = zeros(12);
    r = zeros(12,3);
    
% Eigentliche Spline Koeffizienten
    a = zeros(1,3);
    b = zeros(1,3);
    c = zeros(1,3);
    d = zeros(1,3);
%% Erstelle A-Matrix und r-Vektor
   % Erstelle A-Matrix
    for i = 1:12
        if i == 1
            A(i,i) = 2*h(i);
            A(i,i+1) = h(i);
        elseif i == 12
            A(i,i) = 2*h(i-1);
            A(i,i-1) = h(i-1);
        else
            A(i,i-1) = h(i-1);
            A(i,i) = 2*(h(i)+h(i-1));
            A(i,i+1) = h(i);
        end
    end
    
    % Erstelle r-Vektor
    for i = 1:12
        for koor = 1:3
            if i == 1
                r(i,koor) = 6*(p(koor,i+1)-p(koor,i))/h(i);
            elseif i == 12
                r(i,koor) = -6*(p(koor,i)-p(koor,i-1))/h(i-1);
            else
                r(i,koor) = -6*(p(koor,i)-p(koor,i-1))/h(i-1) + 6*(p(koor,i+1)-p(koor,i))/h(i);
            end
        end
    end

%% Gleichungssystem loesen

    ddot_p = A\r;

%% Koeffizientenberechnung
        for koor = 1:3
            a(koor) = (ddot_p(teil+1,koor)-ddot_p(teil,koor))/(6*h(teil));
            b(koor) = ddot_p(teil,koor)/2;
            c(koor) = (p(koor,teil+1)-p(koor,teil))/h(teil) - h(teil)*(ddot_p(teil+1,koor)+2*ddot_p(teil,koor))/6;
            d(koor) = p(koor,teil);
        end

%% Spline generieren
    s      = zeros(54,3);
    dot_s  = zeros(54,3);
    ddot_s = zeros(54,3);
    
    
    if teil ~=11 % Fuer Teile 1-10 nehmen wir die Koordinaten von 54 Punkte jedes Teils
        for i = 1:54
           t = 0:0.01:0.53;
           for koor = 1:3
               s(i,koor)      = a(koor)*t(i)^3 + b(koor)*t(i)^2 + c(koor)*t(i) + d(koor);
               dot_s(i,koor)  = 3*a(koor)*t(i)^2 + 2*b(koor)*t(i) + c(koor);
               ddot_s(i,koor) = 6*a(koor)*t(i) + 2*b(koor);
           end
        end
    else  % Fuer Teil 11 nehmen wir die Koordinaten von 55 Punkte des Teils, sodass es insgesamt 595 Punkte zu plotten gibt
        for i = 1:55
           t = 0:0.01:0.54;
           for koor = 1:3
               s(i,koor)      = a(koor)*t(i)^3 + b(koor)*t(i)^2 + c(koor)*t(i) + d(koor);
               dot_s(i,koor)  = 3*a(koor)*t(i)^2 + 2*b(koor)*t(i) + c(koor);
               ddot_s(i,koor) = 6*a(koor)*t(i) + 2*b(koor);
           end
        end
    end

% Schleife ueber alle Intervalle

%% --- ENDE ARBEITSBEREICH --------------------------------------------

end



