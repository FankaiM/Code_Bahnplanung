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
%     N_Q = ?;
%     N_I = ?;
%     N_T = ? ;

%T_ges neu setzen, damit T_ges exakt bei N_T erreicht ist
%     T_ges = ?;
%% --- ENDE ARBEITSBEREICH --------------------------------------------

% Aequidistanter h-Vektor
h       = T_ges / ( N_I - 1 ) * ones( 1, N_I - 1 );

%% --- ARBEITSBEREICH: ------------------------------------------------
%     S = zeros( ?);
%     dot_S = zeros( ? );
%     ddot_S = zeros( ? );
%     T = ?;

%     for i = 1:?
%          [?,?,?] = kubischer_spline_skalar( ?);
%     end
%% --- ENDE ARBEITSBEREICH --------------------------------------------

%% Erzeuge Spline fuer jede Komponente

end


%% Hilfsfunktion

%% --- ARBEITSBEREICH: ------------------------------------------------
%function [ ?, ?, ?] = kubischer_spline_skalar( ? )
%% Variablen fuer Gleichungssystem A* ddot_p = r anlegen
%     A = zeros( ?);
%     r = zeros( ? );

% Eigentliche Spline Koeffizienten
%     a = zeros( ? );
%     b = zeros( ? );
%     c = zeros( ? );
%     d = zeros( ? );

%% Erstelle A-Matrix und r-Vektor


%% Gleichungssystem loesen


%% Koeffizientenberechnung
%     for i = 1:?
%         a(i) = ?;
%         b(i) = ?;
%         c(i) = ?;
%         d(i) = ?;
%     end

%% Spline generieren
%     s = zeros(?);
%     dot_s = zeros(?);
%     ddot_s = zeros(?);

% Schleife ueber alle Intervalle

%% --- ENDE ARBEITSBEREICH --------------------------------------------

%end



