function [ S, dot_S, ddot_S, T ] = p2p_quintisch( W_stuetz, T_ges, delta_T)
  % Erzeugt aus N_I Stuetzpunkten in W_stuetz je mit Anfangs- und Endpunkt N_I-1 Trajektorien je in Form eines quintischen Polynoms
  % S         := Trajektorie auf Positionsebene
  % dot_S     := Trajektorie auf Geschwindigkeitsebene
  % ddot_S    := Trajektorie auf Beschleunigungsebene
  % T         := Zeitvektor der Trajektorie

  % W_stuetz  := Stuetzpunkte
  % T_I(end)     := Dauer der Bewegung/Interpolation
  % delta_T   := Taktzeit

  % Anzahl der Freiheitsgrade
  N_Q       = size( W_stuetz,1 );

  % Anzahl der Stuetzpunkte
  N_I       = size( W_stuetz,2 );

  % Zeitintervall fuer Interpolation
  T_I       = 0:delta_T:(T_ges/N_I);  % Zeitintervall fuer ein Teilstueck
  N_T_I     = length(T_I);          % Anzahl der Zeitpunkte eines Teilstuecks

  %% Berechnung der Trajektorie

  % Initialisierung fuer Teilstuecke
  S_I       = zeros( N_Q, N_T_I );
  dot_S_I   = zeros(size(S_I));
  ddot_S_I  = zeros(size(S_I));

  % Initialisierung fuer Gesamttrajektorie
  S         = [];
  dot_S     = [];
  ddot_S    = [];
  T         = [];
  %% --- ARBEITSBEREICH: ------------------------------------------------
% Erzeuge Trajektorie 
%     S = ?;
%     dot_S = ?;
%     ddot_S = ?;

  % Schleife ueber Stuetzpunktepaare
  for i = 1:N_I-1


  end
  %% --- ENDE ARBEITSBEREICH --------------------------------------------
end % function
