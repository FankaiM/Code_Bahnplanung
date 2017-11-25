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
  T_I       = 0:delta_T:T_ges/(N_I-1);  % Zeitintervall fuer ein Teilstueck
  N_T_I     = length(T_I);          % Anzahl der Zeitpunkte eines Teilstuecks

  %% Berechnung der Trajektorie

  % Initialisierung fuer Teilstuecke
  S_I       = zeros( N_Q, N_T_I );
  dot_S_I   = zeros(size(S_I));
  ddot_S_I  = zeros(size(S_I));

  % Initialisierung fuer Gesamttrajektorie
S         = zeros(N_Q,(N_T_I-1)*(N_I-1)+1);%N_Q=3, N_T_I=51, N_I=12
dot_S     = zeros(size(S));%
ddot_S    = zeros(size(S));%
T         = 0:delta_T:T_I(N_T_I)*(N_I-1); % gesamte Zeit wird erneut berechnet
                                         % hier:
                                         % 0.54(Teilzeitintervall)*(12-1)Trajektroien
                                         % =5.94s insgesamt 
  %% --- ARBEITSBEREICH: ------------------------------------------------
% Erzeuge Trajektorie 
tic
% Koeffizientvektor berechnen
ta=0;                                    %Startzeit fuer jede Trajektroie
te=T_I(N_T_I);                           %Endezeit fuer jede Trajektroie
koeffx=zeros(6,11);
koeffy=zeros(6,11);
koeffz=zeros(6,11);
% 5.Ordnung also 6x6 Matrix 
A=[    ta^5,     ta^4,       ta^3,   ta^2, ta, 1;
       te^5,     te^4,       te^3,   te^2, te, 1;
       5*ta^4,   4*ta^3,    3*ta^2,  2*ta, 1,  0;   
       5*te^4,   4*te^3,    3*te^2,  2*te, 1,  0;   
       20*ta^3,  12*ta^2,    6*ta,   2,    0,  0;
       20*te^3,  12*te^2,    6*te,   2,    0,  0 ];

% Schleife ueber Stuetzpunktepaare
for i=1:N_I-1
% koeffx/y/z = inv(A)* Randbedingungsvektor
% Start- und End Geschwindigkeit wird Null eingesetzt.
    koeffx(:,i)=A\[W_stuetz(1,i);W_stuetz(1,i+1);0;0;0;0];
    koeffy(:,i)=A\[W_stuetz(2,i);W_stuetz(2,i+1);0;0;0;0];
    koeffz(:,i)=A\[W_stuetz(3,i);W_stuetz(3,i+1);0;0;0;0];
    for k=1:N_T_I
      % Erzeuge Trajektorie  x
        S(1,(i-1)*(N_T_I-1)+k) = koeffx(1,i)*T_I(k)^5+koeffx(2,i)*T_I(k)^4+koeffx(3,i)*T_I(k)^3+koeffx(4,i)*T_I(k)^2+koeffx(5,i)*T_I(k)+koeffx(6,i);
        dot_S(1,(i-1)*(N_T_I-1)+k) = 5*koeffx(1,i)*T_I(k)^4+4*koeffx(2,i)*T_I(k)^3+3*koeffx(3,i)*T_I(k)^2+2*koeffx(4,i)*T_I(k)+koeffx(5,i);
        ddot_S(1,(i-1)*(N_T_I-1)+k) = 20*koeffx(1,i)*T_I(k)^3+12*koeffx(2,i)*T_I(k)^2+6*koeffx(3,i)*T_I(k)+2*koeffx(4,i); 
       % Erzeuge Trajektorie y
        S(2,(i-1)*(N_T_I-1)+k) = koeffy(1,i)*T_I(k)^5+koeffy(2,i)*T_I(k)^4+koeffy(3,i)*T_I(k)^3+koeffy(4,i)*T_I(k)^2+koeffy(5,i)*T_I(k)+koeffy(6,i);
        dot_S(2,(i-1)*(N_T_I-1)+k) = 5*koeffy(1,i)*T_I(k)^4+4*koeffy(2,i)*T_I(k)^3+3*koeffy(3,i)*T_I(k)^2+2*koeffy(4,i)*T_I(k)+koeffy(5,i);
        ddot_S(2,(i-1)*(N_T_I-1)+k) = 20*koeffy(1,i)*T_I(k)^3+12*koeffy(2,i)*T_I(k)^2+6*koeffy(3,i)*T_I(k)+2*koeffy(4,i); 
      % Erzeuge Trajektorie z   
        S(3,(i-1)*(N_T_I-1)+k) = koeffz(1,i)*T_I(k)^5+koeffz(2,i)*T_I(k)^4+koeffz(3,i)*T_I(k)^3+koeffz(4,i)*T_I(k)^2+koeffz(5,i)*T_I(k)+koeffz(6,i);
        dot_S(3,(i-1)*(N_T_I-1)+k) = 5*koeffz(1,i)*T_I(k)^4+4*koeffz(2,i)*T_I(k)^3+3*koeffz(3,i)*T_I(k)^2+2*koeffz(4,i)*T_I(k)+koeffz(5,i);
        ddot_S(3,(i-1)*(N_T_I-1)+k) = 20*koeffz(1,i)*T_I(k)^3+12*koeffz(2,i)*T_I(k)^2+6*koeffz(3,i)*T_I(k)+2*koeffz(4,i); 
    end  

end
toc
  %% --- ENDE ARBEITSBEREICH --------------------------------------------
end % function
