function [ S, dot_S, ddot_S, T ] = p2p_kubisch( W_stuetz, T_ges, delta_T )
% Erzeugt aus N_I Stuetzpunkten in W_stuetz je mit Anfangs- und Endpunkt N_I-1 Trajektorien je in Form eines kubischen Polynoms
% S         := Trajektorie auf Positionsebene
% dot_S     := Trajektorie auf Geschwindigkeitsebene
% ddot_S    := Trajektorie auf Beschleunigungsebene
% T         := Zeitvektor der Trajektorie

% W_stuetz  := Stuetzpunkte
% T_ges     := Dauer der Bewegung/Interpolation
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
S         = zeros(N_Q,(N_T_I-1)*(N_I-1)+1);%N_Q=3, N_T_I=51, N_I=12
dot_S     = zeros(size(S));%
ddot_S    = zeros(size(S));%
T         = 0:delta_T:delta_T*550;%

%% --- ARBEITSBEREICH: ------------------------------------------------

% Schleife ueber Stuetzpunktepaare
t=0;
koeffx=zeros(4,11);
koeffy=zeros(4,11);
koeffz=zeros(4,11);
A=[t^3,       t^2,     t,  1;
       (t+0.5)^3,(t+0.5)^2,(t+0.5),1;
       3*t^2,   2*t,  1,  0;   
       3*(t+0.5)^2,   2*(t+0.5),  1,  0];


for i=1:N_I-1

    koeffx(:,i)=A\[W_stuetz(1,i);W_stuetz(1,i+1);0;0];
    koeffy(:,i)=A\[W_stuetz(2,i);W_stuetz(2,i+1);0;0];
    koeffz(:,i)=A\[W_stuetz(3,i);W_stuetz(3,i+1);0;0];
    for k=1:N_T_I
      % Erzeuge Trajektorie  x
        S(1,(i-1)*(N_T_I-1)+k) = koeffx(1,i)*T_I(k)^3+koeffx(2,i)*T_I(k)^2+koeffx(3,i)*T_I(k)+koeffx(4,i);
        dot_S(1,(i-1)*(N_T_I-1)+k) = 3*koeffx(1,i)*T_I(k)^2+2*koeffx(2,i)*T_I(k)+koeffx(3,i);
        ddot_S(1,(i-1)*(N_T_I-1)+k) = 6*koeffx(1,i)*T_I(k)+2*koeffx(2,i); 
       % Erzeuge Trajektorie  y
        S(2,(i-1)*(N_T_I-1)+k) = koeffy(1,i)*T_I(k)^3+koeffy(2,i)*T_I(k)^2+koeffy(3,i)*T_I(k)+koeffy(4,i);
        dot_S(2,(i-1)*(N_T_I-1)+k) = 3*koeffy(1,i)*T_I(k)^2+2*koeffy(2,i)*T_I(k)+koeffy(3,i);
        ddot_S(2,(i-1)*(N_T_I-1)+k) = 6*koeffy(1,i)*T_I(k)+2*koeffy(2,i);          
      % Erzeuge Trajektorie z   
        S(3,(i-1)*(N_T_I-1)+k) = koeffz(1,i)*T_I(k)^3+koeffz(2,i)*T_I(k)^2+koeffz(3,i)*T_I(k)+koeffz(4,i);
        dot_S(3,(i-1)*(N_T_I-1)+k) = 3*koeffz(1,i)*T_I(k)^2+2*koeffz(2,i)*T_I(k)+koeffz(3,i);
        ddot_S(3,(i-1)*(N_T_I-1)+k) = 6*koeffz(1,i)*T_I(k)+2*koeffz(2,i); 
    end  

end
%% --- ENDE ARBEITSBEREICH --------------------------------------------
end % function
