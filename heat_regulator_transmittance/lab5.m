clear;
close all ;

%stałe
qgN = 10000; % W
TzewN = -20; % oC
TwewN = 20; % oC
TpN = 10; % oC
a = 0.3; 

Kcw = qgN/((a+1)*TwewN-TzewN-a*TpN); 
Kcp = a*qgN*(TwewN-TpN)/((TpN-TzewN)*((a+1)*TwewN-TzewN-a*TpN)); 
Kcwp = a*Kcw;

%Objętość pomieszczenia
a_w = 5;
b_w = 5;
c_w = 2;
Vw = a_w*b_w*c_w;

a_p = 3;
b_p = 4;
c_p = 2;
Vp = a_p*b_p*c_p;

ro = 1.2;
cp = 1000;

%Obliczenie ciepła
Cv = cp*ro*Vw;
Cvp = cp*ro*Vp;

%Warunki początkowe
qg0 = qgN; 
%Tzew = [-20, -5, 5, 20] ;

%qg = [0, 0.1*qgN, 0.5*qgN, qgN];
Tzew0 = TzewN; 

Tp0 = (Tzew0*(Kcwp*Kcp+Kcwp*Kcw+Kcw*Kcp)+Kcwp*qg0)/(Kcwp*Kcp+Kcwp*Kcw+Kcw*Kcp);%((Twew0*Kcwp)+(Tzew0*Kcp))/(Kcwp+Kcp);
Twew0 = (qg0+(Tzew0*Kcw)+(Tp0*Kcwp))/(Kcwp+Kcw);

%Transmitacja

% M1 = Cvw*s + Kcw + Kcwp;
% M2 = Cvp*s + Kcwp + Kcp;

%Twews = (M2*qgs+(M2*Kcw+Kcwp)*Tzews)/(M2*M1-Kcwp^2)
%Tps = (Kcwp*qgs+(Kcw*Kcwp+M1*Kcp)*Tzews)/(M2*M1-Kcwp^2)

M = [Cvp*Cv, Cvp*Kcw+Cvp*Kcwp+Cv*Kcwp+Cv*Kcp, Kcwp*Kcw + Kcp*Kcw + Kcp*Kcwp];
L_11 = [Cvp, Kcwp + Kcp];
L_12 = [Cvp*Kcw, Kcwp*Kcw + Kcp*Kcw + Kcwp];
L_21 = Kcwp;
L_22 = [Cv*Kcp, Kcw*Kcwp + Kcw*Kcp + Kcwp*Kcp];

% G_11 = L_11/M;
% G_12 = L_12/M;
% G_21 = L_21/M;
% G_22 = L_22/M;


%Dane do symulacji
t0 = 500; %moment wystąpienia skoku
czas_symulacji = 5000;
dTzew = 0; %zakłucenia 
dqg = 0.1*qgN; %0; %zakłucenia
[s] = sim('sim5');

 grid on
 hold on
 subplot(2,1,1);
 plot(s, Twew);
 title("T_w_e_w");
 xlabel("s");
 ylabel("T_w_e_w");
 hold off
 
 
 hold on
 grid on
 subplot(2,1,2);
 plot(s, Tp);
 title("T_p");
 xlabel("s");
 ylabel("T_p");
 hold off
 grid on