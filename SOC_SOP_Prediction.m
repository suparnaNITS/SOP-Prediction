%SOC Estimation and SOP Prediction Algorithm
%% SOC Estimation: Fractional order Extended Kalman Filter  

clear;close all; clc
horizon=195575;% length of input current profile
C_Ah = (6.222*2.7)/3600; Gk = 1/(3600*C_Ah);% Coulombic Efficiency
load('ExpData.mat');% Load experimental data (Current Profile and Terminal Voltage)
T=1;%sampling time in second
Q=1e4*[0.01 0.01;0.01 0.01];% Process noise co-variance
Rk=.0001;% Measurement noise co-variance
%%Initial state and error covariance matrix
vc_capk0=0.001;% Initial voltage across capacitor
soc_capk0=3e-4;% Initial SOC
xcapk0=[vc_capk0;soc_capk0];%[state matrix]=[vc;SOC]
xcapk(:,1)=xcapk0;
vc_capk=xcapk(1,1); soc_capk=xcapk(2,1);
P = 1e-10*[1 0;0 1];%Initial error co-variance  
 for i=1:horizon
u1=Current(1:1:horizon);% Experimental Current profile 
u=u1(i,1);
Vol=Voltage(1:1:horizon);% Experimental Terminal voltage data
Y1(1,i)=Vol(i,1); 
soc_capk=xcapk(2,i);
if i<=54616  % Charging duration
    % Calculate open circuit voltage from estimated SOC
    ocv_capk=-2.069*soc_capk^4 + 4.843*soc_capk^3 - 4.377*soc_capk^2 + 4.245*soc_capk -0.0019;                   % considering leakage
    %Model parameters during charging
    X=[3114.50,13.72,273096.90,0.63,0.81,1];
    ocv_Capk(i)=ocv_capk;
else
    %Model parameters during charging
     X=[3380.75,38.42,223176.09,0.70,0.97,1]; 
     ocv_capk=3.9*soc_capk^5 - 8.493*soc_capk^4 + 5.303*soc_capk^3 - 0.1217*soc_capk^2 + 2.025*soc_capk + 0.02897;     % considering leakage
     ocv_Capk(i)=ocv_capk;
end  
v0=Y1(1,i);
if v0==0
    v0=0.001;
end
%% Time Update 
%x1=R;x2=C;x3=RL;x4=Resr;x5 and x6 are order of differential equation;
x1=X(:,1);x2=X(:,2);x3=X(:,3);x4=X(:,4);x5=X(:,5);x6=X(:,6);
% System matrix and input matrix calculation
A=[(-(1/(x1*x2))-(1/(x3*x2))) 0;-(Gk/x3) 0];
aa=[x5;x6];
aa1=1.^aa;
B=[(1/x2);(Gk)];
Bn=B.*aa1;
fkkest1=-ocv_capk/(x3*x2);
fkkest2=-((Gk*ocv_capk)/x3);
fkest=[fkkest1;fkkest2];
Ad=A.*aa1+(aa.*eye(2));
xcapk(:,i+1)=Ad*xcapk(:,i)+Bn*u+fkest;
% Set memory length for fractional order calculus
if i<=100
               m=i;
   else
        m=100;
end
 sum1=0;
% Calculation of binomial term in fractional order calculus
 for r=2:m
          binom=[(gamma(x5+1)/(gamma(r+1)*gamma(x5-r+1))) 0; 0 (gamma(x6+1)/(gamma(r+1)*gamma(x6-r+1)))];
          sum1=sum1+((-1)^r)*(binom*xcapk(:,i+1-r));
 end
xcapk(:,i+1)=xcapk(:,i+1) - sum1;

vc_capk=xcapk(1,i+1);
soc_capk=xcapk(2,i+1);
% Calculation of Jacobian matrices
if i<=54616  
 Fk=[ x5 - 1/(x1*x2) - 1/(x2*x3)+gamma(x5 + 1)/gamma(x5),((2069*soc_capk^3)/250 - (14529*soc_capk^2)/1000 + (4377*soc_capk)/500 - 849/200)/(x2*x3);
    0,x6 - Gk/x3+ gamma(x6 + 1)/gamma(x6)+(Gk*((2069*soc_capk^3)/250 - (14529*soc_capk^2)/1000 + (4377*soc_capk)/500 - 849/200))/x3];
else
 Fk=[ x5 - 1/(x1*x2) - 1/(x2*x3)+gamma(x5 + 1)/gamma(x5),-((39*soc_capk^4)/2 - (8493*soc_capk^3)/250 + (15909*soc_capk^2)/1000 - (1217*soc_capk)/5000 + 81/40)/(x2*x3);
      0,x6 - Gk/x3+ gamma(x6 + 1)/gamma(x6) - (Gk*((39*soc_capk^4)/2 - (8493*soc_capk^3)/250 + (15909*soc_capk^2)/1000 - (1217*soc_capk)/5000 + 81/40))/x3];
end
if i<=54616  
Hk=[ 1, - (2069*soc_capk^3)/250 + (14529*soc_capk^2)/1000 - (4377*soc_capk)/500 + 849/200];
else
Hk=[ 1, (39*soc_capk^4)/2 - (8493*soc_capk^3)/250 + (15909*soc_capk^2)/1000 - (1217*soc_capk)/5000 + 81/40];
end
F2=Fk+[(gamma(x5+1)/(gamma(1+1)*gamma(x5-1+1))) 0; 0 (gamma(x6+1)/(gamma(1+1)*gamma(x6-1+1)))];
FP=zeros(2,2);
if i<=100
        p=i;
else
        p=100;
end

for m=3:p-1 
    binom2=[(gamma(x5+1)/(gamma(m+1)*gamma(x5-m+1))) 0;0 (gamma(x6+1)/(gamma(m+1)*gamma(x6-m+1)))];
    dbinomi=diag(binom2);
    FP(:,1:2)=FP(:,1:2)+ (binom2*P(:,(2*(i)-(2*m)-1):(2*(i)-2*m))*binom2');
end
%Calculation of a priori error covariance
if i==1
    P(:,(2*(i)-1):2*(i))= P(:,(2*(i)-1):2*(i));
    else
     P(:,(2*(i)-1):2*(i))=P(:,(2*(i)-3):(2*(i)-2));
end
Pi(:,:)= (F2*P(:,(2*(i)-1):2*(i))*F2')+Q+FP;
% Kalman Gain calculation
Kk= (Pi*Hk')*inv(Hk*Pi*Hk'+Rk);
%% Measurement Update
P(:,2*(i)+1:2*(i)+2)=(eye(2)-(Kk*Hk))*Pi*(eye(2)-Kk*Hk)';
yest(1,i)=vc_capk +x4*u + ocv_capk;
xcapk(:,i+1)=xcapk(:,i+1)+Kk*(Y1(1,i)-yest(1,i));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 end
timeVector = (0:T:horizon)/3600;
figure(1)
subplot(2,1,1)
plot(timeVector,xcapk(2,:),'g-','Linewidth',1);
legend('Estimated SOC')
xlabel('Time(Hr)','FontWeight', 'bold')
ylabel('SOC(0-1)','FontWeight', 'bold')
grid on;
 %%%% SOP Prediction

L=1;% Length of rediction horizon
for t1=1:horizon
  %% Current for SOC Constraint (In this study Max_SOC=.9, Min_SOC=.1)
Vc_capk(1,:)=xcapk(1,:);%Estimated voltage across capacitor(Vc)
SOC_capk(1,:)=xcapk(2,:);% Estimated SOC  
Vc=Vc_capk(1,t1);%Estimated Vc at t1 
SOC=SOC_capk(1,t1);%Estimated SOC at t1 
if t1<=54616
      I_SOC1=(.9-SOC)/(Gk*L*T)+(Vc+ocv_Capk(t1))/x3; % Calculation of current for SOC constraint during charging
  if I_SOC1>=0
        I_SOC(:,t1)=I_SOC1;
    else
        I_SOC(:,t1)=0;
   end
else
    I_SOC1=(.1-SOC)/(Gk*L*T)+(Vc+ocv_Capk(t1))/x3; % Calculation of current for SOC constraint during discharging

    if I_SOC1<=0
         I_SOC(:,t1)=I_SOC1;
     else
         I_SOC(:,t1)=0;
     end
end
 
   
    %% Current for Voltage Constraint
if t1<=54616 
     X=[3114.50,13.72,273096.90,0.63,0.81,1];%Model parameters during charging
else
     X=[3380.75,38.42,223176.09,0.70,0.97,1];%Model parameters during discharging
end 
    
sumc=0;
%SOC=xcapk(2,t1);% SOC at t1 when prediction starts taken from FOEKF
OCVK=ocv_Capk(t1);% Open circuit voltage at t1 considering estimated SOC of FOEKF
Ad1=(-(1/(x1*x2))-(1/(x3*x2)))*(T)^x5 +x5; % Discretized system matrix for fractional order state space model
B1=(1/x2)*(T)^x5;% Discretized input matrix for fractional order state space model
sumc0=0;
if t1<=100
        s=t1;       
   else
        s=100;
end
for i1= 0:(L-1)
      for r=2:s
          binom=(gamma(x5+1)/(gamma(r+1)*gamma(x5-r+1)));
          sumc=sumc+((-1)^r)*(binom*Vc_capk(1,t1+i1+1-r));
      end
    sumc0=sumc0+Ad1^(L-1-i1)*sumc;
end
 sumc2(1,t1)=sumc0;  
 sum3c=0;
 for p=0:L-1
          sum3c=sum3c+ ((Ad1^(L-1-p))*B1);
 end
 % Solve the equation current for voltage constraint equation using 'fzero'
if t1<=54616
  f=@(IL_V)(Ad1^L)*Vc+sum3c*IL_V-sumc2(1,t1)-((1/(x3*x2))*(Ad1^(L-1)))* ((-2.069)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(4.843)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(-4.377)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(4.245)*(SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3))-.0019)+x4*IL_V+((-2.069)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(4.843)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(-4.377)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(4.245)*(SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3))-.0019)-2.7;
  I_Voltage(t1)=fzero(f,01.2);
else
  f=@(IL_V)(Ad1^L)*Vc+sum3c*IL_V-sum3c*sumc2(1,t1)-((1/(x3*x2))*(Ad1^(L-1)))* ((3.9)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^5+(- 8.493)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(5.303)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(- 0.1217)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(2.025)*(SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3))+0.02897)+x4*IL_V+((3.9)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^5+(- 8.493)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(5.303)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(- 0.1217)*((SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(2.025)*(SOC+(((Gk*T*L))*IL_V) -((Gk*T*L))*((Vc+OCVK)/x3))+0.02897)-0.001;
  I_Voltage(t1)=fzero(f,01.2);
end
Il_v=I_Voltage(t1);
%% Current for Current Constraint
if t1<=54616
IL_M(:,t1)=1.85;% Calculted maximum current from datasheet of Supercapacitor
else
 IL_M(:,t1)=-1.85;
end
%% Current for Multiconstraint
 Il_v=I_Voltage(t1);
 Il_m=IL_M(:,t1);
 Il_soc=I_SOC(:,t1);
 A1=[Il_v,Il_m,Il_soc];
 if t1<=54616
 Imax(t1)=min(A1);% Finding maximum current
 else
 Imax(t1)=max(A1);% Finding maximum current
 end
 %IL_V=min(A1);
 %% SOP
 IL=Imax(t1);
 % Terminal voltage after applying maximum current for prediction horizon
if t1<=54616
%IL=min(A1);
V_T_Multi(:,t1)=(Ad1^L)*Vc+sum3c*IL-sumc2(1,t1)-((1/(x3*x2))*(Ad1^(L-1)))* ((-2.069)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(4.843)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(-4.377)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(4.245)*(SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3))-.0019)+x4*IL+((-2.069)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^4+(4.843)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^3+(-4.377)*((SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3)))^2+(4.245)*(SOC+(((Gk*T*L))*IL) -((Gk*T*L))*((Vc+OCVK)/x3))-.0019);
else
%IL=-min(A1);
V_T_Multi(:,t1)=(Ad1^L)*Vc+sum3c*IL-sumc2(1,t1)-((1/(x3*x2))*(Ad1^(L-1)))* ((3.9)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^5+(-8.493)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^4+(5.303)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^3+(-0.1217)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^2+(2.025)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))+0.02897)+x4*IL+((3.9)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^5+(-8.493)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^4+(5.303)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^3+(-0.1217)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))^2+(2.025)*(SOC+((Gk*T*L)*IL)-((Gk*T*L)*((Vc+OCVK)/x3)))+0.02897);
end
Power(1,t1)=abs(IL)*V_T_Multi(:,t1);%Power calculation for Multi constraint

end

timeVector1 = (0:T:54616-1)/3600;
timeVector2 = (0:T:140959-1)/3600;
figure(2)
plot(timeVector1,I_Voltage(1,1:54616),'r-' ,timeVector1,I_SOC(1,1:54616),'b-',timeVector1,IL_M(1,1:54616),'g-',timeVector1,Imax(1,1:54616),'k:','Linewidth',2)
legend('I_{Voltage constraint}','I_{SOC constraint}','I_{Current constraint}','I_{Multi constraint}')
xlabel('Time(Hr)','FontWeight', 'bold')
ylabel('Current(A)','FontWeight', 'bold')
title('Current prediction during Charging:L=1')
grid on;
figure(3)
plot(timeVector2,I_Voltage(1,54617:195575),'r-' ,timeVector2,I_SOC(1,54617:195575),'b-',timeVector2,IL_M(1,54617:195575),'g-',timeVector2,Imax(1,54617:195575),'k:','Linewidth',2)
legend('I_{Voltage constraint}','I_{SOC constraint}','I_{Current constraint}','I_{Multi constraint}')
xlabel('Time(Hr)','FontWeight', 'bold')
ylabel('Current(A)','FontWeight', 'bold')
title('Current prediction during disharging:L=1')
grid on;
figure(4)
plot(timeVector1,Power(1,1:54616))
legend('Power')
xlabel('Time(Hr)','FontWeight', 'bold')
ylabel('Power(Watt)','FontWeight', 'bold')
title('SOP Prediction During Charging: L=1')
figure(5)
plot(timeVector2,Power(1,54617:195575))
legend('Power(1,1:54616)')
xlabel('Time(Hr)','FontWeight', 'bold')
ylabel('Power(Watt)','FontWeight', 'bold')
title('SOP Prediction During Disharging: L=1')
