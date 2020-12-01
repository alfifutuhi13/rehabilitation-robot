%PID TEST
% Last edit 2 Mei 2019 05.07 WIB

double i ;%Input
double o ;%Output
double e ;%Error
double mv ;%Measured Value
double int_e ;
double prev_e;
double d_e;
double dt;
double n;
double t;

i=2500;
o=0;
dt=0.1;
t=0;

%Parameter PID Kontrol
Kp=0.4;
Ki=0.3;
Kd=0.01;

%inisiasi
int_e =0;
d_e=0;
prev_e=0;

e=abs(i-o);
n=0;
curve1=animatedline('Color','r');
curve2=animatedline('Color','b');
grid on

while e>0.01 
    t=t+dt;
    n=n+1;
    e=abs(i-o);
    int_e = int_e + e.*dt;
    d_e = (e-prev_e)/dt;
    
    o=Kp.*e + Ki.*int_e + Kd.*d_e;
    
    prev_e=e;
    addpoints(curve1,t,o);
    addpoints(curve2,t,i);
    
    drawnow limitrate
    %plot([0 o],[o 5])
    %set(gca,'Xlim',[0 t],'Ylim',[0 i+0.5]);
    o;
    pause(0.01);
end

