%DTEE6

clear;clc;close all;

hold on
view([40 -270 90])
xlabel('x axis')
axis([-0.5 1 -0.5 1 0 1])
ylabel('y axis')
zlabel('z axis')
grid on
box off

%Inisiasi nilai
L1=0 ;%[m]
L2=0.55 ;%[m]
L3=0.51 ;%[m]

m1=0.38 ;%[kg]
m2=0.38 ;%[kg]
m3=0.50 ;%[kg]

I1=0.2176;
I2=0.0205;
I3=0.0205;

ts=0.1; %[s] TIME SAMPLING
g=9.8; %[m/s^2]

%Kp=0.12 dan KI=0.11 normal, Ki=0.55 overshoot, Kp=0.22 lumayan
Kp=0.15; 
Ki=0.1;

int_e1=0;
int_e2=0;
int_e3=0;


%Step 1
q1(1)=0;
q2(1)=90*pi/180;
q3(1)=-90*pi/180;%

q1deg(1)=q1(1)*180/pi;
q2deg(1)=q2(1)*180/pi;
q3deg(1)=q3(1)*180/pi;

dq1(1)=0;
dq2(1)=0; 
dq3(1)=0;

%INPUT, silakan titik yang diinginkan---------------------------------------------------------    
%Trajectory yang diinginkan


%Pemodelan Jalur Fleksi-Ekstensi
    %Parameter jalur
    Lp=0.33; %Lengan Pasien
    phi=40; %Derajat, sudut fleksi-ekstensi
    
    
%Titik Awal
x0=0.55;                              xo=x0;   %Kalo di arduino, di read lgsg
y0=0;                                 yo=y0;
z0=0.48;                              zo=z0;

%Titik Akhir (fungsi dari titik awal, lengan pasien, dan sudut
%fleksi-ekstensi)
xp=x0+Lp*(1-cosd(phi));
yp=y0;
zp=z0+Lp*sind(phi);

%abc quadratic generator
a=(z0-zp)/(xp-x0).^2;
b=-2*a*xp;
c=zp-a*xp.^2-b*xp;

xtr=x0:0.01:xp;
ydum=ones(1,length(xtr));
ytr=yp*ydum; %Selalu
ztr=a*xtr.^2+b*xtr+c;

hh=plot3(xtr,ytr,ztr,'m');

    %pdesired=[CPfleksi1 CPfleksi2 CPfleksi3 ... CPfleksin dummyvalue dummyvalue];
    xd=[xtr 0 0]; 
    yd=[ytr 0 0];
    zd=[ztr 0 0];
    
n=1;
CPfleksi=length(xd) - 2; %CPfleksi = checkpoint
CPekstensi=0;

%Penyederhanaan agar mempermudah hidup
   s1=sin(q1(1));
   s2=sin(q2(1));
   s23=sin(q2(1)+q3(1));
   s3=sin(q3(1));
   c1=cos(q1(1));
   c2=cos(q2(1));
   c23=cos(q2(1)+q3(1));
   c3=cos(q3(1));

xee(1)= L2*c1*c2+L3*c1*c23;
yee(1)= L2*s1*c2+L3*s1*c23;
zee(1)= L2*s2+L3*s23;



ex(1)=xp-xee;
ey(1)=yp-yee;
ez(1)=zp-zee;
i=0;
l=0;
time=0;

while ((abs(ex(i+1))>0.005) || (abs(ey(i+1))>0.005) || (abs(ez(i+1))>0.005))
   
   i=i+1;
   time=time+ts;
  
   q(i,1)=q1(i);
   q(i,2)=q2(i);
   q(i,3)=q3(i);
   
   %Penyederhanaan agar mempermudah hidup
   s1=sin(q(i,1));
   s2=sin(q(i,2));
   s23=sin(q(i,2)+q(i,3));
   s3=sin(q(i,3));
   c1=cos(q(i,1));
   c2=cos(q(i,2));
   c23=cos(q(i,2)+q(i,3));
   c3=cos(q(i,3));
   
   %Step 2
   %Matriks transformasi T01 T12 T23 = T03
   A1=[c1 0 s1 0; 
        s1 0 -c1 0; 
        0 1 0 0; 
        0 0 0 1];
   A2=[c2 -s2 0 L2*c2; 
       s2 c2 0 L2*s2; 
       0 0 1 0; 
       0 0 0 1];
   A3=[c3 -s3 0 L3*c3; 
       s3 c3 0 L3*s3; 
       0 0 1 0; 
       0 0 0 1];
   
   ea=[0;0;0;1];
   
   xyzlink1=A1*ea;
   xyzlink2=A1*A2*ea;
   xyzlink3=A1*A2*A3*ea;
   
   xlink1(i)=xyzlink1(1);
   ylink1(i)=xyzlink1(2);
   zlink1(i)=xyzlink1(3);

   xlink2(i)=xyzlink2(1);
   ylink2(i)=xyzlink2(2);
   zlink2(i)=xyzlink2(3);

   xlink3(i)=xyzlink3(1);
   ylink3(i)=xyzlink3(2);
   zlink3(i)=xyzlink3(3);
   
   %Step 4
   ex(i+1)=xd(i,1)-xlink3(i);
   ey(i+1)=yd(i,1)-ylink3(i);
   ez(i+1)=zd(i,1)-zlink3(i);
   
   h2=plot3([xlink2(i) xlink3(i)],[ylink2(i) ylink3(i)],[zlink2(i) zlink3(i)],'b');
   h1=plot3([xlink1(i) xlink2(i)],[ylink1(i) ylink2(i)],[zlink1(i) zlink2(i)],'r');
   h0=plot3([0 xlink1(i)],[0 ylink1(i)],[0 zlink1(i)],'k');
   pause(ts);
   h3=plot3(xlink3(i),ylink3(i),zlink3(i),'.g');
   hold on;
        
   set(h2,'Visible','off')
   set(h1,'Visible','off')
   set(h0,'Visible','off')
   
   %Step 3
   %Step pengembangan jalur yang didefinisikan 
   for j=1:(CPfleksi+2)
       xd(i+1,j)=xd(i,j);
       yd(i+1,j)=yd(i,j);
       zd(i+1,j)=zd(i,j);
   end
   
   if n<CPfleksi-1  
        if ((abs(ex(i+1))<0.02) && (abs(ey(i+1))<0.02) && (abs(ez(i+1))<0.02))
        n=n+1;
        %Tambahan biar bisa dituker balik
        xd(i+1,CPfleksi+1)=xd(i+1,1);
        yd(i+1,CPfleksi+1)=yd(i+1,1);
        zd(i+1,CPfleksi+1)=zd(i+1,1);
        
        for v=1:(CPfleksi-1)
            xd(i+1,v)=xd(i+1,v+1);
            yd(i+1,v)=yd(i+1,v+1);
            zd(i+1,v)=zd(i+1,v+1);
            
        end
        xd(i+1,CPfleksi)=xd(i+1,CPfleksi+1);
        yd(i+1,CPfleksi)=yd(i+1,CPfleksi+1);
        zd(i+1,CPfleksi)=zd(i+1,CPfleksi+1);
        
        h4=plot3(xlink3(i),ylink3(i),zlink3(i),'ok');
        
        end
   if n==CPfleksi-1
           if ((abs(ex(i+1))<0.02) && (abs(ey(i+1))<0.02) && (abs(ez(i+1))<0.02))
            n=n+1;
            %Tambahan biar bisa dituker balik
             xd(i+1,CPfleksi+1)=xd(i+1,1);
             yd(i+1,CPfleksi+1)=yd(i+1,1);
             zd(i+1,CPfleksi+1)=zd(i+1,1);
        
        for v=1:(CPfleksi-1)
             xd(i+1,v)=xd(i+1,v+1);
             yd(i+1,v)=yd(i+1,v+1);
             zd(i+1,v)=zd(i+1,v+1);
            
        end
             xd(i+1,CPfleksi)=xd(i+1,CPfleksi+1);
             yd(i+1,CPfleksi)=yd(i+1,CPfleksi+1);
             zd(i+1,CPfleksi)=zd(i+1,CPfleksi+1);
        
        h4=plot3(xlink3(i),ylink3(i),zlink3(i),'ok');
            
            CPekstensi=2*CPfleksi-2; %default 2*CPfleksi-2, 2*CPfleksi-1 itu bisa membuat balik ke (0,0,0)
            end
       end 
        
   end
   
   
   %Step 5
   %Jacobian velocity
   Jv(1,1) = -s1*(c23*L3 + c2*L2);
   Jv(1,2) = -c1*(s23*L3 + s2*L2);
   Jv(1,3) = -c1*(s23*L3); 
   
   Jv(2,1) = c1*(c23*L3 + c2*L2);
   Jv(2,2) = -s1*(s23*L3 + s2*L2);
   Jv(2,3) = -s1*(s23*L3); 

   Jv(3,1) = 0;
   Jv(3,2) = c23*L3+c2*L2;
   Jv(3,3) = c23*L3;
   
   Jinv=inv(Jv);
   
   dx=ex(i); %Biar bisa dimasukkin ke matrix
   dy=ey(i);
   dz=ez(i);
   
   dq=Jinv*[dx;dy;dz]; %Dalam radian
   %dq=dqm'; %Dalam Radian
   
   %Step 6
   %q123desired(i)=q(i)+dq(i);
 %  q123desired(i,1)=q(i,1)+dq(i,1);
 %  q123desired(i,2)=q(i,2)+dq(i,2);
 %  q123desired(i,3)=q(i,3)+dq(i,3);
   
   %Step 7
    eq1(i)=dq(1,1); 
    eq2(i)=dq(2,1);
    eq3(i)=dq(3,1); 
   
   %Step 8
   
   %Kontroler PID
   PID1(i)=Kp*eq1(i) + Ki * int_e1;
   PID2(i)=Kp*eq2(i) + Ki * int_e2;
   PID3(i)=Kp*eq3(i) + Ki * int_e3;
   
   int_e1=int_e1 + eq1(i)*ts;
   int_e2=int_e2 + eq2(i)*ts;
   int_e3=int_e3 + eq3(i)*ts;
   
   %Kompensator
     %Matriks Centrifugal (B)
     B(1,1) = -2*(m1*s1.^2 + m2*c1.^2)*(L2*c2 + L3*c23)*(L2*s2 + L3*s23);%
     B(1,2) = 2*c1*s1*(I2-I1+(m1-m2)*c23*(L2*c2 + L3*c23));
     B(1,3) = -2*(m1*s1.^2 + m2*c1.^2)*s23*L3*(L2*c2 + L3*c23);%
     B(2,1) = -2*c1*s1*(I2-I1 + (m1-m2)*(L2*s2 + L3*s23).^2 );%
     B(2,2) = 2*((m1*c1.^2 + m2*s1.^2)*L3*c23*(L2*s2 + L3*s23)-m3*L3*s23*(L2*c2 + L3*c23)); 
     B(2,3) = -2*c1*s1*(I2-I1 + (m1-m2)*L3*s23*(L2*s2 + L3*s23));%
     B(3,1) = B(2,3);%
     B(3,2) = 2*(m1*c1^2 + m2*s1^2 - m3)*L3.^2*s23*c23;
     B(3,3) = -2*c1*s1*(I2-I1 + (m1-m2)*L3.^2*s23.^2);
     
     %Matriks Coriolis (C) 
     C(1,1) = c1*s1*(m1-m2)*(L2*c2 + L3*c23).^2;
     C(1,2) = c1*s1*(I2-I1 + (m1-m2)*(L2*c2 + L3*c23));
     C(1,3) = c1*s1*(I2-I1 + (m1-m2)*L3*c23*(L2*c2 + L3*c23));
     C(2,1) = (m1*c1.^2 + m2*s1.^2)*(L2*c2 + L3*c23)*(L2*s2 + L3*s23);
     C(2,2) = (m1*c1.^2 + m2*s1.^2 - m3)*(L2*s2 + L3*s23)*(L2*c2 + L3*c23);
     C(2,3) = (m1*c1.^2 + m2*s1.^2)*L3*c23*(L2*s2 + L3*s23)-m3*L3*s23*(L2*c2 + L3*c23);
     C(3,1) = (m1*c1.^2 + m2*s1.^2)*L3*s23*(L2*c2 + L3*c23); 
     C(3,2) = (m1*c1.^2 + m2*s1.^2)*L3*s23*(L2*c2 + L3*c23)- m3*L3*c23*(L2*s2 + L3*s23); 
     C(3,3) = (m1*c1.^2 + m2*s1.^2 - m3)*L3^2*c23*s23;
     
     %Gravitasi
       N(1,1)=0;
       N(2,1)=m2*g*(L2*c2+L3*cos(q(i,2)-q(i,3)))+m3*g*(L2*c2+L3*c23);
       N(3,1)=-m2*g*L2*cos(q(i,2)-q(i,3))+m3*g*L3*c23;
  
   Lagrange=B*[dq(1,1)*dq(2,1);dq(2,1)*dq(3,1);dq(3,1)*dq(1,1)]+C*[(dq(1,1)).^2;(dq(2,1)).^2;(dq(3,1)).^2];
   
   pwm1(i)=PID1(i)+Lagrange(1,1); %Untuk memutar motornya, kalau disimulasi gak kepake
   pwm2(i)=PID2(i)+Lagrange(2,1); %if pwm>(+-)25 (pwm=(+-)25)
   pwm3(i)=PID3(i)+Lagrange(3,1);
   
   %Step 9 PWM
   q1(i+1)=q1(i)+1*pwm1(i);%dq(:,1);
   q2(i+1)=q2(i)+1*pwm2(i);%dq(:,2);
   q3(i+1)=q3(i)+1*pwm3(i);%dq(:,3);
   
   q1deg(i+1)=q1(i+1)*180/pi;
   q2deg(i+1)=q2(i+1)*180/pi;
   q3deg(i+1)=q3(i+1)*180/pi;
   
   %Last step, code pengembalian ke titik awal fleksi
   if n<=CPekstensi 
       
       if ((abs(ex(i+1))<0.02) && (abs(ey(i+1))<0.02) && (abs(ez(i+1))<0.02))
       %Tes mengganti CPfleksi dengan CPekstensi
       n=n+1;
       %CPekstensi=CPekstensi+1;
       xd(i+1,CPfleksi+2)=xd(i+1,CPfleksi+1);
       yd(i+1,CPfleksi+2)=yd(i+1,CPfleksi+1);
       zd(i+1,CPfleksi+2)=zd(i+1,CPfleksi+1);
       l=1;
       for l=l:(CPfleksi)
           xd(i+1,CPfleksi+2-l)=xd(i+1,CPfleksi+1-l);
           yd(i+1,CPfleksi+2-l)=yd(i+1,CPfleksi+1-l);
           zd(i+1,CPfleksi+2-l)=zd(i+1,CPfleksi+1-l);
       end
       xd(i+1,1)=xd(i+1,CPfleksi+2);
       yd(i+1,1)=yd(i+1,CPfleksi+2);
       zd(i+1,1)=zd(i+1,CPfleksi+2);
       
       h5=plot3(xlink3(i),ylink3(i),zlink3(i),'oy');
      
       end
   end %Trajectory kebawah
    
   %img = screencapture(0, 'Position', [0 0 1270 1023]);
   
end %While 1

%time=time+ts; %Merealistiskan waktu

h2=plot3([xlink2(i) xlink3(i)],[ylink2(i) ylink3(i)],[zlink2(i) zlink3(i)],'b');
h1=plot3([xlink1(i) xlink2(i)],[ylink1(i) ylink2(i)],[zlink1(i) zlink2(i)],'r');
h0=plot3([0 xlink1(i)],[0 ylink1(i)],[0 zlink1(i)],'k');
pause(ts); %default 0.05
h3=plot3(xlink3(i),ylink3(i),zlink3(i),'.g');
h5=plot3(xlink3(i),ylink3(i),zlink3(i),'oy');
hold on;

%legend('Lengan 2', 'Lengan 1', 'Trajectory','Trajectory Fleksi')

w1=['Simulation takes ',num2str(time),' seconds with '];
w2=[num2str(i),' iterations.'];
disp(w1)
disp(w2)