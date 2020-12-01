%CP0 to CP1 to CP2 to CP3 to CPn DONE

clear;clc;close all;

dx(1)=0.15;

%xd(1,1)=5;

%Kp=0.45;
%Ki=4.86;

Kp=0.15;
Ki=0.1;

int_e=0;
ts=0.05;

xact(1)=0.1;
i=0;
%desired=[cp1 cp2 cp3 ... cpn dummyvalue dummyvalue];
xd=[0.5 1 0.1 0.8 0.3 0 0]; %Nilai 0 agar stop(?)

n=1;
cp=length(xd) - 2;

while abs(dx)>0.01
    i=i+1;
    time(i)=ts*i;
    
    dx(i+1)= xd(i,1) - xact(i);
    
    xact(i+1)= Kp*dx(i+1) + Ki*int_e;
    int_e= int_e + dx(i+1)*ts;
    
    for j=1:(cp+2)
        xd(i+1,j)=xd(i,j);
    end
    %{
    xd(i+1,2)=xd(i,2);
    xd(i+1,3)=xd(i,3);
    xd(i+1,4)=xd(i,4);
    xd(i+1,5)=xd(i,5);
    
    %Membentuk matriks baru agar bisa ditempati
    xd(i+1,1)=xd(i,1);
    %}
    
   % dx(i+1)=dx(i);%nilai dxnya tetap sama terus
    %x(i+1)=x(i);
    
    if n<cp      
        if abs(dx(i+1))<0.03     
        n=n+1;
        %cp=3 maka:
        %xd(i+1,4)=xd(i+1,1); 
        %xd(i+1,1)=xd(i+1,2); 
        %xd(i+1,2)=xd(i+1,3); 
        %xd(i+1,3)=xd(i+1,5); 
        
        %cp=4, maka
        %xd(i+1,5)=xd(i+1,1); 
        %xd(i+1,1)=xd(i+1,2); 
        %xd(i+1,2)=xd(i+1,3); 
        %xd(i+1,3)=xd(i+1,4);
        %xd(i+1,4)=xd(i+1,6);
        
        %cp=5, maka
        %xd(i+1,6)=xd(i+1,1); 
        %xd(i+1,1)=xd(i+1,2); 
        %xd(i+1,2)=xd(i+1,3); 
        %xd(i+1,3)=xd(i+1,4);
        %xd(i+1,4)=xd(i+1,5);
        %xd(i+1,5)=xd(i+1,7);
        
        %Berhasil
        %cp=p, maka
        xd(i+1,cp+1)=xd(i+1,1); %bagian terpenting 
        for v=1:(cp-1)
            xd(i+1,v)=xd(i+1,v+1);
        end
        xd(i+1,cp)=xd(i+1,cp+2);
        
        %Digeneralkan dengan for
        %xd(i+1,1)=xd*i+1,2);
        %...
        %xd(i+1,p)=xd(i+1,p+2); 
        %xd(i+1,cp+1)=xd(i+1,1);
        %xd(i+1)
            
        
        end 
    end
    
    grid on;
    box off;
    h1=plot(time(i),xact(i),'.g');
    pause(ts) %wajib ada kalo mau gerak
    hold on; %wajib ada kalo mau gerak
    %legend('Value of x')
    xlabel('Time [seconds]')
    ylabel('x')
end

w1=['Simulation takes ',num2str(time(i)),' seconds with '];
w2=[num2str(i),' iterations.'];
disp(w1)
disp(w2)
