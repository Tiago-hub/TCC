%mapeamento de pólos
%%
clc
clear all
close all
s=tf('s');
d=1:1:5;
theta=deg2rad(15:15:75);
c=1;
c2=0.1;
w=0.1:0.1:100;
lambda=9.5;

%%

Kd=c;
for i=1:length(d)
    G=lambda./(j.*w-d(i)).^2;
    F=-(j.*w-d(i))./G;
    Kp=(2.*c.*d(i)*w+imag(F))./(w);
    Ki=(w.*(real(F)+c.*(w.^2+d(i)^2))+d(i)*imag(F))./(w);
    plot(Ki,Kp,'Color','r');
    hold on
end


%%
Kd=c;
for i=1:length(theta)
    G=lambda./(j.*w.*exp(j.*theta(i))).^2;
    E=-(j.*w.*exp(j.*theta(i)))./G;
    Ki=(imag(E)+c.*w.^2.*sin(2*theta(i)))./w.*cos(theta(i));
    Kp=real(E)+tan(theta(i)).*imag(E)+c.*w.^2;
    plot(Ki,Kp,'Color','b');
end
%legend('d1','d2','d3','d4','d5','t15','t30','t45','t60','t75')
%%
Ki=c;
Kp=(0.6447+0.8495)/2;
Kd=(1.854+2.793)/2;
Tn=0.001;

Gc=Kp+Ki/s+Kd*s/(1+Tn*s);
G=lambda/s^2;
figure()
step(feedback(Gc*G,1))
hold on
stepinfo(feedback(Gc*G,1))
