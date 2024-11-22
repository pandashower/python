clearvars;
close all;

radarpower =10*log10(1);
fsrodkowe = 10*log10(38 *10^9);
pasmo = 10*log10(10 *10^6);
zyskmaks = 20;%dBi
f = 10*log10(10);
Tempszumu =10*log10(600);
StratyMocy = 6;%dB
x=-80;
y=110;
vdrona= 16;
RCS = -14;%dBsm
c =10*log10(3*10^8);
lambda = c-fsrodkowe;
k = 10*log10(1.38*10^-23);
polarscatter(0,0)
hold on
for i =1:20
    x_=x+(i-1)*(16/sqrt(2));
    y_=y-(i-1)*(16/sqrt(2));
    [kat,odl] = cart2pol(x_,y_);
    SN = (radarpower+zyskmaks*2+lambda*2+RCS)-(30*log10((4*pi))+(40*log10(sqrt(x_^2+y_^2)))+k+Tempszumu+pasmo+StratyMocy);
    if SN >=3
        polarscatter(kat,odl,"b*")

    end
    if SN < 3
        polarscatter(kat,odl,"r*")
    end

end

hold off
