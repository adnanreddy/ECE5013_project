fp = 1e3; Tp = 1/fp;
Np = 32;
fs = 250e6; Ts = 1/fs; 
t = 0:Ts:Np*Tp-Ts;

k = 1.38e-23; T0 = 293; BW = 125e6; F = 10^(8/10);

sigma_n = sqrt(k*T0*BW*F/2);
noise = sigma_n*(randn(1,length(t))) + 1i*(sigma_n*(randn(1,length(t))));

disp(mean(abs(noise).^2));
disp(k*T0*BW*F);

%%
subplot(2,1,1)
plot(abs(conv(h1,conj(fliplr(h1)))))
subplot(2,1,2)
plot(angle(conv(h1,conj(fliplr(h1))))*180/pi)