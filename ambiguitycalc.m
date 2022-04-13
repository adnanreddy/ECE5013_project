fs=500e6;
Ts=1/fs;
tau=80e-6;
b=120e6;

fp=1e3;
Tp=1/fp;

Np=64;


%%


nugridsize=64;
% 
% M=length(t)-1;
% Mi=floor(M/taugridsize);
% Mi=0.1*Tp*fs;
% 
% taugridsize=floor(M/Mi);
% 
% Mind=-taugridsize*Mi:Mi:taugridsize*Mi; %centered at 0
% Mind=Mind + M+1; %centered at M+1 <=zero lag

maxlag=floor(10/b*fs);
ambiguity=zeros(2*nugridsize+1,2*maxlag+1);


for i=-nugridsize:1:nugridsize;
    
    fd=  i* 2* 1/tau/nugridsize; % fd =-fp/2  to fp/2
    za_fd=za.*exp(j*2*pi*fd*t);
    
[acor,lag]=xcorr(za_fd,za, maxlag);

%2M+1 points in acor, where M=lenght(t)-1, sample taugridsize*2+1 points

tic;
ambiguity(i+nugridsize+1,:)=acor;
toc
end

nugrid=linspace(-2/tau,2/tau,2*nugridsize+1)/1000;
taugrid=(-maxlag:maxlag)*Ts*1e6;

figure(1); mesh(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')
 
figure(2); contour(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')
%%% Creates a rectangular pulse of width tau from 0<t<tau %%%
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end

