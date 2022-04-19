%% 6.

CPI = ceil(2*25*tand(10)*fp/v/64);

initial_x  = 25;
initial_y = 25*sind(theta) + [0:13]*(Np*Tp*v)

%%

for i = 1:CPI
    ECE5013_project;
    
    figure(i);
    
    subplot(2,2,1)
    imagesc(taugrid,1:Np,abs(matcharray1));
    xlabel('Range (m)'); ylabel('Pulse No');

    subplot(2,2,2)
    imagesc(taugrid,nugrid*lambda*fp/2, abs(rangedoppler1))
    xlabel('Range (m)'); ylabel('Velocity (m/sec)');

    subplot(2,2,3)
    imagesc(taugrid,1:Np,abs(matcharray2));
    xlabel('Range (m)'); ylabel('Pulse No');

    subplot(2,2,4)
    imagesc(taugrid,nugrid*lambda*fp/2, abs(rangedoppler2))
    xlabel('Range (m)'); ylabel('Velocity (m/sec)');
    
    filename = strcat('ECE_5013_6_',num2str(i),'.jpg');
    saveas(gcf,filename)
end