%% 6.

initial_x  = 25;
initial_y = 25*sind(theta) + [0:13]*(64*Tp*v);


for i = 1:14
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
    
    angles(i) = angle_CPI; ranges_1(i) = range_CPI_1; velocity_1(i) = velocity_CPI_1;
    ranges_2(i) = range_CPI_2; velocity_2(i) = velocity_CPI_2;
end


figure(15); plot(ranges_1); hold on; plot(ranges_2); title('Range'); xlabel('CPI'); ylabel('Meters'); legend('Tx1','Tx2'); saveas(gcf,'range.jpg')
figure(16); plot(velocity_1); hold on; plot(velocity_2); title('Velocity'); xlabel('CPI'); ylabel('Meters per Second'); legend('Tx1','Tx2'); saveas(gcf,'velocity.jpg')
figure(17); scatter(1:14,angles); title('Azimuth Angle'); xlabel('CPI'); ylabel('Degrees'); saveas(gcf,'angle.jpg')