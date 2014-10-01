
H=figure; hold on;
set(H,'Position',[0 0 1000 300]);
plot(R,'b-o','LineWidth',2);
plot(E,'g-.o','LineWidth',2);
plot(S,'r:o','LineWidth',2);

legend('Rotation','Expansion','Shear');
xlabel('Time Step','FontSize',16);
ylabel('Value','FontSize',16);
ylim([0 4.5]);