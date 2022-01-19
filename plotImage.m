function plotImage(phi,theta,y,plotName)
imagesc(phi, theta, y)
% axis equal
% axis tight
xlabel('\theta (degree)')
xticks([30 90 150])
ylabel('\phi (degree)')
yticks([210 270 330])
set(gca,'FontSize',20)
caxis([-30 40])
% caxis([-150 150])
colorbar
title(plotName)
end

