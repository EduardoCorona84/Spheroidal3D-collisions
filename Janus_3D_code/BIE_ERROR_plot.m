clear;clc;
load('BIE_TEST_DATA');


%compute mean and STD for all things$

for k=1:size(L,3)
 M(k,:)=mean(L(:,:,k));
 Sdev(k,:)=std(L(:,:,k));
end
M=M(:,[3:end])'; Sdev=Sdev(:,[3:end])';
C=linspecer(5);
figure
hold on
for k=1:5
e=errorbar(M(:,k),Sdev(:,k),'s','LineWidth', 2,'color',C(k,:));    
end


set(gca,'FontSize',14);
set(gca, 'FontName', 'Times New Roman')



title('Error Range for BIE test','FontSize',16);
xlabel('k (distance to surface = 10^{-k})','FontSize',16);
ylabel('log-mean error','FontSize',16);

box on
targetdistance= (1:4);
set(gca,'XTick',1:4,'Xticklabel',targetdistance);

legendlabel={'4','8','12','18','24'};
plegend=legend(legendlabel,'Location', 'Eastoutside');
%legend boxoff
title(plegend,'p value:','FontSize',10,'FontName','Times New Roman');
set(plegend, 'FontName','Times New Roman','FontSize',10);


