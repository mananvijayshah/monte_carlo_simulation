close all
clear all
clc
interval=10;
maxtime =60;
dp=520e-6;

figure(1)
name='CVMC_Real_SMP_final-11-Jul-2023-reportA-';
for i=1:maxtime
    data.filename{i}=[name num2str(i) '.mat'];
    load (data.filename{i})
    entities=max(size(reportA));
    part_total=sum(reportA(:,2));
    [a, b]=function_diameter(reportA(:,4));
    d32(i)=a;
    d50(i)=b;
    dmean(i)=mean(reportA(:,4));
    Np_rel=entities/part_total;
    result(i)=Np_rel;
end
x=interval:interval:(maxtime*interval);
plot(x./60,d32./dp,'b*-','linewidth',0.8)
hold on
GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
GR_d32 = (d32(maxtime)-dp)/(maxtime*10)

%  name='CVMC_real_Df_2_-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'go-','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% 
% name='CVMC_real_Df_2.5_-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'rs-','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% 
% name='CVMC_real_Df_3_-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'md-','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% %%%CoP
% name='CVMC_Cop_Df_1.5-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'b*','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% 
% name='CVMC_Cop_Df_2-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'go','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% % 
% % 
% name='CVMC_Cop_Df_2.5-03-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'rs','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% % 
% name='CVMC_Cop_Df_3-04-May-2022-reportA-';
% for i=1:maxtime
%     data.filename{i}=[name num2str(i) '.mat'];
%     load (data.filename{i})
%     entities=max(size(reportA));
%     part_total=sum(reportA(:,2));
%     [a, b]=function_diameter(reportA(:,4));
%     d32(i)=a;
%     d50(i)=b;
%     dmean(i)=mean(reportA(:,4));
%     Np_rel=entities/part_total;
%     result(i)=Np_rel;
% end
% x=interval:interval:(maxtime*interval);
% plot(x./60,d32./dp,'md','linewidth',0.8)
% hold on
% GR_dmean = (dmean(maxtime)-dp)/(maxtime*10);
% GR_d32 = (d32(maxtime)-dp)/(maxtime*10)
% %
% 
% 
% % ylim([1 5])
% % xlim([0 10])
% x=xlabel('{t} [min]');
% % set(x,'Interpreter','latex','FontSize',14);
% y=ylabel('D_{agg,32}/D_{p} [-]');
% % set(y,'Interpreter','latex','FontSize',14);
% h=legend('CVMC-1 ($D_{f}=1.5$)','CVMC-1 ($D_{f}=2$)','CVMC-1 ($D_{f}=2.5$)','CVMC-1 ($D_{f}=3$)','Old model ($D_{f}=1.5$)','Old model ($D_{f}=2$)','Old model ($D_{f}=2.5$)','Old model ($D_{f}=3$)');
% set(h,'Interpreter','latex','FontSize',12,'location', 'Northwest');
% legend boxon
% set(gca,'Box','on','YTick',[0:2:8],'FontSize',12)
% set(gca,'Box','on','XTick',[0:2:10],'FontSize',12)
% grid off
% set(gcf,'Color','white')  %Changing size and color
% set(gcf, 'Position', get(0, 'Screensize'));
% title('Dagg vs Time Cop vs Real CVMC xb=2%')

% % print('CVMC_PTSA_kinetics_Stdef=2Stcoal_Dagg_Asurf','-dpng','-r300')
