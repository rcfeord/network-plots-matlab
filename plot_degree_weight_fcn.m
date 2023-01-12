function plot_degree_weight_fcn(directory, filename)

% 
% script to plot in MEA layout the graph metrics degree and weight
% 
% INPUTS:
% directory i.e. the path to the file
% filename i.e. the name (including extension) of the file
% 
% e.g.    dir =  'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
%         filename =  '200114_FTDOrg_GrpB_3B_Slice11_mSpikes_3_adjM.mat'
%         plot_degree_weight_fcn(dir, filename)
% 
% calculations:
% z - the metric to colour e.g. edge weight
% z2 - the metric to size the circles. e.g. node degree

%{

future updates:
(need to make a proportion of colour thn relabel
colorbar (e.g. mean connection weight)

have the MarkerEdgeColor red if the node is in the rich club
add the edges above the minimum threshold (currently 0.3) plotted with the
LineColor based on the edge weight from grey to black

add in other values e.g. a . marker coloured according to the betweenness
centrality value 

%}
% close all; clear all
% simulate data
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
% load('MPT190403_6B_DIV28_cSpikes_L0.mat','channels')
load(filename);
if ~exist(strcat(filename(1:end-4),'_degreeXweight.png'))
%need to add section that calibrates the z values (node sizes) to values
%between 1 and 25
%need to add section that calbriates the z2 values to values between 0 and
%1; then I need to plot and then do the inverse of this calibration on the
%colorbar and legend e.g. edit cb.TickLabels
if ~exist('channels')
    fprintf(2,'\n  WARNING: channel order not saved in spike matrix \n used default in line 50 of batch_getHeatMaps_fcn \n \n')
    channels = [47,48,46,45,38,37,28,36,27,17,26,16,35,25,15,14,24,34,13,23,12,22,33,21,32,31,44,43,41,42,52,51,53,54,61,62,71,63,72,82,73,83,64,74,84,85,75,65,86,76,87,77,66,78,67,68,55,56,58,57];
end

if ~exist('adjM') & contains(filename,'CTRL')
    adjM = adjM2s(:,:,1);
end
%% get node degree for each channel

count1 = 1; %to track threshold iterations
for cutoff = [0.3 0.5 0.7]
    
    threshold = cutoff;
    
    % % choose either whole mat or only active channels
    %   edges=adjM1;
    edges=adjM;
    edges = edges - eye(size(edges));
    edges(find(isnan(edges))) = 0;
    edges(find(edges < threshold)) = 0;
    edges(find(edges >= threshold))= 1;
    
    DegreeVec(:,count1)       = sum(edges);
    
    count1 = count1 + 1;
end

z = round(mean(DegreeVec,2));

%% acquire edge weights

weights = adjM;
weights = weights - eye(size(weights));
weights(find(isnan(weights))) = 0;
weights(find(weights < 0)) = 0;
z2 = mean(weights)';

% % simulate node degree data
% z = round(abs(rand(60,1)*25));
% % z(find(z==0))=1;
% z(find(channels == 15)) = 0;

% % simulate edge weight data
% z2 = abs(rand(60,1));
% z2(find(z2 < 0)) = 0;

%% plot

% hist
 hd1 = round(z);
 F2 = figure;
 F2.OuterPosition = [750   533   400   400];
 histogram(hd1,'BinWidth',1);
 set(gca,'color','none');
 box off
 xlabel('node degree')
 ylabel('Frequency')
 set(gca,'FontName','Arial');
 set(gca,'FontSize',14);
 saveas(F2,strcat(filename(1:end-4),'_hist_degree.png'));
close(F2)

 hd1 = z2;
 F3 = figure;
 F3.OuterPosition = [750   533   400   400];
 histogram(hd1,'BinWidth',0.1);
 set(gca,'color','none');
 box off
 xlabel('edge weight')
 ylabel('Frequency')
 set(gca,'FontName','Arial');
 set(gca,'FontSize',14);
 saveas(F3,strcat(filename(1:end-4),'_hist_edgeweight.png'));
close(F3)

%option 1: create colour vector that colours red or blue depending if in
%rich club; requires adding a legend
%option 2: colour according to another metric; requires adding a colorbar

%get colormap

coords = num2str(channels);
x = str2num(coords(1:end,1));
y = 9 - str2num(coords(1:end,2));

F1 = figure;
F1.OuterPosition = [750   533   456   467];
hold on

mycolours = colormap;
legdata = ['01'; '10'; '20'];
l1 = plot(1,-1,'o','MarkerSize',round(str2num(legdata(1,:))/2),'MarkerEdgeColor','k');
l2 = plot(2,-1,'o','MarkerSize',round(str2num(legdata(2,:))/2),'MarkerEdgeColor','k');
l3 = plot(3,-1,'o','MarkerSize',round(str2num(legdata(3,:))/2),'MarkerEdgeColor','k');

for i = 1:length(channels)
    if z(i)>0
        plot(x(i),y(i),'o','MarkerSize',round(z(i)/2),'MarkerFaceColor',...
            mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
            'MarkerEdgeColor','k') ;
        % where x and y are MEA coords and z is graph metric value
    end
end

set(gca,'color','none')
ylim([0 9])
xlim([0 9])

a = gca;
xticks(1:8)
yticks(1:8)
yvalues = a.YTickLabels ;
yvalues = sort(str2num(cell2mat(yvalues)),'descend');
a.YTickLabels = yvalues;


% % set(a,'xcolor','none')
% a.XAxis.Label.Color=[0 0 0];
% a.XAxis.Label.Visible='on';
%
set(a, 'XAxisLocation', 'top')
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;
% a.XAxis.Color = [0.5 0.7 0.8]
% a.YAxis.Color = [0.5 0.7 0.8]

cb = colorbar;
cb.FontName = 'Arial';
cb.FontSize = 14;
cb.FontName = 'Arial';
cb.FontSize = 14;
cb.Label.String = 'Mean edge weight';
cb.Label.FontSize = 14;
cb.Label.FontName = 'Arial';
cb.Label.FontWeight = 'bold';
cb.Label.FontSize = 14;

% need to add colorbar labelling for node color

%may require having colorbar at the bottom
% cb.Location = 'Southoutside'

% need to add legend for node sizes

l4 = legend([l1 l2 l3],'Orientation','horizontal','Position',...
    [0.3 0.07 0 0],'Box','off','String',legdata,...
    'FontSize',14);
l4.Title.String = 'Mean node degree';
l4.Title.FontSize = 14;
l4.Title.FontName = 'Arial';
l4.Title.FontWeight = 'bold';

for j=1:length(cb.TickLabels)
    cbtlabels(j) = str2num(cb.TickLabels{j});
end

saveas(F1,strcat(filename(1:end-4),'_degreeXweight.png'));
close(F1)
end
end



