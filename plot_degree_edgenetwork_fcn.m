function plot_degree_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, NetMet)

%  script to plot in MEA layout the graph network with edges above a threshold
% 
% INPUTS:
% directory i.e. the path to the file
% filename i.e. the name (including extension) of the file
% edge_thresh i.e. a value between 0 and 1 for the minimum correlation to
% plot
% 
% NOTE please ensure adjacency matrix is absolute values only if you are
% using a directional correlation method if not, load your adjacency matrix
% file, use 'abs(adjacencymatrix)' and save
% 
% e.g.    dir =  'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
%         filename =  '200114_FTDOrg_GrpB_3B_Slice11_mSpikes_3_adjM.mat'
%         edge_thresh = 0.803468468435468454354168465435468454686840001
%         plot_degree_edgenetwork_fcn(directory, filename, edge_thresh)
% 
% calculations:
% z - the metric to colour e.g. edge weight
% z2 - the metric to size the circles. e.g. node degree

%%

if ~exist('channels')
    fprintf(2,'\n  WARNING: channel order not saved in spike matrix \n used default in line 50 of batch_getHeatMaps_fcn \n \n')
    channels = [47,48,46,45,38,37,28,36,27,17,26,16,35,25,15,14,24,34,13,23,12,22,33,21,32,31,44,43,41,42,52,51,53,54,61,62,71,63,72,82,73,83,64,74,84,85,75,65,86,76,87,77,66,78,67,68,55,56,58,57];
end

%% get node degree for each channel
% 
% count1 = 1; %to track threshold iterations
% 
% if strcmp(node_thresh,'average')
%     for cutoff = [0.3 0.5 0.7]
%         threshold = cutoff;
%         
%         % % choose either whole mat or only active channels
%         %   edges=adjM1;
%         edges=adjM;
%         edges = edges - eye(size(edges));
%         edges(find(isnan(edges))) = 0;
%         edges(find(edges < threshold)) = 0;
%         edges(find(edges >= threshold))= 1;
%         
%         DegreeVec(:,count1)       = sum(edges);
%         
%         count1 = count1 + 1;
%     end
%     
%     z = round(mean(DegreeVec,2));
% else
%     for cutoff = node_thresh
%         threshold = cutoff;
%         
%         % % choose either whole mat or only active channels
%         %   edges=adjM1;
%         edges=adjM;
%         edges = edges - eye(size(edges));
%         edges(find(isnan(edges))) = 0;
%         edges(find(edges < threshold)) = 0;
%         edges(find(edges >= threshold))= 1;
%         
%         DegreeVec(:,count1)       = sum(edges);
%         
%         count1 = count1 + 1;
%     end
%     
%     z = round(mean(DegreeVec,2));
% end
% 
% %% acquire edge weights
% 
% weights = adjM;
% weights = weights - eye(size(weights));
% weights(find(isnan(weights))) = 0;
% weights(find(weights < 0)) = 0;
% z2 = mean(weights)';
% 
% % % simulate node degree data
% % z = round(abs(rand(60,1)*25));
% % % z(find(z==0))=1;
% % z(find(channels == 15)) = 0;
% 
% % % simulate edge weight data
% % z2 = abs(rand(60,1));
% % z2(find(z2 < 0)) = 0;

%% plot
%get colormap

coords = num2str(channels);
x = str2num(coords(1:end,1));
y = 9 - str2num(coords(1:end,2));

F1 = figure;
F1.OuterPosition = [100   100   600   580];
hold on

%% add edges

max_ew = 1.5; % maximum edge width for plotting
min_ew = 0.5;     % min edge width

for elecA = 1:length(channels)
    for elecB = 1:length(channels)
        if adjM(elecA,elecB) >= edge_thresh & elecA ~= elecB && ~isnan(adjM(elecA,elecB));
            plot([x(elecA),x(elecB)],[y(elecA),y(elecB)],'LineWidth',...
                min_ew + (max_ew-min_ew)*((adjM(elecA,elecB)-edge_thresh)/(1-edge_thresh)),'Color',[0.8 0.8 0.8]);
        end
    end
end

%% add nodes

mycolours = colormap;
legdata = ['01'; '10'; '20']; % data for the legend

for i = 1:length(channels)
    if z(i)>0
        plot(x(i),y(i),'o','MarkerSize',round(z(i)/2),'MarkerFaceColor',...
            [1 0.7 0],'LineStyle','none','LineWidth',1.2,...
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
a.FontSize = 14;



%%


% % % set(a,'xcolor','none')
% % a.XAxis.Label.Color=[0 0 0];
% % a.XAxis.Label.Visible='on';
% %
set(a, 'XAxisLocation', 'top')
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;
% % a.XAxis.Color = [0.5 0.7 0.8]
% % a.YAxis.Color = [0.5 0.7 0.8]
% 
% cb = colorbar;
% cb.FontName = 'Arial';
% cb.FontSize = 14;
% cb.FontName = 'Arial';
% cb.FontSize = 14;
% cb.Label.String = 'Mean edge weight';
% cb.Label.FontSize = 14;
% cb.Label.FontName = 'Arial';
% cb.Label.FontWeight = 'bold';
% cb.Label.FontSize = 14;
% 
% % need to add colorbar labelling for node color
% 
% %may require having colorbar at the bottom
% % cb.Location = 'Southoutside'
% 
% % need to add legend for node sizes

% ax1 = axes(F1);
% ax2 = copyobj(ax1, gcf);
% delete( get(b,'Children') ) 

l1 = plot(1,-1,'o','MarkerSize',round(str2num(legdata(1,:))/2),'MarkerEdgeColor','k','MarkerFaceColor',[1 0.7 0],'LineWidth',1.2);
l2 = plot(2,-1,'o','MarkerSize',round(str2num(legdata(2,:))/2),'MarkerEdgeColor','k','MarkerFaceColor',[1 0.7 0],'LineWidth',1.2);
l3 = plot(3,-1,'o','MarkerSize',round(str2num(legdata(3,:))/2),'MarkerEdgeColor','k','MarkerFaceColor',[1 0.7 0],'LineWidth',1.2);
LX = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
MOLE = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w');
Emin = edge_thresh; %min value an edge can be (for legend)
Emed = edge_thresh+((1-edge_thresh)/2); %middle of min and max
Emax = 1; %max value an edge weight can be
l4 = plot([1 2],[-1 -1],'LineWidth',min_ew + (max_ew-min_ew)*((Emin -edge_thresh)/(1-edge_thresh)),'Color',[0.7 0.7 0.7]);
l5 = plot([1 2],[-1 -1],'LineWidth',min_ew + (max_ew-min_ew)*((Emed -edge_thresh)/(1-edge_thresh)),'Color',[0.7 0.7 0.7]);
l6 = plot([1 2],[-1 -1],'LineWidth',min_ew + (max_ew-min_ew)*((Emax -edge_thresh)/(1-edge_thresh)),'Color',[0.7 0.7 0.7]);

% l7 = legend([l1 l2 l3],'Orientation','vertical','Position',...
%     [0.92 0.8 0 0],'Box','off','String',legdata,...
%     'FontSize',14);
% % l7.Title.String = 'Degree';
% l7.Title.FontSize = 14;
% l7.Title.FontName = 'Arial';
% l7.Title.FontWeight = 'normal';
% % hold off
% l8 = legend([l4 l5 l6],'Orientation','vertical','Position',...
%     [0.92 0.4 0 0],'Box','off','String',legdata,...
%     'FontSize',14);
% % l8.Title.String = 'Weight';
% l8.Title.FontSize = 14;
% l8.Title.FontName = 'Arial';
% l8.Title.FontWeight = 'normal';

colorbar('Visible','Off') %create space for legend
legdata_string = ['  1 ' ; ' 10 ' ; ' 20 ';'    ';'    '; num2str(Emin,'%.2f') ; num2str(Emed,'%.2f') ; num2str(Emax,'%.2f')];

l9 = legend([l1 l2 l3 LX MOLE l4 l5 l6],'Orientation','vertical','Position',...
    [0.9 0.6 0 0],'Box','off','String',legdata_string,...
    'FontSize',12);
text(9.5,8.5,'Degree','FontSize',14);
text(9.5,5.25,'Weight','FontSize',14);
% l9.Title.String = 'Weight';
% l9.Title.FontSize = 14;
% l9.Title.FontName = 'Arial';
% l9.Title.FontWeight = 'normal';

% for j=1:length(cb.TickLabels)
%     cbtlabels(j) = str2num(cb.TickLabels{j});
% end

saveas(F1,strcat(filename(1:end-4),'_degreeXedgenetwork.png'));
close(F1)
end




