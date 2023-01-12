function plot_FR_edgenetwork_fcn(directory, filename, edge_thresh, node_thresh, scale, cblim)



% 
%  script to plot in MEA layout the graph network with edges above a threshold
% 
% INPUTS:
% directory i.e. the path to the file
% filename i.e. the name (including extension) of the file
% edge_thresh i.e. a value between 0 and 1 for the minimum correlation to
% plot
% scale: optional input argument, omit or set as blank vector [] for actual
% values or choose 'logc' 'logr' 'rate' 'count'
% cblim: colorbar limits two item vector with min and max for colorbar
% scale; if using logc still enter raw limits e.g. [1 100000] 
% if using logr, enter raw rate limits in Hz e.g.  [0.001 10]
% type: choose either 'rate' or 'count'
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










%% set params
textxpos = 9.3;
%%
% close all; clear all
% simulate data
% cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
% load('MPT190403_6B_DIV28_cSpikes_L0_adjM_0.05.mat')
load(filename);
if ~isempty(strfind(filename,'CTRL'))
    adjM = adjM2s(:,:,1); % take first matrix
end
% get spikes
if ~isempty(strfind(filename,'CTRL'))
    spikefile = strcat(filename(1:strfind(filename,'_adjM')),'spikeMat.mat');
    load(spikefile);
    spikeMat = full(rand_spikeMat_sparse_all{1,1});
    if ~exist('fs')
        fs = 1000; 
    end
else
    spikefile = strcat(filename(1:strfind(filename,'_adjM')-1),'.mat');
    load(spikefile);
    try
        spikeMat = full(mSpikes);
    catch
        spikeMat = full(cSpikes);
    end
    if ~exist('fs')
        fs = 25000; 
    end
end

% set minimum and maximum values to plot
Vmin = cblim(1); Vmax = cblim(2);
Dmin = 0; Dmax = 1;
z3 = sum(spikeMat);
if strcmp(scale,'logr') | strcmp(scale,'rate')
    z3 = z3 ./ (length(spikeMat) / fs) ;
end
z3(find(z3>Vmax)) = Vmax;
z3(find(z3<Vmin)) = Vmin;
if ~exist('scale')
    scale = 'count';
end

% z3t = Dmin + (Dmax-Dmin)*((z3-Vmin)/(Vmax-Vmin)) % 

if strcmp(scale,'logc') | strcmp(scale,'logr')
    z3 = log10(z3);
    z3(find(z3 == -Inf)) = log10(Vmin);
% % % % %     take range of values between min and max in cblim
% % % % %     and evenly distribute values to between 0 and 1
    V  = z3;
    z3 = Dmin + (Dmax-Dmin)  *  ((V-log10(Vmin)) / (log10(Vmax)-log10(Vmin))); 
    z3(find(z3>1)) = 1; % set max value to 1 (i.e. anything greater than max colorbar value appears as max colorbar value)
elseif strcmp(scale,'count') | strcmp(scale,'rate')
    V  = z3;
    z3 = Dmin + (Dmax-Dmin)  *  ((V-Vmin) / (Vmax-Vmin)); 
    z3(find(z3>1)) = 1; 
end

% make 0s 0.001 so that when plotting, this is the minimum colour on the
% colorbar; note, an alternative is to keep them as 0s, then when plotting
% use try loop and then catch and plot as grey to indicate that they are
% less than the minimum
z3(find(z3 <= 0)) = 0.0001;
   
if ~exist(strcat(filename(1:end-4),'_FireRateXedgenetwork.png'))
%need to add section that calibrates the z values (node sizes) to values
%between 1 and 25
%need to add section that calbriates the z2 values to values between 0 and
%1; then I need to plot and then do the inverse of this calibration on the
%colorbar and legend e.g. edit cb.TickLabels
if ~exist('channels')
    fprintf(2,'\n  WARNING: channel order not saved in spike matrix \n used default in line 50 of batch_getHeatMaps_fcn \n \n')
    channels = [47,48,46,45,38,37,28,36,27,17,26,16,35,25,15,14,24,34,13,23,12,22,33,21,32,31,44,43,41,42,52,51,53,54,61,62,71,63,72,82,73,83,64,74,84,85,75,65,86,76,87,77,66,78,67,68,55,56,58,57];
end
%% get node degree for each channel
% edge_thresh = 0.7;
count1 = 1; %to track threshold iterations

if strcmp(node_thresh,'average')
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
else
    for cutoff = node_thresh
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
end

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
%get colormap

coords = num2str(channels);
x = str2num(coords(1:end,1));
y = 9 - str2num(coords(1:end,2));

F1 = figure;
F1.OuterPosition = [750   533   456   440];
hold on

%% add edges

max_ew = 1.5; % maximum edge width for plotting
min_ew = 0.5;     % min edge width
linealpha = 0.2; % 0 means completely transparent

for elecA = 1:length(channels)
    for elecB = 1:length(channels)
        if adjM(elecA,elecB) >= edge_thresh & elecA ~= elecB & ~isnan(adjM(elecA,elecB));
            p = plot([x(elecA),x(elecB)],[y(elecA),y(elecB)],'LineWidth',...
                min_ew + (max_ew-min_ew)*((adjM(elecA,elecB)-edge_thresh)/(1-edge_thresh)),'Color',0.5*[1 1 1]);
            p.Color(4) = linealpha;
        end
    end
end

%% add nodes

mycolours = colormap;
legdata = ['01'; '10'; '20']; % data for the legend

% 'MarkerFaceColor',...
%             mycolours(ceil(length(mycolours)*z2(i)),1:3),
sizeScalar = 1.3;
for i = 1:length(channels)
    if z(i)>0
        plot(x(i),y(i),'o','MarkerSize',sizeScalar * (round(z(i)/2)),'MarkerFaceColor',...
            mycolours(ceil(length(mycolours)*z3(i)),1:3),'LineStyle','none','LineWidth',1.2,...
            'MarkerEdgeColor','k') ;
        % where x and y are MEA coords and z is graph metric value
    end
end

set(gca,'color','none');
ylim([0 9]);
xlim([0 9]);

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
set(a, 'XAxisLocation', 'top');
a.XAxis.FontName = 'Arial';
a.XAxis.FontSize = 14;
a.YAxis.FontName = 'Arial';
a.YAxis.FontSize = 14;
% % a.XAxis.Color = [0.5 0.7 0.8]
% % a.YAxis.Color = [0.5 0.7 0.8]
% 

cb = colorbar;
cb.Location = 'eastoutside';
cb.FontName = 'Arial';
cb.FontSize = 12;
cb.FontName = 'Arial';
cb.Label.String = '';
cb.Label.FontSize = 14;
cb.Label.FontName = 'Arial';
cb.Label.FontWeight = 'Normal';
cb.Label.FontSize = 14;
cb.TickLabels = cb.Ticks * cblim(2);
cb.Position = [0.82    0.0595    0.0285    0.2];

if strcmp(scale,'logc') | strcmp(scale,'logr')
    yminim_cbar = floor(log10(cblim(1)));
    ylimit_cbar = log10(cblim(2));    
    cb.Ticks = linspace(0,ylimit_cbar,ylimit_cbar-yminim_cbar+1) ./ ylimit_cbar;%(start,end,number of numbers)
%     numticks = 4; cb.Ticks = linspace(0,1,numticks);%(start,end,number of numbers);
    
    if yminim_cbar >= 0
        roundvalue = - (length(num2str(10^yminim_cbar)) - 1);
    else
        roundvalue = abs(yminim_cbar);
    end
    
    cb.TickLabels = round( 10.^(linspace(yminim_cbar,ylimit_cbar,ylimit_cbar-yminim_cbar+1)) , roundvalue);    
%     cb.TickLabels = round( 10.^(linspace(yminim_cbar,ylimit_cbar,numticks)) , roundvalue);
    
elseif strcmp(scale,'count') | strcmp(scale,'rate')
    ylimit_cbar = cblim(2); yminim_cbar = cblim(1);
    numticks = 4; cb.Ticks = linspace(0,1,numticks);%(start,end,number of numbers);
    ticklabelsvec = linspace(yminim_cbar,ylimit_cbar,numticks);
    ticklabelsvec(find(ticklabelsvec >  1                             )) = round(ticklabelsvec(find(ticklabelsvec >  1                             ))  );
    ticklabelsvec(find(ticklabelsvec <= 1     &  ticklabelsvec > 0.1  )) = round(ticklabelsvec(find(ticklabelsvec <= 1     &  ticklabelsvec > 0.1  )),1);
    ticklabelsvec(find(ticklabelsvec <= 0.1   &  ticklabelsvec > 0.01 )) = round(ticklabelsvec(find(ticklabelsvec <= 0.1   &  ticklabelsvec > 0.01 )),2);
    ticklabelsvec(find(ticklabelsvec <= 0.01                          )) = round(ticklabelsvec(find(ticklabelsvec <= 0.01                          )),3);
    cb.TickLabels = ticklabelsvec;
end

if strcmp(scale,'logr') | strcmp(scale,'rate')
    text(textxpos,2.3,'Spikes/s','FontSize',14);
else
    text(textxpos,2.3,'# Spikes','FontSize',14);
end
cb.TickDirection = 'out';
cb.TickLength = 0.05;

% 
% % need to add colorbar labelling for node color
% 
% %may require having colorbar at the bottom
% % cb.Location = 'Southoutside'
% 
% % need to add legend for node sizes

% ax1 = axes(F1);cbl
% ax2 = copyobj(ax1, gcf);
% delete( get(b,'Children') ) 

l1 = plot(1,-1,'o','MarkerSize',sizeScalar*(round(str2num(legdata(1,:))/2)),'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',1.2);
l2 = plot(2,-1,'o','MarkerSize',sizeScalar*(round(str2num(legdata(2,:))/2)),'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',1.2);
l3 = plot(3,-1,'o','MarkerSize',sizeScalar*(round(str2num(legdata(3,:))/2)),'MarkerEdgeColor','k','MarkerFaceColor',[0.9 0.9 0.9],'LineWidth',1.2);
LX = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w','Visible','off');
MOLE = plot(1,-1,'o','MarkerSize',0.000000001,'MarkerEdgeColor','w','MarkerFaceColor','w','Visible','off');
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

colorbar('Visible','Off'); %create space for legend
legdata_string = ['  1 ' ; ' 10 ' ; ' 20 ';'    ';'    '; num2str(Emin,'%.2f') ; num2str(Emed,'%.2f') ; num2str(Emax,'%.2f')];

l9 = legend([l1 l2 l3 LX MOLE l4 l5 l6],'Orientation','vertical','Position',...
    [0.9 0.6 0 0],'Box','off','String',legdata_string,...
    'FontSize',12);
text(textxpos,8.5,'Degree','FontSize',14);
text(textxpos,5.25,'Weight','FontSize',14);
% l9.Title.String = 'Weight';
% l9.Title.FontSize = 14;
% l9.Title.FontName = 'Arial';
% l9.Title.FontWeight = 'normal';

% for j=1:length(cb.TickLabels)
%     cbtlabels(j) = str2num(cb.TickLabels{j});
% end

% turn axes off 
axis off

saveas(F1,strcat(filename(1:end-4),'_FireRateXedgenetwork_',scale,'.png'));
saveas(F1,strcat(filename(1:end-4),'_FireRateXedgenetwork_',scale,'1','.epsc'));
saveas(F1,strcat(filename(1:end-10),'_FireRateXedgenetwork_',scale,'2'),'epsc');
saveas(F1,strcat(filename(1:end-10),'_FireRateXedgenetwork_',scale,'3'),'epsc2');
saveas(F1,strcat(filename(1:end-10),'_FireRateXedgenetwork_',scale,'4'),'svg');

close(F1)
end
end



