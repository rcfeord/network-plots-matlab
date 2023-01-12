function plot_rc_bc_fcn(directory, filename,thr_type,thresholds)

%
% script to plot in MEA layout
% rich clubs and betweenness centrality
%
% INPUTS:
% directory i.e. the path to the file
% filename i.e. the name (including extension) of the file
% thr_type i.e. 'absolute' or 'proportional'
% thresholds i.e. the percentiles or the absolute thresholds to average
% over e.g. [60 75 90] averages over 40, 25 and 10 % density
%
% e.g.    dir =  'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
%         filename =  '200114_FTDOrg_GrpB_3B_Slice11_mSpikes_3_adjM.mat'
%         plot_rc_bc_fcn(dir, filename)
%
% calculations:
% z - the metric to colour -- rich club
% z2 - the metric to size the circles. -- betw. centrality

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
if ~exist(strcat(filename(1:end-4),'_RCxBC.png'))
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
    
    %% get RC and BC scores
    count2 = 1;
    
    RC_score = zeros(length(channels),1);
    
    if ~exist('thr_type')
        thr_type = 'proportional';
    end
    if ~exist('thresholds')
        thresholds = [60 75 90];
    end
    
    for cutoff = thresholds;
        
        if strcmp(thr_type,'proportional')
            threshold = prctile(adjM(:), cutoff);
        elseif strcmp(thr_type,'absolute')
            threshold = cutoff;
        end
        
        edges=adjM;
        edges = edges - eye(size(edges));
        edges(find(isnan(edges))) = 0;
        edges(find(edges < threshold)) = 0;
        edges(find(edges >= threshold))= 1;
        DegreeVec(:,count2)       = sum(edges);
        
        
        
        if round(max(sum(edges))/5)*5>0 % check for enough nodes
            % betweenness centrality
            BC_vec(:,count2) = betweenness_bin(edges);
            
            % rich club
            k               =   max(sum(edges));
            [Rc,Nk,Ek]      =   rich_club_bu(edges,k);
            RC              =   max(Rc);
            maxKrand        =   min(find(Rc==RC));
            RC_nodes_vec    =   find(sum(edges) >= maxKrand);
            RC_node_IDs     =   channels(RC_nodes_vec);
            
            RC_score(RC_nodes_vec)      =    RC_score(RC_nodes_vec)+1;
            
        else
            BC_vec(:,count2) = betweenness_bin(edges);            
        end
        
        count2 = count2 + 1;
        
    end
    
    if ~~exist('BC_vec')
        dg = round(mean(DegreeVec,2));
        z2 = 10*(mean(BC_vec,2)/(59*58)); % z2 is betweenness centrality as a proportion of shortest paths passing through
        z = RC_score/(count2-1); % z is the proportion of times node was in rich club
        
        %% plot
        
        % hist
        hd1 = round(mean(BC_vec,2));
        F2 = figure;
        F2.OuterPosition = [750   533   400   400];
        histogram(hd1(find(RC_score==0)),'BinWidth',20);
        hold on
        histogram(hd1(find(RC_score>0)),'BinWidth',20);
        set(gca,'color','none');
        box off
        xlabel('count of shortest paths')
        ylabel('frequency')
        legend('Not in rich club','In rich club')
        set(gca,'FontName','Arial');
        set(gca,'FontSize',14);
        saveas(F2,strcat(filename(1:end-4),'_histBC.png'));
        close(F2)
        
        %option 1: create colour vector that colours red or blue depending if in
        %rich club; requires adding a legend
        %option 2: colour according to another metric; requires adding a colorbar
        
        %get colormap
        
        coords = num2str(channels);
        x = str2num(coords(1:end,1));
        y = 9 - str2num(coords(1:end,2));
        z2(find(z2==0))=0.000001;
        z2(find(z2>1))=1;
        
        F1 = figure;
        F1.OuterPosition = [750   533   456   467];
        hold on
        
        mycolours = colormap;
        % legdata = ['01'; '10'; '20'];
        l1 = plot(1,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[.333 0 0]);
        l2 = plot(2,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[.666 0 0]);
        l3 = plot(3,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
        
        for i = 1:length(channels)
            if dg(i) > 0
                if z(i)>0
                    plot(x(i),y(i),'o','MarkerSize',14,'MarkerFaceColor',...
                        mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
                        'MarkerEdgeColor',[z(i) 0 0],'LineWidth',4) ;
                    % where x and y are MEA coords and z is graph metric value
                else
                    plot(x(i),y(i),'o','MarkerSize',10,'MarkerFaceColor',...
                        mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
                        'MarkerEdgeColor',[0 0 0],'LineWidth',1) ;
                end
            elseif channels(i) == 15
                % dont plot ref
            else
                plot(x(i),y(i),'o','MarkerSize',1,'MarkerFaceColor',...
                    'k','LineStyle','none',...
                    'MarkerEdgeColor',[0 0 0],'LineWidth',1) ;
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
        cb.Label.String = '% shortest paths';
        cb.Label.FontSize = 14;
        cb.Label.FontName = 'Arial';
        cb.Label.FontWeight = 'bold';
        cb.Label.FontSize = 14;
        
        % need to add colorbar labelling for node color
        
        %may require having colorbar at the bottom
        % cb.Location = 'Southoutside'
        
        % need to add legend for node sizes
        
        l4 = legend([l1 l2 l3],'Orientation','horizontal','Position',...
            [0.4 0.07 0 0],'Box','off','String',['.33';'.66';'1.0'],...
            'FontSize',14);
        l4.Title.String = 'rich club status (proportion)';
        l4.Title.FontSize = 12;
        l4.Title.FontName = 'Arial';
        l4.Title.FontWeight = 'bold';
        
        for j=1:length(cb.TickLabels)
            cbtlabels(j) = str2num(cb.TickLabels{j});
        end
        cbtlabels = cbtlabels*10;
        cb.TickLabels = cbtlabels;
        %% save
        saveas(F1,strcat(filename(1:end-4),'_RCxBC.png'));
        close(F1)
    else
        disp(strcat('max # edges = ',num2str(max(sum(edges)))))
        
        dg = round(mean(DegreeVec,2));
        z2 = zeros(length(channels),1); % z2 is betweenness centrality as a proportion of shortest paths passing through
        z = zeros(length(channels),1); % z is the proportion of times node was in rich club
        
        %% plot
        
        % hist
        hd1 = round(mean(BC_vec,2));
        F2 = figure;
        F2.OuterPosition = [750   533   400   400];
        histogram(hd1(find(RC_score==0)),'BinWidth',20);
        hold on
        histogram(hd1(find(RC_score>0)),'BinWidth',20);
        set(gca,'color','none');
        box off
        xlabel('count of shortest paths')
        ylabel('frequency')
        legend('Not in rich club','In rich club')
        set(gca,'FontName','Arial');
        set(gca,'FontSize',14);
        saveas(F2,strcat(filename(1:end-4),'_histBC.png'));
        close(F2)
        
        %option 1: create colour vector that colours red or blue depending if in
        %rich club; requires adding a legend
        %option 2: colour according to another metric; requires adding a colorbar
        
        %get colormap
        
        coords = num2str(channels);
        x = str2num(coords(1:end,1));
        y = 9 - str2num(coords(1:end,2));
        z2(find(z2==0))=0.000001;
        z2(find(z2>1))=1;

        F1.OuterPosition = [750   533   456   467];
        hold on
        
        mycolours = colormap;
        % legdata = ['01'; '10'; '20'];
        l1 = plot(1,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[.333 0 0]);
        l2 = plot(2,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[.666 0 0]);
        l3 = plot(3,-1,'o','MarkerSize',14,'LineWidth',4,'MarkerEdgeColor',[1 0 0]);
        
        for i = 1:length(channels)
            if dg(i) > 0
                if z(i)>0
                    plot(x(i),y(i),'o','MarkerSize',14,'MarkerFaceColor',...
                        mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
                        'MarkerEdgeColor',[z(i) 0 0],'LineWidth',4) ;
                    % where x and y are MEA coords and z is graph metric value
                else
                    plot(x(i),y(i),'o','MarkerSize',10,'MarkerFaceColor',...
                        mycolours(ceil(length(mycolours)*z2(i)),1:3),'LineStyle','none',...
                        'MarkerEdgeColor',[0 0 0],'LineWidth',1) ;
                end
            elseif channels(i) == 15
                % dont plot ref
            else
                plot(x(i),y(i),'o','MarkerSize',1,'MarkerFaceColor',...
                    'k','LineStyle','none',...
                    'MarkerEdgeColor',[0 0 0],'LineWidth',1) ;
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
        cb.Label.String = '% shortest paths';
        cb.Label.FontSize = 14;
        cb.Label.FontName = 'Arial';
        cb.Label.FontWeight = 'bold';
        cb.Label.FontSize = 14;
        
        % need to add colorbar labelling for node color
        
        %may require having colorbar at the bottom
        % cb.Location = 'Southoutside'
        
        % need to add legend for node sizes
        
        l4 = legend([l1 l2 l3],'Orientation','horizontal','Position',...
            [0.4 0.07 0 0],'Box','off','String',['.33';'.66';'1.0'],...
            'FontSize',14);
        l4.Title.String = 'rich club status (proportion)';
        l4.Title.FontSize = 12;
        l4.Title.FontName = 'Arial';
        l4.Title.FontWeight = 'bold';
        
        for j=1:length(cb.TickLabels)
            cbtlabels(j) = str2num(cb.TickLabels{j});
        end
        cbtlabels = cbtlabels*10;
        cb.TickLabels = cbtlabels;
        %% save
        saveas(F1,strcat(filename(1:end-4),'_RCxBC.png'));
        close(F1)

    end
end
end



