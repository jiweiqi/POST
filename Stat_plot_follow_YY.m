%% PLOT follow YY
%%  Script for read stat-1Dy data file from NGA Code
%% Weiqi Ji jan 2rd. 2016

clear; clc;close all;
filepath = pwd;
linestyle_list = { '-', '--', ':', '-.' };
fontname = 'Times New Roman';
x_loc_list = {'15', '30', '45' };

% x_loc_list = {'000' ,'001', '002', '003', '015', '030', '045', '060' };
%% Plot Mean Profile
%% Page set
%# centimeters units
X = 21.0;                  %# A3 paper size
Y = 29.7;                  %# A3 paper size
xMargin = 0; %3.175;               %# left/right margins from page borders
yMargin = 0; %2.54;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)

%# create figure/axis
hFig = figure();
%# figure size displayed on screen (50% scaled, but same aspect ratio)
set(hFig, 'Units','centimeters', 'Position',[0 0 X Y])
movegui(hFig, 'center')
% plot([0 1 nan 0 1], [0 1 nan 1 0]), axis tight
% set(gca, 'XTickLabel',[], 'YTickLabel',[], ...
%     'Units','normalized', 'Position',[0 0 1 1])

% %# figure size printed on paper
% set(hFig, 'PaperUnits','centimeters')
% set(hFig, 'PaperSize',[X Y])
% set(hFig, 'PaperPosition',[xMargin yMargin xSize ySize])
% set(hFig, 'PaperOrientation','portrait')
% 
% %# export to PDF and open file
% print -dpdf -r0 out.pdf
% winopen out.pdf

for ix = 1:3
    x_loc = x_loc_list(ix);
    %% EXP
    exp_dir = '\pmD.stat';
    filename_exp = strcat('D', cellstr(x_loc), '.Yfav');
    tmp_exp = importdata( fullfile(filepath,exp_dir,filename_exp{1}), ' ',4);
    data_exp = tmp_exp.data;
    title_list_exp = tmp_exp.colheaders;
%     N_Column_exp = size(title_list_exp,2);
    C1_exp = data_exp( :, 1 ); %array of r/D
    list_variables_exp = [ 2, 4, 16, 18 ];%For mean value
    
    %% LES/PDF
    filename = strcat('stat-1Dy-0', cellstr(x_loc), '.dat');
    tmp = importdata( fullfile(filepath,filename{1}), ' ',1);
    data = tmp.data;
    title = tmp.textdata{1};
    title = strrep(title, '"', []);
    title_list = regexp(title, '\s+', 'split');
    title_list(1) = [];
    N_Column = size(title_list);
    C1 = data( :, 1 ); %array of r/D
    list_variables = [ 13, 15, 18, 19 ];%For mean value
    list_x_range = [5, 6, 8];
    list_y_range = [1, 2200, 0.06, 0.12];
    list_name = {'<\xi>','<T>','<Y_{CO}>','<Y_{CO2}>'};
    for iy=1:4
        subplot( 4,3,(iy-1)*3+ix );
        C2_exp = data_exp(:, list_variables_exp(iy) );
        C2 = data( :, list_variables(iy) );
        h_exp = plot( C1_exp, C2_exp ,'or','LineWidth',2, 'MarkerSize',8,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','r' );
        hold all;
        h_les = plot( C1, C2 ,'-b','LineWidth',2 );
        set(gca,'XMinorTick','on','YMinorTick','on');
        xlim( [0,list_x_range(ix)] );
        set(gca,'XTick',[0:1:list_x_range(ix)])
        ylim( [0,list_y_range(iy)] );
        xlabel( '\rm{r/D}' );
        ylabel( list_name{iy} );
        if iy == 4
            set(gca,'YTick',[0: 0.04: list_y_range(iy)])
            switch(ix)
                case(1)
                    xlabel( {'{r/D}'; '(a) \it{x/D}\rm = 15'}, 'FontName', 'Times New Roman' );
                case 2
                    xlabel( {'{r/D}'; '(b) \it{x/D}\rm = 30'}, 'FontName', 'Times New Roman' );
                case 3
                    xlabel( {'{r/D}'; '(c) \it{x/D}\rm = 45'}, 'FontName', 'Times New Roman');
            end
        end
        if ix == 1 && iy==1
            h_leg = legend( 'Exp', 'LES/PDF');
            legend boxoff;
            set(h_leg, 'FontSize',8, 'FontName', 'Times New Roman')
        end
        hold all;
    end
end

%% Plot RMS Profile
%% Page set
%# centimeters units
X = 21.0;                  %# A3 paper size
Y = 29.7;                  %# A3 paper size
xMargin = 0; %3.175;               %# left/right margins from page borders
yMargin = 0; %2.54;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)

%# create figure/axis
hFig = figure();
%# figure size displayed on screen (50% scaled, but same aspect ratio)
set(hFig, 'Units','centimeters', 'Position',[0 0 X Y])
movegui(hFig, 'center')
for ix = 1:3
    x_loc = x_loc_list(ix);
    %% EXP
    exp_dir = '\pmF.stat'; %chenges according the LFame D/E/F
    filename_exp = strcat('F', cellstr(x_loc), '.Yfav');
    tmp_exp = importdata( fullfile(filepath,exp_dir,filename_exp{1}), ' ',4);
    data_exp = tmp_exp.data;
    title_list_exp = tmp_exp.colheaders;
%     N_Column_exp = size(title_list_exp,2);
    C1_exp = data_exp( :, 1 ); %array of r/D
    list_variables_exp = [ 2, 4, 16, 18 ] + 1;%For rms value
    
    %% LES/PDF
    filename = strcat('stat-1Dy-0', cellstr(x_loc), '.dat');
    tmp = importdata( fullfile(filepath,filename{1}), ' ',1);
    data = tmp.data;
    title = tmp.textdata{1};
    title = strrep(title, '"', []);
    title_list = regexp(title, '\s+', 'split');
    title_list(1) = [];
    N_Column = size(title_list);
    C1 = data( :, 1 ); %array of r/D
    list_variables = [ 34, 36, 39, 40 ];%For rms value
    list_x_range = [5, 6, 8];
    list_y_range = [0.22, 600, 0.03, 0.04];
    list_name = {'<\xi''>','<T''>','<Y_{CO}''>','<Y_{CO2}''>'};
    for iy=1:4
        subplot( 4,3,(iy-1)*3+ix );
        C2_exp = data_exp(:, list_variables_exp(iy) );
        C2 = data( :, list_variables(iy) );
        h_exp = plot( C1_exp, C2_exp ,'or','LineWidth',2, 'MarkerSize',8,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','r' );
        hold all;
        h_les = plot( C1, C2 ,'-b','LineWidth',2 );
        set(gca,'XMinorTick','on','YMinorTick','on');
        xlim( [0,list_x_range(ix)] );
        set(gca,'XTick',[0:1:list_x_range(ix)])
        ylim( [0,list_y_range(iy)] );
        xlabel( '\rm{r/D}' );
        ylabel( list_name{iy} );
        if iy == 4
            set(gca,'YTick',[0: 0.01: list_y_range(iy)])
            switch(ix)
                case(1)
                    xlabel( {'{r/D}'; '(a) \it{x/D}\rm = 15'}, 'FontName', 'Times New Roman' );
                case 2
                    xlabel( {'{r/D}'; '(b) \it{x/D}\rm = 30'}, 'FontName', 'Times New Roman' );
                case 3
                    xlabel( {'{r/D}'; '(c) \it{x/D}\rm = 45'}, 'FontName', 'Times New Roman');
            end
        end
        if ix == 1 && iy==1
            h_leg = legend( 'Exp', 'LES/PDF');
            legend boxoff;
            set(h_leg, 'FontSize',8, 'FontName', 'Times New Roman')
        end
        hold all;
    end
end
