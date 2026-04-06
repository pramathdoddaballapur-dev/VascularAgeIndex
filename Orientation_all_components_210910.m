%% Primary orientation and alignment estimation from 2PM images
%  Aggrion: 10 sept 2021 / Separation adv and medial components
clear all; clc; close all;

% Path
% foldername = 'G:\Ericf5_Projects\Training\Training_rotation\Orientation\';
foldername = 'F:\_F24mo\20200918_RPA_F23_1_0\15_0/';
% code will only process one stack at a time
load([foldername 'Thickness\15_0' '_thickness.mat'], 'thickness');

stack = [15];
for j =1:length(stack)
    
    eval(['seriename = ''OJ_', num2str(stack(j)),'_C3'';']); %OJ_#stack#_C#channel# ch1 and 3 not 2
%     seriename = 'OJ_1_C1'; %OJ_#stack#_C#channel#

    % Stack info (optional)
    p = 'XX'; %pressure mmHg
    lambda = 'XX'; % 2 = in vivo, 1 = -5%, 3 = +5%
    region = 'XX';
    age = 'XX';
    sample = '1';

    ch = str2double(seriename(end)); %change channel name in both

    if ch == 1  subj = 'Collagen'; elseif ch == 2  subj = 'Elastin'; elseif ch == 3  subj = 'Cell nuclei'; end

    filename1 = [foldername 'Orientation\' seriename '.txt'];
    delimiter = ' ';
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename1,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);
    fclose(fileID);

    nslice = size(dataArray{:,1},1)/180;


    x = dataArray{:, 1}(1:180);
    Y = zeros(180,nslice); Y_gr = zeros(180,10);


    %run from here when you clean Y. CLOSE THE FIGURES
    close all
    k = 1; gr = ceil(nslice/10);
    lgd = cell(gr,1);
    f1 = figure(1); cc = jet(gr); hold on
    for i = 1 : nslice
        pointer = 180*(i-1);
        Y(:,i) = dataArray{:, 2}((1+pointer):(180+pointer));  %comment if you  clean Y
        eval(['Y_gr(:,i-(10*(k-1))) = Y(:,i);']);
        % groups of 10 slices to get suborientation
        if mod(i,10) == 0
            eval(['y_mean_' num2str(k) ' = mean(Y_gr,2);']);
            eval(['plot(x,y_mean_' num2str(k) ', ''LineWidth'',2, ''color'',cc(k,:));']);
            lgd{k}=strcat(num2str(i-10),'-', num2str(i));
            k = k + 1;
            Y_gr = zeros(180,10); %titl
        elseif i == nslice
            eval(['y_mean_' num2str(k) ' = mean(Y_gr(:,1:nslice-(10*(k-1))),2);']);
            eval(['plot(x,y_mean_' num2str(k) ', ''LineWidth'',2, ''color'',cc(k,:));']);
            lgd{k}=strcat(num2str(k - 1),'0-', num2str(i));
        end
    end

    legend(lgd, 'Location', 'eastoutside'), eval(['title(''Layer based orientation (groups of 10 slides) - ' subj ''')']);
    xlabel('Orientation in degrees'), ylabel('Incidence (px)')


    % all stack orientation
    y_sum = sum(Y,2);

    % normalized stack orientation (rad)
    alpha = x*pi/180;
    y_normalized = y_sum ./ trapz(alpha,y_sum);


    f3 = figure(3), imagesc(Y), colormap('jet'); c = colorbar;
    yticks([1 90 180])
    yticklabels({'-90', '0','90'})
    eval(['title(''Layer based orientation - ' subj ''')']);
    xlabel('Depth (um)'), ylabel('Orientation in degrees')


    answer = questdlg('Do you want to delete external slices?','Yes','No');
    switch answer
        case 'Yes'
            try
                fun = @(x) thickness(x).stack == stack(j) % useful for complicated fields
                index = find(arrayfun(fun, 1:numel(thickness)))
                bottom = num2str(ceil(mean(thickness(index).fin)))
                top = num2str(floor(mean(thickness(index).in)))
                interface = num2str(floor(thickness(index).z_adv_med))

                definput = {top,interface,bottom};
            catch
                definput = {'1','50','100'};
            end
            prompt = {'Enter min limit:','Enter interface:','Enter max limit:'};
            dlgtitle = 'Input';
            numlines = 2;
            answer = inputdlg(prompt,dlgtitle,numlines,definput);

            % answer = inputdlg(prompt,dlgtitle,[1 40],definput);
            Y_compl = Y;
            if ch == 1
                Y(:,1:str2num(answer{1})) = 0;  Y(:,str2num(answer{2}):end) = 0;
                Y_compl(:,1:(str2num(answer{2})+3)) = 0;  Y_compl(:,str2num(answer{3}):end) = 0; % complementary
            else
                Y(:,1:str2num(answer{2})) = 0;  Y(:,str2num(answer{3}):end) = 0;
                Y_compl(:,1:str2num(answer{1})) = 0;  Y_compl(:,(str2num(answer{2})-3):end) = 0; % complementary
            end
            %run again
            close all
            k = 1;
            f1 = figure(1); cc = jet(gr); hold on
            %         if ch == 2,    set(gca, 'YScale', 'log'), end
            for i = 1 : nslice
                pointer = 180*(i-1);
                eval(['Y_gr(:,i-(10*(k-1))) = Y(:,i);']);
                eval(['Y_compl_gr(:,i-(10*(k-1))) = Y_compl(:,i);']);
                % groups of 10 slices to get suborientation
                if mod(i,10) == 0
                    eval(['y_compl_mean_' num2str(k) ' = mean(Y_compl_gr,2);']);
                    eval(['y_mean_' num2str(k) ' = mean(Y_gr,2);']);
                    eval(['plot(x,y_mean_' num2str(k) ', ''LineWidth'',2, ''color'',cc(k,:));']);
                    eval(['plot(x,y_compl_mean_' num2str(k) ', ''LineStyle'','':'', ''LineWidth'',2, ''color'',cc(k,:));']);
                    lgd{k}=strcat(num2str(i-10),'-', num2str(i));
                    k = k + 1;
                    Y_gr = zeros(180,10); %titl
                    Y_compl_gr = zeros(180,10);
                elseif i == nslice
                    eval(['y_mean_' num2str(k) ' = mean(Y_gr(:,1:nslice-(10*(k-1))),2);']);
                    eval(['plot(x,y_mean_' num2str(k) ', ''LineWidth'',2, ''color'',cc(k,:));']);
                    eval(['y_compl_mean_' num2str(k) ' = mean(Y_compl_gr(:,1:nslice-(10*(k-1))),2);']);
                    eval(['plot(x,y_compl_mean_' num2str(k) ', ''LineStyle'','':'', ''LineWidth'',2, ''color'',cc(k,:));']);
                    lgd{k}=strcat(num2str(k - 1),'0-', num2str(i));
                end
            end
            lgd2= repelem(lgd,2);
            legend(lgd2, 'Location', 'eastoutside'), eval(['title(''Layer based orientation (groups of 10 slides) - ' subj ''')']);
            xlabel('Orientation in degrees'), ylabel('Incidence (px)')


            % normalized stack orientation (degrees)
            y_sum = sum(Y,2);
            y_norm = y_sum./sum(y_sum);
            y_compl_sum = sum(Y_compl,2);
            y_compl_norm = y_compl_sum./sum(y_compl_sum);
            % normalized stack orientation (rad)
            alpha = x * pi /180;
            y_normalized = y_sum ./ trapz(alpha,y_sum);
            y_compl_normalized = y_compl_sum ./ trapz(alpha,y_compl_sum);

            f3 = figure(3),
            imagesc(Y), colormap('jet');
            yticks([1 90 180])
            yticklabels({'-90', '0','90'})
            eval(['title(''' subj '- ' age ' #' sample ', p = ' p ' mmHg, \lambda_z^{iv}, Dorsal'')']);
            xlabel('Depth (\mum)'), ylabel('Orientation (Degrees)')
            line([str2num(answer{1}) str2num(answer{1})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{3}) str2num(answer{3})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{2}) str2num(answer{2})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle',':'); %plot for each distribution
            %
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            set(findall(gcf,'-property','FontWeight'),'FontWeight','b') 
            
            f4 = figure(4),
            imagesc(Y_compl), colormap('jet');
            yticks([1 90 180])
            yticklabels({'-90', '0','90'})
            eval(['title(''' subj 'complem. - ' age ' #' sample ', p = ' p ' mmHg, \lambda_z^{iv}, Dorsal'')']);
            xlabel('Depth (\mum)'), ylabel('Orientation (Degrees)')
            line([str2num(answer{1}) str2num(answer{1})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{3}) str2num(answer{3})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{2}) str2num(answer{2})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle',':'); %plot for each distribution
            %
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            set(findall(gcf,'-property','FontWeight'),'FontWeight','b') 

            %each layer normalized
            f5 = figure(10),
            Y_normalized = Y ./ trapz(x,Y);

            imagesc(Y_normalized), colormap('jet');
            yticks([1 90 180])
            yticklabels({'-90', '0','90'})
            eval(['title(''' subj '- ' age ' #' sample ', p = ' p ' mmHg, \lambda_z^{iv}, Dorsal'')']);
            xlabel('Depth (\mum)'), ylabel('Orientation (Degrees)')
            %%line([27 27],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle',':'); %plot for each distribution
            line([str2num(answer{1}) str2num(answer{1})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{3}) str2num(answer{3})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution
            line([str2num(answer{2}) str2num(answer{2})],  [0,180],'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle',':'); %plot for each distribution
            %
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            set(findall(gcf,'-property','FontWeight'),'FontWeight','b')


        case 'No'
            %nothing
 
    end

    % Fit mean orientation with Von Mises distribution
%     if ch==3
%         tmp = y_normalized;
%         tmp_compl = y_compl_normalized;
%         y_normalized(1:90) = tmp(91:180); y_normalized(91:180) = tmp(1:90);
%         y_compl_normalized(1:90) = tmp_compl(91:180); y_compl_normalized(91:180) = tmp_compl(1:90);
%         x_ = 0:179; alpha = x_*pi/180;
% 
%     
%     [thetahat kappa] = circ_vmpar(alpha,y_normalized);
%     [thetahat_compl kappa_compl] = circ_vmpar(alpha,y_compl_normalized);
%     
%     [p alpha] = circ_vmpdf(alpha, thetahat, kappa);
%     [p_compl alpha] = circ_vmpdf(alpha, thetahat_compl, kappa_compl);
%     
%     y_normalized = y_normalized/trapz(alpha,y_normalized);
%     y_compl_normalized = y_compl_normalized/trapz(alpha,y_compl_normalized);
%     end
   
    tmp = y_normalized;
    tmp_compl = y_compl_normalized;
    for rot = 0:89
        y_normalized = circshift(tmp,rot);
        y_compl_normalized = circshift(tmp_compl,rot);
        [thetahat kappa] = circ_vmpar(alpha,y_normalized);
        [thetahat_compl kappa_compl] = circ_vmpar(alpha,y_compl_normalized);
        [p alpha] = circ_vmpdf(alpha, thetahat, kappa);
        [p_compl alpha] = circ_vmpdf(alpha, thetahat_compl, kappa_compl);
        err(rot+1) = sqrt(sum(y_normalized-p).^2);
        err_compl(rot+1) = sqrt(sum(y_compl_normalized-p_compl).^2);
    end
    [min_err min_id] = min(err);
    [min_err_compl min_id_compl] = min(err_compl);
    y_normalized = circshift(tmp,min_id-1);
    y_compl_normalized = circshift(tmp_compl,min_id_compl-1);
    [thetahat kappa] = circ_vmpar(alpha,y_normalized);
    kappa = kappa+2
    [thetahat_compl kappa_compl] = circ_vmpar(alpha,y_compl_normalized);
%     thetahat = thetahat - (min_id-1)*pi/180;
%     thetahat_compl = thetahat_compl - (min_id_compl-1)*pi/180;
    [p alpha] = circ_vmpdf(alpha, thetahat, kappa);
    [p_compl alpha] = circ_vmpdf(alpha, thetahat_compl, kappa_compl);
    y_normalized = tmp/trapz(alpha,tmp);
    y_compl_normalized = tmp_compl/trapz(alpha,tmp_compl);
       
    thetahat = thetahat - (min_id-1)*pi/180;
    p = circshift(p,-min_id+1);
    thetahat_compl = thetahat_compl - (min_id_compl-1)*pi/180;
    p_compl = circshift(p_compl,-min_id_compl+1);
    
    if thetahat < -pi/2 
        thetahat = thetahat + pi; 
    end
    if thetahat_compl < -pi/2 
        thetahat_compl = thetahat_compl + pi; 
    end 

    % Plot result
    f2 = figure(2);
    if ch == 2
        set(gca, 'YScale', 'log'), plot(x,y_normalized,x,p,'LineWidth',2)
        hold on, plot(x,y_compl_normalized,x,p_compl,'LineWidth',2,'LineStyle',':')
        % elseif ch == 3
        %     plot(x,tmp,x,p,'LineWidth',2)
    else
        plot(x,y_normalized,x,p,'LineWidth',2)
        hold on, plot(x,y_compl_normalized,x,p_compl,'LineWidth',2,'LineStyle',':')
    end
    eval(['title(''Mean orientation - ' subj ''')']);
    xlabel('Orientation in degrees'), ylabel('Incidence (%)')
    legend('Exp data','Fitted Von Mises PDF','Exp data (complem)','Fitted Von Mises PDF (complem)')
%     annotation('textbox',[0.587 0.7 0.3 0.1],'FontSize',9,...
%         'String',{'Von Mises PDF: ',['Theta =' num2str(thetahat*180/pi, '%4.2f°'), ', Kappa =' num2str(kappa, '%4.2f')]});
     annotation('textbox',[0.587 0.7 0.3 0.1],'FontSize',9,...
        'String',{'Von Mises PDF: ',['Theta =' num2str(thetahat*180/pi, '%4.2f°'), ', Kappa =' num2str(kappa, '%4.2f'), 'Theta'' =' num2str(thetahat_compl*180/pi, '%4.2f°'), ', Kappa'' =' num2str(kappa_compl, '%4.2f')]});


    %reduce matrix Y
    Y2 = Y;
    if ch==1
        Y2(:,str2num(answer{2}):end) = [];   Y2(:,1:str2num(answer{1})) = [];
    else
         Y2(:,str2num(answer{3}):end) = [];   Y2(:,1:str2num(answer{2})) = [];
    end
    Y2 = Y2';
    f4=figure(8), hold on,
    options.handle=figure(8);
    options.color_area = [180 180 180]./255;    % Red theme
    options.color_line = [142 142 142]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    options.x_axis = x';
    options.x_axis = options.x_axis(:)';

    % Computing the mean and standard deviation of the Y2 matrix
    Y2_mean = mean(Y2,1); a2= trapz(x,Y2_mean); Y2_mean = Y2_mean./a2;
    Y2_std  = std(Y2,0,1)./a2;
    error = (Y2_std./sqrt(size(Y2,1)));
    p = p./trapz(x,p);

    % Plotting the result
    figure(options.handle);
    x_vector = [options.x_axis, fliplr(options.x_axis)];
    patch = fill(x_vector, [Y2_mean+error,fliplr(Y2_mean-error)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, Y2_mean, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    hold on;
    if ch == 2
        set(gca, 'YScale', 'log'), plot(x,p,'LineWidth',2,'color', [1 1 1]./255,'LineStyle',':')
        % elseif ch == 3
        %     plot(x,tmp,x,p,'LineWidth',2)
%     elseif ch == 3
%         x2=[0.5:1:89.5, -89.5:-0.5]';
%         plot(x2(1:90),p(1:90),'LineWidth',2,'color', [1 1 1]./255,'LineStyle',':')
%         plot(x2(91:180),p(91:180),'LineWidth',2,'color', [1 1 1]./255,'LineStyle',':')
    else 
        x2=[-89.5:-0.5, 0.5:89.5,]';
        plot(x(1:90),p(1:90),'LineWidth',2,'color', [1 1 1]./255,'LineStyle',':')
        plot(x(91:180),p(91:180),'LineWidth',2,'color', [1 1 1]./255,'LineStyle',':')
    end
    eval(['title(''Mean orientation - ' subj ''')']);
    xlabel('Orientation (Degrees)'), ylabel('Normalized Incidence')
    legend('Exp data','Exp data','Fitted VM PDF')
    annotation('textbox',[0.587 0.67 0.3 0.1],'FontSize',10,...
        'String',{'VM PDF: ',['\theta =' num2str(thetahat*180/pi, '%4.2f°'), ', \kappa =' num2str(kappa, '%4.2f')]});

    hold off;


    % line([15.27,15.27],[0 0.05],  'Color',[255 255 255]./255,'LineWidth',1.5,'LineStyle','-'); %plot for each distribution


    save([foldername 'Orientation\Orientation_highcontrast' seriename ], 'Y', 'y_norm', 'x', 'f1', 'f2', 'f3','f4','f5','kappa','thetahat', 'kappa_compl', 'thetahat_compl', 'p', 'p_compl', 'lambda', 'region');

    results(length(stack),:) = [thetahat*180/pi kappa, thetahat_compl*180/pi kappa_compl;];

end