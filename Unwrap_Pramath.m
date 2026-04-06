
clear all, close all, clc

%% IMPORT IMAGES AND DEFINE POSITIONS FOR CARTESIAN-POLAR TRANSFORMATION

 
%% input (only part that change)
stacks = [15]; %stack names

full_path = 'F:\_M18mo\251213_16-47-17_M18mo_M3_RPA_15mmHg/';
%%
%'D:\Dr. Cavinato\Pramath\20240515_Dar_Ao_arches\20240515_Dar_Ao_arches\Arch1\240515_11-05-24_Arch1_1a\'
for j = 1 : length(stacks)
    ImageJ
    stack_name = stacks(j);
    stack_path = fullfile(full_path, [num2str(stack_name),'_edit\']);
    
    % from image to matrix
    eval(['ij.IJ.run(''Open...'', ''path=[',full_path, num2str(stack_name),'_edit\C1-',num2str(stack_name),'.tif]'');']);
    IJM.getDatasetAs('C1') 
    xy_1 = im2double(uint16(C1)); clear C1; ij.IJ.run("Close",""); 
    eval(['ij.IJ.run(''Open...'', ''path=[',full_path, num2str(stack_name),'_edit\C2-',num2str(stack_name),'.tif]'');']);
    IJM.getDatasetAs('C2')
    xy_2 = im2double(uint16(C2)); clear C2; ij.IJ.run("Close","");
    eval(['ij.IJ.run(''Open...'', ''path=[',full_path, num2str(stack_name),'_edit\C3-',num2str(stack_name),'.tif]'');']);
    IJM.getDatasetAs('C3')   
    xy_3 = im2double(uint16(C3)); clear C3; ij.IJ.run("Close","");    
   
 
 
    n = size(xy_3,1); m = size(xy_3,2); l = size(xy_3,3); 
    filename = [num2str(stack_name),'.mat'];
    save([stack_path,filename], 'm', 'l', 'n', 'filename');

    yz = zeros(l,m,n); yz_1 = zeros(l,m,n);yz_2 = zeros(l,m,n); yz_3 = zeros(l,m,n);

    for nn = 1:n,    yz_1(:,:,nn) = imrotate(reshape(xy_1(nn,:,:),[m,l]),-90); end 
    save([stack_path,filename], 'yz_1', 'xy_1','-append'); clear xy_1;
    for nn = 1:n,    yz_2(:,:,nn) = imrotate(reshape(xy_2(nn,:,:),[m,l]),-90); end
    save([stack_path,filename], 'yz_2', 'xy_2','-append'); clear xy_2;
    for nn = 1:n,    yz_3(:,:,nn) = imrotate(reshape(xy_3(nn,:,:),[m,l]),-90); end
    save([stack_path,filename], 'xy_3', 'yz_3','-append'); clear xy_3;
    yz = yz_1 + yz_2 + yz_3;  clear yz_1 yz_2 yz_3; 
    disp((floor(n/2)+floor(n/4)));
    pt = [floor(n/2)-floor(n/4), floor(n/2), (floor(n/2)+floor(n/4))];
    % pt = [1030, 1035, 1042];
    if mod(m,2)==0    v = 7000; else    v = 7001; end

    ze_b = zeros(v,v);
    ze = zeros(v,v); 
    se = strel('sphere',2); se2 = strel('sphere',2);

    for i = 1:3
        % draw 3 points defining the profile
        yz_c = yz(:,:,pt(i));
        counts = imhist(yz_c, 16);   T = otsuthresh(counts); 
        yz_cb = imbinarize(yz_c);  %one central yz profile, binarized 

        yz_cb = imdilate(yz_cb, se);
        yz_cb = imerode(yz_cb, se); yz_cb = imclose(yz_cb, se2); % cleaner shape
        figure(), imshow(yz_c) 
        eval(['set(gcf,''Name'',''Zoom in if needed, then press Enter'');']);
        pause(), hold on

        eval(['set(gcf,''Name'',''Select THREE POINTS along the profile'');']);

        %pick the first point on the inner surface
        [xR1,yR1] = ginput(1);
        hold on
        plot(xR1,yR1,'+b','LineWidth',2)
        %pick the second point on the inner surface
        [xR2,yR2] = ginput(1);
        hold on
        plot(xR2,yR2,'+b','LineWidth',2)

        %pick the third point on the inner surface
        [xR3,yR3] = ginput(1);
        hold on
        plot(xR3,yR3,'+b','LineWidth',2)

        AR = [xR1 yR1 1 ; xR2 yR2 1; xR3 yR3 1];
        BR = [-xR1^2-yR1^2; -xR2^2-yR2^2; -xR3^2-yR3^2 ];
        CR = AR\BR;

        a = CR(1);
        b = CR(2);

        xcR = -a/2; xcR = xcR + floor(v/2)-floor(m/2);
        ycR = -b/2;
        centerpt(i,:) = [xcR ycR];  %center for the polar coordinates


        close all
    end

    %apply transformation
    xcR = mean(centerpt(:,1)); ycR = mean(centerpt(:,2)); %average of 3 measurements
    ze(1:l,floor(v/2)+1-floor(m/2):(floor(v/2)+floor(m/2)))=yz_c;
    figure(), imshow(ze), hold on, plot(xcR, ycR,'+r','LineWidth',2)

    clearvars a AR b BR C1 C2 C3 CR xb yz_cb xR1 yR1 xR2 yR2 xR3 yR3
    [h,w,~] = size(ze);
    [X,Y] = meshgrid((1:w)-floor(xcR), (1:0.5:h)-floor(ycR));
    [theta,rho] = cart2pol(X(1:2*l,floor(v/2)+1-floor(m/2):floor(v/2)+floor(m/2)), Y(1:2*l,floor(v/2)+1-floor(m/2):floor(v/2)+floor(m/2))); 
    theta = theta/pi*180;

    yz_c_2 = imresize((yz_c), [2*l m]);
    pcoor = [theta(:) rho(:) yz_c_2(:)];
    rho_theta90 = pcoor(find(pcoor(:,1)==-90),2);
    d_rho = (max(rho_theta90)- min(rho_theta90)+1)/l;
    rhospace = min(pcoor(:,2)):d_rho:max(pcoor(:,2));
    k = length(rhospace);

    close all

    z = zeros(2*l,m);
    fig = figure, h = warp(theta, -rho, z, yz_c_2), view(2), axis off, xlim([min(pcoor(:,1)) max(pcoor(:,1))]), ylim([-max(pcoor(:,2)) -min(pcoor(:,2))]),
    caxis([0 1]);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 m/150 k/150],'color','k');
    set(gca, 'position', [0 0 1 1],'color',[0 0 0]);
    fig.InvertHardcopy = 'off';

    eval(['set(gcf,''Name'',''Press Enter: Apply trasformation to all images'');']);

   
%     eval(['saveas(fig,[stack_path,''check_', num2str(stack_name),'.tiff''])']);
    save([stack_path,filename], 'theta', 'rho', 'xcR', 'ycR', 'z', 'pcoor','k', 'stack_name','-append')
    clearvars -except stacks full_path
    close()
    ij.IJ.run("Quit","");

end



for j = 1 : length(stacks)
%apply to all channels and save YZ images
stack_name = stacks(j);
stack_path = fullfile(full_path, [num2str(stack_name),'_edit\']);
load([stack_path,[num2str(stack_name),'.mat']], 'l', 'm', 'n', 'k', 'theta', 'rho', 'z', 'pcoor','yz_1','yz_2','yz_3')
yz_s_1 = zeros(2*l,m);
yz_2x = zeros(2*l,m,n);

for c = 1:3
eval(['yz_2x = imresize(yz_',num2str(c),', [2*l,m]);']);
eval(['w = waitbar(0,''Process ',num2str(c),'/3...'');']);
for nn = 1:n 
    yz_s = yz_2x(:,:,nn);
    fig = figure('visible','off'); 
    h = warp(theta, -rho, z, yz_s), view(2), axis off, 
    xlim([min(pcoor(:,1)) max(pcoor(:,1))]), ylim([-max(pcoor(:,2)) -min(pcoor(:,2))]),
    caxis([0 1]);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 m/150 2*k/150],'color','k');
    set(gca, 'position', [0 0 1 1],'color',[0 0 0]);
    fig.InvertHardcopy = 'off';
    saveas(fig,[stack_path,num2str(stack_name),'_C',num2str(c),'_',num2str(nn),'.tiff']);
    close all

    waitbar(nn / n)
end
close(w)
end


clearvars -except stacks full_path


end
