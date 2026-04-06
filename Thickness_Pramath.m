%% Thickness estimation from 2PM images
clear all, close all, clc


%% INPUTS

sampleName = '5'; %edit sample name
folderName = 'F:\_F18mo\20201207_RPA_F_18m_0_3\unloaded/5\Thickness\'; %folder including sample + \Thickness
vein = 0;  %% A vein or aneurysm? Do not change 
%  stacks = [7 13 1 10 16 4]; 
stacks = [5]; 
%%

ns = length(stacks); 
for i=1:ns    fileNamesC(i) = {num2str(stacks(i))}; end

[thickness, load_info] = profile2thickness(sampleName, fullfile(folderName, fileNamesC(:)),vein);


%% extraction and threshold profile for thickness 

function [thickness, load_info] = profile2thickness(sampleName, stacks, vein)
   
% stacks=fullfile(folderName, fileNamesC(:)); % when need to skip function
%   nregions = input('How many regions acquired?')
  
  [fp,~,~] = fileparts(stacks{1});
  outFileNameBase = fullfile(fp,[sampleName, '_thickness']);
  a = [];
  ns = length(stacks);
  
  % Import raw data from selected stack
  for i = 1:ns
      
    [~,fp2,~] = fileparts(stacks{i}); fp2 = str2double(fp2);
    
    
    load_info{1,i} = 'XX'; load_info{2,i} = 'XX';
    
    new_1 = importdata([stacks{i},'_1.csv']);  
    new_2 = importdata([stacks{i},'_2.csv']); 
    new_3 = importdata([stacks{i},'_3.csv']); 
   
    index_1 = 1:size(new_1.data,1); 
    n_mes = size(new_1.data,2); pos_ch2 = n_mes+2; pos_ch3 = n_mes*2+2;
    
    data_spec = [index_1' new_1.data new_2.data new_3.data];
    
    %if data==0 means that the slice was empty, sometimes for the first and last slices
    %background elastin values are high, put background to 0,
    data_spec(data_spec==0) = NaN; 
    a = min(data_spec(:,pos_ch2:pos_ch3-1)); 
    b = min(data_spec(:,pos_ch3:end)); 
    data_spec(:,pos_ch2:pos_ch3-1) = data_spec(:,pos_ch2:pos_ch3-1) - a; 
    data_spec(:,pos_ch3:end) = data_spec(:,pos_ch3:end) - b;
    data_spec(data_spec<0) = NaN; data_spec(isnan(data_spec))=0;  
    data_spec(:,end+1)=sum(data_spec(:,2:end),2);
    %normalized distributions
    ndata = data_spec; ndata(:,2:end) = data_spec(:,2:end)./trapz((data_spec(:,2:end))).*100;
    
    avg_data(:,1) = mean(data_spec(:,2:pos_ch2-1),2);
    avg_data(:,2) = mean(data_spec(:,pos_ch2:pos_ch3-1),2);
    avg_data(:,3) = mean(data_spec(:,pos_ch3:19),2);
    avg_data = avg_data./trapz(avg_data).*100;


    %plot normalized distributions
    figure(3),
    subplot(1,ns,i)    %subplot(ceil(ns/3),3,i)
    ax(i) = gca; 
    plot(ndata(:,1),avg_data(:,1),'r',ndata(:,1),avg_data(:,2),'g',ndata(:,1),avg_data(:,3),'b')
    title(['St:', num2str(fp2) ,', p:',load_info{1,i},', ',load_info{2,i}]),    hold on
    
    for k = 1:n_mes
        
        %initial and final point for each channel, defined ad pin%<=P<=pfin%
        n1 = ndata(:,1+k); n2 = ndata(:,1+n_mes+k); n3 = ndata(:,1+2*n_mes+k); n_tot = ndata(:,end);
        zin1(k) = find(cumtrapz(n1)<=5,1,'last'); zfin1(k) = find(cumtrapz(n1)<90,1,'last'); th1(k) = zfin1(k)-zin1(k);
        zin2(k) = find(cumtrapz(n2)<=5,1,'last'); zfin2(k) = find(cumtrapz(n2)<95,1,'last'); th2(k) = zfin2(k)-zin2(k);
        zin3(k) = find(cumtrapz(n3)<=15,1,'last'); zfin3(k) = find(cumtrapz(n3)<95,1,'last'); th3(k) = zfin3(k)-zin3(k);
        %      zin_tot = find(cumtrapz(n_tot)<=8,1,'last'); zfin_tot = find(cumtrapz(n_tot)<92,1,'last'); th_tot=zfin_tot-zin_tot;
        th_tot(k)=zfin2(k)-zin1(k);
        
        
        % interface adventitia-media considered the change in ch3 profile
        dydz = gradient(n3) ./ gradient(ndata(:,1)); dif = diff(dydz); % derivative ch3 distribution, then differentiale
        [~,z_max_1]=max(n1); [~,z_max_3]=max(n3); [~,z_max_2]=max(n2); % limits within which to search (max of ch 1,2,3)
        [diff_inter,index] = max(dif(max(z_max_1-5,zin1(k)):min(z_max_1+3,z_max_3-4)));
        try
            if vein == 0
                if length(find(dif==diff_inter))>1
                    tmp = find(dif==diff_inter);
                    z_interface(k) = tmp(1)+1;
                else
                z_interface(k) = find(dif==diff_inter)+1;
                end
            end
        % z_interface = find(dif==max(dif(z_max_1:(zfin1-1))));  % interface point (second choice, end of ch1)
        catch
            z_interface(k) = ceil(mean([z_max_1 max([z_max_2 z_max_3])])+2);
        end
        z_interface(k) = max([z_max_1 z_interface(k)]);
    end 
    
    %structure with all thickness
    thickness(i).stack = fp2;
    thickness(i).thickness_ch1 = mean(th1); thickness(i).thickness_ch2 = mean(th2); thickness(i).thickness_ch3 = mean(th3);
    thickness(i).thickness_wall = mean(th_tot);
    thickness(i).in = mean(zin1);
    thickness(i).fin = mean(zfin2);
    
    local(i).thickness_ch1 = th1; local(i).thickness_ch2 = th2; local(i).thickness_ch3 = th3;
    local(i).thickness_wall = th_tot;
    local(i).in = zin1;
    local(i).fin = zfin2;
    
    
    Zin1 = floor(mean(zin1)); Zfin1 = floor(mean(zfin1)); Zin2 = floor(mean(zin2)); Zfin2 = floor(mean(zfin2));
    if vein == 0
    	Z_interface = floor(mean(z_interface));
    end
    f3 = figure(3); 
    subplot(1,ns,i)    %subplot(ceil(ns/3),3,i)
    plot(Zin1,avg_data(Zin1,1),'or',Zin2,avg_data(Zin2,2),'og')
    plot(Zfin1,avg_data(Zfin1,1),'or',Zfin2,avg_data(Zfin2,2),'og')
    hold on
    
    if vein == 0
        line([Z_interface Z_interface], get(gca, 'ylim')); %plot for each distribution

        th_a = z_interface - zin1; Th_a = Z_interface - Zin1;
        th_m = zfin2 - z_interface; Th_m = Zfin2 - Z_interface;
        thickness(i).z_adv_med = Z_interface;   local(i).z_adv_med = Z_interface;  
    else
        th_a = zfin2 - zin1; Th_a = Zfin2 - Zin1;
        th_m = zfin2 - zin1; Th_m = Zfin2 - Zin1;
        thickness(i).z_adv_med = NaN;   local(i).z_adv_med = NaN;
    end
    
    thickness(i).thickness_adv = Th_a;      local(i).thickness_adv = Th_a;
    thickness(i).thickness_med = Th_m;      local(i).thickness_med = Th_m;
    
    txt = {'th. adv = ' num2str(Th_a), 'th. media = ' num2str(Th_m)};
    text(80,2.5,txt,'FontSize',8)

%   legend('Collagen','Elastin','Nucleus')
  linkaxes(ax)
    

    f4 = figure(4), hold on, 
    subplot(1,ns,i)   %subplot(ceil(ns/3),3,i)
    options.handle=figure(4);
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem'; 
    options.x_axis = index_1;
    options.x_axis = options.x_axis(:)';

    options.color_area = [1 0 0];    % Red theme
    options.color_line = [1 0 0];
    % Computing the mean and standard deviation of the Y2 matrix
    area1= trapz(new_1.data); new_1_norm = new_1.data./area1*100; avg_new_1_norm=mean(new_1_norm,2);
    new_1_std  = std(new_1_norm');
    error = (new_1_std./sqrt(size(new_1.data,2)))';
    % Plotting the results    
    figure(options.handle);
    x_vector = [options.x_axis, fliplr(options.x_axis)];
    patch = fill(x_vector, [avg_new_1_norm'+error',fliplr(avg_new_1_norm'-error')], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, avg_new_1_norm, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    
    hold on
    options.color_area = [0 1 0];    % Green theme
    options.color_line = [0 1 0];
    % Computing the mean and standard deviation of the Y2 matrix
    area2= trapz(new_2.data); new_2_norm = new_2.data./area2*100; avg_new_2_norm=mean(new_2_norm,2);
    new_2_std  = std(new_2_norm');
    error = (new_2_std./sqrt(size(new_2.data,2)))';
    % Plotting the results    
    figure(options.handle);
    x_vector = [options.x_axis, fliplr(options.x_axis)];
    patch = fill(x_vector, [avg_new_2_norm'+error',fliplr(avg_new_2_norm'-error')], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    plot(options.x_axis, avg_new_2_norm, 'color', options.color_line, ...
        'LineWidth', options.line_width); 
    
    hold on
    options.color_area = [0 0 1];    % Blue theme
    options.color_line = [0 0 1];
    % Computing the mean and standard deviation of the Y2 matrix
    area3= trapz(new_3.data); new_3_norm = new_3.data./area3*100; avg_new_3_norm=mean(new_3_norm,2);
    new_3_std  = std(new_3_norm');
    error = (new_3_std./sqrt(size(new_3.data,2)))';
    % Plotting the results    
    figure(options.handle);
    x_vector = [options.x_axis, fliplr(options.x_axis)];
    patch = fill(x_vector, [avg_new_3_norm'+error',fliplr(avg_new_3_norm'-error')], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, avg_new_3_norm, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    
    title('Intensity profile')
    xlabel('Depth (um)'), ylabel('Norm. intensity profile (%)')
    if i==ns
      legend('Collagen SEM','Collagen AVG','Elastin SEM','Elastin AVG','Nuclei SEM','Nuclei AVG')
    end
	hold on;
    save([outFileNameBase '.mat'], 'thickness', 'local' ,'f4', 'f3')
    
    clear data_spec avg_data


  end
  
if ns==1
    tmp = [thickness.in thickness.fin thickness.z_adv_med thickness.thickness_adv thickness.thickness_med];
    thickness.line = tmp;
    save([outFileNameBase '.mat'], 'thickness', 'local' ,'f4', 'f3')
end

end  

  


