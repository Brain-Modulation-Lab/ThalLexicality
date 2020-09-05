temp = dir('contact_*.mat');

fs = 1000;

for i_order = 1:length(temp)
    
    if strcmp('contact_64_session1_highgamma_ref.mat',temp(i_order).name)...
            || strcmp('contact_93_beta_ref.mat',temp(i_order).name)
    
    
    
        % decide band type according to temp(i_order).name
        if strcmp(temp(i_order).name(end-11:end-4),'beta_ref')
            band_id = 1; % beta_ref
        else
            band_id = 2; % highgamma_ref
        end

        clearvars -except i_order temp fs band_id

        % load in timing table and data
        load([temp(i_order).folder filesep temp(i_order).name]);
        load(temp(i_order).name);


        % sort z_oi and epoch_oi according to the trialorder stored in timing
        % table
        st_z = z_oi(ephy_timing.TrialOrder); % each trial is 2s precue to 2s post speech offset
        st_epoch = epoch_oi(ephy_timing.TrialOrder,:);

        tr2del = find(isnan(ephy_timing.ActOn));

        st_z(tr2del) = [];
        st_epoch(tr2del,:) = [];
        ephy_timing(tr2del,:) = [];

        st_z_smooth = cellfun(@(x) smooth(x',200)',st_z,'UniformOutput',0);   

        %baseline (which is 1s pre-cue)
        Bt = cellfun(@(x) x(1001:2000), st_z_smooth, 'UniformOutput',0);

        % normalize to baseline
        zz = cellfun(@(x,y) (x-mean(y))./std(y), st_z_smooth, Bt,'UniformOutput',0);

        % define region to be plotted

        roi_starts = num2cell((st_epoch.stimulus_starts - st_epoch.starts - 0.5) * fs)'; 

        if min(cell2mat(cellfun(@length, zz, 'UniformOutput',0))) < 5000   
           length_used = floor(min(cell2mat(cellfun(@length, zz, 'UniformOutput',0)))/100) * 100;

           roi_ends = num2cell(repmat(length_used,[1 length(zz)]));

        else
           roi_ends = num2cell((st_epoch.stimulus_starts - st_epoch.starts + 3) * fs)';
        end

        zz_crop = cellfun(@(x,y,z) z(x:y), roi_starts, roi_ends, zz,'UniformOutput',0);

        plot_mat = cell2mat(zz_crop');

        figure;

        tp = linspace(-0.5, length(plot_mat)/1000 -0.5 - 0.001, length(plot_mat));
        imagesc(tp,1:size(plot_mat,1),plot_mat); set(gca, 'YDir', 'Normal');


        % adjust some plotting parameters according to different bands plotted
        switch band_id
            case 1 % beta
                hold on; time_spon =plot(ephy_timing.ReactionT,(1:length(ephy_timing.ReactionT)),'-k','LineWidth',3); %,'Color',[11/255 102/255 35/255]
                hold on; time_cue = plot([0,0],ylim,'k--');
                hold on; plot(ephy_timing.ActOn,(1:length(ephy_timing.ActOn)),'r*','MarkerSize',5); % red for lb
            %   hold on; plot(ephy_timing.ActMax,(1:length(ephy_timing.ActMax)),'k*','MarkerSize',5);
            %   hold on; plot(ephy_timing.ActOff,(1:length(ephy_timing.ActOff)),'m*','MarkerSize',5);
            case 2 % highgamma
                hold on; time_spon =plot(ephy_timing.ReactionT,(1:length(ephy_timing.ReactionT)),'-k','LineWidth',3);
                hold on; time_cue = plot([0,0],ylim,'k--');
                hold on; plot(ephy_timing.ActOn,(1:length(ephy_timing.ActOn)),'b*','MarkerSize',5); % blue for hg
            %   hold on; plot(ephy_timing.ActMax,(1:length(ephy_timing.ActMax)),'k*','MarkerSize',5);
            %   hold on; plot(ephy_timing.ActOff,(1:length(ephy_timing.ActOff)),'m*','MarkerSize',5);            
        end

        hh = gca;
        set(hh,'box','on');
        set(hh,'TickLength',[0.005,0.005]);
        set(hh,'XTickLabel',-0.5:0.5:3);
        colormap parula; caxis([-6,6]);
        pcolorbar = colorbar;
        set(pcolorbar,'TickLength',0.005);
        pcolorbar.Label.String = "z-scored power";

        xlabel("Time Relative to Cue Presentation (s)");
        ylabel("Trial Number");
        
        if band_id == 2
            saveas(gcf, 'Figure5a.fig');
        else
            saveas(gcf, 'Figure5b.fig');
        end
        close all;

    else
    end
    
end   