function  Optotagging_hj()
% PSTH Computes the peri-stimulus time histogram from spike times.
% Script to quickly look at the PSTH of a cell
% M Haberl 1/28/2023
% adapted from RP
% RP 2/2023
% adapted from MH 23
clear all;
close all;
clc;
restoredefaultpath;
rehash toolboxcache;
 % group1 = [12:15, 22:25, 32:37, 44:45, 47, 54:56]; % Haseeb cFos/tTA, only groups, mice, sessions with final clusters and opto1 folders
%delay = [0, 3, 5];
delay = [0];
 
do_plots=0;

 % group1 = [165, 184:187, 203:204, 205:216];
 %group1 = [241:247, 141:153, 161:165, 222:227, 235:239];
   group1 = [11:17 21:25, 31:38, 41:47, 52:56, 81:82, 84:88, 101:110, ...
     121, 141:142, 145, 148, 151, 161:162, 165, 168, 171, 181:182, 185, 188, 191, ...
   194, 197, 201:202, 205, 208, 211, 214, 217, 222, 225, 228, 231, 234, ...
    237, 241:242, 245, 248, 251, 254, ...
    143:144, 146:147, 149:150, 152:153, 163:164, 166:167, 169:170, ...
    183:184, 186:187, 189:190, 192:193, 195:196, 198:199, 203:204, 206:207, ...
    209:210, 212:213, 215:216, 218:219, 223:224, 226:227, 229:230, 232:233, 235:236, ...
    238:239, 243:244, 246:247, 249:250, 252:253, 255:256];
    % deleted cells in 225, 171, 170
  

%group2 = [71 85 100 104 112]; %control
combined_groups = {group1}; %
sFolders = {'opto1'};
load('W:\Haseeb\Analysis_scripts\DataOrganization\sessionInfo.mat');

%addpath('W:\LABanalysis\SilviaProjectCode\wavelet');
addpath('W:\Haseeb\Matlab_code');
addpath('W:\LABanalysis\SilviaProjectCode\AnalysisPerTrial\RunAnalysis\Chris\binFinderRedux\Parameters\');
addpath('W:\LABanalysis\SilviaProjectCode\MousePath\MatlabImportExport_v6.0.0');
addpath('W:\LABanalysis\SilviaProjectCode\Tools'); % just needed for write_csvfile
%rmpath(genpath('C:\Users\rinap\Dropbox\SilviaProjectCode\AnalysisPerTrial\RunAnalysis\IO'));
addpath(genpath('W:\LABanalysis\cluster_cutting\Newest MClust (microvolts, history logging)')); %required for reading the LFPs
%addpath('W:\LABanalysis\SilviaProjectCode');

%addpath(genpath('W:\LABanalysis\cluster_cutting\NlxVideoFix_V2_NotBroke'));
addpath('W:\Haseeb\Analysis_scripts\DataOrganization');
addpath('W:\LABanalysis\SilviaProjectCode\AnalysisPerTrial\+sleep');
%addpath('W:\LABanalysis\SilviaProjectCode\+COX_Project\+cox_functions');
ticker = 0;
SpikesPerGroup = [];
BeforeAfterStim_AllCellsPerG = [];
SpikeTimeStamps_AllCellsPerG = [];
for del = 1:1
delay_time = delay(del) / 1000; % delay time in ms
for gg = 1:1 %1:
    all_i_thisgroup = combined_groups{gg}; %sessionsToRun;
    tt = 0;
    SpikesPerSess = [];
    for i_iter = 1:numel(all_i_thisgroup)
        ticker = ticker +1;
        tt = tt+1;
        iii = all_i_thisgroup(tt);
        aaa = sessInfo(iii).animal;
        animal_matrix(ticker) = aaa;
        iii_matrix(ticker) = iii;
        group_matrix(ticker) = gg;

        analyze_dir = sFolders{1};
        %% load Opto TTL timestamps
         optoDir = fullfile(sessInfo(iii).mainDir, analyze_dir);
         [TTLData.ON_ts, TTLData.OFF_ts, TTLData.lengthONStim] = GetLaserTTL_hj(optoDir);
%        TTL_timestamps = fullfile(sessInfo(iii).mainDir, 'TTL_Timestamps', 'TTLData_opto.mat');
%        load(TTL_timestamps);
        %% Read TT List
        tt_files = {};
        cell_no = 0;
        fid=fopen(fullfile(sessInfo(iii).mainDir, sessInfo(iii).tList));
        display(iii)
        display(fid)
        while 1
            
            cell_no = cell_no +1;
         
            tline = fgetl(fid);
           %tline = fopen(fid);
            %display(tline)
            if ~ischar(tline), break, end
            tt_files{cell_no} = tline;
            %display(tt_files)
        end
        fclose(fid);
        numCells = numel(tt_files);
        %%display(tt_files)
        
        %% Load spikeData
        spikeData = readSpikeDataOnly(fullfile(sessInfo(iii).mainDir, analyze_dir) ,tt_files);
        
        %spikeTimes is a cell of arrays of spike time stamps of all the cells 
        spikeTimes = fixSpikes(spikeData);
        %display(spikeTimes)
        
        %% find spikes within TTLs
        % TO store all spiking data for every TTL
        SpikesWithinTTLs = {};
        SpikeTS_withinTTL = {};
        PerStimSpikeRate = [];

        day_TTL_spikes = {};
        sum_TTL_spikes = {};

        for cell = 1:numel(spikeTimes) % analyze cell by cell
            cellSpikes = spikeTimes{cell};
           % display(cellSpikes);
            [~, cellname]  = fileparts(tt_files{cell});
            if isempty(cellSpikes)
                fprintf('Skipping: No spikes in cell nr. %s\n', num2str(cell));
                continue
            end
            c_PSTH = [];
          
            max_number_trials = size(TTLData.ON_ts,2);
            total_cell_TTL_spikes = [];

           % testCell = {}
            
            for pulse = 1:max_number_trials % to go through all TTLs in the TTL Data file
                ttlON = TTLData.ON_ts(pulse) + delay_time; %0.000;  % adding a laser 0 ms delay
                %ttlON = TTLData.ON_ts(pulse);
                %ttlOFF = TTLData.OFF_ts(pulse);
                ttlOFF = ttlON + 0.015;  % 15 ms pulse duration
                %x = length(ttlON:ttlOFF);  %stimulation length (this is just the number of values, not the real time)
                %min_x = ttlON - x;
                min_PSTH = ttlON - 1;  % Start peri-stim histogram 100ms before and finish 100ms afterwards
                % max_x_plot = ttlOFF + 3*x;
                max_PSTH = ttlON + 1;
                %max_x_data = ttlOFF + 5*x;
                %max_x_data = ttlON + 1;
                StimLength(pulse) = ttlOFF-ttlON; %gives length of stim in seconds, Haseeb> 15ms
                % to store spike info together per cell
               % StimLength(pulse) = 0.015
               

          % spikesInTTL gives the index number of the spikes during TTLs out of the total
          % spikes in the session. For instance for cell 16 in mouse 520,
          % Day 2, TTL 3 has one spike with index 219 out of a total 9311
          % spikes.

          % cellSpikes gives the timestamp of each spike

                spikesInTTL = find(ttlON <= cellSpikes & cellSpikes <= ttlOFF);
                
                %testCell{pulse} = spikesInTTL
                
                TTL_array = size(spikesInTTL);
                total_cell_TTL_spikes(pulse) =  TTL_array(1);

                spikesPSTH=  find(min_PSTH <= cellSpikes & cellSpikes <= max_PSTH);
                % now get the real-time of the spikes in the PSTH and then
                % the time relative to the onset of the Stimulation
                spt_PSTH = [cellSpikes(spikesPSTH)]' - ttlON;
                bins_sec = 0.05; % 50ms size bins
                edges = [-1 : bins_sec : 1];  % define bins of the PSTH
                [N_sp, edges] = histcounts(spt_PSTH, edges); % how many spikes occur at each timepoint / bin
                sp_per_sec(pulse,:) = (N_sp / bins_sec) ; %how many spikes per second occur at each bin ; divide by time of bin in sec. (equivalent to e.g. for 50ms bins multiply by 20 to get to Hz)
                average_firing = numel(cellSpikes) / (cellSpikes(end) -cellSpikes(1) );
                [TTL_on_off , ~] = histcounts([0, 0+StimLength(pulse)], edges);
                c_PSTH = [c_PSTH; N_sp]; % accumulate the number of spikes of all stimulations of this cells for PSTH plots
                % make a raster plot for this cell with every pulse
                
                if do_plots == 1
                if pulse==1
                    figure
                end
                scatter(spt_PSTH, pulse, '.k'), hold on  % plot all spikes occuring during PSTH and use max_number_trials-pulse to invert the axis
                clearvars spt_PSTH
                end
                
            end % going through every TTL

            

            day_TTL_spikes{1,cell} = total_cell_TTL_spikes;
            day_array = day_TTL_spikes{1,cell};
            sum_TTL_spikes{cell} = sum(day_array);
            %disp(day_array);

            if do_plots == 1
            set(gca, 'ydir', 'reverse' )
            plot([0, 0], [0, 100] ,'--bs','LineWidth',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',1)
            plot([StimLength(pulse), StimLength(pulse)], [0,100] ,'--bs','LineWidth',0.5,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',1)
           xlabel('time from laser pulse in sec')
           ylabel('trials')
           filename = fullfile(optoDir ,sprintf('Cell_%s_%s_PST-spikes', num2str(cell),cellname));
           set(gcf, 'renderer', 'painters');
           print(filename,'-dpng')
           print(filename,'-dpdf')
          % exportgraphics(gcf, strcat(filename, '.pdf'))
           %saveas(gcf, strcat(filename, '.svg'))
           
           figure
           % plot(edges(2:end)-0.05, mean_PSTH), hold on
           % plot_delin = edges(2:end);
           histogram('BinEdges', edges, 'BinCounts', mean(sp_per_sec,1), ...
               'EdgeColor','none', 'FaceColor',[0 0 0.3]) % Plot the peri-stimulus histogram, as the
           ylabel('spikes / sec')
           title(sprintf('Mean firing rate: %s', num2str(average_firing)));
           filename = fullfile(optoDir ,sprintf('Cell_%s_%s_PSTH', num2str(cell), cellname));
           print(filename,'-dpng')
           print(filename,'-dpdf')
           close all
            end
           
           %{
           mean_PSTH = mean(c_PSTH,1);
            scatter(plot_delin(find(TTL_on_off)),[0, 0]);
            on_off_nmbrs = find(TTL_on_off);
            on_off = plot_delin(find(TTL_on_off));
            plot([on_off(1), on_off(1)], [0, max(mean_PSTH)] ,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
            plot([on_off(2), on_off(2)], [0,max(mean_PSTH)] ,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)
           %}
            all_c_spt_PSTH{cell} = mean(sp_per_sec,1); % combine all single-units to write the analysis to a file
%}
        end % looping through cells

        %display(cell)
        mkdir(fullfile(sessInfo(iii).mainDir, 'processedData'));
        cells_spikes_intervals_file = fullfile(sessInfo(iii).mainDir, 'processedData', sprintf('Cells_Spikes_in_Intervals_delay%sms.mat', num2str(delay(del))) );
        %total_TTLs_file = fullfile(sessInfo(iii).mainDir, 'processedData', 'total_daily_TTL_spikes.mat');
        save(cells_spikes_intervals_file, 'day_TTL_spikes', 'sum_TTL_spikes'); 

    end  %  looping through every animal
end  % looping through every group

end % end going through 0,3, 5 ms delays
end  % end of function
