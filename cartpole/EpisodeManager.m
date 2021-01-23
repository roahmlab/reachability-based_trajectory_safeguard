classdef EpisodeManager < matlab.apps.AppBase
    % EPISODEMANAGER creates a simple UI to observe the progress of
    % training and to allow manually stopping trining without losing
    % progress
    
    % Copyright 2018 The MathWorks, Inc.
        
    properties (SetObservable,AbortSet)
        RequestToStopTraining (1,1) logical = false
    end
    
    properties (SetAccess = private, Hidden)
        % Info panel
        MainPanel 
        
        % Legend Info Panel
        LegendTextLabel
        LegendLine = []
        LegendAxes 
        
        
        % Progress Bar
        ProgressBar
        
        % Episode Info Panel
        EpisodeNumberLabel
        EpisodeRewardLabel
        EpisodeStepsLabel
        EpisodeQ0Label
        GlobalStepsLabel
        
        % Average Info Panel
        AverageRewardLabel
        AverageStepsLabel
        AverageWindowLengthLabel
        
        % Training Options Info Panel
        DeviceLabel
        LearnRateLabel
        ActorDeviceLabel
        CriticDeviceLabel
        ActorLearnRateLabel
        CriticLearnRateLabel
        MaxEpisodesLabel
        MaxStepsPerEpisodeLabel
        
        % Final Info Panel
        FinalReasonLabel
        FinalReasonResultLabel
        FinalValueLabel
        FinalValueResultLabel
        ElapsedTimeLabel
        
        
        % Message for Parallel
        MessageLabel
        
        % Plot Options Panel
        EnableQ0CheckBox
        ShowLastNEpisodesCheckBox
        ShowLastNEpisodesTextBox
    end
    
    properties (SetAccess = private,Hidden,Transient)
        Figure
        Axes
        StopButton
        PlotObjects = struct()
        
        EpisodeNum = 1
        
        % for plot control
        ShowEpisodeQ0 (1,1) logical = false
        HistoryLength = inf
        RawXData = []
        RawYData = []
    end
    
    properties (SetAccess = private, Hidden)
        StopTrainingFromButtonFlag = false
        MaxEpisodes = 1000
    end
    
    methods
        function this = EpisodeManager(trainingOptions,devices,learnRates,hasQ0)
            % set MaxEpisodes
            this.MaxEpisodes = trainingOptions.MaxEpisodes;
            %% build the app components
            % main figure
            createMainFigure(this);
            % create grid layout for uifigure
            createLayoutForFigure(this,hasQ0);
            % create grid layout for main panel
            createLayoutForInfoPanel(this,trainingOptions,devices,learnRates,hasQ0);
            this.Figure.Visible = 'on'; % set to off
            % intialize the action message
            setActionMessage(this,'');
        end
        
        %% set action message
        function setActionMessage(this,msg)
            this.MessageLabel.Text = msg;
        end
        
        %% get action message
        function msg = getActionMessage(this)
            msg = this.MessageLabel.Text;
        end
        
        %% add title, x/y labels for plot
        function setTextOnAxis(this,ttl,xlbl,ylbl) %#ok<*ADAPPREF>
            title(this.Axes,ttl);
            xlabel(this.Axes,xlbl);
            ylabel(this.Axes,ylbl);
        end
        
        %% reset 
        function reset(this)
            this.EpisodeNum = 1;
            this.RequestToStopTraining = false;
            cla(this.Axes);
        end
        
        %% delete
        function delete(this)
            delete(this.Figure);
        end
        
        %% update with show Episode Q0 or not
        function setShowEpisodeQ0(this,val)
            if this.ShowEpisodeQ0 ~= val
                this.ShowEpisodeQ0 = val;
                this.EnableQ0CheckBox.Value = val;
                updateLegendInfo(this,val);
                updatePlot(this);
            end
        end
        
        %% set hisotry length (latest N episodes)
        function toggleHistoryLength(this,val)
            if val
                setHistoryLength(this,this.ShowLastNEpisodesTextBox.Value);
            else
                setHistoryLength(this,inf);
            end
        end
        
        %% show latest N episode stats
        function setHistoryLength(this,val)
            if this.HistoryLength ~= val
                this.HistoryLength = val;
                if isinf(val)
                    this.ShowLastNEpisodesTextBox.Enable = 'off';
                else
                    this.ShowLastNEpisodesTextBox.Enable = 'on';
                    this.ShowLastNEpisodesTextBox.Value = val;
                end
                updatePlot(this);
            end
        end
        
        %% set start time
        function setStartTime(this,value)
            this.MainPanel.Title = getString(message('rl:general:MainPanelTitle',value));
        end
        
        %% set elapsed time
        function setElapsedTime(this,value)
            this.ElapsedTimeLabel.Text = value;
        end
        
        %% step one episode
        function stepEpisode(this,computedInfo)
            % call at the end of the episode to update the episode manager
            
            %% update plot
            plotfieldnames = plotEpisodeInfo(this,computedInfo,this.EpisodeNum);
            %% update labels and progress
            averageSteps = computedInfo.AverageSteps;
            globalStepCount = computedInfo.GlobalStepCount;
            maxEpisodes = this.MaxEpisodes;
            updateEpisodeInfo(this,this.EpisodeNum,computedInfo,maxEpisodes,plotfieldnames,averageSteps,globalStepCount);
            %% update the step count
            this.EpisodeNum = this.EpisodeNum + 1;
            
            % flush the event queue
            drawnow;
        end
                
        %% update with training result structure
        function updateWithTrainingResult(this,trainingResult,postTraining)
            % called post-training to recreate plot of the episode manager
            
            %% get information
            info = trainingResult.Information;
            opts = info.TrainingOpts;
            if postTraining
                stats = trainingResult;% postTraining: stats expanded 
            else
                stats = trainingResult.TrainingStats; % savedAgent: stats as struct
            end
            lastEpisodeIncomplete = isfield(info,'StopTrainingCriteria') && strcmp(info.StopTrainingCriteria,'Stop Training Button');
            if lastEpisodeIncomplete
                hasCritic = isfield(stats,'EpisodeQ0');
                stats = removeIncompleteEpisodeFromStats(stats,hasCritic);
            end
            episodeNum = stats.EpisodeIndex(end);
            globalStepCount = stats.TotalAgentSteps(end);
            maxEpisodes = opts.MaxEpisodes;
            averageSteps = stats.AverageSteps(end);
            startTime = info.TraningStartTime;
            elapsedTime = info.ElapsedTime;
            
            
            %% update starting time
            this.MainPanel.Title = getString(message('rl:general:MainPanelTitle',startTime));
            %% update stop training button
            this.StopButton.Text = getString(message('rl:general:TextTrainingStopped'));
            this.StopButton.Enable = false;
            %% update plot
            plotfieldnames = plotEpisodeInfo(this,stats,stats.EpisodeIndex);
            %% update labels and progress
            updateEpisodeInfo(this,episodeNum,stats,maxEpisodes,plotfieldnames,averageSteps,globalStepCount);
            %% update final info
            if  postTraining % plot with final training result
                reason = info.StopTrainingCriteria;
                value = info.StopTrainingValue;
                if this.StopTrainingFromButtonFlag
                    updateFinalInfo(this,reason,value,episodeNum);
                else
                    updateFinalInfo(this,reason,value);
                end
            else % plot with save agent result
                reason = opts.SaveAgentCriteria;
                value = opts.SaveAgentValue;
                updateSavedAgentInfo(this,reason,value);
            end
            %% update elpsed time
            this.ElapsedTimeLabel.Text = elapsedTime;
        end
        
        
        function stepLearn(this,data)
        end
        
        %% stop training
        function stopTraining(this,reason,value)
            % stop training without providing episodeIdx signifies that the
            % stop training button was hit (manually stopped). Otherwise
            % stop training will be called if the training reaches a
            % stopping condition (e.g. maximum number of episodes reached)
            if nargin < 2
                % training stopped manually
                episodeIdx = this.EpisodeNum;
                updateFinalInfo(this,getString(message('rl:general:TextStopButton')),[],episodeIdx);
                this.StopTrainingFromButtonFlag = true;
            elseif ~this.RequestToStopTraining
                % training stopped from "train" command
                updateFinalInfo(this,reason,value);
            end
            % stop training button
            this.RequestToStopTraining = true;
            this.StopButton.Text = getString(message('rl:general:TextTrainingStopped'));
            this.StopButton.Enable = false;
            
            % update the action message
            setActionMessage(this,'');
        end
        
        %% update with final info
        function updateFinalInfo(this,reason,value,varargin)
            % training stopped reason
            this.FinalReasonLabel.Text = getString(message('rl:general:TextTrainingStoppedBy'));
            this.FinalReasonResultLabel.Text = num2str(reason);
            % training stopped value
            if nargin<4 % training stopped from "train" command                
                this.FinalValueLabel.Text = getString(message('rl:general:TextTrainingStoppedAtValue'));
                this.FinalValueResultLabel.Text = getStopTrainingValueStr(value);
            else % training stopped manully
                this.FinalValueLabel.Text = getString(message('rl:general:TextTrainingStoppedAtEpisode'));
                this.FinalValueResultLabel.Text = num2str(varargin{1});
            end
        end
        
        %% update with saved agent info
        function updateSavedAgentInfo(this,reason,value)
            % saved agent reason
            this.FinalReasonLabel.Text = getString(message('rl:general:TextAgentSavedBy'));
            this.FinalReasonResultLabel.Text = num2str(reason);
            % saved agent value
            this.FinalValueLabel.Text = getString(message('rl:general:TextAgentSavedAtValue'));
            this.FinalValueResultLabel.Text = getStopTrainingValueStr(value);
        end
    end
    
    methods (Access = private)
        % layout of the Episode manager
        
        %% main figure
        function createMainFigure(this)
            % main figure
            this.Figure = uifigure;
            % this.Figure.Visible = 'off';
            this.Figure.Position = [195 120 1565 940]; % [left bottom width height]
            this.Figure.Name = getString(message('rl:general:FigureName'));
            % this.Figure.Resize = 'off';
            this.Figure.CloseRequestFcn = @(src,ed) delete(this);
        end
        
        %% uifigure layout       
        function createLayoutForFigure(this,hasQ0)
            % two columns (plot and info)
            mainGrid = uigridlayout(this.Figure,[1 2]);
            mainGrid.ColumnWidth = {'1x',450};
            % grid on the left (plot and options)
            leftGrid = uigridlayout(mainGrid,[2 1]);
            leftGrid.RowHeight = {'8x','1x'};
            % grid on the left: main axes for plot
            this.Axes = uiaxes(leftGrid);
            hold(this.Axes,'on');
            setTextOnAxis(this,getString(message('rl:general:PlotTitleDefault')),...
                getString(message('rl:general:PlotXLabel')),getString(message('rl:general:PlotYLabel'))); %#ok<*ADMTHDINV>
            % grid on the left: options for plot
            createPlotOptionsPanel(this,leftGrid,hasQ0);
            % grid on the right: main panel
            this.MainPanel = uipanel(mainGrid);
            this.MainPanel.FontSize = 14;
            this.MainPanel.FontWeight = 'bold';
            this.MainPanel.BorderType = 'line';
        end
        
        %% main panel layout
        function createLayoutForInfoPanel(this,trainingOptions,devices,learnRates,hasQ0)
            % grid for info panel          
            grid = uigridlayout(this.MainPanel,[7 1]);
            if numel(learnRates)>1
                grid.RowHeight = {'1.5x', '5.5x', '3.5x','6.5x', '3.5x', '3x','0.5x'};                
            else
                grid.RowHeight = {'1.5x', '6x',  '4x',    '5x',   '4x',  '3x','0.5x'};
            end
            % add 7 panels: progress, episode info, average info, training
            % otps, final info and legend info.
            createProgressBarPanel(this,grid);
            createEpisodeInfoPanel(this,grid,hasQ0);
            createAverageInfoPanel(this,grid,trainingOptions);
            createOptionInfoPanel(this,grid,trainingOptions,devices,learnRates);
            createFinalInfoPanel(this,grid);
            createLegendInfoPanel(this,grid,hasQ0);
            createMessage(this,grid);
        end
        
        %% progress panel
        function createProgressBarPanel(this,grid)
            % Creates axes with no border and no tick marks.
            % Limits are set to [0,1] for x and y.
            panel = uipanel('Parent', grid, 'BorderType', 'none');
            % local grid
            grid2 = uigridlayout(panel,[1 2]);
            grid2.ColumnWidth = {'3x','1x'};
            grid2.RowHeight = {30};
            progressAxes = uiaxes(...
                'Parent', grid2, ...
                'XGrid', 'off', 'YGrid', 'off', ...
                'XTick', [], 'YTick', [], ...
                'XLim', [0,1], 'YLim', [0,1], ...
                'Box', 'on');
            progressAxes.XAxis.Visible = 'off';
            progressAxes.YAxis.Visible = 'off';
            progressAxes.Interactions = [];
            progressAxes.Toolbar.Visible = 'off';
            % background color for progress bar
            patch(progressAxes, ...
                'XData', [0,0,1,1], ...
                'YData', [0,1,1,0], ...
                'FaceColor', [230,230,230]/256, ...
                'EdgeColor', [188,188,188]/256, ...
                'LineWidth', 0.1);
            % progress color
            this.ProgressBar =  patch(progressAxes, ...
                [0,0,0,0], ...
                [0,1,1,0], ...
                [0, 114, 189]/256, ...
                'EdgeColor', 'none');
            % stop button on the right
            this.StopButton = uibutton(grid2,'push');
            this.StopButton.Text = getString(message('rl:general:TextStopTraining'));
            this.StopButton.ButtonPushedFcn = createCallbackFcn(this,@stopTrainingPushed,false);            
        end
        
        %% legend info panel
        function createLegendInfoPanel(this,grid,hasQ0)
            % Creates axes with no border and no tick marks.
            % Limits are set to [0,1] for x and y.
            panel = uipanel('Parent', grid, 'BorderType', 'none');
            panel.BackgroundColor = [1 1 1];
            % local grid
            grid2 = uigridlayout(panel,[1 2]);
            grid2.ColumnWidth = {'1x','2x'};           
            % legend axes
            legendAxes = uiaxes(...
                'Parent', grid2, ...
                'XGrid', 'off', 'YGrid', 'off', ...
                'XTick', [], 'YTick', [], ...
                'XLim', [0,1], 'YLim', [0,1], ...
                'Box', 'on');
            legendAxes.XAxis.Visible = 'off';
            legendAxes.YAxis.Visible = 'off';
            legendAxes.Interactions = [];
            legendAxes.BackgroundColor = [1 1 1];
            legendAxes.Toolbar.Visible = 'off';  
            this.LegendAxes = legendAxes;
            % Episode reward
            line(legendAxes, 'XData', [0,0.5,1], 'YData', [0.82,0.82,0.82], ...
                'LineStyle', '-', ...
                'LineWidth', 0.5, ...
                'Color', '#0072BD', ...
                'Marker', 'o',...
                'MarkerIndices',2);
            % Average reward 
            line(legendAxes, 'XData', [0,0.5,1], 'YData', [0.5,0.5,0.5], ...
                'LineStyle', '-', ...
                'LineWidth', 0.5, ...
                'Color', '#D95319', ...
                'Marker', '*',...
                'MarkerIndices',2);
            % local grid
            grid3 = uigridlayout(grid2,[3 1]);
            uilabel(grid3,'Text','Episode Reward');
            uilabel(grid3,'Text','Average Reward');
            this.LegendTextLabel = uilabel(grid3,'Text','');
            % has Q0
            updateLegendInfo(this,hasQ0);
        end
        
        
        %% episode info panel
        function createEpisodeInfoPanel(this,grid,hasQ0)
            % panel
            panel = uipanel('Parent',grid);
            panel.Title = getString(message('rl:general:EpisodePanelTitle'));
            panel.FontSize = 13;
            panel.FontWeight = 'bold';
            panel.BorderType = 'none';
            % local grid
            if hasQ0
                numRows = 5;   
                rowHeight = {20,20,20,20,20};
            else
               numRows = 4;
               rowHeight = {20,20,20,20};
            end
            grid2 = uigridlayout(panel,[numRows 2]);
            grid2.ColumnWidth = {200,'1x'};
            grid2.RowHeight = rowHeight;
            % episode number
            uilabel(grid2,'text',getString(message('rl:general:TextEpisodeNumber')));
            this.EpisodeNumberLabel = uilabel(grid2,'Text','...');
            % episode reward
            uilabel(grid2,'text',getString(message('rl:general:TextEpisodeReward')));
            this.EpisodeRewardLabel = uilabel(grid2,'Text','...');
            % episode steps
            uilabel(grid2,'text',getString(message('rl:general:TextEpisodeSteps')));
            this.EpisodeStepsLabel = uilabel(grid2,'Text','...');
            % episode Q0
            if hasQ0
                 uilabel(grid2,'text',getString(message('rl:general:TextEpisodeQ0')));
                 this.EpisodeQ0Label = uilabel(grid2,'Text','...');
            end
            % global (total) step count
            uilabel(grid2,'text',getString(message('rl:general:TextTotalNumSteps')));
            this.GlobalStepsLabel = uilabel(grid2,'Text','...');
        end
        
        %% average info panel
        function createAverageInfoPanel(this,grid,trainingOptions)
            % panel
            panel = uipanel('Parent',grid);
            panel.Title = getString(message('rl:general:AveragePanelTtile'));
            panel.FontSize = 13;
            panel.FontWeight = 'bold';
            panel.BorderType = 'none';
            % local grid
            grid2 = uigridlayout(panel,[3 2]);
            grid2.ColumnWidth = {200,'1x'};
            grid2.RowHeight = {20,20,20};
            % average reward
            uilabel(grid2,'text',getString(message('rl:general:TextAverageReward')));
            this.AverageRewardLabel = uilabel(grid2,'Text','...');
            % average steps
            uilabel(grid2,'text',getString(message('rl:general:TextAverageSteps')));
            this.AverageStepsLabel = uilabel(grid2,'Text','...');
            % average window length
            windowLength = trainingOptions.ScoreAveragingWindowLength;
            uilabel(grid2,'text',getString(message('rl:general:TextAverageWindow')));
            this.AverageWindowLengthLabel = uilabel(grid2,'text',num2str(windowLength));
        end
        
        %% option info panel
        function createOptionInfoPanel(this,grid,trainingOptions,devices,learnRates)
            % panel
            panel = uipanel('Parent',grid);
            panel.Title = getString(message('rl:general:OptionsPanelTitle'));
            panel.FontSize = 13;
            panel.FontWeight = 'bold';
            panel.BorderType = 'none';
            % local grid
            if numel(learnRates)>1  
                numRows = 6;
                rowHeight = {20,20,20,20,20,20};
            else
                numRows = 4;
                rowHeight = {20,20,20,20};
            end
            grid2 = uigridlayout(panel,[numRows 2]);
            grid2.ColumnWidth = {200,'1x'};
            grid2.RowHeight = rowHeight;
            % agent uses device and learn rate
            fnames = fields(devices);
            if numel(learnRates)>1  
                % devices
%                 uilabel(grid2,'text',getString(message('rl:general:TextHardwareResourceForActor')));
%                 this.ActorDeviceLabel = uilabel(grid2,'Text',num2str(devices.(fnames{1})));
%                 uilabel(grid2,'text',getString(message('rl:general:TextHardwareResourceForCritic')));
%                 this.CriticDeviceLabel = uilabel(grid2,'Text',num2str(devices.(fnames{2})));
%                 learn rates
%                 uilabel(grid2,'text',getString(message('rl:general:TextLearnRateForActor')));
%                 this.ActorLearnRateLabel = uilabel(grid2,'Text',num2str(learnRates(1)));
%                 uilabel(grid2,'text',getString(message('rl:general:TextLearnRateForCritic')));
%                 this.CriticLearnRateLabel = uilabel(grid2,'Text',num2str(learnRates(2)));
            else
                % device
                uilabel(grid2,'text',getString(message('rl:general:TextHardwareResource')));
                this.DeviceLabel = uilabel(grid2,'Text',num2str(devices.(fnames{1})));
                % learn rate
                uilabel(grid2,'text',getString(message('rl:general:TextLearnRate')));
                this.LearnRateLabel = uilabel(grid2,'Text',num2str(learnRates));
            end
            % maximum number of episodes
            maxEpisodes = trainingOptions.MaxEpisodes;
            uilabel(grid2,'text',getString(message('rl:general:TextOptMaxEpisodes')));
            this.MaxEpisodesLabel = uilabel(grid2,'text',num2str(maxEpisodes));
            % maximum steps per episode
            maxSteps = trainingOptions.MaxStepsPerEpisode;
            uilabel(grid2,'text',getString(message('rl:general:TextOptMaxStepsPerEpisode')));
            this.MaxStepsPerEpisodeLabel = uilabel(grid2,'text',num2str(maxSteps));
        end
        
        %% final info panel
        function createFinalInfoPanel(this,grid)
            % panel
            panel = uipanel('Parent',grid);
            panel.Title = getString(message('rl:general:FinalPanelTitle'));
            panel.FontSize = 13;
            panel.FontWeight = 'bold';
            panel.BorderType = 'none';
            % local grid 
            grid2 = uigridlayout(panel,[3 2]);
            grid2.RowHeight = {20,20,20};
            grid2.ColumnWidth = {200,'1x'};
            % training stopped reason
            this.FinalReasonLabel = uilabel(grid2,'text',getString(message('rl:general:TextTrainingStoppedBy')));
            this.FinalReasonResultLabel = uilabel(grid2,'text','...');
            % training stopped value
%             this.FinalValueLabel = uilabel(grid2,'text',getString(message('rl:general:TextTrainingStoppedAt')));
%             this.FinalValueResultLabel = uilabel(grid2,'text','...');
            % elapsed time
            uilabel(grid2,'text',getString(message('rl:general:TextElapsedTime')));
            this.ElapsedTimeLabel = uilabel(grid2,'text','...');
        end
        
        %% message for parallel
        function createMessage(this,grid)
            % message
            this.MessageLabel = uilabel(grid,'text',' ');
            this.MessageLabel.HorizontalAlignment = 'right';
            this.MessageLabel.FontSize = 13;
        end
        
        %% plot options panel
        function createPlotOptionsPanel(this,grid,hasQ0)
            % panel
            panel = uipanel(grid);
            panel.BorderType = 'none'; %'line' for debugging
            panel.FontWeight = 'bold';
            panel.FontSize = 13;
            panel.Title = 'Plot Options';
            % main grid (3 columns)
            grid = uigridlayout(panel,[1 3]);
            grid.ColumnWidth = {'5x','1x','1x'};
            grid.RowHeight = {25};
            % show Episode Q0?
            this.EnableQ0CheckBox = uicheckbox('parent',grid,'text','Show Episode Q0','Value',hasQ0);
            if hasQ0
                this.ShowEpisodeQ0 = hasQ0;
            else                
               this.EnableQ0CheckBox.Enable = 'off';
            end              
            % show last N episode? 
            this.ShowLastNEpisodesCheckBox = uicheckbox('parent',grid,'text','Show Last N Episodes');
            % N value
            this.ShowLastNEpisodesTextBox = uieditfield('numeric','parent',grid,'Value',100,'Enable','off');
            this.ShowLastNEpisodesTextBox.Limits = [1 Inf];
            this.ShowLastNEpisodesTextBox.RoundFractionalValues = 'on';
            % callbacks
            this.EnableQ0CheckBox.ValueChangedFcn = ...
                @(src,ed) setShowEpisodeQ0(this,ed.Value);
            this.ShowLastNEpisodesCheckBox.ValueChangedFcn = ...
                @(src,ed) toggleHistoryLength(this,ed.Value);
            this.ShowLastNEpisodesTextBox.ValueChangedFcn = ...
                @(src,ed) setHistoryLength(this,ed.Value);
            
        end
        
        %% utility function: plot episode information
        function plotfieldnames = plotEpisodeInfo(this,computedInfo,episodeNum)
            % update the data
            this.RawXData = [this.RawXData,episodeNum];
            this.RawYData = [this.RawYData,computedInfo];
            % update plot
            plotfieldnames = updatePlot(this);
        end
        
        %% utility function: update plot
        function plotfieldnames = updatePlot(this)
            % get the data
            yDataStruct = this.RawYData;
            xData_ = this.RawXData;
            if ~isempty(xData_)
                if ~isfield(yDataStruct,'EpisodeQ0')
                    setShowEpisodeQ0(this,false);
                end
                N = this.HistoryLength;
                % Add error checking
                if ~isinf(N) && numel(xData_) > N
                    xData_ = xData_(end-N+1:end);
                    if isfield(yDataStruct,'SimulationInfo') % when plotting with saved agent/training result.
                        yDataStruct = getLastNTrainingStats(this,yDataStruct);
                    else
                        yDataStruct = yDataStruct(end-N+1:end);
                    end
                end
                lineStyle = {'-','-','--'}; % EpisodeReward, AverageReward and EpisodeQ0
                marker = {'o','*','x'};
                color = {'#0072BD','#D95319','#EDB120'}; % Use default colors
                if this.ShowEpisodeQ0
                    expectedfieldnames = {'EpisodeReward','AverageReward','EpisodeQ0'};
                else
                    expectedfieldnames = {'EpisodeReward','AverageReward'};
                end
                fieldnames = fields(yDataStruct);
                % the same order as expectedfieldnames
                plotfieldnames = intersect(expectedfieldnames,fieldnames,'stable');
                % iterate over the fields names
                for i = 1:numel(plotfieldnames)
                    % extract data from the struct
                    fname = plotfieldnames{i};
                    yData_ = [yDataStruct.(fname)];
                    % if we have already plotted this object, go ahead and
                    % append xdata, ydata, otherwise create a new plot
                    % object
                    if isfield(this.PlotObjects,fname) && isvalid(this.PlotObjects.(fname))
                        gobj = this.PlotObjects.(fname);
                        gobj.XData = xData_;
                        gobj.YData = yData_;
                    else
                        gobj = plot(this.Axes,xData_,yData_,'LineStyle',lineStyle{i},'Marker',marker{i},'Color',color{i});
                        this.PlotObjects.(fname) = gobj;
                    end
                end
                % discard plots that are no longer used (like Q0)
                discardplots = setdiff(fieldnames,expectedfieldnames);
                for i = 1:numel(discardplots)
                    fname = discardplots{i};
                    if isfield(this.PlotObjects,fname)
                        delete(this.PlotObjects.(fname));
                        this.PlotObjects = rmfield(this.PlotObjects,fname);
                    end
                end
            end           
        end
        
        %% utility function: get last N epsidoe Training stats
         function result = getLastNTrainingStats(this,stats)
            N = this.HistoryLength; 
            % Clean up unused training statistics            
            result.EpisodeIndex     = stats.EpisodeIndex   (end-N+1:end);
            result.EpisodeReward    = stats.EpisodeReward  (end-N+1:end);
            result.EpisodeSteps     = stats.EpisodeSteps   (end-N+1:end);
            result.AverageReward    = stats.AverageReward  (end-N+1:end);
            result.TotalAgentSteps  = stats.TotalAgentSteps(end-N+1:end);
            result.AverageSteps     = stats.AverageSteps   (end-N+1:end);
            if isfield(stats,'EpisodeQ0')
                result.EpisodeQ0 = stats.EpisodeQ0(end-N+1:end);
            end
            if ~isempty(stats.SimulationInfo)
                result.SimulationInfo = stats.SimulationInfo(end-N+1:end);
            end
         end
        
        %% utility function: update episode information display
        function updateEpisodeInfo(this,episodeNum,info,maxEpisodes,plotfieldnames,averageSteps,globalStepCount)
            % update labels for episode info panel
            this.EpisodeNumberLabel.Text = num2str(episodeNum);
            this.EpisodeRewardLabel.Text = num2str(info.EpisodeReward(end));
            this.EpisodeStepsLabel.Text = num2str(info.EpisodeSteps(end));
            if ismember('EpisodeQ0',plotfieldnames)
                this.EpisodeQ0Label.Text = num2str(info.EpisodeQ0(end));
            end
            this.GlobalStepsLabel.Text = num2str(globalStepCount);
            % update labels for average info panel
            this.AverageRewardLabel.Text = num2str(info.AverageReward(end));
            this.AverageStepsLabel.Text = num2str(averageSteps);
            % update progress bar
            progressRatio = max(0.01,episodeNum/maxEpisodes);            
            this.ProgressBar.XData = [0,0,progressRatio,progressRatio];
        end
        
        %% utility function: update legend info
        function updateLegendInfo(this,showQ0)
            if isempty(this.LegendLine)
                this.LegendLine = line(this.LegendAxes, 'XData', [0,0.5,1], 'YData', [0.18,0.18,0.18], ...
                    'LineStyle', '--', ...
                    'LineWidth', 0.5, ...
                    'Color', '#D95319', ...
                    'Marker', 'x',...
                    'MarkerIndices',2);
            end
            if showQ0
                % Episode Q0                
                this.LegendLine.Visible = 'on';
                this.LegendTextLabel.Text = 'Episode Q0';
            else
                this.LegendLine.Visible = 'off';
                this.LegendTextLabel.Text = '';
            end
        end
    end
    
    methods (Access = private)
        % call back
        function stopTrainingPushed(this)
            stopTraining(this);
        end
    end
    
end

%% local utility functions
function str = getStopTrainingValueStr(value)
if isa(value,'function_handle')
    str = ['@' func2str(value)];
else
    str = num2str(value);
end
end

function stats = removeIncompleteEpisodeFromStats(stats,hasCritic)
stats.EpisodeIndex   (end) = [];
stats.EpisodeReward  (end) = [];
stats.EpisodeSteps   (end) = [];
stats.AverageReward  (end) = [];
stats.TotalAgentSteps(end) = [];
stats.AverageSteps   (end) = [];
if hasCritic
    stats.EpisodeQ0  (end) = [];
end
stats.SimulationInfo(end) = [];
end


