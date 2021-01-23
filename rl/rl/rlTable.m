classdef rlTable
    % rlTable - Create a value or Q table that represents the value for a
    % given state or state/action pair
    %
    % V(s)
    %    V = rlTable(ObservationInfo)
    %
    % Q(s,a)
    %    Q = rlTable(ObservationInfo,ActionInfo)
    %
    %    where ObservationInfo and ActionInfo are discrete data
    %    specficiations rlFiniteSetSpec.
    %
    %    See also: rlQAgent, rlSARSAAgent, rlTableRepresentation
    
    % Copyright 2017-2018 The MathWorks Inc.
    
   
    %% Public Properties
    properties (Dependent)
        % Value or Q table matrix
        Table
    end
    
    %% Protected properties
    properties (Access = private)
        % Value or Q table matrix
        Table_   
        % Defines Observation space
        ObservationInfo rl.util.rlFiniteSetSpec
        % Defines Action space
        ActionInfo rl.util.rlFiniteSetSpec  
        % Type of table "Value" or "Q"
        Type string
    end
    
    %% Public methods
    methods                
        % constructor, creating a q table from stategrid and actiongrid
        function obj = rlTable(ObservationInfo,ActionInfo)
            % Revisit for case when ObservationInfo/ActionInfo is an array
            narginchk(1, 2)
            if nargin == 1
                % Value Table
                validateInput(obj, ObservationInfo)
                obj.Type = "Value";
                obj.ObservationInfo = ObservationInfo;
                obj.Table_ = zeros(numel(obj.ObservationInfo.Elements),1);
            else
                % Q Table
                validateInput(obj, ObservationInfo)
                validateInput(obj, ActionInfo)
                
                obj.Type = "Q";
                obj.ObservationInfo = ObservationInfo;
                obj.ActionInfo  = ActionInfo;
                obj.Table_ = zeros(numel(obj.ObservationInfo.Elements), numel(obj.ActionInfo.Elements));
            end
        end  
        
        function Value = evaluate(obj,Data)
            % Data is size [observation] or [observation action]
            % REVISIT: Data could be a mini-batch
            if strcmp(obj.Type,"Q")
                Value = evaluateQ(obj,Data);
            else
                Value = evaluateV(obj,Data);
            end
        end
        function [obsinfo,actinfo] = getInfo(obj)
            obsinfo = obj.ObservationInfo;
            actinfo = obj.ActionInfo;
        end
        % Compute gradient of output with respect to parameters
        function grad = gradientWithRespectToParameters(obj,Data)
            state_idx = getElementIndex(obj.ObservationInfo,Data{1});
            if strcmp(obj.Type,"Q")
                action_idx = getElementIndex(obj.ActionInfo,Data{2});
            else
                action_idx = ones(size(state_idx));
            end
            [m,n] = size(obj.Table);
            grad = zeros(m,n);
            for ct = 1:numel(state_idx)
                grad(state_idx(ct), action_idx(ct))=grad(state_idx(ct), action_idx(ct))+1;
            end
            grad = grad/numel(state_idx);
        end
        
        
        % Compute grad of loss with respect to parameters
        function grad = gradientOfLossWithRespectToParameters(obj,Data,Target)
            state_idx = getElementIndex(obj.ObservationInfo,Data{1});
            if strcmp(obj.Type,"Q")
                action_idx = getElementIndex(obj.ActionInfo,Data{2});
            else
                action_idx = ones(size(state_idx));
            end
            [m,n] = size(obj.Table);
            grad = zeros(m,n);
            for ct = 1:numel(Target)
                grad(state_idx(ct), action_idx(ct))=grad(state_idx(ct), action_idx(ct))-2*(Target(ct)-obj.Table(state_idx(ct), action_idx(ct)));
            end
            grad = grad/numel(Target);
        end
        
    end
    
    %% Set/Get Methods
    methods
        function obj = set.Table(obj,Value)
            % Set the values for the table matrix
            if isnumeric(Value) && isequal(size(Value),size(obj.Table_))
                obj.Table_ = Value;
            else
                error(message('rl:agent:errSetRLTable'))
            end
        end
        
        function Value = get.Table(obj)
            % Get table matrix
            Value = obj.Table_;
        end
    end
    
    %% Protected methods
    methods (Access = private)
        function validateInput(obj, Input)
            % Validate inputs to the constructor
            if ~(isa(Input,'rl.util.rlFiniteSetSpec')  && isscalar(Input))
                error(message('rl:agent:errValidateRLTableInput'))
            end
        end
        
        function Value = evaluateQ(obj,Data)
            % Evaluate Q-Table where Data is [s,a]
            % state_idx = strcmp(Data{1},obj.ObservationInfo.Elements);
            state_idx = getElementIndex(obj.ObservationInfo,Data{1});
            % action_idx = strcmp(Data{2},obj.ActionInfo.Elements);
            action_idx = getElementIndex(obj.ActionInfo,Data{2});
            for ct = 1: numel(state_idx)
                Value(ct,1) = obj.Table_(state_idx(ct),action_idx(ct));
            end
        end
        
        function Value = evaluateV(obj,Data)
            % Evaluate V-Table where Data is [s]
            % state_idx = strcmp(Data{1},obj.ObservationInfo.Elements);
            state_idx = getElementIndex(obj.ObservationInfo,Data{1});
            for ct = 1: numel(state_idx)
                Value(ct,1) = obj.Table_(state_idx(ct),1);
            end
        end
    end
end
 