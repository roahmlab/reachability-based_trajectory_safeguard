classdef PolicyInstance < handle
% POLICYINSTANCE
% set and get a policy instance in this class workspace. Useful for
% resolving policies before simulation

% Revised: 7-1-2019
% Copyright 2019 The MathWorks, Inc.

    properties (Access = private)
        Policy
    end
    methods (Access = private)
        function this = PolicyInstance()
        end
    end
    methods (Static)
        function set(policy)
            if ~isempty(policy)
                validateattributes(policy,"rl.agent.AbstractPolicy",{},"setgetPolicyInstance","policy");
            end
            this = rl.util.PolicyInstance.getInstance();
            this.Policy = policy;
        end
        function policy = get()
            this = rl.util.PolicyInstance.getInstance();
            policy = this.Policy;
        end
    end
    methods (Static,Access = private)
        function this = getInstance()
            persistent this_
            if isempty(this_) || ~isvalid(this_)
                this_ = rl.util.PolicyInstance;
            end
            this = this_;
        end
    end
end