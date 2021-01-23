classdef QuadraticLayer < nnet.layer.Layer
% QUADRATIC
% y = [x^2 x*u u^2]'

% Copyright 2018 The MathWorks, Inc.

% Revised: 10-30-2018
% Copyright 2018 The MathWorks, Inc.

% check with:
%   quad = rl.layer.QuadraticLayer('Name','quad')
%   checkLayer(quad,[15 1],'ObservationDimension',4)
%   checkLayer(quad,[1 1 15],'ObservationDimension',4)

    methods
        function this = QuadraticLayer(varargin)
            
            parser = inputParser;
            addParameter(parser,'Name',"quadratic",...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Name'));
            addParameter(parser,'Description',getString(message('rl:general:QuadraticLayerDescription')),...
                @(val)validateattributes(val,{'string','char'},{'scalartext'},'','Description'));
            parse(parser,varargin{:});
            r = parser.Results;
            
            this.Type = "QuadraticLayer";
            this.Name = r.Name;
            this.Description = r.Description;
        end
        function y = predict(~,u)
            %forward : the result of evaluating this node
            %
            % evaluate this node in the graph
            
            % validate that u is "vector based" and get the permute index
            [pidx,sz,vdim,nbatch] = localGetPermuteIdx(u);
            
            u = permute(u,pidx);
            
            % get the number of samples
            N  = size(u,2);
            % get the number of z elements (per batch)
            Nu = size(u,1);
            Nq = sum(1:Nu);
            TRI = triu(true(Nu));
            
            % the basis can be expressed as follows
            %
            %   basis = [x u]'*[x u] = [x*x x*u
            %                           x*u u*u]
            %
            % if we extract the upper traingular portion
            %
            % basis = [x*x x*u u*u]'
            
            u  = reshape(u,Nu,1,N);
            ut = reshape(u,1,Nu,N);
            uhat = u.*ut;
            y = uhat(repmat(TRI,1,1,N));
            y = reshape(y,Nq,N);
            
            % might need to update if batch dim changes
            Ny = size(y,1);
            rdim = ones(1,max(numel(sz),4));
            rdim(4) = nbatch;
            rdim(vdim) = Ny;
            y = reshape(y,rdim);
        end
        function dydu = backward(this,u,y,dY,memory) %#ok<INUSL,INUSD>
            %backward: backward pass
            %
            % [dYdu1,dYdu2,...dYdun,dYdq] = backward(this,u1,u2,...,un,y,dY,memory)
            
            [pidx,sz,vdim,nbatch] = localGetPermuteIdx(u);
            dY = permute(dY,pidx);
            u  = permute(u ,pidx);
            
            % UT is defined as the upper triangular operation that extracts
            % the upper triangular elements of the input and puts it in a
            % vector
            %
            % y = UT(a) = UT(z*z')
            %
            % dyda = upper triangular matrix(dY)
            % dydz = (dyda + dyda')*z (see below)
            
            [Nu,N] = size(u);
            TRI = triu(true(Nu));
            
            % push the gradients onto the basis inputs
            dydu = zeros(Nu,N,'like',dY);
            dyda = zeros(Nu  ,'like',dY);
            for i = 1:N
                
                % a = u*u'
                % dadu = dyda*u + (u'*dydu)' = (dydu + dydu')*z
                %   where dyda is the chain gradient
                
                dyda(TRI) = dY(:,i);
                dydu(:,i) = (dyda + dyda')*u(:,i);  
            end
            
            % reshape the gradients
            rdim = ones(1,max(numel(sz),4));
            rdim(4) = nbatch;
            rdim(vdim) = Nu;
            dydu = reshape(dydu,rdim);
        end
    end
end
function [pidx,sz,vdim,nbatch] = localGetPermuteIdx(u)
% find a permute index so we can leverage matrix multiplication y = a*u
sz = size(u);
sz1 = sz;
if numel(sz1) > 3
    nbatch = sz1(4);
    sz1(4) = [];
    [~,vdim] = max(sz1);
    pidx = [vdim,4,setdiff(1:numel(sz),[vdim,4])];
else
    nbatch = 1;
    [~,vdim] = max(sz1);
    pidx = [vdim,setdiff(1:numel(sz),vdim)];
end
if nnz(sz1 > 1) > 1
    error(message('rl:general:QuadraticLayerExpectsVector',mat2str(sz1)));
end
end