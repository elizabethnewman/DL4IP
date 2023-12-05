classdef dOperatorDictionaryRepresentation
    % classdef dOperator
    %
    % Authors:
    %   (c) Matthias Chung (matthias.chung@emory.edu) and Rosemary Renaut in October 2021
    %
    % MATLAB Version: 9.11.0.1769968 (R2021b)
    %
    % Description:
    %
    %  D = dOperator(fctn_handle, inDim,outDim)
    %
    % Properties:
    %     A          - function handle @(x,transFlag) ...
    %     sizes      - provides the input and output dimensions [outDim, inDim]
    %     transposed - flag if operator is transposed or not
    %
    % Example:
    %     D = dOperatorFctn(A,inDim,outDim)
    %
    % References:


    properties
        A       % matrix or operator
        D       % dictionary
        patchSize
        imgSize
        sizes
        lambda  = 0.0       % Tikhonov regularization parameter
        L                   % Tikhonov regularization matrix or operator
        transposed
    end

    methods

        function obj = dOperatorDictionaryRepresentation(A,D,patchSize,imgSize,lambda,L)
            obj.A          = A;
            obj.D          = D;
            obj.patchSize  = patchSize;
            obj.imgSize    = imgSize;
            
            obj.transposed = false;

            if exist('lambda','var'), obj.lambda = lambda; end
            if ~exist('L','var') || isempty(L)
                obj.L = zeros(0,size(A,2));
            else
                obj.L = L;
            end

            obj.sizes = [size(A,1) + size(obj.L,1),size(D,2) * prod(floor(imgSize ./ patchSize))];
        end

        %% transpose method
        function obj = ctranspose(obj)
            obj.transposed = not(obj.transposed);
            obj.sizes      = flip(obj.sizes);
        end

        %% size method
        function s = size(obj,dim)
            if nargin < 2
                s = obj.sizes;
            else
                s = obj.sizes(dim);
            end
        end

        %% mtimes method
        function x = mtimes(obj, x)
            if obj.sizes(2) ~= size(x,1), error('size mismatch'); end 

            if obj.transposed
                y = obj.A' * x(1:size(obj.A,1)) + sqrt(obj.lambda) * (obj.L' * x(size(obj.A,1)+1:end));
                x = dictionaryOperator(y,1,obj.D,obj.patchSize,obj.imgSize);
            else
                y = dictionaryOperator(x,0,obj.D,obj.patchSize,obj.imgSize);
                x = [obj.A * y; sqrt(obj.lambda) * (obj.L * y)];
            end
        end

    end
end
