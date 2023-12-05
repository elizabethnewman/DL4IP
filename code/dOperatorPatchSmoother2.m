classdef dOperatorPatchSmoother2
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
        L       
        patchSize
        imgSize
        idx         % permutation
        sizes
        transposed
    end

    methods

        %% initialize D operator
        function obj = dOperatorPatchSmoother2(patchSize, imgSize)
            obj.L          = patchSmoother2(patchSize, imgSize);
            obj.patchSize  = patchSize;
            obj.imgSize    = imgSize;
            obj.sizes      = size(obj.L);
            obj.transposed = false;

            % % permutation index

            

            % This is the transposed permutation, 
            % Let X = patchify(img).  Then, vec(X(idx)) = vec(img)
            % tmp     = reshape(1:prod(imgSize),prod(patchSize),[]);
            % tmp     = col2im(tmp,patchSize,imgSize,'distinct');
            % obj.idx = tmp(:)';

            % test
            % img = rand(imgSize);
            % X = im2col(img,patchSize,"distinct");
            % t = norm(X(obj.idx(:)) - img(:),'fro');

            
            % This is the forward permutation: vec(patchify(img)) = vec(img(idx))

            % tmp     = reshape(1:prod(imgSize),imgSize);
            % tmp     = im2col(tmp,patchSize,"distinct");
            % obj.idx = tmp(:)';

            % % test
            % img = rand(imgSize);
            % X = im2col(img,patchSize,"distinct");
            % t = norm(X(:) - img(obj.idx(:)),'fro');

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
                x = obj.L' * x;
            else
                x = obj.L * x;
            end
        end

    end
end
