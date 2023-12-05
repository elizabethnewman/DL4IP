classdef dOperatorPatchSmoother
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
        function obj = dOperatorPatchSmoother(patchSize, imgSize)
            obj.L          = patchSmoother(patchSize, imgSize);
            obj.patchSize  = patchSize;
            obj.imgSize    = imgSize;
            obj.sizes      = [size(obj.L,2) * prod(patchSize),prod(imgSize)];
            obj.transposed = false;

            % % permutation index

            

            % This is the transposed permutation, 
            % Let X = patchify(img).  Then, vec(X(idx)) = vec(img)
            tmp     = reshape(1:prod(imgSize),prod(patchSize),[]);
            tmp     = col2im(tmp,patchSize,imgSize,'distinct');
            obj.idx = tmp(:)';

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
                % kron(L,eye(p * q)) * P' * vec(x)

                % unvectorized
                x = reshape(x,prod(obj.patchSize),[]);

                % apply Kronecker product transpose
                x = x * obj.L';

                % apply transpose of permutation
                x = col2im(x,obj.patchSize,obj.imgSize,"distinct");

                % vectorize image in proper order
                x = x(obj.idx(:));
            else
                % kron(L',eye(p * q)) * P * vec(x)

                % this is effectively a permutation of vec(x)
                x = im2col(reshape(x,obj.imgSize),obj.patchSize,"distinct");
                
                % this this is applying the Kronecker product
                x = x * obj.L;

                % we want to vectorize this result, no permutation
                x = x(:);
            end
        end

    end
end
