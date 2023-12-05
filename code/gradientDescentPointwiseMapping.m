function[x,info,trajectory] = gradientDescentPointwiseMapping(A,b,w,dw,options)
%

if nargin == 0, runMinimalExample; return; end

% ----------------------------------------------------------------------- %
% setup options
defaultOptions = struct('z0','none', 'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, ...
    'absTol', 1e-8, 'relTol', 1e-8,...
    'NoStop', 'off','verbose',false);

% overwrite
if exist('options','var') && ~isempty(options)
    optFields = fields(defaultOptions);
    for k = 1:length(optFields)
        if ~isfield(options,optFields{k})
            options.(optFields{k}) = defaultOptions.(optFields{k});
        end
    end
else
    options = defaultOptions;
end

optOptions = optimset('TolX',1e-7);

% ----------------------------------------------------------------------- %
% setup printouts
info.headers = {'iter', 'f', '|x1-x0|', '|r|', '|r|/|b|','|g|', '|g|/|g0|', '|s|', '|s|/|s0|','nnz(x)','nnz(x)/numel(x)','alpha','stopCrit'};
info.frmt    = {'%-15d','%-15.2e','%-15.2e','%-15.2e','%-15.2e','%-15.2e','%-15.2e','%-15.2e','%-15.2e','%-15d','%-15.4f','%-15.2e','%d%d%d%d%d%d'};
info.values  = zeros(0,length(info.headers)+5);
trajectory.z = zeros(size(A,2),size(b,2),0);

% ----------------------------------------------------------------------- %
% initial evaluation

% function to evaluate
f = @(A,b,z) 0.5 * norm(A * w(z) - b,'fro').^2;

% initialization (TODO: find good initialization)
if isempty(options.z0) || strcmp(options.z0,'none')
    % z = abs(A \ b);
    z = rand(size(A,2),size(b,2));
    % coeffx0 = A * ones(size(A,2),1); 
    % alpha   = (coeffx0' * b) / norm(coeffx0)^2;
    % z       = log(alpha) * ones(size(A,2),1);

    % MRNSD initial guess
    % coeffx0 = A * ones(size(A,2),1); 
    % alpha   = (coeffx0' * b) / norm(coeffx0)^2;
    % 
    % tmp = fzero(@(z) w(z) - alpha,1);
    % % tmp = fzero(@(z) exp(z) - alpha,1);
    % 
    % z   = tmp * ones(size(A,2),1);

else
    z = options.z0;
end

if nargout > 2, trajectory.z = cat(3,trajectory.z,z); end

% setup stopping criteria
nrm_b  = norm(b,'fro');
nrm_r0 = norm(A * w(z) - b,'fro');
nrm_g0 = norm(A' * (A * w(z) - b),'fro');
nrm_s0 = norm(dw(z) .* (A' * (A * w(z) - b)),'fro');
numelX = size(A.A,2);

% store values
info.values = cat(1,info.values,[0,f(A,b,z),0,nrm_r0,nrm_r0 / nrm_b,nrm_g0,1,nrm_s0,1,nnz(w(z)),nnz(w(z))/numelX,0,zeros(1,6)]);

% print
if options.verbose
    fprintf([repmat('%-15s',1,length(info.headers)), '\n'], info.headers{:});
    fprintf([info.frmt{:}, '\n'], info.values)
end

% ----------------------------------------------------------------------- %
% main iteration
for k = 1:options.MaxIter 
    % clear stopping criteria
    stoppingCritera = zeros(1,6);
    
    % compute residual
    r     = A * w(z) - b;
    nrm_r = norm(r,'fro');
    stoppingCritera(1) = (nrm_r         < options.absTol);
    stoppingCritera(2) = (nrm_r / nrm_b < options.relTol);

    % compute gradient w.r.t. x
    g     = A' * r;
    nrm_g = norm(g,'fro');
    stoppingCritera(3) = (nrm_g          < options.absTol);
    stoppingCritera(4) = (nrm_g / nrm_g0 < options.relTol);

    % compute search direction
    s     = -dw(z) .* g;
    nrm_s = norm(s,'fro');
    stoppingCritera(5) = (nrm_s          < options.absTol);
    stoppingCritera(6) = (nrm_s / nrm_s0 < options.relTol);

    % stopping criteria
    if all(stoppingCritera > 0)
        fprintf('stopping criteria met!\n')
        break
    end

    % compute optimal step size
    h     = @(alpha) 0.5 * norm(A * w(z + alpha * s) - b,'fro').^2;
    % alpha = fminbnd(h,0,1e8,optOptions);

    alphas = [0,logspace(-8,8,200)];
    tmp = zeros(size(alphas));
    for i = 1:length(alphas)
        tmp(i) = h(alphas(i));
    end
    [~,idx] = min(tmp);
    alpha   = alphas(idx);

    % h     = @(alpha) 0.5 * norm(A * (w(z) + alpha * s) - b,'fro').^2;
    % alpha = fminbnd(h,0,1e2);
    % u = s(:) .* dw(z(:));
    % alpha = -u(:)' * g(:) / (u(:)' * u(:));
    
    % if ~exist('superres_gdnn','dir'), mkdir('superres_gdnn'); end
    % % 
    % tmp3 = reshape(w(z),400,[]);
    % tmp2 = col2im(A.D * tmp3,[16,16],[512,512],'distinct');
    % fig = figure(1); clf;
    % imagesc(tmp2); colormap parula; % colorbar; 
    % clim([0,1])
    % axis('off'); axis('image');
    % set(gca,'FontSize',24)
    % exportgraphics(fig,['superres_gdnn/img_',num2str(k),'.png'])
    % title(k)
    % drawnow

    % update
    zOld = z;
    z    = z + alpha * s;
    

    if nargout > 2, trajectory.z = cat(3,trajectory.z,z); end


    % update
    info.values = cat(1,info.values,[k,f(A,b,z),norm(z - zOld,'fro'),nrm_r, nrm_r / nrm_b, nrm_g,nrm_g / nrm_g0,nrm_s,nrm_s / nrm_s0,nnz(w(z)),nnz(w(z))/numelX,alpha,stoppingCritera]);
    if options.verbose
        fprintf([info.frmt{:}, '\n'], info.values(end,:))
    end

    if alpha < 1e-16, disp('alpha < 0'); break; end


end

% return x
x = w(z);

end


function runMinimalExample

rng(42);





% setup simple problem
A = [1,0;0,1];
% A     = [1,0;0,3];
% A     = [2,1;1,4];
xTrue = [0.05;0.05];
b     = A * xTrue;


[cc,aa,pp] = deal(0,1,1);
f  = @(x) (x - cc > 0) .* (exp(aa * (x - cc)) - aa * (x - cc).^pp - 1);
df = @(x) (x - cc > 0) .* (aa * exp(aa * (x - cc)) - pp * aa * (x - cc).^(pp - 1));


c = 1;
myMaps.w     = {@(x) exp(x),   @(x) c * x.^2,    @(x) c * abs(x),   @(x) max(x,0), @(x) f(x)};
myMaps.dw    = {@(x) exp(x),   @(x) 2 * c * x,   @(x) c * sign(x),  @(x) (x > 0),  @(x) df(x)};
myMaps.names = {'e^x', 'x^2', '|x|', 'max(x,0)', 'f(x)'};


options = struct('z0',[1;2], 'MaxIter',100, 'x_true','none', ...
    'NoiseLevel','none', 'eta',1.01, 'NE_Rtol',1e-12, ...
    'absTol', 1e-8, 'relTol', 1e-8,...
    'NoStop', 'off','verbose',false);

results = cell(1,length(myMaps.w));
for i = 1:length(myMaps.w)
    [x,info,trajectory] = gradientDescentPointwiseMapping(A,b,myMaps.w{i},myMaps.dw{i},options);
    results{i}.x = x;
    results{i}.info = info;
    results{i}.trajectory = trajectory;
end


% convergence in x-space
[aLow,aHigh] = deal(-0.25,3); 
N            = 100;
[xx,yy]      = meshgrid(linspace(aLow,aHigh,N), linspace(aLow,aHigh,N));
zz           = 0.5 * sum((A * [xx(:)'; yy(:)'] - b(:)).^2,1);

fig = figure(1); clf; fig.Name = 'x-Space';

contourf(xx,yy,reshape(zz,size(xx)),'HandleVisibility','off'); colormap gray
hold on;

myColors = get(gca,'ColorOrder');
for i = 1:length(myMaps.w)
    xTraj = squeeze(myMaps.w{i}(results{i}.trajectory.z));
    plot(xTraj(1,:),xTraj(2,:), 'Color',myColors(i,:),'LineWidth',6,'DisplayName',myMaps.names{i})
    plot(xTraj(1,1),xTraj(2,1),'o', 'MarkerFaceColor',myColors(i,:),'MarkerSize',20,'MarkerEdgeColor','k','LineWidth',1,'HandleVisibility','off')
end

yline(0,'w','Linewidth',2)
xline(0,'w','Linewidth',2)

plot(xTrue(1),xTrue(2),'kp','LineWidth',2,'MarkerFaceColor','y','MarkerSize',20,'HandleVisibility','off')

hold off;
xlabel('x1')
ylabel('x2')
legend()
set(gca,'FontSize',18)


% convergence in z-space
[aLow,aHigh] = deal(-3,3); 
N       = 100;
[xx,yy] = meshgrid(linspace(aLow,aHigh,N), linspace(aLow,aHigh,N));

fig = figure(2); clf; fig.Name = 'z-Space';

myColors = get(gca,'ColorOrder');
for i = 1:length(myMaps.w)
    subplot(1,length(myMaps.w),i);
    zz = 0.5 * sum((A * myMaps.w{i}([xx(:)'; yy(:)']) - b(:)).^2,1);
    contourf(xx,yy,reshape(zz,size(xx)),'HandleVisibility','off'); colormap gray
    hold on;
    zTraj = squeeze(results{i}.trajectory.z);
    plot(zTraj(1,:),zTraj(2,:), 'Color',myColors(i,:),'LineWidth',6,'DisplayName',myMaps.names{i})
    plot(zTraj(1,1),zTraj(2,1),'o', 'MarkerFaceColor',myColors(i,:),'MarkerSize',20,'MarkerEdgeColor','k','LineWidth',1,'HandleVisibility','off')
    title(myMaps.names{i})
    xlabel('z1')
    ylabel('z2')
    hold off;
    set(gca,'FontSize',18)
end




fig = figure(3); clf; fig.Name = 'Convergence';
for i = 1:length(myMaps.w)
    semilogy(results{i}.info.values(:,2), 'LineWidth',6,'DisplayName',myMaps.names{i})
    hold on;
end
xlabel('iteration')
ylabel('$\frac{1}{2} * \|A * w(z) - b\|_2^2$','Interpreter','latex')
legend()
grid
set(gca,'FontSize',18)

end

