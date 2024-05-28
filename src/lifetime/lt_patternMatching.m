function [gridMLE_tau,gridMLE_b] = lt_patternMatching(tail_decay,tcspc_t, tcspc_bin, max_lt)
%LT_PATTERNMATCHING Determain flourecent lifetime with pattern matching. Based on
%MLE grid search by C. Thile

% gridMLE_tau lifetimes corresponding to lable
% gridMLE_b   background values corresponding to label

% tail_decay matrix of decays to fit. 1st dim is label, 2nd time
% tcspc_t    edges or time points for TCSPC histogram in ns
% tcspc_bin  width of individual TCSPC bin in ns
% max_lt     Max lifetime for pattern generation

%% Simulate TCSPC decays

% Normalised monoexponetial decay with background
pfun_monoexp    = @(tau) tcspc_bin(:).*exp(-tcspc_t(:)./tau);
pfun_monoexpBG  = @(tau,b) b./numel(tcspc_t(:))+(1-b).*pfun_monoexp(tau)./sum(pfun_monoexp(tau),1);

%% Calculate patterns

% define grid
grid_lts = linspace(0.01,max_lt,250); % lifetime values
% grid_lts = [1.4, 2.4, 3];
grid_bs  = linspace(0.0,0.6,30);  % background values

% get all combinations
[grid_ltmat,grid_bmat] = meshgrid(grid_lts,grid_bs);
grid_ltmat = grid_ltmat(:)';
grid_bmat  = grid_bmat(:)';

grid_logdecays = log(pfun_monoexpBG(grid_ltmat,grid_bmat));

%% Find best matching pattern
% Calculate Lambda and directly find the index of the minimum value
[~,grid_ind] = min(-tail_decay(:,1:end) * grid_logdecays,[],2);

gridMLE_tau = grid_ltmat(grid_ind).';
gridMLE_b   = grid_bmat(grid_ind).';

end

