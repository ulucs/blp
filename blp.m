load("nevofiles/ps2.mat"); load("nevofiles/iv.mat");
ns = 20; nmkt = 94; nbrn = 24; np = 4; n_inst = 20;
% keeping data always in the long format creates unnecessary extra work
% these functions aim to widen the long data
widen = @(x) permute(reshape(x, nbrn, nmkt, []), [2 1 3]);
widens = @(x) permute(reshape(x, nmkt, 1, ns, np), [1 4 2 3]);

% dimensions: market * parameter (if exists) * brand * sample
% parameter is second since multipliaction will flatten it
ws = widen(s_jt); y = ws./sum(ws, 2); 
x2w = permute(widen(x2), [1 3 2]);
vs = widens(v); ds = widens(demogr);
IV = [iv(:,2:(n_inst+1)) x1(:,2:(nbrn+1))];

% define the initial values and the enforced shape: zeros are by design
th_init = [0.3302   5.4819   0       0.2037   0;
       2.4526  15.8935  -1.2000  0        2.6342;
       0.0163  -0.2506   0       0.0511   0;
       0.2441   1.2650   0      -0.8091   0];
[it, jt, th0] = find(th_init); shp = @(th) sparse(it, jt, th);

% now minimize loss with given initial values, 150 iterations is quite low
% but I use the example code's optimset to obtain exact results
opts = optimset('GradObj','off','MaxIter',150,'Display','iter','TolFun',0.1,'TolX',0.01);
th2 = fminsearch(@(th) blp_loss(shp(th), ws, x1, x2w, IV, vs, ds), th0, opts);
shp(th2) % print the result, maybe the linear coefficients too?

function loss = blp_loss(th, ws, x1, x2w, IV, vs, ds)
nmkt = size(ws, 1); nbrn = size(ws, 2); 
pt = @pagemtimes; smax = @(x) exp(x)./(1+sum(exp(x),2));
% calculate nonlinear part of utility; replaces mufunc.m
% pt multiplies the matrices across each brand and sample
mu = squeeze(pt(x2w.*vs, th(:, 1)) + sum(x2w.*pt(ds, th(:, 2:end)'), 2));
% get market shares given delta; replaces mktsh.m, ind_sh.m
psh = @(d) mean(smax(d+mu), 3); % flatten samples by taking the mean (MC Integration)
% use the fixed point iteration to obatin delta, replaces meanval.m
delta = converge(@(d) d + log(ws) - log(psh(d)), zeros(size(ws)), 1e-6);
% reshape delta, and run IV-GMM to obtain the linear part
dl = reshape(delta', nbrn*nmkt, 1); % this is just inverse of widen
% IV estimation with GMM, the bundled function is optimized
[~, ~, loss] = lscov(IV'*x1, IV'*dl, IV'*IV); % maybe save the linear vars?
end
