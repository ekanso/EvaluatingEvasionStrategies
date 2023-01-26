function est_fun = BeyesInferenceEstimator(prior_pdf,noise_pdf)

% phi_grid = 0:pi/1000:2*pi;
est = zeros(201,1);
for i = 0:200
    obs = i*pi/100;
    x_grid = 0:pi/500:2*pi;
    pr = prior_pdf(x_grid);
    likelihood = noise_pdf(x_grid,obs);
    post = pr.*likelihood;
    [~,index] = max(post);
    est(i+1) = (index-1)*pi/500;
end
est_fun = griddedInterpolant(0:pi/100:2*pi,unwrap(est),'cubic');
end