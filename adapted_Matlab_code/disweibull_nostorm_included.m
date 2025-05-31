function dgrid = disweibull_nostorm_included(shape, scale, support, par, numbins)
    %% Distribution of cyclone windspeed during a short period of duration f
    % CDF of cyclone windspeed (Weibull). Shape and scale already account
    % for no storm years
    %fcdfcond = 1 - exp(-(support / scale).^(1 / shape)); 
    fcdf = 1 - exp(-(support / scale).^shape);
    %unconditional cdf of cyclone windspeed (including zero) i.e. no cyclone occurs
    fpdf = [diff([0,fcdf]),1-sum(diff([0,fcdf]))]; %calculate pdf
    %disp(['Printing length of fpdf: ', num2str(length(fpdf))]);
    
    %% Convoluted PDF of windspeed over full model period (5 years)
    pdf = fpdf;

    num_periods = par.period / par.f; 

    for k = 1:num_periods
        pdf = conv(pdf, fpdf);  
        pdf = pdf / sum(pdf); % ensure pdf is a valid pdf

    end

    cdf = cumsum(pdf);
    
    %% Truncate the distribution (limit to 99.9% cumulative probability)
    maxind = find(cdf > 0.999, 1);
    support = 0:maxind - 1;
    pdf = pdf(1:maxind) / sum(pdf(1:maxind));
    cdf = cumsum(pdf);
    dgrid.pdfwspeed = pdf; dgrid.cdfwspeed = cdf; dgrid.dwspeed = support;
    
    % Discretize cdfcond
    %foodata = interp1(cdfcond, support(2:end), rand(1, 1e6));
    foodata = interp1(cdf, support, rand(1, 1e6));
    [counts, edges] = histcounts(foodata, numbins);
    midedges = (edges(2:end) + edges(1:end-1)) / 2;
    frequency = counts / sum(counts);
    dgrid.pdfcond = frequency;
    %disp(['Printing length of dgrid.pdfcond: ', num2str(length(dgrid.pdfcond))]);
    dgrid.cdfcond = cumsum(dgrid.pdfcond);
    %disp(['Printing length of dgrid.cdfcond: ', num2str(length(dgrid.cdfcond))]);
    
    %% Create Final Unconditional CDF
    dgrid.d = [0, midedges * par.marginaldamage];
    %dgrid.cdf = dgrid.cdfcond; % Directly use cdfcond as the unconditional CDF
    dgrid.cdf = [0, dgrid.cdfcond];
    %disp(['Printing length of dgrid.cdf: ', num2str(length(dgrid.cdf))]);
    dgrid.pdf = diff([0, dgrid.cdf]);
    %disp(['Printing length of dgrid.pdf: ', num2str(length(dgrid.pdf))]);
    %dgrid.pdf = diff(dgrid.cdf);

end