# Latent Gaussian models for multivariate and high-dimensional discrete-valued time series

Latent Gaussian models are designed for modeling count time series in multivariate, possibly high-dimensional settings. The key idea is to link the second-order properties of the discrete-valued observations to those of the latent Gaussian process through so-called link functions. This repository primarily reproduces the results in Kim, Düker, Fisher, and Pipiras ([2025](#ref-latent_dfm)), where the Gaussian dynamic factor model is assumed as the latent process. It also includes code for related models, such as those in Kim ([2023](#ref-latent_thesis)) and Düker, Lund, Pipiras ([2024](#ref-latent_var)), where the latent Gaussian processes are assumed to be vector autoregressive models. We are developing a fully functioning package, which will be released soon.

## References

<div id="ref-latent_var" class="references" style="margin-bottom: 1em;">
Marie-Christine Düker, Robert Lund, and Vladas Pipiras. 2024. "High-dimensional latent Gaussian count time series: Concentration results for autocovariances and applications." *Electronic Journal of Statistics* 18(2) pp. 5484--5562 
<a href="https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-18/issue-2/High-dimensional-latent-Gaussian-count-time-series--Concentration-results/10.1214/24-EJS2292.full" target="_blank">
    https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-18/issue-2/High-dimensional-latent-Gaussian-count-time-series--Concentration-results/10.1214/24-EJS2292.full
  </a>
</div>

<div id="ref-latent_dfm" class="references" style="margin-bottom: 1em;">
Younghoon Kim, Marie-Christine Düker, Zachary F. Fisher, and Vladas Pipiras. 2025. "Latent Gaussian dynamic factor modeling and forecasting for multivariate count time series." *Journal of Time Series Analysis* pp. 1--16 
 <a href="https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.70016" target="_blank">
    https://onlinelibrary.wiley.com/doi/abs/10.1111/jtsa.70016
  </a>
</div>

<div id="ref-latent_thesis" class="references" style="margin-bottom: 1em;">
Younghoon Kim. 2023. <em>Modeling Multiple-Subject and Discrete-Valued High-Dimensional Time Series.</em> PhD thesis, The University of North Carolina at Chapel Hill. <a href="https://www.proquest.com/openview/b7d1eae2131e518bc3af4ca4f2816513/1?pq-origsite=gscholar&cbl=18750&diss=y">https://www.proquest.com/openview/b7d1eae2131e518bc3af4ca4f2816513/1?pq-origsite=gscholar&cbl=18750&diss=y</a>
</div>
