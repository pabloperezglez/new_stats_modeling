# Plots manifest -- Arctic shipping GAM experiments
_49 figures in `/plots`, named `<exp_id>_<type>.png`. Generated 2026-05-26._

Legend: M0a/M0b = spatial bake-off (hex random-effect vs tensor te); M1-M5 = nested Tweedie mean build-up; D_* = distribution comparison at the final predictor; S_twlss = Tweedie location-scale (dispersion). Stage 04 (ziplss) did not converge -- no plots.

 1. `M0a_spatial_re_smooths.png` -- M0a baseline, hex random-effect spatial -- smooth-effect panels, one per term (free y-axes, 95% CI)
 2. `M0a_spatial_re_resid.png` -- M0a baseline, hex random-effect spatial -- QQ of deviance residuals + residuals vs linear predictor
 3. `M0a_spatial_re_residmap.png` -- M0a baseline, hex random-effect spatial -- map of per-hex mean deviance residual (spatial autocorrelation)
 4. `M0a_spatial_re_acf.png` -- M0a baseline, hex random-effect spatial -- ACF of monthly-aggregated residuals (temporal autocorrelation)
 5. `M0b_spatial_te_smooths.png` -- M0b baseline, te() tensor spatial (frozen winner) -- smooth-effect panels, one per term (free y-axes, 95% CI)
 6. `M0b_spatial_te_resid.png` -- M0b baseline, te() tensor spatial (frozen winner) -- QQ of deviance residuals + residuals vs linear predictor
 7. `M0b_spatial_te_residmap.png` -- M0b baseline, te() tensor spatial (frozen winner) -- map of per-hex mean deviance residual (spatial autocorrelation)
 8. `M0b_spatial_te_acf.png` -- M0b baseline, te() tensor spatial (frozen winner) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
 9. `M1_ice_smooths.png` -- M1 + s(ice_conc) -- smooth-effect panels, one per term (free y-axes, 95% CI)
10. `M1_ice_resid.png` -- M1 + s(ice_conc) -- QQ of deviance residuals + residuals vs linear predictor
11. `M1_ice_residmap.png` -- M1 + s(ice_conc) -- map of per-hex mean deviance residual (spatial autocorrelation)
12. `M1_ice_acf.png` -- M1 + s(ice_conc) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
13. `M2a_air_temp_smooths.png` -- M2a + s(air_temp) raw -- smooth-effect panels, one per term (free y-axes, 95% CI)
14. `M2a_air_temp_resid.png` -- M2a + s(air_temp) raw -- QQ of deviance residuals + residuals vs linear predictor
15. `M2a_air_temp_residmap.png` -- M2a + s(air_temp) raw -- map of per-hex mean deviance residual (spatial autocorrelation)
16. `M2a_air_temp_acf.png` -- M2a + s(air_temp) raw -- ACF of monthly-aggregated residuals (temporal autocorrelation)
17. `M2b_air_temp_resid_smooths.png` -- M2b + s(air_temp_resid) ice-orthogonalised (carried fwd) -- smooth-effect panels, one per term (free y-axes, 95% CI)
18. `M2b_air_temp_resid_resid.png` -- M2b + s(air_temp_resid) ice-orthogonalised (carried fwd) -- QQ of deviance residuals + residuals vs linear predictor
19. `M2b_air_temp_resid_residmap.png` -- M2b + s(air_temp_resid) ice-orthogonalised (carried fwd) -- map of per-hex mean deviance residual (spatial autocorrelation)
20. `M2b_air_temp_resid_acf.png` -- M2b + s(air_temp_resid) ice-orthogonalised (carried fwd) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
21. `M3_wind_smooths.png` -- M3 + s(wind_speed) -- smooth-effect panels, one per term (free y-axes, 95% CI)
22. `M3_wind_resid.png` -- M3 + s(wind_speed) -- QQ of deviance residuals + residuals vs linear predictor
23. `M3_wind_residmap.png` -- M3 + s(wind_speed) -- map of per-hex mean deviance residual (spatial autocorrelation)
24. `M3_wind_acf.png` -- M3 + s(wind_speed) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
25. `M4_port_smooths.png` -- M4 + s(port_dist) -- smooth-effect panels, one per term (free y-axes, 95% CI)
26. `M4_port_resid.png` -- M4 + s(port_dist) -- QQ of deviance residuals + residuals vs linear predictor
27. `M4_port_residmap.png` -- M4 + s(port_dist) -- map of per-hex mean deviance residual (spatial autocorrelation)
28. `M4_port_acf.png` -- M4 + s(port_dist) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
29. `M5_post2017_smooths.png` -- M5 + post2017 (FINAL mean predictor) -- smooth-effect panels, one per term (free y-axes, 95% CI)
30. `M5_post2017_resid.png` -- M5 + post2017 (FINAL mean predictor) -- QQ of deviance residuals + residuals vs linear predictor
31. `M5_post2017_residmap.png` -- M5 + post2017 (FINAL mean predictor) -- map of per-hex mean deviance residual (spatial autocorrelation)
32. `M5_post2017_acf.png` -- M5 + post2017 (FINAL mean predictor) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
33. `D_tw_smooths.png` -- Distribution: Tweedie -- smooth-effect panels, one per term (free y-axes, 95% CI)
34. `D_tw_resid.png` -- Distribution: Tweedie -- QQ of deviance residuals + residuals vs linear predictor
35. `D_tw_residmap.png` -- Distribution: Tweedie -- map of per-hex mean deviance residual (spatial autocorrelation)
36. `D_tw_acf.png` -- Distribution: Tweedie -- ACF of monthly-aggregated residuals (temporal autocorrelation)
37. `D_nb_smooths.png` -- Distribution: negative binomial -- smooth-effect panels, one per term (free y-axes, 95% CI)
38. `D_nb_resid.png` -- Distribution: negative binomial -- QQ of deviance residuals + residuals vs linear predictor
39. `D_nb_residmap.png` -- Distribution: negative binomial -- map of per-hex mean deviance residual (spatial autocorrelation)
40. `D_nb_acf.png` -- Distribution: negative binomial -- ACF of monthly-aggregated residuals (temporal autocorrelation)
41. `D_gaussian_smooths.png` -- Distribution: gaussian (naive ref) -- smooth-effect panels, one per term (free y-axes, 95% CI)
42. `D_gaussian_resid.png` -- Distribution: gaussian (naive ref) -- QQ of deviance residuals + residuals vs linear predictor
43. `D_gaussian_residmap.png` -- Distribution: gaussian (naive ref) -- map of per-hex mean deviance residual (spatial autocorrelation)
44. `D_gaussian_acf.png` -- Distribution: gaussian (naive ref) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
45. `S_twlss_scale_ice_smooths.png` -- Tweedie location-scale (dispersion ~ ice) -- smooth-effect panels, one per term (free y-axes, 95% CI)
46. `S_twlss_scale_ice_resid.png` -- Tweedie location-scale (dispersion ~ ice) -- QQ of deviance residuals + residuals vs linear predictor
47. `S_twlss_scale_ice_residmap.png` -- Tweedie location-scale (dispersion ~ ice) -- map of per-hex mean deviance residual (spatial autocorrelation)
48. `S_twlss_scale_ice_acf.png` -- Tweedie location-scale (dispersion ~ ice) -- ACF of monthly-aggregated residuals (temporal autocorrelation)
49. `S_twlss_scale_ice_scale_vs_ice.png` -- Tweedie location-scale (dispersion ~ ice) -- estimated Tweedie dispersion (scale) vs sea-ice concentration

