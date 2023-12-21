function res = tracking_estimate_outliers(sim,flag)
switch flag
    case 1
        %"the traditional kalman filter with input processment and measurement noise covariance matrices"
        res = tracking_trueKF(sim);
    case 2
        %2017 student_based 2 KF "A novel robust Student’s t based Kalman filter"
        res = tracking_Stu2KF(sim);
    case 3
        %2017 "Robust Kalman Filters Based on Gaussian Scale Mixture Distributions With Application to Target Tracking"
        res = tracking_GSTMKF(sim);
    case 4
        %2021 "A Novel Robust Kalman Filtering Framework Based on Normal-Skew Mixture Distribution"
        res = tracking_NVGMNGHVGMKF(sim);
    case 5
        %"A novel outlier-robust Kalman filtering framework based on statistical similarity measure"
        res = tracking_SSMKF(sim);
    case 6
        %2022 "Statistical Similarity Measure-Based Adaptive Outlier-Robust State Estimator With Applications"
        res = tracking_AORKF(sim);
    case 7
        %2021 "A Novel Robust Nonlinear Kalman Filter Based on Multivariate Laplace Distribution"
        res = tracking_MLKF(sim);
    case 8
        res = tracking_REAKF2_4(sim);
    case 9
        res = tracking_REAKF2_4_student(sim);
end