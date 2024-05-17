function res = tracking_estimate(sim,flag)
switch flag
    case 1
        %"the traditional kalman filter with input processment and measurement noise covariance matrices"
        res = tracking_trueKF(sim);
    case 2
        %2009, "Embedded cubature Kalman filter with adaptive setting of  free parameter"
        res = tracking_CIRKF(sim);
    case 3
        %2017 "A Novel Adaptive Kalman Filter With Inaccurate Process and Measurement Noise Covariance Matrices"
        res = tracking_VBAKF(sim);
    case 4
        %2019 "A novel adaptive Kalman filter with unknown probability of measurement loss"
        res = tracking_UPAKF(sim);
    case 5
        %2019 "Sequential Bayesian estimation of state and input in dynamical systems using output-only measurements"
        res = tracking_SBEKF(sim);
    case 6
        %2022 "A Bayesian Expectation-Maximization (BEM) methodology for joint input-state estimation and virtual sensing of structures"
        res = tracking_EMKF(sim);
    case 7
        %2020 "Force localization and reconstruction based on a novel sparse Kalman filter"
        res = tracking_SparseKF(sim);
    case 8
        %2021 "Multiple fading factors-based strong tracking variational Bayesian adaptive Kalman filter"
        res = tracking_MSTVBAKF(sim);
    case 9
        res = tracking_EAKF2_4(sim);
end
end