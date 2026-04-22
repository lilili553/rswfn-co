pub fn rks_c_lyp(ideriv: i32, npt: usize, rhoa1: &[f64], sigmaaa1: &[f64], zk: &mut [f64], vrhoa: &mut [f64], vsigmaaa: &mut [f64], v2rhoa2: &mut [f64], v2rhoasigmaaa: &mut [f64], v2sigmaaa2: &mut [f64]) {
    const TOL: f64 = 1.0e-20;
    if ideriv == 0 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = 1.0 / t2;
                let t6 = 1.0 / (1.0 + 0.349 * t3);
                let t9 = 0.2533 * t3;
                let t10 = (-t9).exp();
                let t12 = rho * rho;
                let t14 = t2 * t2;
                let t20 = t3 * t6;
                zk[i] = -0.04918 * t6 * rho - 0.00649176 * t10 * t6 / t14 / t12 / rho * (0.25 * t12 * (0.1148493600075277e2 * t14 * t12 + (0.2611111111111111e1 - 0.09850555555555556 * t3 - 0.1357222222222222 * t20) * sigma - 0.5 * (0.25e1 - 0.014072222222222222 * t3 - 0.01938888888888889 * t20) * sigma - 0.027777777777777776 * (t9 + 0.349 * t20 - 11.0) * sigma) - 0.4583333333333333 * t12 * sigma);
            } else {
                zk[i] = 0.0;
            }
        }
    } else if ideriv == 1 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = 1.0 / t2;
                let t5 = 1.0 + 0.349 * t3;
                let t6 = 1.0 / t5;
                let t9 = 0.2533 * t3;
                let t10 = (-t9).exp();
                let t11 = t10 * t6;
                let t12 = rho * rho;
                let t14 = t2 * t2;
                let t16 = 1.0 / t14 / t12 / rho;
                let t20 = t3 * t6;
                let t22 = 0.2611111111111111e1 - 0.09850555555555556 * t3 - 0.1357222222222222 * t20;
                let t30 = t9 + 0.349 * t20 - 11.0;
                let t33 = 0.1148493600075277e2 * t14 * t12 + t22 * sigma - 0.5 * (0.25e1 - 0.014072222222222222 * t3 - 0.01938888888888889 * t20) * sigma - 0.027777777777777776 * t30 * sigma;
                let t38 = 0.25 * t12 * t33 - 0.4583333333333333 * t12 * sigma;
                zk[i] = -0.04918 * t6 * rho - 0.00649176 * t11 * t16 * t38;
                let t45 = t14 * rho;
                let t49 = t30 / rho * sigma;
                let t54 = rho * sigma;
                let t60 = t5 * t5;
                let t61 = 1.0 / t60;
                let t64 = t12 * t12;
                let t66 = 1.0 / t64 / rho;
                let t81 = 1.0 / t2 / rho;
                let t83 = t81 * t6;
                let t85 = 1.0 / t45;
                let t86 = t85 * t61;
                vrhoa[i] = -0.04918 * t6 - 0.00649176 * t11 * t16 * (0.5 * rho * t33 + 0.25 * t12 * (0.3062649600200738e2 * t45 - 0.027777777777777776 * t49) - 0.25 * t54) - 0.005721273333333333 * t61 * t3 - 0.000548120936 * t66 * t10 * t6 * t38 - 0.00075520808 * t10 * t61 * t66 * t38 + 0.02380312 * t11 / t14 / t64 * t38
                    - 0.00649176 * t11 * t16 * (0.25 * t12 * ((0.03283518518518519 * t81 + 0.04524074074074074 * t83 - 0.015789018518518518 * t86) * sigma - 0.5 * (0.004690740740740741 * t81 + 0.006462962962962963 * t83 - 0.002255574074074074 * t86) * sigma - 0.027777777777777776 * (-0.08443333333333333 * t81 - 0.11633333333333333 * t83 + 0.04060033333333333 * t86) * sigma + 0.027777777777777776 * t49) - 0.6666666666666666 * t54);
                vsigmaaa[i] = 0.0007213066666666667 * t11 * t85 - 0.02596704 * t11 * t16 * (0.25 * t12 * t22 - 0.6666666666666666 * t12);
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vsigmaaa[i] = 0.0;
            }
        }
    } else if ideriv == 2 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = 1.0 / t2;
                let t5 = 1.0 + 0.349 * t3;
                let t6 = 1.0 / t5;
                let t9 = 0.2533 * t3;
                let t10 = (-t9).exp();
                let t11 = t10 * t6;
                let t12 = rho * rho;
                let t13 = t12 * rho;
                let t14 = t2 * t2;
                let t16 = 1.0 / t14 / t13;
                let t17 = t14 * t12;
                let t18 = 0.1148493600075277e2 * t17;
                let t20 = t3 * t6;
                let t22 = 0.2611111111111111e1 - 0.09850555555555556 * t3 - 0.1357222222222222 * t20;
                let t23 = t22 * sigma;
                let t28 = 0.5 * (0.25e1 - 0.014072222222222222 * t3 - 0.01938888888888889 * t20) * sigma;
                let t30 = t9 + 0.349 * t20 - 11.0;
                let t32 = 0.027777777777777776 * t30 * sigma;
                let t33 = t18 + t23 - t28 - t32;
                let t38 = 0.25 * t12 * t33 - 0.4583333333333333 * t12 * sigma;
                zk[i] = -0.04918 * t6 * rho - 0.00649176 * t11 * t16 * t38;
                let t45 = t14 * rho;
                let t47 = 1.0 / rho;
                let t48 = t30 * t47;
                let t49 = t48 * sigma;
                let t51 = 0.3062649600200738e2 * t45 - 0.027777777777777776 * t49;
                let t54 = rho * sigma;
                let t56 = 0.5 * rho * t33 + 0.25 * t12 * t51 - 0.25 * t54;
                let t60 = t5 * t5;
                let t61 = 1.0 / t60;
                let t64 = t12 * t12;
                let t65 = t64 * rho;
                let t66 = 1.0 / t65;
                let t67 = t66 * t10;
                let t68 = t6 * t38;
                let t71 = t10 * t61;
                let t76 = 1.0 / t14 / t64;
                let t81 = 1.0 / t2 / rho;
                let t83 = t81 * t6;
                let t85 = 1.0 / t45;
                let t86 = t85 * t61;
                let t88 = 0.03283518518518519 * t81 + 0.04524074074074074 * t83 - 0.015789018518518518 * t86;
                let t99 = -0.08443333333333333 * t81 - 0.11633333333333333 * t83 + 0.04060033333333333 * t86;
                let t103 = t88 * sigma - 0.5 * (0.004690740740740741 * t81 + 0.006462962962962963 * t83 - 0.002255574074074074 * t86) * sigma - 0.027777777777777776 * t99 * sigma + 0.027777777777777776 * t49;
                let t107 = 0.25 * t12 * t103 - 0.6666666666666666 * t54;
                vrhoa[i] = -0.04918 * t6 - 0.00649176 * t11 * t16 * t56 - 0.005721273333333333 * t61 * t3 - 0.000548120936 * t67 * t68 - 0.00075520808 * t71 * t66 * t38 + 0.02380312 * t11 * t76 * t38 - 0.00649176 * t11 * t16 * t107;
                let t116 = 0.25 * t12 * t22 - 0.6666666666666666 * t12;
                vsigmaaa[i] = 0.0007213066666666667 * t11 * t85 - 0.02596704 * t11 * t16 * t116;
                let t121 = 1.0 / t60 / t5;
                let t126 = t64 * t12;
                let t127 = 1.0 / t126;
                let t137 = t99 * t47 * sigma;
                let t141 = t30 / t12 * sigma;
                let t150 = rho * t51;
                let t165 = 1.0 / t2 / t126;
                let t166 = t165 * t10;
                let t188 = 1.0 / t2 / t12;
                let t190 = t188 * t6;
                let t192 = 1.0 / t17;
                let t193 = t192 * t61;
                let t195 = 1.0 / t13;
                let t196 = t195 * t121;
                let s1 = -0.002662299191111111 * t121 * t85 - 0.007628364444444444 * t61 * t81 + 0.009500762890666667 * t127 * t10 * t68 + 0.09521248 * t11 * t76 * t56 - 0.02596704 * t11 * t16 * (0.5 * rho * t103 + 0.25 * t12 * (-0.027777777777777776 * t137 + 0.027777777777777776 * t141)) - 0.00649176 * t11 * t16 * (t18 + t23 - t28 - t32 + t150) + 0.09521248 * t11 * t76 * t107 - 0.002192483744 * t67 * t6 * t56 - 0.00302083232 * t71 * t66 * t56;
                let s2 = s1 - 0.0002550589422186667 * t166 * t61 * t38 - 0.00302083232 * t71 * t66 * t107 - 0.0003514234932266667 * t10 * t121 * t165 * t38 - 0.002192483744 * t67 * t6 * t107;
                v2rhoa2[i] = s2
                    - 0.00009255935539253333 * t166 * t68
                    - 0.2221624533333333 * t11 / t14 / t65 * t38
                    - 0.01298352
                        * t11
                        * t16
                        * (0.25
                            * t12
                            * ((-0.04378024691358025 * t188 - 0.06032098765432099 * t190 + 0.03157803703703704 * t193 - 0.003673578308641975 * t196) * sigma - 0.5 * (-0.006254320987654321 * t188 - 0.008617283950617284 * t190 + 0.004511148148148148 * t193 - 0.0005247969012345679 * t196) * sigma - 0.027777777777777776 * (0.11257777777777778 * t188 + 0.15511111111111112 * t190 - 0.08120066666666667 * t193 + 0.009446344222222222 * t196) * sigma + 0.05555555555555555 * t137
                                - 0.05555555555555555 * t141)
                            - 0.6666666666666666 * sigma)
                    - 0.00649176 * t11 * t16 * (1.0 * t150 + 0.2552208000167282e2 * t17 - 0.5 * sigma)
                    + 0.013090273386666667 * t71 * t127 * t38;
                v2rhoasigmaaa[i] = -0.00649176 * t11 * t16 * (-0.9444444444444444 * rho - 0.027777777777777776 * rho * t30) + 0.00006090232622222222 * t195 * t10 * t6 + 0.00008391200888888889 * t71 * t195 + 0.009978075555555556 * t11 * t192 - 0.01298352 * t11 * t16 * (0.25 * t12 * (-3.0e-23 * t83 + 1.0e-23 * t86 + 0.05555555555555555 * t48) + 1.3333333333333333 * rho) - 0.01298352 * t11 * t192 * t22 - 0.002192483744 * t67 * t6 * t116 - 0.00302083232 * t71 * t66 * t116
                    + 0.09521248 * t11 * t76 * t116
                    - 0.02596704 * t11 * t16 * (0.25 * t12 * t88 - 1.3333333333333333 * rho);
                v2sigmaaa2[i] = 0.0;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                v2rhoa2[i] = 0.0;
                vsigmaaa[i] = 0.0;
                v2rhoasigmaaa[i] = 0.0;
                v2sigmaaa2[i] = 0.0;
            }
        }
    }
}

pub fn rks_c_p86(ideriv: i32, npt: usize, rhoa1: &[f64], sigmaaa1: &[f64], zk: &mut [f64], vrhoa: &mut [f64], vsigmaaa: &mut [f64], v2rhoa2: &mut [f64], v2rhoasigmaaa: &mut [f64], v2sigmaaa2: &mut [f64]) {
    const TOL: f64 = 1.0e-20;
    if ideriv == 0 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = 1.0 / rho;
                let t3 = t2.powf(1.0 / 3.0);
                let t4 = 0.6203504908994 * t3;
                let t6 = t2.powf(1.0 / 6.0);
                let t12 = t4.ln();
                let t18 = if 1.0 <= t4 { -0.1423 / (1.0 + 0.8292885914166397 * t6 + 0.20682485366586 * t3) } else { 0.0311 * t12 - 0.048 + 0.0012407009817988 * t3 * t12 - 0.00719606569443304 * t3 };
                let t20 = sigma.sqrt();
                let t22 = t3 * t3;
                let t24 = 0.002568 + 0.01443307452126544 * t3 + 0.000002843543831490386 * t22;
                let t28 = 1.0 + 5.411317332115466 * t3 + 0.1816419932959077 * t22 + 0.01763993811759022 * t2;
                let t31 = 0.001667 + t24 / t28;
                let t34 = rho.powf(1.0 / 6.0);
                let t39 = (-0.00081290825 * t20 / t31 / t34 / rho).exp();
                let t41 = rho.powf(1.0 / 3.0);
                zk[i] = rho * t18 + 1.0 * t39 * t31 * sigma / t41 / rho;
            } else {
                zk[i] = 0.0;
            }
        }
    } else if ideriv == 1 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = 1.0 / rho;
                let t3 = t2.powf(1.0 / 3.0);
                let t4 = 0.6203504908994 * t3;
                let t5 = 1.0 <= t4;
                let t6 = t2.powf(1.0 / 6.0);
                let t9 = 1.0 + 0.8292885914166397 * t6 + 0.20682485366586 * t3;
                let t12 = t4.ln();
                let t18 = if t5 { -0.1423 / t9 } else { 0.0311 * t12 - 0.048 + 0.0012407009817988 * t3 * t12 - 0.00719606569443304 * t3 };
                let t20 = sigma.sqrt();
                let t22 = t3 * t3;
                let t24 = 0.002568 + 0.01443307452126544 * t3 + 0.000002843543831490386 * t22;
                let t28 = 1.0 + 5.411317332115466 * t3 + 0.1816419932959077 * t22 + 0.01763993811759022 * t2;
                let t29 = 1.0 / t28;
                let t31 = 0.001667 + t24 * t29;
                let t33 = t20 / t31;
                let t34 = rho.powf(1.0 / 6.0);
                let t36 = 1.0 / t34 / rho;
                let t39 = (-0.00081290825 * t33 * t36).exp();
                let t40 = t39 * t31;
                let t41 = rho.powf(1.0 / 3.0);
                let t43 = 1.0 / t41 / rho;
                let t44 = sigma * t43;
                zk[i] = rho * t18 + 1.0 * t40 * t44;
                let t47 = if t5 { 0.0 } else { 0.0 };
                let t49 = t9 * t9;
                let t51 = t6 * t6;
                let t52 = t51 * t51;
                let t55 = rho * rho;
                let t56 = 1.0 / t55;
                let t59 = 1.0 / t22;
                let t60 = t59 * t56;
                let t62 = -0.1382147652361066 / (t52 * t6) * t56 - 0.06894161788661999 * t60;
                let t66 = t59 * t12;
                let t69 = t3 * t2;
                let t73 = if t5 { 0.1423 / t49 * t62 } else { -0.010366666666666667 * t2 - 0.0004135669939329333 * t66 * t56 - 0.0004135669939329333 * t69 + 0.002398688564811013 * t60 };
                let t75 = t31 * t31;
                let t76 = 1.0 / t75;
                let t77 = t20 * t76;
                let t79 = 1.0 / t3;
                let t80 = t79 * t56;
                let t82 = -0.004811024840421814 * t60 - 0.000001895695887660258 * t80;
                let t84 = t28 * t28;
                let t85 = 1.0 / t84;
                let t86 = t24 * t85;
                let t90 = -0.1803772444038489 * t60 - 0.1210946621972718 * t80 - 0.01763993811759022 * t56;
                let t93 = t82 * t29 - 1.0 * t86 * t90;
                let t94 = t36 * t93;
                let t98 = 1.0 / t34 / t55;
                let t101 = 0.00081290825 * t77 * t94 + 0.0009483929583333333 * t33 * t98;
                let t102 = t101 * t39;
                let t103 = t31 * sigma;
                let t104 = t103 * t43;
                let t107 = t39 * t93;
                let t111 = 1.0 / t41 / t55;
                let t112 = sigma * t111;
                vrhoa[i] = rho * t47 + t18 + rho * t73 + 1.0 * t102 * t104 + 1.0 * t107 * t44 - 1.3333333333333333 * t40 * t112;
                let t115 = rho.sqrt();
                let t117 = 1.0 / t115 / t55;
                vsigmaaa[i] = -0.0016258165 * t20 * t117 * t39 + 4.0 * t40 * t43;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vsigmaaa[i] = 0.0;
            }
        }
    } else if ideriv == 2 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa1[i].max(0.0);
                let t2 = 1.0 / rho;
                let t3 = t2.powf(1.0 / 3.0);
                let t4 = 0.6203504908994 * t3;
                let t5 = 1.0 <= t4;
                let t6 = t2.powf(1.0 / 6.0);
                let t9 = 1.0 + 0.8292885914166397 * t6 + 0.20682485366586 * t3;
                let t11 = 0.1423 / t9;
                let t12 = t4.ln();
                let t14 = t3 * t12;
                let t18 = if t5 { -t11 } else { 0.0311 * t12 - 0.048 + 0.0012407009817988 * t14 - 0.00719606569443304 * t3 };
                let t20 = sigma.sqrt();
                let t22 = t3 * t3;
                let t24 = 0.002568 + 0.01443307452126544 * t3 + 0.000002843543831490386 * t22;
                let t28 = 1.0 + 5.411317332115466 * t3 + 0.1816419932959077 * t22 + 0.01763993811759022 * t2;
                let t29 = 1.0 / t28;
                let t31 = 0.001667 + t24 * t29;
                let t32 = 1.0 / t31;
                let t33 = t20 * t32;
                let t34 = rho.powf(1.0 / 6.0);
                let t36 = 1.0 / t34 / rho;
                let t39 = (-0.00081290825 * t33 * t36).exp();
                let t40 = t39 * t31;
                let t41 = rho.powf(1.0 / 3.0);
                let t43 = 1.0 / t41 / rho;
                let t44 = sigma * t43;
                zk[i] = rho * t18 + 1.0 * t40 * t44;
                let t47 = if t5 { 0.0 } else { 0.0 };
                let t48 = rho * t47;
                let t49 = t9 * t9;
                let t50 = 1.0 / t49;
                let t51 = t6 * t6;
                let t52 = t51 * t51;
                let t53 = t52 * t6;
                let t54 = 1.0 / t53;
                let t55 = rho * rho;
                let t56 = 1.0 / t55;
                let t59 = 1.0 / t22;
                let t60 = t59 * t56;
                let t62 = -0.1382147652361066 * t54 * t56 - 0.06894161788661999 * t60;
                let t66 = t59 * t12;
                let t69 = t3 * t2;
                let t73 = if t5 { 0.1423 * t50 * t62 } else { -0.010366666666666667 * t2 - 0.0004135669939329333 * t66 * t56 - 0.0004135669939329333 * t69 + 0.002398688564811013 * t60 };
                let t75 = t31 * t31;
                let t76 = 1.0 / t75;
                let t77 = t20 * t76;
                let t79 = 1.0 / t3;
                let t80 = t79 * t56;
                let t82 = -0.004811024840421814 * t60 - 0.000001895695887660258 * t80;
                let t84 = t28 * t28;
                let t85 = 1.0 / t84;
                let t86 = t24 * t85;
                let t90 = -0.1803772444038489 * t60 - 0.1210946621972718 * t80 - 0.01763993811759022 * t56;
                let t93 = t82 * t29 - 1.0 * t86 * t90;
                let t94 = t36 * t93;
                let t98 = 1.0 / t34 / t55;
                let t101 = 0.00081290825 * t77 * t94 + 0.0009483929583333333 * t33 * t98;
                let t102 = t101 * t39;
                let t103 = t31 * sigma;
                let t104 = t103 * t43;
                let t107 = t39 * t93;
                let t111 = 1.0 / t41 / t55;
                let t112 = sigma * t111;
                vrhoa[i] = t48 + t18 + rho * t73 + 1.0 * t102 * t104 + 1.0 * t107 * t44 - 1.3333333333333333 * t40 * t112;
                let t115 = rho.sqrt();
                let t117 = 1.0 / t115 / t55;
                vsigmaaa[i] = -0.0016258165 * t20 * t117 * t39 + 4.0 * t40 * t43;
                let t129 = (-0.0843 / (1.0 + 0.1101176160755631 * t6 + 0.1619735131738333 * t3) + t11) * t56;
                let t135 = (-0.01555 * t12 + 0.0211 - 0.00080645563816922 * t14 + 0.00421838333811592 * t3) * t56;
                let t137 = if t5 { 1.709920934161366 * t129 } else { 1.709920934161366 * t135 };
                let t139 = t55 * rho;
                let t147 = t62 * t62;
                let t152 = t55 * t55;
                let t153 = 1.0 / t152;
                let t156 = 1.0 / t139;
                let t160 = 1.0 / t22 / t2;
                let t161 = t160 * t153;
                let t163 = t59 * t156;
                let t180 = if t5 {
                    -0.2846 / t49 / t9 * t147 + 0.1423 * t50 * (-0.1151789710300888 / t53 / t2 * t153 + 0.2764295304722132 * t54 * t156 - 0.04596107859241333 * t161 + 0.13788323577724 * t163)
                } else {
                    0.010366666666666667 * t56 - 0.0002757113292886222 * t160 * t12 * t153 - 0.004521665800333405 * t163 + 0.0008271339878658667 * t66 * t156 + 0.0004135669939329333 * t3 * t56 + 0.001599125709874009 * t161
                };
                let t185 = if t5 { -1.709920934161366 * t129 } else { -1.709920934161366 * t135 };
                let t193 = t93 * t93;
                let t203 = 1.0 / t69 * t153;
                let t205 = t79 * t156;
                let t215 = t90 * t90;
                let t226 = (-0.003207349893614542 * t161 + 0.009622049680843627 * t163 - 0.0000006318986292200858 * t203 + 0.000003791391775320515 * t205) * t29 - 2.0 * t82 * t85 * t90 + 2.0 * t24 / t84 / t28 * t215 - 1.0 * t86 * (-1.202514962692326 * t161 + 3.607544888076978 * t163 - 0.04036488739909061 * t203 + 0.2421893243945437 * t205 + 0.03527987623518044 * t156);
                let t247 = t101 * t101;
                v2rhoa2[i] = rho * t137 + 6.222222222222222 * t40 * sigma / t41 / t139 + 2.0 * rho * t180 + rho * t185 - 5.333333333333333 * t102 * t103 * t111 + 2.0 * (-0.0016258165 * t20 / t75 / t31 * t36 * t193 - 0.001896785916666667 * t77 * t98 * t93 + 0.00081290825 * t77 * t36 * t226 - 0.002054851409722222 * t33 / t34 / t139) * t39 * t104 + 4.0 * t102 * t93 * sigma * t43 + 2.0 * t39 * t226 * t44 - 5.333333333333333 * t107 * t112 + 2.0 * t247 * t39 * t104 + 4.0 * t48 + 4.0 * t47 + 4.0 * t73;
                let t254 = 1.0 / t20;
                let t266 = t117 * t39;
                v2rhoasigmaaa[i] = 4.0 * (0.000406454125 * t254 * t76 * t94 + 0.0004741964791666667 * t254 * t32 * t98) * t39 * t104 - 0.0016258165 * t101 * t20 * t266 + 4.0 * t102 * t31 * t43 - 0.0016258165 * t33 * t266 * t93 + 4.0 * t107 * t43 + 0.002167755333333333 * t20 / t115 / t139 * t39 - 5.333333333333333 * t40 * t111;
                let t287 = t41 * t41;
                v2sigmaaa2[i] = -0.009754899 * t254 * t117 * t39 + 0.00000264327929167225 / t287 / t139 * t32 * t39;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                v2rhoa2[i] = 0.0;
                vsigmaaa[i] = 0.0;
                v2rhoasigmaaa[i] = 0.0;
                v2sigmaaa2[i] = 0.0;
            }
        }
    }
}

/// Becke 88 exchange functional
///
/// # Arguments
/// * `ideriv` - 导数阶数 (0, 1, 2)
/// * `npt`    - 网格点数
/// * `rhoa`   - 电荷密度 (长度 >= npt)
/// * `sigmaaa`- 密度梯度平方 (长度 >= npt)
/// * `zk`     - 泛函值输出 (长度 >= npt)
/// * `vrhoa`  - 一阶导数 wrt rho (长度 >= npt)
/// * `vsigmaaa`- 一阶导数 wrt sigma (长度 >= npt)
/// * `v2rhoa2`- 二阶导数 wrt rho,rho (长度 >= npt)
/// * `v2rhoasigmaaa`- 二阶导数 wrt rho,sigma (长度 >= npt)
/// * `v2sigmaaa2`- 二阶导数 wrt sigma,sigma (长度 >= npt)
pub fn rks_x_b88(ideriv: i32, npt: usize, rhoa: &[f64], sigmaaa: &[f64], zk: &mut [f64], vrhoa: &mut [f64], vsigmaaa: &mut [f64], v2rhoa2: &mut [f64], v2rhoasigmaaa: &mut [f64], v2sigmaaa2: &mut [f64]) {
    const TOL: f64 = 1.0e-20;
    if ideriv == 0 {
        for i in 0..npt {
            let rho = rhoa[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = t2 * rho;
                let t5 = 1.0 / t3;
                let t7 = sigma.sqrt();
                let t8 = t7 * t5;
                let t10 = (1.259921049894873 * t8 + (1.0 + 1.587401051968199 * t8 * t8).sqrt()).ln();
                let t13 = 1.0 + 0.0317500104573508 * t8 * t10;
                let t14 = 1.0 / t13;
                zk[i] = -0.7385587663820224 * t3 - 0.005291668409558467 * t5 * sigma * t14;
            } else {
                zk[i] = 0.0;
            }
        }
    } else if ideriv == 1 {
        for i in 0..npt {
            let rho = rhoa[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = t2 * rho;
                let t5 = 1.0 / t3;
                let t6 = t5 * sigma;
                let t7 = sigma.sqrt();
                let t8 = t7 * t5;
                let t10 = (1.259921049894873 * t8 + (1.0 + 1.587401051968199 * t8 * t8).sqrt()).ln();
                let t13 = 1.0 + 0.0317500104573508 * t8 * t10;
                let t14 = 1.0 / t13;
                zk[i] = -0.7385587663820224 * t3 - 0.005291668409558467 * t6 * t14;

                let t18 = rho * rho;
                let t20 = 1.0 / t2 / t18;
                let t24 = t13 * t13;
                let t25 = 1.0 / t24;
                let t30 = t2 * t2;
                let t35 = 1.0 / t30 / t18;
                let t38 = 1.0 + 1.587401051968199 * sigma * t35;
                let t39 = t38.sqrt();
                let t40 = 1.0 / t39;

                let t43 = -0.08466669455293548 * t7 * t20 * t10 - 0.106673350692263 * sigma / t30 / (rho * rho * rho) * t40;
                let t44 = t25 * t43;

                vrhoa[i] = -0.9847450218426965 * t2 + 0.00705555787941129 * t20 * sigma * t14 + 0.002645834204779234 * t6 * t44;

                let t49 = 1.0 / t7;
                let t55 = 0.06350002091470161 * t49 * t5 * t10 + 0.08000501301919725 * t35 * t40;
                let t56 = t25 * t55;

                vsigmaaa[i] = -0.02116667363823387 * t5 * t14 + 0.005291668409558467 * t6 * t56;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vsigmaaa[i] = 0.0;
            }
        }
    } else if ideriv == 2 {
        for i in 0..npt {
            let rho = rhoa[i].max(0.0);
            if rho > TOL {
                let sigma = sigmaaa[i].max(0.0);
                let t2 = rho.powf(1.0 / 3.0);
                let t3 = t2 * rho;
                let t5 = 1.0 / t3;
                let t6 = t5 * sigma;
                let t7 = sigma.sqrt();
                let t8 = t7 * t5;
                let t10 = (1.259921049894873 * t8 + (1.0 + 1.587401051968199 * t8 * t8).sqrt()).ln();
                let t13 = 1.0 + 0.0317500104573508 * t8 * t10;
                let t14 = 1.0 / t13;
                zk[i] = -0.7385587663820224 * t3 - 0.005291668409558467 * t6 * t14;

                let t18 = rho * rho;
                let t20 = 1.0 / t2 / t18;
                let t21 = t20 * sigma;
                let t24 = t13 * t13;
                let t25 = 1.0 / t24;
                let t29 = t18 * rho;
                let t30 = t2 * t2;
                let t32 = 1.0 / t30 / t29;
                let t35 = 1.0 / t30 / t18;
                let t38 = 1.0 + 1.587401051968199 * sigma * t35;
                let t39 = t38.sqrt();
                let t40 = 1.0 / t39;

                let t43 = -0.08466669455293548 * t7 * t20 * t10 - 0.106673350692263 * sigma * t32 * t40;
                let t44 = t25 * t43;
                vrhoa[i] = -0.9847450218426965 * t2 + 0.00705555787941129 * t21 * t14 + 0.002645834204779234 * t6 * t44;

                let t49 = 1.0 / t7;
                let t55 = 0.06350002091470161 * t49 * t5 * t10 + 0.08000501301919725 * t35 * t40;
                let t56 = t25 * t55;
                vsigmaaa[i] = -0.02116667363823387 * t5 * t14 + 0.005291668409558467 * t6 * t56;

                let t62 = 1.0 / t2 / t29;
                let t69 = 1.0 / t24 / t13;
                let t70 = t43 * t43;
                let t77 = t18 * t18;
                let t83 = sigma * sigma;
                let t89 = 1.0 / t39 / t38;

                v2rhoa2[i] = -0.6564966812284644 / t30 - 0.03292593677058602 * t62 * sigma * t14 - 0.01411111575882258 * t21 * t44 - 0.005291668409558467 * t6 * t69 * t70 + 0.002645834204779234 * t6 * t25 * (0.3951112412470322 * t7 * t62 * t10 + 1.06673350692263 * sigma / t30 / t77 * t40 - 0.4515557042823225 * t83 / t2 / t77 / t29 * t89);

                v2rhoasigmaaa[i] = 0.02822223151764516 * t20 * t14 - 0.00705555787941129 * t21 * t56 + 0.01058333681911693 * t5 * t25 * t43 - 0.005291668409558467 * t6 * t69 * t43 * t55 + 0.002645834204779234 * t6 * t25 * (-0.169333389105871 * t49 * t20 * t10 - 0.640040104153578 * t32 * t40 + 0.3386667782117419 * sigma / t2 / t77 / t18 * t89);

                let t124 = t55 * t55;
                v2sigmaaa2[i] = 0.04233334727647774 * t5 * t25 * t55 - 0.01058333681911693 * t6 * t69 * t124 + 0.005291668409558467 * t6 * t25 * (-0.1270000418294032 / t7 / sigma * t5 * t10 + 0.1600100260383945 / sigma * t35 * t40 - 0.2540000836588064 / t2 / t77 / rho * t89);
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vsigmaaa[i] = 0.0;
                v2rhoa2[i] = 0.0;
                v2rhoasigmaaa[i] = 0.0;
                v2sigmaaa2[i] = 0.0;
            }
        }
    }
}

pub fn rks_x_lda(
    ideriv: i32,               // 阶数
    npt: usize,                // 格点数量
    rhoa1: &[f64],             // 电子密度
    sigmaaa1: &[f64],          // 未用到
    zk: &mut [f64],            // 泛函的数值
    vrhoa: &mut [f64],         // 泛函对总密度的一阶导数，即df/dρ
    vsigmaaa: &mut [f64],      // 未用到
    v2rhoa2: &mut [f64],       // 泛函对总密度的二阶导数
    v2rhoasigmaaa: &mut [f64], // 未用到
    v2sigmaaa2: &mut [f64],    // 未用到
) {
    const TOL: f64 = 1.0e-20;

    if ideriv == 0 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let t1 = rho.powf(1.0 / 3.0);
                zk[i] = -0.7385587663820224e0 * t1 * rho;
            } else {
                zk[i] = 0.0;
            }
        }
    } else if ideriv == 1 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let t1 = rho.powf(1.0 / 3.0);
                zk[i] = -0.7385587663820224e0 * t1 * rho;
                vrhoa[i] = -0.9847450218426965e0 * t1;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
            }
        }
    } else if ideriv == 2 {
        for i in 0..npt {
            let rho = rhoa1[i].max(0.0);
            if rho > TOL {
                let t1 = rho.powf(1.0 / 3.0);
                zk[i] = -0.7385587663820224e0 * t1 * rho;
                vrhoa[i] = -0.9847450218426965e0 * t1;
                let t5 = t1 * t1;
                v2rhoa2[i] = -0.6564966812284644e0 / t5;
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                v2rhoa2[i] = 0.0;
            }
        }
    }
}

pub fn uks_x_lda(
    ideriv: i32,
    rhoa1: &[f64],
    rhob1: &[f64],
    sigmaaa1: &[f64],
    sigmabb1: &[f64],
    sigmaab1: &[f64],
    zk: &mut [f64],
    vrhoa: &mut [f64],
    vrhob: &mut [f64],
    vsigmaaa: &mut [f64],
    vsigmabb: &mut [f64],
    vsigmaab: &mut [f64],
    v2rhoa2: &mut [f64],
    v2rhob2: &mut [f64],
    v2rhoab: &mut [f64],
    v2rhoasigmaaa: &mut [f64],
    v2rhoasigmaab: &mut [f64],
    v2rhoasigmabb: &mut [f64],
    v2rhobsigmabb: &mut [f64],
    v2rhobsigmaab: &mut [f64],
    v2rhobsigmaaa: &mut [f64],
    v2sigmaaa2: &mut [f64],
    v2sigmaaaab: &mut [f64],
    v2sigmaaabb: &mut [f64],
    v2sigmaab2: &mut [f64],
    v2sigmaabbb: &mut [f64],
    v2sigmabb2: &mut [f64],
    npt: usize,
) {
    const TOL: f64 = 1.0e-20;

    if ideriv == 0 {
        for i in 0..npt {
            let rhoa = rhoa1[i].max(0.0);
            let rhob = rhob1[i].max(0.0);
            let rho = rhoa + rhob;
            if rho > TOL {
                if rhoa < TOL {
                    let t1 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhob;
                } else if rhob < TOL {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa;
                } else {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    let t4 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa - 0.9305257363491e0 * t4 * rhob;
                }
            } else {
                zk[i] = 0.0;
            }
        }
    } else if ideriv == 1 {
        for i in 0..npt {
            let rhoa = rhoa1[i].max(0.0);
            let rhob = rhob1[i].max(0.0);
            let rho = rhoa + rhob;
            if rho > TOL {
                if rhoa < TOL {
                    let t1 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhob;
                    vrhoa[i] = 0.0;
                    vrhob[i] = -0.12407009817988e1 * t1;
                } else if rhob < TOL {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa;
                    vrhoa[i] = -0.12407009817988e1 * t1;
                    vrhob[i] = 0.0;
                } else {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    let t4 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa - 0.9305257363491e0 * t4 * rhob;
                    vrhoa[i] = -0.12407009817988e1 * t1;
                    vrhob[i] = -0.12407009817988e1 * t4;
                }
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vrhob[i] = 0.0;
            }
        }
    } else if ideriv == 2 {
        for i in 0..npt {
            let rhoa = rhoa1[i].max(0.0);
            let rhob = rhob1[i].max(0.0);
            let rho = rhoa + rhob;
            if rho > TOL {
                if rhoa < TOL {
                    let t1 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhob;
                    vrhoa[i] = 0.0;
                    vrhob[i] = -0.12407009817988e1 * t1;
                    v2rhoa2[i] = 0.0;
                    v2rhoab[i] = 0.0;
                    let t5 = t1 * t1;
                    v2rhob2[i] = -0.4135669939329333e0 / t5;
                } else if rhob < TOL {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa;
                    vrhoa[i] = -0.12407009817988e1 * t1;
                    vrhob[i] = 0.0;
                    let t5 = t1 * t1;
                    v2rhoa2[i] = -0.4135669939329333e0 / t5;
                    v2rhob2[i] = 0.0;
                    v2rhoab[i] = 0.0;
                } else {
                    let t1 = rhoa.powf(1.0 / 3.0);
                    let t4 = rhob.powf(1.0 / 3.0);
                    zk[i] = -0.9305257363491e0 * t1 * rhoa - 0.9305257363491e0 * t4 * rhob;
                    vrhoa[i] = -0.12407009817988e1 * t1;
                    vrhob[i] = -0.12407009817988e1 * t4;
                    let t9 = t1 * t1;
                    v2rhoa2[i] = -0.4135669939329333e0 / t9;
                    let t12 = t4 * t4;
                    v2rhob2[i] = -0.4135669939329333e0 / t12;
                    v2rhoab[i] = 0.0;
                }
            } else {
                zk[i] = 0.0;
                vrhoa[i] = 0.0;
                vrhob[i] = 0.0;
                v2rhoa2[i] = 0.0;
                v2rhob2[i] = 0.0;
                v2rhoab[i] = 0.0;
            }
        }
    }
}
