# Code-for-Ischemic-Preconditioning-on-Post-Exercise-Injury-Biomarkers-DOMS-and-Endothelial-Function.
This repository contains the raw data, code, and all analysis results from the included studies.
 #  Protective Effects of Ischemic Preconditioning combined Acute Exercise On human body
# === need package ===
install.packages("pacman")
pacman::p_load(
  smplot2,
  meta,
  readxl,
  metafor,
  dmetar,
  magrittr,
  tidyverse,
  tidybayes,
  devtools,
  clubSandwich,
  psych,
  broom,
  knitr,
  kableExtra,
  scales,
  rms,
  gridExtra,
  ggExtra,
  rjags,
  igraph,
  splines,
  mgcv,
  robvis,
  orchaRd,
  bayestestR,
  multcomp,
  ggplotify,
  patchwork,
  stringr,
  forcats
)

# need function 
three_level_results <- function(res_final, res_three_level){
  
  res_three_level_i2=i2_ml(res_three_level)
  
  I1 <- paste(round(as.numeric(res_three_level_i2[2]),2), "%", sep = "")
  I2 <- paste(round(as.numeric(res_three_level_i2[3]),2), "%", sep = "")
  
  results <- data.frame(
    K = res_final$s.nlevels[1], 
    N = res_final$s.nlevels[2], 
    Estimate = res_final$b,     
    t_value = res_final$zval,   
    df = res_final$df,          
    p_value = res_final$pval,   
    CI_Lower = res_final$ci.lb, 
    CI_Upper = res_final$ci.ub, 
    PI_Lower = predict(res_final)$pi.lb,
  )
  
  n1<- data.frame(round(res_three_level$pval,4))
  n2<- data.frame(round(res_three_level$ci.lb,4))
  n3<- data.frame(round(res_three_level$ci.ub,4))
  n4<- cbind(n1, n2,n3) 
  colnames(n4)<- c("pval", "ci_lb", "ci_ub")
  
  x3<- n4 %>%   mutate(across(c(pval), ~ round(., 4)),
                       across(c(ci_lb, ci_ub), ~ sprintf("%.4f", .))) %>%
    mutate(p.value_adj = ifelse(pval == 0, "0.000", sprintf("%.4f", pval)),
           conf_interval_adj = paste("[", ci_lb, "; ", ci_ub, "]", sep = "")) %>%
    select(p.value_adj,conf_interval_adj)  
  
  results <- 
    cbind(results,x3)%>%
    mutate(
      Estimate = round(Estimate, 4), 
      t_value = round(t_value, 4),
      df = round(df, 2),
      p_value = ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 4)),
      CI_Lower = round(CI_Lower, 4),
      CI_Upper = round(CI_Upper, 4),
      PI_Lower = round(PI_Lower, 4),
      PI_Upper = round(PI_Upper, 4),
    ) %>%
    mutate(conf_interval = paste("[", CI_Lower, "; ", CI_Upper, "]", sep = ""),
           PI_interval = paste("[", PI_Lower, "; ", PI_Upper, "]", sep = "")) %>%
    mutate( 
      Qtest = ifelse(row_number() == 1, round(res_three_level$QE[1], 3), ""),I1,I2
      
    ) %>%
    select(K, N, Estimate, conf_interval, conf_interval_adj, p_value,p.value_adj,
           t_value, Qtest, I1,I2, PI_interval)
  
  x1<-  kable(results, caption = "Multivariate Meta-Analysis Model Results", format = "html") %>%
    kable_styling()
  
  return(x1)
} # three-level model 
display_res <- function(res) {
  for (i in 1:length(res)) {
    formatted_value <- format(res[[i]], scientific = FALSE, digits = 5)   
    cat("res_four_level_i2[", i, "] = ", formatted_value, "\n", sep = "")  
  }
}
egger_test <- function(data){
  
  data$sei <- sqrt(data$vi)
  
  model <- rma.mv(
    yi = yi,
    V = vi,
    random = ~ 1 | study/id,
    mods = ~ sei,
    data = data,
    method = "REML"
  )
  
  result <- coef_test(model, vcov = "CR2")
  
  return(result)
}
subgroup_analysis <- function(model, data2) {
  res_i2 <- i2_ml(model)
  I1 <- paste0(round(as.numeric(res_i2[2]), 2), "%")  
  I2 <- paste0(round(as.numeric(res_i2[3]), 2), "%")  
  overall_text <- paste0("F (", model$QMdf[1], " , ", model$QEdf[1], ") = ", round(model$QM[1], 3))
  x2 <- broom::tidy(model, conf.int = TRUE, conf.level = 0.95) %>%
    dplyr::mutate(
      conf_interval = paste0("[", round(conf.low, 4), "; ", round(conf.high, 4), "]"),
      Overall      = ifelse(row_number() == 1, overall_text, ""),
      Qtest        = ifelse(row_number() == 1, round(model$QE[1], 3), ""),
      I1           = ifelse(row_number() == 1, I1, ""),
      I2           = ifelse(row_number() == 1, I2, "")
    ) %>%
    dplyr::mutate(
      across(c(std.error, p.value, statistic), ~ round(., 4)),
      p.value = ifelse(p.value == 0, "0.000", sprintf("%.4f", p.value))
    ) %>%
    dplyr::select(term, estimate, conf_interval,
                  statistic, p.value,
                  Overall, Qtest, I1, I2) 
  res2 <- metafor::robust(model, cluster = data2$study, clubSandwich = TRUE)
  x3 <- data.frame(
    p.value_adj      = ifelse(res2$pval == 0, "0.000", sprintf("%.4f", round(res2$pval, 4))),
    conf_interval_adj = paste0("[",
                               sprintf("%.4f", round(res2$ci.lb, 4)), "; ",
                               sprintf("%.4f", round(res2$ci.ub, 4)), "]")
  )
  
  K_list  <- vapply(seq_len(model$QMdf[1]),
                    function(i) length(unique(data2$study[model$X.f[,i]==1])),
                    integer(1))
  ES_list <- vapply(seq_len(model$QMdf[1]),
                    function(i) sum(model$X.f[,i]),
                    numeric(1))
  
  result <- cbind(
    K = K_list,                                  
    `Effect size(n)` = ES_list,                  
    `Effect size` = round(x2$estimate, 4),       
    `95%CI` = x2$conf_interval,                  
    `95%CI adj` = x3$conf_interval_adj,          
    `pval` = x2$p.value,                         
    `p-value adj` = x3$p.value_adj,              
    `t-value` = x2$statistic,                    
    Overall = x2$Overall,                        
    Qtest = x2$Qtest,                            
    `I1` = x2$I1,                                
    `I2` = x2$I2                                 
  )
  return(as.data.frame(result))
}
three_level_results1 <- function(res_three_level){
  
  res_three_level_i2=i2_ml(res_three_level)
  
  I1 <- paste(round(as.numeric(res_three_level_i2[2]),2), "%", sep = "")
  I2 <- paste(round(as.numeric(res_three_level_i2[3]),2), "%", sep = "")
  
  results <- data.frame(
    K = res_three_level$s.nlevels[1],
    N = res_three_level$s.nlevels[2],
    Estimate = res_three_level$b,
    t_value = res_three_level$zval,
    df = res_three_level$df,
    p_value = res_three_level$pval,
    CI_Lower = res_three_level$ci.lb,
    CI_Upper = res_three_level$ci.ub,
    PI_Lower = predict(res_three_level)$pi.lb,
    PI_Upper = predict(res_three_level)$pi.ub
  )
  n1<- data.frame(round(res_three_level$pval,4))
  n2<- data.frame(round(res_three_level$ci.lb,4))
  n3<- data.frame(round(res_three_level$ci.ub,4))
  n4<- cbind(n1, n2,n3) 
  colnames(n4)<- c("pval", "ci_lb", "ci_ub")
  x3<- n4 %>%   mutate(across(c(pval), ~ round(., 4)),
                       across(c(ci_lb, ci_ub), ~ sprintf("%.4f", .))) %>%
    mutate(p.value_adj = ifelse(pval == 0, "0.000", sprintf("%.4f", pval)),
           conf_interval_adj = paste("[", ci_lb, "; ", ci_ub, "]", sep = "")) %>%
    dplyr::select(p.value_adj,conf_interval_adj)  
  
  results <- 
    cbind(results,x3)%>%
    mutate(
      Estimate = round(Estimate, 4), 
      t_value = round(t_value, 4),
      df = df,
      p_value = ifelse(p_value < 0.0001, "< 0.0001", round(p_value, 4)),
      CI_Lower = round(CI_Lower, 4),
      CI_Upper = round(CI_Upper, 4),
      PI_Lower = round(PI_Lower, 4),
      PI_Upper = round(PI_Upper, 4),
    ) %>%
    mutate(conf_interval = paste("[", CI_Lower, "; ", CI_Upper, "]", sep = ""),
           PI_interval = paste("[", PI_Lower, "; ", PI_Upper, "]", sep = "")) %>%
    mutate( 
      Qtest = ifelse(row_number() == 1, round(res_three_level$QE[1], 3), ""),I1,I2
      
    ) %>%
    dplyr::select(K, N, Estimate, conf_interval, conf_interval_adj, p_value,p.value_adj,
           t_value, Qtest, I1,I2, PI_interval)
  
  
  x1<-  kable(results, caption = "Multivariate Meta-Analysis Model Results", format = "html") %>%
    kable_styling()
  
  return(  results)
} # three-level model 
leave1out_threelevel <- function(data, r = 0.6){
  
  studies <- unique(data$study)
  
  results <- lapply(studies, function(s){
    
    data_sub <- data[data$study != s, ]
    
    V <- with(data_sub,
              impute_covariance_matrix(
                vi = vi,
                cluster = study,
                r = r
              ))
    
    model <- rma.mv(
      yi = yi,
      V = V,
      data = data_sub,
      random = ~ 1 | study/id,
      control = list(optimizer = "optim"),
      method = "REML",
      test = "t",
      slab = study
    )
    res_three_level_i2=i2_ml( model)
    I1 <- paste(round(as.numeric(res_three_level_i2[2]),2), "%", sep = "")
    I2 <- paste(round(as.numeric(res_three_level_i2[3]),2), "%", sep = "")
    model_robust <- robust(
      model,
      cluster = data_sub$study,
      clubSandwich = TRUE
    )
    est <- model_robust$beta
    se  <- model_robust$se
    
    data.frame(
      study_removed = s,
      estimate = est,
      SE = se,
      CI_low = est - 1.96*se,
      CI_high = est + 1.96*se,
      p = model_robust$pval,
      I1 = I1,
      I2=I2
    )
    
  })
  
  dplyr::bind_rows(results)
}
r_sensitivity <- function(data){
  
  r_values <- seq(0.5, 0.9, by = 0.1)
  
  results <- lapply(r_values, function(r){
    
    V <- with(data,
              impute_covariance_matrix(
                vi = vi,
                cluster = study,
                r = r
              ))
    
    model1 <- rma.mv(
      yi = yi,
      V = V,
      data = data,
      random = ~ 1 | study/id,
      control = list(optimizer = "optim"),
      method = "REML",
      test = "t",
      slab = study
    )
    
    model <- robust(model1,cluster = study,clubSandwich=TRUE)
    
    res <- three_level_results1(model)
    
    res$r <- r
    
    return(res)
    
  })
  
  result_table <- dplyr::bind_rows(results)
  
  return(result_table)
}

#====  Effects of CK ====
data_CK<- read_excel("C:/Users/Administrator/Desktop/Data and R code.xlsx", sheet = 1)
data_CK_SMD <- escalc(measure = "SMD",m1i = IPCmeanchange, m2i =  CONmeanchange,sd1i = IPCsdchange, 
                sd2i = CONsdchange,n1i = IPCn, n2i = CONn,data = data_CK,slab=study)
#`impute_covariance_matrix()` was deprecated in clubSandwich 0.5.11. We use `metafor::vcalc()` instead.
V <- metafor::vcalc(
  vi = vi,
  cluster = studyid,
  obs = 1:nrow(data_CK_SMD),
  data = data_CK_SMD,
  rho = 0.7
)
res_CK <- rma.mv(yi = yi, V = V,data = data_CK_SMD,random = ~ 1 |  study/id,
              method = "REML",  test = "t", slab = study )
res_CK_final<- robust(res_CK,cluster = study,clubSandwich=TRUE)
res_CK_final
i2_ml(res_CK_final)
egger_test(data_CK_SMD) 
funnel(res_CK_final) # Code to generate the funnel in Figure 3(a)
weights <- weights(res_CK_final)
weights <- round(weights / sum(weights) * 100, 1)
predict(res_CK_final)

#  CK forest plot(code to generate the forest plot in Figure 4)
par(tck=-0.01, mgp=c(1,0,0), mar=c(3,4,3,2))
dd <- c(0,diff(data_CK$studyid))
rows <- (1:res_CK_final$k) + cumsum(dd)
x2 <- forest(res_CK_final, rows = rows, ylim = c(-6, max(rows) + 3), xlim = c(-22, 14), cex = 1,
             ilab = cbind(
               IPCn,
               sprintf("%.2f", IPCmeanchange),
               sprintf("%.2f", IPCsdchange),
               sprintf("%.0f", CONn),
               sprintf("%.2f", CONmeanchange),
               sprintf("%.2f", CONsdchange),
               sprintf("%.2f%%", weights)
               ),
             ilab.xpos = c(-13.5, -12, -10, -8.5, -7, -5,7), 
             efac = c(0, 0.5), header = T, mlab = "Pooled Estimate", shade = TRUE, addpred = TRUE)
text(c(-13.5, -12, -10, -8.5, -7, -5,7), max(rows)+2,cex=1, c("N1", "Mean1", "SD1", "N2","Mean2", "SD2","weight"), font = 2) 
text(c(-10.5,-5.5),    max(rows)+3.5, cex=1,c("EXP", "CON"), font = 2) 
text(c(0),max(rows)+2,cex=1, c("Standardized mean difference, III,\nREML (95% CI)"), font = 3) 
abline(h = rows[c(1,diff(rows)) == 3] - 1, lty="dotted") 
text(x2$xlim[1], -1.75, pos=4, cex=1,bquote(paste("Test for heterogeneity: ",
                                                    tau^2 , "2.1=", .(fmtx( res_CK_final$sigma2[1], digits=4)), "; ",
                                                    tau^2 , "2.2=", .(fmtx( res_CK_final$sigma2[2], digits=4)), "; ",
                                                    chi^2, "=", .(fmtx( res_CK_final$QE, digits=2)),
                                                    ", df=", .( res_CK_final$k -  res_CK_final$p), ", ",
                                                    .(fmtp( res_CK_final$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)))))
text(x2$xlim[1], -2.5, pos=4, cex=1,bquote(paste("Test for overall effect: ",
                                                   "Z=", .(fmtx( res_CK_final$zval, digits=2)), ", ",
                                                   .(fmtp( res_CK_final$pval, digits=3, pname="P", add0=TRUE, equal=TRUE)))))



text(x2$xlim[1], -3.25, pos = 4, cex = 1, "Predict interval:")
                                               
pred_final <- predict(res_CK_final)
segments(x0 = pred_final$pi.lb,
         x1 = pred_final$pi.ub,
         y0 = -3.25,  
         y1 = -3.25,
         col = "red", lwd = 2)
text(x = 12,  
     y = -3.25,  
     labels = paste0("[", 
                     sprintf("%.2f", pred_final$pi.lb), 
                     "; ", 
                     sprintf("%.2f", pred_final$pi.ub), 
                     "]"),
     cex = 1, font = 2)

#subgroup_analysis (code to generate the values in Table 5)
res_CK_window <- rma.mv(yi = yi, V = V,data = data_CK_SMD,random = ~ 1 |  study/id,
                 mods = ~ Window-1,
                 method = "REML",  test = "t", slab = study )
summary(res_CK_window)
subgroup_analysis(res_CK_window,data_CK_SMD)

L <- rbind(
  "0-2 vs 48"  = c( 1,  0, -1,  0),  # 0-2 vs 48
  "0-2 vs 24"  = c( 1, -1,  0,  0),  # 0-2 vs 24
  "0-2 vs 72"  = c( 1,  0,  0, -1),  # 0-2 vs 72
  "24 vs 48"   = c( 0,  1, -1,  0),  # 24 vs 48
  "24 vs 72"   = c( 0,  1,  0, -1),  # 24 vs 72
  "48 vs 72"   = c( 0,  0,  1, -1)   # 48 vs 72
)
res_contrast <- glht(res_CK_window, linfct = L)
summary(res_contrast)
confint(res_contrast)         

#Regression_linear and non-linear
res_CK_details <- rma.mv(yi = yi, V = V,data = data_CK_SMD,random = ~ 1 |  study/id,
                        mods = ~ Details,
                        method = "REML",  test = "t", slab = study )

res_CK_details_poly2 <- rma.mv(yi = yi,
                 V = V,
                 data = data_CK_SMD,
                 random = ~ 1 |  study/id,
                 method = "REML", 
                 test = "t",
                 slab = study,
                 mods= ~  poly(Details, degree=2, raw=TRUE))

res_CK_details_poly3 <- rma.mv(yi = yi,
                V = V,
                data = data_CK_SMD,
                random = ~ 1 |  study/id,
                method = "REML", 
                test = "t",
                slab = study,
                mods= ~  poly(Details, degree=3, raw=TRUE))

#Regression_RCS
res_CK_details_rcs3 <- rma.mv(yi = yi,
                               V = V,
                               data = data_CK_SMD,
                               random = ~ 1 |  study/id,
                               method = "REML", 
                               test = "t",
                               slab = study,
                               mods= ~  rcs(Details, 3))

res_CK_details_rcs4 <- rma.mv(yi = yi,
                              V = V,
                              data = data_CK_SMD,
                              random = ~ 1 |  study/id,
                              method = "REML", 
                              test = "t",
                              slab = study,
                              mods= ~  rcs(Details, 4))

# Code to generate the Supplementary Material in Table 5
fitstats(res_CK_details,
         res_CK_details_poly2,
         res_CK_details_poly3,
         res_CK_details_rcs3, 
         res_CK_details_rcs4 
         )

#choose and found best to Plot, here is lieanr-noly is best: AIC,combing RVE, then plotting.
# CK Meta-regression(code to generate the Meta-regression in Figure 8(a))
res_CK_details_final<- robust(res_CK_details ,cluster = study,clubSandwich=TRUE)
res_CK_details_final
regplot(res_CK_details_final) 

# ====  CK sensitivity ==== 、
x1<- r_sensitivity(data_CK_SMD)

kable(x1, format = "html") %>%
  kable_styling()

# CK_r forest plot(code to generate the Supplementary in Figure 1)
data <- data.frame(
  r = c(0.5, 0.6, 0.7, 0.8, 0.9),
  K = rep(6, 5),
  N = rep(18, 5),
  
  estimate = c(-0.3856, -0.3815, -0.3774, -0.3732, -0.3688),
  ci_lower = c(-1.2751, -1.2659, -1.2566, -1.2473, -1.2379),
  ci_upper = c(0.5040, 0.5028, 0.5018, 0.5010, 0.5002),
  
  p_value = c(0.3154, 0.3174, 0.3195, 0.3218, 0.3242),
  t_value = c(-1.1159, -1.1111, -1.1058, -1.1001, -1.0941),
  
  Qtest = c(43.749, 47.081, 53.480, 67.275, 109.812),
  I1 = c(67.37, 63.62, 59.88, 56.14, 52.38),
  I2 = c(7.59, 11.20, 14.80, 18.41, 22.02),
  
  PI_lower = c(-2.5841, -2.5719, -2.5596, -2.5472, -2.5348),
  PI_upper = c(1.8130, 1.8088, 1.8048, 1.8008, 1.7971)
)
data <- data %>%
  mutate(
    label_R = sprintf("R = %.1f", r),
    
    label_stats = sprintf("K=%d  N=%d\nQ=%.2f  I²=%.2f%%  I¹=%.2f%%", 
                          K, N, Qtest, I2, I1),
    
    label_effect = sprintf("%.2f (%.2f, %.2f)", 
                           estimate, ci_lower, ci_upper),
    
    label_p = sprintf("p = %.3f", p_value)
  )
ggplot(data, aes(x = estimate, y = fct_rev(label_R))) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = 0.15, size = 0.7, color = "grey") +
  geom_point(size = 3, shape = 18, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
  geom_text(aes(x = -2.69, label = label_stats), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.575, label = label_effect), size = 3.5) +
  geom_text(aes(x = 1.29, label = label_p), size = 3.5) +
  annotate("text", x = -2.69, y = 5.53, label = "Heterogeneity Statistics", 
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 1.195, y = 5.53, label = "Effect Size (95% CI)", 
           hjust = 1.42, fontface = "bold", size = 3.5) +
  annotate("text", x = 1.545, y = 5.53, label = "p-value", 
           hjust = 1.5, fontface = "bold", size = 3.5) +
  scale_x_continuous(limits = c(-3, 1.7), 
                     breaks = seq(-1.5, 0, 0.5),
                     expand = c(0, 0)) +
  labs(title = "MP——Sensitivity Multivariate Analysis with Different r Values", 
       x = "Effect Size", 
       y = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.3, size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_line(color = "grey90", size = 0.2),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

#====  CK sensitivity ====
# sensitive_leave1out(code to generate the Value in Supplementary Material Table 8)
leave1out_threelevel(data_CK_SMD, r= 0.7)

# sensitive_cook(code to generate the Value in Supplementary Material Table 11)
cooks <- cooks.distance(res_CK_final) 
# View outliers with Cooks > 3 * mean(code to generate the Value in Supplementary Material Table 11)
outliers_cooks <- cooks %>% 
  cbind(data_CK_SMD$esid) %>%         
  subset(cooks > 3.0*mean(cooks)) %>% 
  View()

data_CK_leave_cooks <- data_CK_SMD[-18, ]

V_cooks  <- with(data_CK_leave_cooks, impute_covariance_matrix(vi = vi,cluster = studyid,r = 0.7))

res_CK_cooks  <- rma.mv(yi = yi, V = V_cooks,data = data_CK_leave_cooks,random = ~ 1 |  study/id,
                 method = "REML",  test = "t", slab = study )

res_CK_cooks_final<- robust(res_CK_cooks,cluster = study,clubSandwich=TRUE)
res_CK_cooks_final

#code to generate the Value in Supplementary Material Table 11
#Based on reviewer's valuable selective suggestions, we did the cluster of study and then run the analysis again.
x6<- aggregate(data_CK_SMD, cluster=study, rho=0.7)
res<- rma(yi,vi, data =x6)
forest(res)
summary(res)

#sensitive_residuals(code to generate the Value in Supplementary Material Table 11)
resid <- residuals(res_CK_final) %>%
  scale(center = F, scale = T)
#View outliears with residuals > 3 * mean
outliers_resid <- resid %>%
  cbind(data_CK_SMD$esid) %>%               
  subset(resid > 3.0 | resid < - 3.0) %>%  
  View()
#Non-outlier

#sensitive_hat(code to generate the Value in Supplementary Material Table 11)
high_hat <- data.frame(hatvalues.rma.mv(res_CK_final)) %>% 
  rename(hat = hatvalues.rma.mv.res_CK_final.) %>% 
  filter(hat > 3 * mean(hat))
# Non-outlier

# trim-and-fill 
res_tf <- trimfill(res_simple)
funnel(res_tf, main = "Funnel Plot with Trim-and-Fill")
summary(res_tf)
regtest(res_simple, model = "lm", predictor = "sei") 

# 1. PET regression (code to generate the Forest Plot in Supplementary Material Figure 5)
pet <- rma(yi = yi, sei = sei, mods = ~ sei, method = "REML", data = data2)
summary(pet)
forest(pet)

# 2. PEESE regression(code to generate the Forest Plot in Supplementary Material Figure 5)
peese <- rma(yi = yi, sei = sei, mods = ~ I(sei^2), method = "REML", data = data2)
summary(peese)
forest(peese)


#====  Effects of MP ====
data_MP <- read_excel("C:/Users/Administrator/Desktop/Data and R code.xlsx", sheet = 2)

data_MP_SMD <- escalc(measure = "SMD", m1i = IPCmeanchange, m2i = CONmeanchange,
                      sd1i = IPCsdchange, sd2i = CONsdchange,
                      n1i = IPCn, n2i = CONn, data = data_MP, slab = study)
# `impute_covariance_matrix()` was deprecated in clubSandwich 0.5.11. We use `metafor::vcalc()` instead.
V <- metafor::vcalc(
  vi = vi,
  cluster = studyid,
  obs = 1:nrow(data_MP_SMD),
  data = data_MP_SMD,
  rho = 0.7
)
res_MP <- rma.mv(yi = yi, V = V, data = data_MP_SMD,
                 random = ~ 1 | study/id,
                 method = "REML", test = "t", slab = study)

res_MP_final <- robust(res_MP, cluster = study, clubSandwich = TRUE)

res_MP_final
i2_ml(res_MP_final)                     # i2
egger_test(data_MP_SMD)                  # 
funnel(res_MP_final)                     # Code to generate the funnel in Figure 3(b)
weights <- weights(res_MP_final)
weights <- round(weights / sum(weights) * 100, 1)
predict(res_MP_final)

#  MP forest plot (code to generate the forest plot in Figure 5)
par(tck = -0.01, mgp = c(1, 0, 0), mar = c(3, 4, 3, 2))
dd <- c(0, diff(data_MP$studyid))
rows <- (1:res_MP_final$k) + cumsum(dd)

x2 <- forest(res_MP_final,
             rows = rows,
             ylim = c(-6, max(rows) + 3),
             xlim = c(-22, 14),
             cex = 1,
             ilab = cbind(
               IPCn,
               sprintf("%.2f", IPCmeanchange),
               sprintf("%.2f", IPCsdchange),
               sprintf("%.0f", CONn),
               sprintf("%.2f", CONmeanchange),
               sprintf("%.2f", CONsdchange),
               sprintf("%.2f%%", weights)
             ),
             ilab.xpos = c(-13.5, -12, -10, -8.5, -7, -5, 7),
             efac = c(0, 0.5), header = T, mlab = "Pooled Estimate", shade = TRUE, addpred = TRUE)
text(c(-13.5, -12, -10, -8.5, -7, -5, 7),
     max(rows) + 2, cex = 1,
     c("N1", "Mean1", "SD1", "N2", "Mean2", "SD2", "weight"), font = 2)
text(c(-10.5, -5.5), max(rows) + 3.5, cex = 1,
     c("EXP", "CON"), font = 2)
text(c(0), max(rows) + 2, cex = 1,
     c("Standardized mean difference, III,\nREML (95% CI)"), font = 3)
abline(h = rows[c(1, diff(rows)) == 3] - 1, lty = "dotted")
text(x2$xlim[1], -2, pos = 4, cex = 1,
     bquote(paste("Test for heterogeneity: ",
                  tau^2 , "2.1=", .(fmtx( res_MP_final$sigma2[1], digits = 4)), "; ",
                  tau^2 , "2.2=", .(fmtx( res_MP_final$sigma2[2], digits = 4)), "; ",
                  chi^2, "=", .(fmtx( res_MP_final$QE, digits = 2)),
                  ", df=", .( res_MP_final$k -  res_MP_final$p), ", ",
                  .(fmtp( res_MP_final$QEp, digits = 2, pname = "P", add0 = TRUE, equal = TRUE)))))
text(x2$xlim[1], -3, pos = 4, cex = 1,
     bquote(paste("Test for overall effect: ",
                  "Z=", .(fmtx( res_MP_final$zval, digits = 2)), ", ",
                  .(fmtp( res_MP_final$pval, digits = 3, pname = "P", add0 = TRUE, equal = TRUE)))))

text(x2$xlim[1], -4, pos = 4, cex = 1, "Predict interval:")

pred_final <- predict(res_MP_final)
segments(x0 = pred_final$pi.lb,
         x1 = pred_final$pi.ub,
         y0 = -4,
         y1 = -4,
         col = "red", lwd = 2)
text(x = 12,
     y = -4,
     labels = paste0("[", sprintf("%.2f", pred_final$pi.lb),
                     "; ", sprintf("%.2f", pred_final$pi.ub), "]"),
     cex = 1, font = 2)

# subgroup_analysis（code to generate the values in Table 5）
res_MP_window <- rma.mv(yi = yi, V = V, data = data_MP_SMD,
                        random = ~ 1 | study/id,
                        mods = ~ Window - 1,
                        method = "REML", test = "t", slab = study)
summary(res_MP_window)
subgroup_analysis(res_MP_window, data_MP_SMD)
L <- rbind(
  "0-2 vs 48"  = c(1,  0, -1,  0),   # 0-2 vs 48
  "0-2 vs 24"  = c(1, -1,  0,  0),   # 0-2 vs 24
  "0-2 vs 72"  = c(1,  0,  0, -1),   # 0-2 vs 72
  "24 vs 48"   = c(0,  1, -1,  0),   # 24 vs 48
  "24 vs 72"   = c(0,  1,  0, -1),   # 24 vs 72
  "48 vs 72"   = c(0,  0,  1, -1)    # 48 vs 72
)
res_contrast <- glht(res_MP_window, linfct = L)
summary(res_contrast)          
confint(res_contrast)          

# Regression_linear and non-linear
res_MP_details <- rma.mv(yi = yi, V = V, data = data_MP_SMD,
                         random = ~ 1 | study/id,
                         mods = ~ Details,
                         method = "REML", test = "t", slab = study)

res_MP_details_poly2 <- rma.mv(yi = yi,
                               V = V,
                               data = data_MP_SMD,
                               random = ~ 1 | study/id,
                               method = "REML",
                               test = "t",
                               slab = study,
                               mods = ~ poly(Details, degree = 2, raw = TRUE))

res_MP_details_poly3 <- rma.mv(yi = yi,
                               V = V,
                               data = data_MP_SMD,
                               random = ~ 1 | study/id,
                               method = "REML",
                               test = "t",
                               slab = study,
                               mods = ~ poly(Details, degree = 3, raw = TRUE))

# Regression_RCS
res_MP_details_rcs3 <- rma.mv(yi = yi,
                              V = V,
                              data = data_MP_SMD,
                              random = ~ 1 | study/id,
                              method = "REML",
                              test = "t",
                              slab = study,
                              mods = ~ rcs(Details, 3))

res_MP_details_rcs4 <- rma.mv(yi = yi,
                              V = V,
                              data = data_MP_SMD,
                              random = ~ 1 | study/id,
                              method = "REML",
                              test = "t",
                              slab = study,
                              mods = ~ rcs(Details, 4))

# Code to generate the Supplementary Material in Table 6
fitstats(res_MP_details,
         res_MP_details_poly2,
         res_MP_details_poly3,
         res_MP_details_rcs3,
         res_MP_details_rcs4)

# choose and found best to Plot, here is linear-only is best: AIC, combing RVE, then plotting.
# MP Meta-regression (code to generate the Meta-regression in Figure 8(b))
res_MP_details_final <- robust(res_MP_details, cluster = study, clubSandwich = TRUE)
res_MP_details_final
regplot(res_MP_details_final)

#====  MP sensitivity ====
x1 <- r_sensitivity(data_MP_SMD)

kable(x1, format = "html") %>%
  kable_styling()

# MP_r forest plot (code to generate the Supplementary Material in Figure 2)
data <- data.frame(
  r = c(0.5, 0.6, 0.7, 0.8, 0.9),
  estimate = c(0.866, 0.975, 1.133, 1.390, 1.938),
  ci_lower = c(0.083, 0.108, 0.141, 0.187, 0.258),
  ci_upper = c(1.648, 1.843, 2.125, 2.592, 3.618),
  p_value = c(0.0302, 0.0275, 0.0252, 0.0235, 0.0238),
  PI_lower = c(-1.966, -2.268, -2.700, -3.411, -4.992),
  PI_upper = c(3.698, 4.219, 4.966, 6.191, 8.868),
  I2 = c(58.6, 66.0, 73.3, 80.7, 88.2),
  Q_p = c(0.0894, 0.053, 0.0235, 0.0056, 0.0002)
)
data <- data %>%
  mutate(
    label_R = sprintf("R = %.1f", r),
    label_effect = sprintf("%.3f [%.3f, %.3f]", estimate, ci_lower, ci_upper),
    label_p = sprintf("p = %.4f", p_value),
    label_hetero = sprintf("I²=%.1f%%  Q p=%.4f", I2, Q_p)
  )
x_hetero <- min(data$ci_lower) - 1.2
x_effect <- max(data$ci_upper) + 0.2
x_p <- max(data$ci_upper) + 1.2
ggplot(data, aes(x = estimate, y = fct_rev(label_R))) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, size = 0.7, color = "grey") +
  geom_point(size = 3, shape = 18, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
  geom_text(aes(x = x_hetero, label = label_hetero), hjust = 0, size = 3.5) +
  geom_text(aes(x = x_effect, label = label_effect), hjust = 0, size = 3.5) +
  geom_text(aes(x = x_p, label = label_p), hjust = 0, size = 3.5) +
  annotate("text", x = x_hetero, y = 5.5, label = "Heterogeneity Stats",
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = x_effect, y = 5.5, label = "Effect Size (95% CI)",
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = x_p, y = 5.5, label = "p-value",
           hjust = 0, fontface = "bold", size = 3.5) +
  scale_x_continuous(limits = c(x_hetero - 0.5, x_p + 0.5), expand = c(0, 0)) +
  labs(title = "MP——Sensitivity Multivariate Analysis with Different r Values (Random Effects)",
       x = "Effect Size (SMD)",
       y = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_line(color = "grey90", size = 0.2),
    plot.margin = margin(1, 2, 1, 1, "cm")
  )



# sensitive_leave1out (Code to generate the Supplementary Material in Table 9)
leave1out_threelevel(data_MP_SMD, r = 0.7)

# sensitive_cook(Code to generate the Supplementary Material in Table 12)
cooks <- cooks.distance(res_MP_final)
# View outliers with Cooks > 3 * mean
outliers_cooks <- cooks %>%
  cbind(data_MP_SMD$esid) %>%            # bind study names for reference
  subset(cooks > 3.0 * mean(cooks)) %>%  # subset outliers
  View()
rows_to_remove <- c(9, 13, 38)   #outliers
data_MP_leave_cooks <- data_MP_SMD[-rows_to_remove, ]
V_cooks <- with(data_MP_leave_cooks,
                impute_covariance_matrix(vi = vi, cluster = studyid, r = 0.7))
res_MP_cooks <- rma.mv(yi = yi, V = V_cooks, data = data_MP_leave_cooks,
                       random = ~ 1 | study/id,
                       method = "REML", test = "t", slab = study)

res_MP_cooks_final <- robust(res_MP_cooks, cluster = study, clubSandwich = TRUE)
res_MP_cooks_final
#(Code to generate the Supplementary Material in Table 12)
# Based on reviewer's valuable selective suggestions, we did the cluster of study and then run the analysis again.
x6 <- aggregate(data_MP_SMD, cluster = study, rho = 0.7)
res <- rma(yi, vi, data = x6)
forest(res)
summary(res)
# sensitive_residuals(Code to generate the Supplementary Material in Table 12)
resid <- residuals(res_MP_final) %>%
  scale(center = FALSE, scale = TRUE)
# View outliers with residuals > 3 * mean
outliers_resid <- resid %>%
  cbind(data_MP_SMD$esid) %>%
  subset(resid > 3.0 | resid < -3.0) %>%
  View()
# Non-outlier

# sensitive_hat(Code to generate the Supplementary Material in Table 12)
high_hat <- data.frame(hatvalues.rma.mv(res_MP_final)) %>%
  rename(hat = hatvalues.rma.mv.res_MP_final.) %>%
  filter(hat > 3 * mean(hat))
# Non-outlier

#====  Effects of cTn ====
data_cTn <- read_excel("C:/Users/Administrator/Desktop/Data and R code.xlsx", sheet = 3)

data_cTn_SMD <- escalc(measure = "SMD", m1i = IPCmeanchange, m2i = CONmeanchange,
                       sd1i = IPCsdchange, sd2i = CONsdchange,
                       n1i = IPCn, n2i = CONn, data = data_cTn, slab = study)
# `impute_covariance_matrix()` was deprecated in clubSandwich 0.5.11. We use `metafor::vcalc()` instead.
V <- metafor::vcalc(
  vi = vi,
  cluster = studyid,
  obs = 1:nrow(data_cTn_SMD),
  data = data_cTn_SMD,
  rho = 0.6
)

res_cTn <- rma.mv(yi = yi, V = V, data = data_cTn_SMD,
                  random = ~ 1 | study/id,
                  method = "REML", test = "t", slab = study)

res_cTn_final <- robust(res_cTn, cluster = study, clubSandwich = TRUE)

res_cTn_final
i2_ml(res_cTn_final)                     
egger_test(data_cTn_SMD)                  
funnel(res_cTn_final)                     # Code to generate the funnel in Figure 3(c)
weights <- weights(res_cTn_final)
weights <- round(weights / sum(weights) * 100, 1)
predict(res_cTn_final)

#  cTn forest plot (code to generate the forest plot in Figure 7)
par(tck = -0.01, mgp = c(1, 0, 0), mar = c(3, 4, 3, 2))
dd <- c(0, diff(data_cTn$studyid))
rows <- (1:res_cTn_final$k) + cumsum(dd)

x2 <- forest(res_cTn_final,
             rows = rows,
             ylim = c(-6, max(rows) + 3),
             xlim = c(-22, 14),
             cex = 1,
             ilab = cbind(
               IPCn,
               sprintf("%.2f", IPCmeanchange),
               sprintf("%.2f", IPCsdchange),
               sprintf("%.0f", CONn),
               sprintf("%.2f", CONmeanchange),
               sprintf("%.2f", CONsdchange),
               sprintf("%.2f%%", weights)
             ),
             ilab.xpos = c(-13.5, -12, -10, -8.5, -7, -5, 7),
             efac = c(0, 0.5), header = T, mlab = "Pooled Estimate", shade = TRUE, addpred = TRUE)
text(c(-13.5, -12, -10, -8.5, -7, -5, 7),
     max(rows) + 2, cex = 1,
     c("N1", "Mean1", "SD1", "N2", "Mean2", "SD2", "weight"), font = 2)
text(c(-10.5, -5.5), max(rows) + 3.5, cex = 1,
     c("EXP", "CON"), font = 2)
text(c(0), max(rows) + 2, cex = 1,
     c("Standardized mean difference, III,\nREML (95% CI)"), font = 3)
abline(h = rows[c(1, diff(rows)) == 3] - 1, lty = "dotted")

text(x2$xlim[1], -1.75, pos = 4, cex = 1,
     bquote(paste("Test for heterogeneity: ",
                  tau^2 , "2.1=", .(fmtx( res_cTn_final$sigma2[1], digits = 4)), "; ",
                  tau^2 , "2.2=", .(fmtx( res_cTn_final$sigma2[2], digits = 4)), "; ",
                  chi^2, "=", .(fmtx( res_cTn_final$QE, digits = 2)),
                  ", df=", .( res_cTn_final$k -  res_cTn_final$p), ", ",
                  .(fmtp( res_cTn_final$QEp, digits = 2, pname = "P", add0 = TRUE, equal = TRUE)))))

text(x2$xlim[1], -2.5, pos = 4, cex = 1,
     bquote(paste("Test for overall effect: ",
                  "Z=", .(fmtx( res_cTn_final$zval, digits = 2)), ", ",
                  .(fmtp( res_cTn_final$pval, digits = 3, pname = "P", add0 = TRUE, equal = TRUE)))))

text(x2$xlim[1], -3.25, pos = 4, cex = 1, "Predict interval:")

pred_final <- predict(res_cTn_final)
segments(x0 = pred_final$pi.lb,
         x1 = pred_final$pi.ub,
         y0 = -3.25,
         y1 = -3.25,
         col = "red", lwd = 2)
text(x = 12,
     y = -3.25,
     labels = paste0("[", sprintf("%.2f", pred_final$pi.lb),
                     "; ", sprintf("%.2f", pred_final$pi.ub), "]"),
     cex = 1, font = 2)

# Regression_linear and non-linear
res_cTn_details <- rma.mv(yi = yi, V = V, data = data_cTn_SMD,
                          random = ~ 1 | study/id,
                          mods = ~ Details,
                          method = "REML", test = "t", slab = study)

res_cTn_details_poly2 <- rma.mv(yi = yi,
                                V = V,
                                data = data_cTn_SMD,
                                random = ~ 1 | study/id,
                                method = "REML",
                                test = "t",
                                slab = study,
                                mods = ~ poly(Details, degree = 2, raw = TRUE))

res_cTn_details_poly3 <- rma.mv(yi = yi,
                                V = V,
                                data = data_cTn_SMD,
                                random = ~ 1 | study/id,
                                method = "REML",
                                test = "t",
                                slab = study,
                                mods = ~ poly(Details, degree = 3, raw = TRUE))

# Regression_RCS
res_cTn_details_rcs3 <- rma.mv(yi = yi,
                               V = V,
                               data = data_cTn_SMD,
                               random = ~ 1 | study/id,
                               method = "REML",
                               test = "t",
                               slab = study,
                               mods = ~ rcs(Details, 3))

res_cTn_details_rcs4 <- rma.mv(yi = yi,
                               V = V,
                               data = data_cTn_SMD,
                               random = ~ 1 | study/id,
                               method = "REML",
                               test = "t",
                               slab = study,
                               mods = ~ rcs(Details, 4))

# Code to generate the Supplementary Material in Table 7
fitstats(res_cTn_details,
         res_cTn_details_poly2,
         res_cTn_details_poly3,
         res_cTn_details_rcs3,
         res_cTn_details_rcs4)

# choose and found best to Plot, here is linear-only is best: AIC, combing RVE, then plotting.
# cTn Meta-regression (code to generate the Meta-regression in Figure 8(c))
res_cTn_details_final <- robust(res_cTn_details, cluster = study, clubSandwich = TRUE)
res_cTn_details_final
regplot(res_cTn_details_final)

#====  cTn sensitivity ====
x1 <- r_sensitivity(data_cTn_SMD)

kable(x1, format = "html") %>%
  kable_styling()

# cTn_r forest plot (code to generate the Supplementary in Figure 4)
data <- data.frame(
  r = c(0.5, 0.6, 0.7, 0.8, 0.9),
  K = rep(3, 5),
  N = rep(12, 5),
  estimate = c(0.0799, 0.0847, 0.0897, 0.0948, 0.0999),
  ci_lower = c(-1.7918, -1.788, -1.7851, -1.7834, -1.7828),
  ci_upper = c(1.9516, 1.9574, 1.9645, 1.9729, 1.9826),
  p_value = c(0.8704, 0.8627, 0.8548, 0.8468, 0.8388),
  t_value = c(0.185, 0.1963, 0.2078, 0.2196, 0.2314),
  Qtest = c(28.791, 31.687, 37.626, 50.723, 91.805),
  I1 = c(67.84, 65.20, 62.43, 59.53, 56.55),
  I2 = c(10.21, 12.79, 15.56, 18.51, 21.61),
  PI_lower = c(-1.7918, -1.788, -1.7851, -1.7834, -1.7828),
  PI_upper = c(1.9516, 1.9574, 1.9645, 1.9729, 1.9826)
)

data <- data %>%
  mutate(
    label_R = sprintf("R = %.1f", r),
    
    label_stats = sprintf("K=%d  N=%d\nQ=%.2f  I²=%.2f%%  I¹=%.2f%%", 
                          K, N, Qtest, I2, I1),
    label_effect = sprintf("%.2f (%.2f, %.2f)", 
                           estimate, ci_lower, ci_upper),
    label_p = sprintf("p = %.3f", p_value)
  )
ggplot(data, aes(x = estimate, y = fct_rev(label_R))) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
                 height = 0.15, size = 0.7, color = "grey") +
  geom_point(size = 3, shape = 18, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
  geom_text(aes(x = -2.69, label = label_stats), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.7, label = label_effect), size = 3.5) +           
  geom_text(aes(x = 1.5, label = label_p), size = 3.5) +                
  annotate("text", x = -2.69, y = 5.53, label = "Heterogeneity Statistics",
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 0.7, y = 5.53, label = "Effect Size (95% CI)",
           hjust = 0.5, fontface = "bold", size = 3.5) +                
  annotate("text", x = 1.5, y = 5.53, label = "p-value",
           hjust = 0.5, fontface = "bold", size = 3.5) +                
  scale_x_continuous(limits = c(-3, 2.5),
                     breaks = seq(-0.5, 1.5, 0.5),
                     expand = c(0, 0)) +
  labs(title = "MP——Sensitivity Multivariate Analysis with Different r Values",
       x = "Effect Size",
       y = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.3, size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_line(color = "grey90", size = 0.2),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

# sensitive_leave1out (Code to generate the Supplementary Material in Table 10)
leave1out_threelevel(data_cTn_SMD, r = 0.6)

# sensitive_cook(Code to generate the Supplementary Material in Table 13)
cooks <- cooks.distance(res_cTn_final)
# View outliers with Cooks > 3 * mean
outliers_cooks <- cooks %>%
  cbind(data_cTn_SMD$esid) %>%            # bind study names for reference
  subset(cooks > 3.0 * mean(cooks)) %>%  # subset outliers
  View()

data_cTn_leave_cooks <- data_cTn_SMD[-3,]  # outliers

V_cooks <- with(data_cTn_leave_cooks,
                impute_covariance_matrix(vi = vi, cluster = studyid, r = 0.6))

res_cTn_cooks <- rma.mv(yi = yi, V = V_cooks, data = data_cTn_leave_cooks,
                        random = ~ 1 | study/id,
                        method = "REML", test = "t", slab = study)

res_cTn_cooks_final <- robust(res_cTn_cooks, cluster = study, clubSandwich = TRUE)
res_cTn_cooks_final

#(Code to generate the Supplementary Material in Table 13)
# Based on reviewer's valuable selective suggestions, we did the cluster of study and then run the analysis again.
x6 <- aggregate(data_cTn_SMD, cluster = study, rho = 0.6)
res <- rma(yi, vi, data = x6)
forest(res)
summary(res)
# sensitive_residuals(Code to generate the Supplementary Material in Table 13)
resid <- residuals(res_cTn_final) %>%
  scale(center = FALSE, scale = TRUE)
# View outliers with residuals > 3 * mean
outliers_resid <- resid %>%
  cbind(data_cTn_SMD$esid) %>%
  subset(resid > 3.0 | resid < -3.0) %>%
  View()
# Non-outlier

# sensitive_hat(Code to generate the Supplementary Material in Table 13)
high_hat <- data.frame(hatvalues.rma.mv(res_cTn_final)) %>%
  rename(hat = hatvalues.rma.mv.res_cTn_final.) %>%
  filter(hat > 3 * mean(hat))
# Non-outlier


#====  Effects and sensitive of FMD ==== 

r_values <- seq(0.5, 0.9, by = 0.1)

meta_results <- list()

for(i in seq_along(r_values)){
  
  r <- r_values[i]
  
  data_FMD <- read_excel("C:/Users/Administrator/Desktop/Data and R code.xlsx", sheet = 4)
  
  # ===== Step 1: calculate SDchange =====
  data_FMD$IPCsdchange <- sqrt(data_FMD$IPCsdpre^2 + data_FMD$IPCsdpost^2 - 
                             2 * r * data_FMD$IPCsdpre * data_FMD$IPCsdpost)
  
  data_FMD$CONsdchange <- sqrt(data_FMD$CONsdpre^2 + data_FMD$CONsdpost^2 - 
                             2 * r * data_FMD$CONsdpre * data_FMD$CONsdpost)
  
  # ===== Step 2: calculate mean change =====
  data_FMD$IPCmeanchange <- data_FMD$IPCmeanpost - data_FMD$IPCmeanpre
  data_FMD$CONmeanchange <- data_FMD$CONmeanpost - data_FMD$CONmeanpre
  
  # ===== Step 3: meta analysis =====
  meta_results[[i]] <- metacont(
    IPCn,
    IPCmeanchange,
    IPCsdchange,
    CONn,
    CONmeanchange,
    CONsdchange,
    data = data_FMD,
    sm = "SMD",
    random = TRUE,
    common = TRUE,
    studlab = paste(data_FMD$study),
    label.e = "Experiment",
    label.c = "Control",
    prediction = TRUE
  )
  
  cat("\n====================\n")
  cat("Result for r =", r, "\n")
  print(summary(meta_results[[i]]))
}


results_summary <- data.frame()

for(i in seq_along(meta_results)){
  
  m <- meta_results[[i]]
  
  results_summary <- rbind(results_summary, data.frame(
    r = r_values[i],
    TE = m$TE.random,
    lower = m$lower.random,
    upper = m$upper.random,
    p = m$pval.random,
    I2 = m$I2
  ))
}

print(results_summary)


results_table <- data.frame(
  r = sprintf("%.1f", r_values),
  
  Fixed_Effect = sprintf(
    "%.2f (%.2f to %.2f)",
    sapply(meta_results, function(x) x$TE.common),
    sapply(meta_results, function(x) x$lower.common),
    sapply(meta_results, function(x) x$upper.common)
  ),
  Fixed_p = round(sapply(meta_results, function(x) x$pval.common), 3),
  
  Random_Effect = sprintf(
    "%.2f (%.2f to %.2f)",
    sapply(meta_results, function(x) x$TE.random),
    sapply(meta_results, function(x) x$lower.random),
    sapply(meta_results, function(x) x$upper.random)
  ),
  Random_p = round(sapply(meta_results, function(x) x$pval.random), 3),
  
  Q = round(sapply(meta_results, function(x) x$Q), 3),
  Q_df = sapply(meta_results, function(x) x$df.Q),
  Q_p = round(sapply(meta_results, function(x) x$pval.Q), 3),
  
  I2 = round(sapply(meta_results, function(x) summary(x)$I2), 3),
  tau2 = round(sapply(meta_results, function(x) (x$tau)^2), 3),
  
  PI = sprintf(
    "%.2f to %.2f",
    sapply(meta_results, function(x) x$lower.predict),
    sapply(meta_results, function(x) x$upper.predict)
  )
)

print(results_table)


# Code to generate the forest in Figure 6
forest(meta_results[[2]])

# FMD_r forest plot (code to generate the Supplementary in Figure 3)
data <- data.frame(
  R = c(0.5, 0.6, 0.7, 0.8, 0.9),
  K = rep(3,5),
  n = rep(3,5),
  Q = c(0.46,0.49,0.53,0.60,0.87), 
  I2 = c(0,0,0,0,0),
  p_het = c(0.7952,0.7840,0.7682,0.7404,0.6487),
  estimate = c(1.7451,1.8367,1.9531,2.1356,2.4897),
  ci_lower = c(1.2041,1.2785,1.3828,1.5453,1.8572),
  ci_upper = c(2.3040,2.3949,2.5235,2.7260,3.1221),
  p_effect = c("<0.0001", "<0.0001", "<0.0001", "<0.0001", "<0.0001")
)
data <- data %>%
  mutate(
    label_R = sprintf("R = %.1f", R),
    label_stats = sprintf("K=%d  n=%d\nQ=%.2f  I²=(%.2f)%%\np(het)=%s", 
                          K, n, Q, I2, p_het),
    label_effect = sprintf("%.2f (%.2f, %.2f)", estimate, ci_lower, ci_upper),
    label_p = sprintf("%s", p_effect)
  )
ggplot(data, aes(x = estimate, y = fct_rev(label_R))) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), 
                 height = 0.15, size = 0.7, color = "grey") +
  geom_point(size = 3, shape = 18, color = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 0.5) +
  geom_text(aes(x = 0.38, label = label_stats), hjust = 0, size = 3.5) +
  geom_text(aes(x = 3.80, label = label_effect), size = 3.5) +
  geom_text(aes(x = 4.4, label = label_p), size = 3.5) +
  annotate("text", x = 0.38, y = 5.53, label = "Heterogeneity Statistics", 
           hjust = 0, fontface = "bold", size = 3.5) +
  annotate("text", x = 4.42, y = 5.53, label = "Effect Size (95% CI)", 
           hjust = 1.4, fontface = "bold", size = 3.5) +
  annotate("text", x = 4.64, y = 5.53, label = "p-value", 
           hjust = 1.5, fontface = "bold", size = 3.5) +
  scale_x_continuous(limits = c(0.2, 4.7), 
                     breaks = seq(1,3, 0.5),
                     expand = c(0, 0)) +
  labs(title = "FMD——Sensitivity Analysis with Different r Values", 
       x = "Effect Size", 
       y = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.3, size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_line(color = "grey90", size = 0.2),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(1, 1, 1, 1, "cm")
  )

