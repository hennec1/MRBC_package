#' Mean red blood cell age 
#' 
#' Estimates the mean red blood cell age from an estimated survival curve
#' @param time A vector of time points at which the relative APE is collected
#' @param APE A vector of relative APE for each time measurement
#' @return The mean red blood cell age over the time measured
#' @examples
#' time <- c(0,1,7,14,21,28,42,56,70,84,98,105,112,119,126,140,161,182); 
#' APE <- c(NA,2.3, 56.62, 87.15, 86.56, 90.60, 98.46, 89.64, 100, 81.19, 86.95,83.69, 83.88, 63.53, 64.30, 33.97, 22.84, 18.43);
#' mrbc <- meanRBC.function(time, APE)
#' @export
meanRBC.function <- function(time, APE) {
  
  # Step 0: Create data frame
    sub1 <- data.frame(time, APE)
  
  # Step 1:Define the starting point for all calculations: the time at which the %APE is equal to 50% of the maximum %APE value 
    # Do this because the %APE increases immediately after labeling, so we need to get a consistent starting point
    for(i in 1:nrow(sub1)){
      sub1$start_time[i] <- ifelse(is.na(sub1$time[i]), NA, 
                                 ifelse(sub1$APE[i] < 50 & sub1$APE[i+1] > 50,
                                        (sub1$time[i+1]-sub1$time[i])*((50-sub1$APE[i])/(sub1$APE[i+1]-sub1$APE[i])) 
                                        + sub1$time[i], 0))
    }
  
    # To account for any non-monotonically decreasing points after 100% APE
      for(i in 1:nrow(sub1)){
        if(any(is.na(sub1$APE[(i-i+1):(i-1)]))){
          sub1$start_time[i] <- sub1$start_time[i]
        } else if(any(sub1$APE[(i-i+1):(i-1)] == 100)) {
          sub1$start_time[i] <- 0
        } else {
          sub1$start_time[i] <- sub1$start_time[i]
        }
      }
  
    # Reset the time based on the new start time. 
  
      for(i in 1:nrow(sub1)){
        sub1$new_start_time[i] <- ifelse(is.na(sub1$time[i]), NA,
                                     ifelse((sub1$time[i]-max(sub1$start_time, na.rm = T))<0, 0, 
                                            (sub1$time[i]-max(sub1$start_time, na.rm = T))))
    
      }
  
  # Step 2: Set the early time points %APE to 100% 
    # For the time points that happen before the max %APE, set the %APE to 100% in order for the curve to be monotonically decreasing. 
    for(i in 1:nrow(sub1)){
      if(i < nrow(sub1)){
        sub1$reset_APE[i] <- ifelse(is.na(sub1$time[i]), NA,
                                  ifelse(max(sub1$APE[i+1:nrow(sub1)], na.rm = T) == 100, 100, 
                                         ifelse(sub1$new_start_time[i] == 0, 100, 
                                                ifelse(sub1$APE[i] < sub1$APE[i+1], 
                                                       (sub1$APE[i-1]-((sub1$APE[i-1]-sub1$APE[i+1])*
                                                                         ((sub1$new_start_time[i]-sub1$new_start_time[i-1])/(sub1$new_start_time[i+1]-sub1$new_start_time[i-1])))),
                                                       sub1$APE[i]))))
      } else {
        sub1$reset_APE[i] <- ifelse(is.na(sub1$time[i]), NA,
                                  ifelse(sub1$APE[i] == 100, 100, 
                                         ifelse(sub1$new_start_time[i] == 0, 100, 
                                                ifelse(sub1$APE[i] < 0, 
                                                       (sub1$APE[i-1]-((sub1$APE[i-1]-0)*
                                                                         ((sub1$new_start_time[i]-sub1$new_start_time[i-1])/(0-sub1$new_start_time[i-1])))),
                                                       sub1$APE[i]))))
      }
    }
  
  
  
  # Step 3: Set the rate of change between two consecutive time points
    for(i in 1:nrow(sub1)){
      sub1$rate_change[i] <- ifelse(is.na(sub1$time[i]), NA, 
                                  ifelse((sub1$new_start_time[i] == 0 & sub1$new_start_time[i-1] == 0), 0, 
                                         (sub1$reset_APE[i-1]-sub1$reset_APE[i])/(sub1$new_start_time[i]-sub1$new_start_time[i-1])))
    }
  
  
  
  # Step 4: Define the final end time as the time when the rate of change of %APE has decreased to less than 0.5% per day 
    # When new %APE is less than 40
    for(i in 1:nrow(sub1)){
      if(i == 1){
        sub1$end_point[i] <- 0
      } else {
        sub1$end_point[i] <- ifelse(sub1$end_point[i-1] != 0, NA, 
                                  ifelse(sub1$reset_APE[i] < 40 & sub1$rate_change[i] < 0.5, 
                                         (((1-(.5-sub1$rate_change[i])/(sub1$rate_change[i-1]-sub1$rate_change[i])))*
                                            (sub1$new_start_time[i]-sub1$new_start_time[i-1])+sub1$new_start_time[i-1]),0))
      }
    }
  
  
  # Step 5: For any time points after this final end time -> change these to the final end time
    for(i in 1:nrow(sub1)){
      sub1$curve_time[i] <- ifelse(is.na(sub1$end_point[i]), NA, 
                                 ifelse(sub1$end_point[i] == 0, sub1$new_start_time[i], sub1$end_point[i]))
    }
  
  
  # Step 6: Define points for curves unadjusted for ending percent 
    # This adjusts any %APE that occurs at time points after the new end time
    for(i in 1:nrow(sub1)){
      if(i == 1){
        sub1$curve_pt_unaj[i] <- 100
      } else {
        sub1$curve_pt_unaj[i] <- ifelse(is.na(sub1$rate_change[i]), NA, 
                                      ifelse(sub1$curve_time[i] > max(sub1$end_point, na.rm = T), NA, 
                                             ifelse(sub1$new_start_time[i] != sub1$new_start_time[i-1],
                                                    (sub1$reset_APE[i-1]-((sub1$reset_APE[i-1]-sub1$reset_APE[i])*
                                                                            (sub1$curve_time[i]-sub1$curve_time[i-1])/(sub1$new_start_time[i]-sub1$new_start_time[i-1]))),
                                                    100)))
      }
    }
  
  # Step 7: Adjust the APE for ending percent
    # Correct the %APE using the last %APE to avoid overestimating the survival curve 
    for(i in 1:nrow(sub1)){
      sub1$adj_APE[i] <- ifelse(is.na(sub1$end_point[i]), NA, 
                              ifelse(sub1$curve_pt_unaj[i]-(min(sub1$curve_pt_unaj, na.rm = T)/max(sub1$curve_time, na.rm = T)*sub1$curve_time[i]) > 100, 
                                     100, sub1$curve_pt_unaj[i]-(min(sub1$curve_pt_unaj, na.rm = T)/max(sub1$curve_time, na.rm = T)*sub1$curve_time[i])))
    }
    
    
  
  # Step 8: Calculate daily survival for 160 days
    # Fit the survival curve with a 5th order polynomial 
    mod1 <- lm(formula = adj_APE ~ poly(curve_time,5, raw = TRUE), data = sub1)
  
    surv_day <- data.frame(x5 = c(NA, mod1$coefficients[6], NA, NA, NA), x4 = c(NA, mod1$coefficients[5], NA, NA, NA), 
                           x3 = c(NA, mod1$coefficients[4], NA, NA, NA), x2 = c(NA, mod1$coefficients[3], NA, NA, NA), 
                           x = c(NA, mod1$coefficients[2], NA, NA, NA), cf = c(NA, mod1$coefficients[1], NA, NA, NA),
                           row.names = c("Day","%APE","L(x)", "c(x)","c(x)2"))
    
    # Survival for each day 
      for(i in c((0+7):(160+7))){
        surv_day[i] <- c(i-7, NA, NA, NA, NA)
        temp <- (surv_day$x5[2]*surv_day[1,i]^5 + surv_day$x4[2]*surv_day[1,i]^4 +
                 surv_day$x3[2]*surv_day[1,i]^3 + surv_day$x2[2]*surv_day[1,i]^2 +
                 surv_day$x[2]*surv_day[1,i]^1 + surv_day$cf[2])
        surv_day[2,i] <- ifelse(temp < 0 | surv_day[2,(i-2)] == 0, 0, ifelse(temp > surv_day[2,c(i-1)], surv_day[2,c(i-1)], 
                                                                           ifelse(temp < 0 | surv_day[2,c(i-1)] == 0, 0, temp)))
        cname <-paste0("Day_", c(i-7))
      
        names(surv_day)[names(surv_day) == paste0("V", i)] <- cname
      }
    
    # Survival: 
      for(i in c((0+7):(160+7))){
        surv_day[3,i] <- surv_day[2,i]/100
      }
    
    # Sum of survival: 
      sum_surv <- rowSums(surv_day[3,], na.rm = TRUE)
    
    # Compute c(x) for each day: survival/sum survival
      for(i in c((0+7):(160+7))){
        surv_day[4,i] <- surv_day[3,i]/sum_surv
      }
      
      for(i in c((0+7):(160+7))) {
        surv_day[5,i] <- surv_day[4,i]*(2*surv_day[1,i]+1)/2
      }
      
      mean_rbc <- rowSums(surv_day[5,], na.rm = TRUE)
      
      unname(mean_rbc)
  
}




