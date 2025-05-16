package com.quant.finance.utils;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;

/**
 * Utility class for date-related calculations commonly used in finance
 */
public class DateUtils {

    /**
     * Calculate the fraction of a year between two dates using Actual/365 day count convention
     * 
     * @param startDate the start date
     * @param endDate the end date
     * @return fraction of year between start and end dates
     */
    public static double yearFraction(LocalDate startDate, LocalDate endDate) {
        // Using Actual/365 day count convention
        long daysBetween = ChronoUnit.DAYS.between(startDate, endDate);
        return daysBetween / 365.0;
    }
    
    /**
     * Calculate the fraction of a year using the Actual/360 day count convention
     * 
     * @param startDate the start date
     * @param endDate the end date
     * @return fraction of year using Actual/360 convention
     */
    public static double yearFractionAct360(LocalDate startDate, LocalDate endDate) {
        long daysBetween = ChronoUnit.DAYS.between(startDate, endDate);
        return daysBetween / 360.0;
    }
    
    /**
     * Calculate the fraction of a year using the 30/360 day count convention
     * 
     * @param startDate the start date
     * @param endDate the end date
     * @return fraction of year using 30/360 convention
     */
    public static double yearFraction30360(LocalDate startDate, LocalDate endDate) {
        int d1 = startDate.getDayOfMonth();
        int m1 = startDate.getMonthValue();
        int y1 = startDate.getYear();
        
        int d2 = endDate.getDayOfMonth();
        int m2 = endDate.getMonthValue();
        int y2 = endDate.getYear();
        
        // Adjust day counts according to 30/360 convention
        if (d1 == 31) d1 = 30;
        if (d2 == 31 && d1 == 30) d2 = 30;
        
        // Calculate 30/360 day count
        return (360 * (y2 - y1) + 30 * (m2 - m1) + (d2 - d1)) / 360.0;
    }
    
    /**
     * Calculate the number of business days between two dates
     * This is a simplified version that only excludes weekends
     * 
     * @param startDate the start date (inclusive)
     * @param endDate the end date (inclusive)
     * @return number of business days
     */
    public static int businessDaysBetween(LocalDate startDate, LocalDate endDate) {
        if (startDate.isAfter(endDate)) {
            throw new IllegalArgumentException("Start date cannot be after end date");
        }
        
        int businessDays = 0;
        LocalDate currentDate = startDate;
        
        while (!currentDate.isAfter(endDate)) {
            int dayOfWeek = currentDate.getDayOfWeek().getValue();
            // If it's not a weekend (Saturday=6, Sunday=7)
            if (dayOfWeek < 6) {
                businessDays++;
            }
            currentDate = currentDate.plusDays(1);
        }
        
        return businessDays;
    }
    
    /**
     * Add a specified number of business days to a date
     * 
     * @param date the starting date
     * @param businessDays number of business days to add
     * @return the resulting date
     */
    public static LocalDate addBusinessDays(LocalDate date, int businessDays) {
        LocalDate result = date;
        int addedBusinessDays = 0;
        
        while (addedBusinessDays < businessDays) {
            result = result.plusDays(1);
            if (result.getDayOfWeek().getValue() < 6) { // not a weekend
                addedBusinessDays++;
            }
        }
        
        return result;
    }
} 