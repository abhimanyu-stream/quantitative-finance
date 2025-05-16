package com.quant.finance.models;

import com.quant.finance.derivatives.Option;
import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.util.HashMap;
import java.util.Map;

/**
 * Abstract base class for option pricing models
 */
public abstract class OptionPricingModel {
    
    /**
     * Calculate option price
     * 
     * @param option Option to price
     * @param valuationDate Valuation date
     * @param marketData Market data for pricing
     * @return Option price
     */
    public abstract double calculatePrice(Option option, LocalDate valuationDate, Map<String, Object> marketData);
    
    /**
     * Calculate option Greeks
     * 
     * @param option Option to analyze
     * @param valuationDate Valuation date
     * @param marketData Market data for pricing
     * @return Map of Greeks
     */
    public abstract Map<String, Double> calculateGreeks(Option option, LocalDate valuationDate, Map<String, Object> marketData);
    
    /**
     * Factory method to create an appropriate pricing model based on option type and style
     * 
     * @param modelType Type of pricing model to create
     * @return OptionPricingModel instance
     */
    public static OptionPricingModel createModel(String modelType) {
        switch (modelType.toLowerCase()) {
            case "black-scholes":
                return new BlackScholesModel();
            case "binomial":
                return new BinomialTreeModel();
            case "montecarlo":
                return new MonteCarloOptionModel();
            default:
                throw new IllegalArgumentException("Unknown model type: " + modelType);
        }
    }
}

/**
 * Black-Scholes-Merton option pricing model
 */
class BlackScholesModel extends OptionPricingModel {
    
    @Override
    public double calculatePrice(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        double dividendYield = (double) marketData.getOrDefault("dividendYield", 0.0);
        
        // Calculate time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        
        if (timeToExpiry <= 0) {
            // Option is expired, return intrinsic value
            return calculateIntrinsicValue(option, spotPrice);
        }
        
        // Calculate d1 and d2
        double d1 = calculateD1(spotPrice, option.getStrikePrice(), riskFreeRate, 
                               dividendYield, volatility, timeToExpiry);
        double d2 = d1 - volatility * Math.sqrt(timeToExpiry);
        
        // Calculate option price based on type
        if (option.getOptionType() == Option.OptionType.CALL) {
            return spotPrice * Math.exp(-dividendYield * timeToExpiry) * standardNormalCDF(d1) - 
                   option.getStrikePrice() * Math.exp(-riskFreeRate * timeToExpiry) * standardNormalCDF(d2);
        } else {
            return option.getStrikePrice() * Math.exp(-riskFreeRate * timeToExpiry) * standardNormalCDF(-d2) - 
                   spotPrice * Math.exp(-dividendYield * timeToExpiry) * standardNormalCDF(-d1);
        }
    }
    
    @Override
    public Map<String, Double> calculateGreeks(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> greeks = new HashMap<>();
        
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        double dividendYield = (double) marketData.getOrDefault("dividendYield", 0.0);
        
        // Calculate time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        
        if (timeToExpiry <= 0) {
            // Option is expired, return zero Greeks
            greeks.put("delta", 0.0);
            greeks.put("gamma", 0.0);
            greeks.put("theta", 0.0);
            greeks.put("vega", 0.0);
            greeks.put("rho", 0.0);
            return greeks;
        }
        
        // Calculate d1 and d2
        double d1 = calculateD1(spotPrice, option.getStrikePrice(), riskFreeRate, 
                               dividendYield, volatility, timeToExpiry);
        double d2 = d1 - volatility * Math.sqrt(timeToExpiry);
        
        // Calculate Greek values based on option type
        if (option.getOptionType() == Option.OptionType.CALL) {
            // Delta
            greeks.put("delta", Math.exp(-dividendYield * timeToExpiry) * standardNormalCDF(d1));
            
            // Gamma (same for call and put)
            double gamma = Math.exp(-dividendYield * timeToExpiry) * standardNormalPDF(d1) / 
                          (spotPrice * volatility * Math.sqrt(timeToExpiry));
            greeks.put("gamma", gamma);
            
            // Theta
            double theta1 = -spotPrice * volatility * Math.exp(-dividendYield * timeToExpiry) * 
                           standardNormalPDF(d1) / (2 * Math.sqrt(timeToExpiry));
            double theta2 = -riskFreeRate * option.getStrikePrice() * Math.exp(-riskFreeRate * timeToExpiry) * 
                           standardNormalCDF(d2);
            double theta3 = dividendYield * spotPrice * Math.exp(-dividendYield * timeToExpiry) * 
                           standardNormalCDF(d1);
            greeks.put("theta", theta1 + theta2 + theta3);
            
            // Vega (same for call and put)
            double vega = spotPrice * Math.exp(-dividendYield * timeToExpiry) * 
                         standardNormalPDF(d1) * Math.sqrt(timeToExpiry);
            greeks.put("vega", vega);
            
            // Rho
            greeks.put("rho", option.getStrikePrice() * timeToExpiry * Math.exp(-riskFreeRate * timeToExpiry) * 
                             standardNormalCDF(d2));
        } else {
            // Delta
            greeks.put("delta", -Math.exp(-dividendYield * timeToExpiry) * standardNormalCDF(-d1));
            
            // Gamma (same for call and put)
            double gamma = Math.exp(-dividendYield * timeToExpiry) * standardNormalPDF(d1) / 
                          (spotPrice * volatility * Math.sqrt(timeToExpiry));
            greeks.put("gamma", gamma);
            
            // Theta
            double theta1 = -spotPrice * volatility * Math.exp(-dividendYield * timeToExpiry) * 
                           standardNormalPDF(d1) / (2 * Math.sqrt(timeToExpiry));
            double theta2 = riskFreeRate * option.getStrikePrice() * Math.exp(-riskFreeRate * timeToExpiry) * 
                           standardNormalCDF(-d2);
            double theta3 = -dividendYield * spotPrice * Math.exp(-dividendYield * timeToExpiry) * 
                           standardNormalCDF(-d1);
            greeks.put("theta", theta1 + theta2 + theta3);
            
            // Vega (same for call and put)
            double vega = spotPrice * Math.exp(-dividendYield * timeToExpiry) * 
                         standardNormalPDF(d1) * Math.sqrt(timeToExpiry);
            greeks.put("vega", vega);
            
            // Rho
            greeks.put("rho", -option.getStrikePrice() * timeToExpiry * Math.exp(-riskFreeRate * timeToExpiry) * 
                             standardNormalCDF(-d2));
        }
        
        return greeks;
    }
    
    /**
     * Calculate the d1 term in the Black-Scholes formula
     */
    private double calculateD1(double spotPrice, double strikePrice, double riskFreeRate, 
                              double dividendYield, double volatility, double timeToExpiry) {
        return (Math.log(spotPrice / strikePrice) + 
               (riskFreeRate - dividendYield + 0.5 * volatility * volatility) * timeToExpiry) / 
               (volatility * Math.sqrt(timeToExpiry));
    }
    
    /**
     * Calculate the intrinsic value of an option
     */
    private double calculateIntrinsicValue(Option option, double spotPrice) {
        if (option.getOptionType() == Option.OptionType.CALL) {
            return Math.max(0, spotPrice - option.getStrikePrice());
        } else {
            return Math.max(0, option.getStrikePrice() - spotPrice);
        }
    }
    
    /**
     * Calculate the standard normal cumulative distribution function
     */
    private double standardNormalCDF(double x) {
        // Use an approximation of the standard normal CDF
        // Abramowitz and Stegun approximation 7.1.26
        double b1 = 0.319381530;
        double b2 = -0.356563782;
        double b3 = 1.781477937;
        double b4 = -1.821255978;
        double b5 = 1.330274429;
        double p = 0.2316419;
        
        if (x >= 0.0) {
            double t = 1.0 / (1.0 + p * x);
            return 1.0 - standardNormalPDF(x) * 
                   (b1 * t + 
                    b2 * t * t + 
                    b3 * t * t * t + 
                    b4 * t * t * t * t + 
                    b5 * t * t * t * t * t);
        } else {
            return 1.0 - standardNormalCDF(-x);
        }
    }
    
    /**
     * Calculate the standard normal probability density function
     */
    private double standardNormalPDF(double x) {
        return Math.exp(-x * x / 2.0) / Math.sqrt(2.0 * Math.PI);
    }
}

/**
 * Binomial tree option pricing model
 */
class BinomialTreeModel extends OptionPricingModel {
    private int steps = 100; // Default number of steps
    
    public BinomialTreeModel() {
        this(100);
    }
    
    public BinomialTreeModel(int steps) {
        this.steps = steps;
    }
    
    @Override
    public double calculatePrice(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        double dividendYield = (double) marketData.getOrDefault("dividendYield", 0.0);
        
        // Calculate time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        
        if (timeToExpiry <= 0) {
            // Option is expired, return intrinsic value
            return calculateIntrinsicValue(option, spotPrice);
        }
        
        // Calculate parameters for the binomial model
        double dt = timeToExpiry / steps;
        double u = Math.exp(volatility * Math.sqrt(dt));
        double d = 1.0 / u;
        double p = (Math.exp((riskFreeRate - dividendYield) * dt) - d) / (u - d);
        double discount = Math.exp(-riskFreeRate * dt);
        
        // Initialize arrays for stock prices and option values at expiration
        double[] stockPrices = new double[steps + 1];
        double[] optionValues = new double[steps + 1];
        
        // Calculate stock prices at expiration
        for (int i = 0; i <= steps; i++) {
            stockPrices[i] = spotPrice * Math.pow(u, steps - i) * Math.pow(d, i);
            optionValues[i] = calculateIntrinsicValue(option, stockPrices[i]);
        }
        
        // Work backwards through the tree
        for (int j = steps - 1; j >= 0; j--) {
            for (int i = 0; i <= j; i++) {
                // Calculate stock price at node (j,i)
                stockPrices[i] = spotPrice * Math.pow(u, j - i) * Math.pow(d, i);
                
                // Calculate option value at node (j,i)
                double exerciseValue = calculateIntrinsicValue(option, stockPrices[i]);
                double continuationValue = discount * (p * optionValues[i] + (1 - p) * optionValues[i + 1]);
                
                // For American options, we need to check if early exercise is beneficial
                if (option.getExerciseStyle() == Option.ExerciseStyle.AMERICAN) {
                    optionValues[i] = Math.max(exerciseValue, continuationValue);
                } else {
                    optionValues[i] = continuationValue;
                }
            }
        }
        
        // Value at root node is the option price
        return optionValues[0];
    }
    
    @Override
    public Map<String, Double> calculateGreeks(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> greeks = new HashMap<>();
        
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        
        // Calculate time to expiry in years
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        
        if (timeToExpiry <= 0) {
            // Option is expired, return zero Greeks
            greeks.put("delta", 0.0);
            greeks.put("gamma", 0.0);
            greeks.put("theta", 0.0);
            greeks.put("vega", 0.0);
            greeks.put("rho", 0.0);
            return greeks;
        }
        
        // Calculate delta using finite difference
        double deltaShift = 0.01 * spotPrice;
        Map<String, Object> upMarketData = new HashMap<>(marketData);
        upMarketData.put("spotPrice", spotPrice + deltaShift);
        
        Map<String, Object> downMarketData = new HashMap<>(marketData);
        downMarketData.put("spotPrice", spotPrice - deltaShift);
        
        double priceUp = calculatePrice(option, valuationDate, upMarketData);
        double priceDown = calculatePrice(option, valuationDate, downMarketData);
        
        // Delta
        double delta = (priceUp - priceDown) / (2 * deltaShift);
        greeks.put("delta", delta);
        
        // Gamma
        double price = calculatePrice(option, valuationDate, marketData);
        double gamma = (priceUp - 2 * price + priceDown) / (deltaShift * deltaShift);
        greeks.put("gamma", gamma);
        
        // Theta (1 day change)
        LocalDate nextDay = valuationDate.plusDays(1);
        double priceNextDay = calculatePrice(option, nextDay, marketData);
        double theta = (priceNextDay - price) / (1.0 / 365.0);
        greeks.put("theta", theta);
        
        // Vega (1% volatility change)
        double volShift = 0.01;
        Map<String, Object> upVolMarketData = new HashMap<>(marketData);
        upVolMarketData.put("volatility", volatility + volShift);
        
        double priceUpVol = calculatePrice(option, valuationDate, upVolMarketData);
        double vega = (priceUpVol - price) / volShift;
        greeks.put("vega", vega);
        
        // Rho (1% rate change)
        double rateShift = 0.0001;
        Map<String, Object> upRateMarketData = new HashMap<>(marketData);
        upRateMarketData.put("riskFreeRate", riskFreeRate + rateShift);
        
        double priceUpRate = calculatePrice(option, valuationDate, upRateMarketData);
        double rho = (priceUpRate - price) / rateShift;
        greeks.put("rho", rho);
        
        return greeks;
    }
    
    /**
     * Calculate the intrinsic value of an option
     */
    private double calculateIntrinsicValue(Option option, double spotPrice) {
        if (option.getOptionType() == Option.OptionType.CALL) {
            return Math.max(0, spotPrice - option.getStrikePrice());
        } else {
            return Math.max(0, option.getStrikePrice() - spotPrice);
        }
    }
}

/**
 * Monte Carlo simulation-based option pricing model
 */
class MonteCarloOptionModel extends OptionPricingModel {
    private int numSimulations = 10000;
    private int stepsPerYear = 252;
    private boolean antithetic = true;
    private boolean controlVariate = true;
    
    public MonteCarloOptionModel() {
        this(10000, 252, true, true);
    }
    
    public MonteCarloOptionModel(int numSimulations, int stepsPerYear, 
                                boolean antithetic, boolean controlVariate) {
        this.numSimulations = numSimulations;
        this.stepsPerYear = stepsPerYear;
        this.antithetic = antithetic;
        this.controlVariate = controlVariate;
    }
    
    @Override
    public double calculatePrice(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        // Use MonteCarloPricer to implement this
        MonteCarloPricer pricer = new MonteCarloPricer(
            numSimulations,          // Number of simulations
            stepsPerYear,            // Number of time steps
            antithetic,              // Use antithetic variance reduction
            controlVariate,          // Use parallel execution (repurposing controlVariate flag)
            0                        // Random seed (0 for random seed)
        );
        
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        double dividendYield = (double) marketData.getOrDefault("dividendYield", 0.0);
        
        Map<String, Double> results = pricer.priceOption(
            option, valuationDate, spotPrice, volatility, riskFreeRate, dividendYield
        );
        
        return results.get("price");
    }
    
    @Override
    public Map<String, Double> calculateGreeks(Option option, LocalDate valuationDate, Map<String, Object> marketData) {
        Map<String, Double> greeks = new HashMap<>();
        
        // Extract market data
        double spotPrice = (double) marketData.getOrDefault("spotPrice", 0.0);
        double riskFreeRate = (double) marketData.getOrDefault("riskFreeRate", 0.0);
        double volatility = (double) marketData.getOrDefault("volatility", 0.0);
        double dividendYield = (double) marketData.getOrDefault("dividendYield", 0.0);
        
        // Time to expiry
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        
        if (timeToExpiry <= 0) {
            // Option is expired, return zero Greeks
            greeks.put("delta", 0.0);
            greeks.put("gamma", 0.0);
            greeks.put("theta", 0.0);
            greeks.put("vega", 0.0);
            greeks.put("rho", 0.0);
            return greeks;
        }
        
        // Calculate delta using finite difference
        double deltaShift = 0.01 * spotPrice;
        Map<String, Object> upMarketData = new HashMap<>(marketData);
        upMarketData.put("spotPrice", spotPrice + deltaShift);
        
        Map<String, Object> downMarketData = new HashMap<>(marketData);
        downMarketData.put("spotPrice", spotPrice - deltaShift);
        
        double priceUp = calculatePrice(option, valuationDate, upMarketData);
        double priceDown = calculatePrice(option, valuationDate, downMarketData);
        
        // Delta
        double delta = (priceUp - priceDown) / (2 * deltaShift);
        greeks.put("delta", delta);
        
        // Current price
        double price = calculatePrice(option, valuationDate, marketData);
        
        // Gamma
        double gamma = (priceUp - 2 * price + priceDown) / (deltaShift * deltaShift);
        greeks.put("gamma", gamma);
        
        // Theta (1 day change)
        LocalDate nextDay = valuationDate.plusDays(1);
        double priceNextDay = calculatePrice(option, nextDay, marketData);
        double theta = (priceNextDay - price) / (1.0 / 365.0);
        greeks.put("theta", theta);
        
        // Vega (1% volatility change)
        double volShift = 0.01;
        Map<String, Object> upVolMarketData = new HashMap<>(marketData);
        upVolMarketData.put("volatility", volatility + volShift);
        
        double priceUpVol = calculatePrice(option, valuationDate, upVolMarketData);
        double vega = (priceUpVol - price) / volShift;
        greeks.put("vega", vega);
        
        // Rho (1% rate change)
        double rateShift = 0.0001;
        Map<String, Object> upRateMarketData = new HashMap<>(marketData);
        upRateMarketData.put("riskFreeRate", riskFreeRate + rateShift);
        
        double priceUpRate = calculatePrice(option, valuationDate, upRateMarketData);
        double rho = (priceUpRate - price) / rateShift;
        greeks.put("rho", rho);
        
        return greeks;
    }
} 