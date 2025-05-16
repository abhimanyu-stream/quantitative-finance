package com.quant.finance.models;

import com.quant.finance.derivatives.Option;
import com.quant.finance.utils.DateUtils;

import java.time.LocalDate;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Monte Carlo simulation model for option pricing
 */
public class MonteCarloPricer {
    private final int numSimulations;
    private final int timeSteps;
    private final boolean antithetic;
    private final boolean parallel;
    private final Random random;
    
    /**
     * Constructor for the Monte Carlo pricer
     * 
     * @param numSimulations Number of simulations to run
     * @param timeSteps Number of time steps in each simulation
     * @param antithetic Whether to use antithetic variance reduction
     * @param parallel Whether to run simulations in parallel
     * @param seed Random seed for reproducibility (0 for random seed)
     */
    public MonteCarloPricer(int numSimulations, int timeSteps, boolean antithetic, boolean parallel, long seed) {
        this.numSimulations = numSimulations;
        this.timeSteps = timeSteps;
        this.antithetic = antithetic;
        this.parallel = parallel;
        this.random = seed == 0 ? new Random() : new Random(seed);
    }
    
    /**
     * Price an option using Monte Carlo simulation
     * 
     * @param option The option to price
     * @param valuationDate The date on which to value the option
     * @param spotPrice Current price of the underlying
     * @param volatility Volatility of the underlying
     * @param riskFreeRate Risk-free interest rate
     * @param dividendYield Dividend yield of the underlying
     * @return Map containing price and error estimates
     */
    public Map<String, Double> priceOption(Option option, LocalDate valuationDate, 
                                          double spotPrice, double volatility, 
                                          double riskFreeRate, double dividendYield) {
        // Calculate time to expiry
        double timeToExpiry = DateUtils.yearFraction(valuationDate, option.getExpiryDate());
        if (timeToExpiry <= 0) {
            throw new IllegalArgumentException("Option has already expired");
        }
        
        // Determine the effective number of simulations
        int effectiveSimulations = antithetic ? numSimulations / 2 : numSimulations;
        
        // Run the simulations
        List<Double> payoffs = parallel ? 
                runParallelSimulations(option, spotPrice, volatility, riskFreeRate, dividendYield, timeToExpiry, effectiveSimulations) :
                runSequentialSimulations(option, spotPrice, volatility, riskFreeRate, dividendYield, timeToExpiry, effectiveSimulations);
        
        // Calculate price and error estimates
        double dt = timeToExpiry / timeSteps;
        double discountFactor = Math.exp(-riskFreeRate * timeToExpiry);
        
        double sum = 0.0;
        double sumSquared = 0.0;
        
        for (double payoff : payoffs) {
            double discountedPayoff = payoff * discountFactor;
            sum += discountedPayoff;
            sumSquared += discountedPayoff * discountedPayoff;
        }
        
        double price = sum / payoffs.size();
        double variance = (sumSquared - sum * sum / payoffs.size()) / (payoffs.size() - 1);
        double standardError = Math.sqrt(variance / payoffs.size());
        
        Map<String, Double> result = new HashMap<>();
        result.put("price", price);
        result.put("standardError", standardError);
        result.put("confidenceLower95", price - 1.96 * standardError);
        result.put("confidenceUpper95", price + 1.96 * standardError);
        
        return result;
    }
    
    private List<Double> runSequentialSimulations(Option option, double spotPrice, double volatility,
                                                 double riskFreeRate, double dividendYield, 
                                                 double timeToExpiry, int effectiveSimulations) {
        List<Double> payoffs = new ArrayList<>(numSimulations);
        double dt = timeToExpiry / timeSteps;
        double drift = (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * dt;
        double diffusion = volatility * Math.sqrt(dt);
        
        for (int i = 0; i < effectiveSimulations; i++) {
            double[] paths = simulatePaths(spotPrice, drift, diffusion);
            
            // Calculate payoff based on final price
            double payoff = calculatePayoff(option, paths[paths.length - 1]);
            payoffs.add(payoff);
            
            // If using antithetic variance reduction, add the antithetic path
            if (antithetic) {
                double[] antitheticPaths = simulateAntitheticPaths(spotPrice, drift, diffusion);
                double antitheticPayoff = calculatePayoff(option, antitheticPaths[antitheticPaths.length - 1]);
                payoffs.add(antitheticPayoff);
            }
        }
        
        return payoffs;
    }
    
    private List<Double> runParallelSimulations(Option option, double spotPrice, double volatility,
                                              double riskFreeRate, double dividendYield, 
                                              double timeToExpiry, int effectiveSimulations) {
        double dt = timeToExpiry / timeSteps;
        double drift = (riskFreeRate - dividendYield - 0.5 * volatility * volatility) * dt;
        double diffusion = volatility * Math.sqrt(dt);
        
        ExecutorService executor = Executors.newFixedThreadPool(
                Math.min(Runtime.getRuntime().availableProcessors(), effectiveSimulations));
        
        try {
            List<CompletableFuture<List<Double>>> futures = IntStream.range(0, effectiveSimulations)
                .mapToObj(i -> CompletableFuture.supplyAsync(() -> {
                    List<Double> results = new ArrayList<>();
                    
                    // Generate random seed for this thread
                    Random threadRandom = new Random(random.nextLong());
                    
                    // Simulate paths using thread-local random generator
                    double[] paths = simulatePaths(spotPrice, drift, diffusion, threadRandom);
                    double payoff = calculatePayoff(option, paths[paths.length - 1]);
                    results.add(payoff);
                    
                    // If using antithetic variance reduction
                    if (antithetic) {
                        double[] antitheticPaths = simulateAntitheticPaths(spotPrice, drift, diffusion, threadRandom);
                        double antitheticPayoff = calculatePayoff(option, antitheticPaths[antitheticPaths.length - 1]);
                        results.add(antitheticPayoff);
                    }
                    
                    return results;
                }, executor))
                .collect(Collectors.toList());
            
            return futures.stream()
                    .map(CompletableFuture::join)
                    .flatMap(List::stream)
                    .collect(Collectors.toList());
        } finally {
            executor.shutdown();
        }
    }
    
    private double[] simulatePaths(double spotPrice, double drift, double diffusion) {
        return simulatePaths(spotPrice, drift, diffusion, random);
    }
    
    private double[] simulatePaths(double spotPrice, double drift, double diffusion, Random random) {
        double[] prices = new double[timeSteps + 1];
        prices[0] = spotPrice;
        
        for (int t = 0; t < timeSteps; t++) {
            double z = random.nextGaussian();
            prices[t + 1] = prices[t] * Math.exp(drift + diffusion * z);
        }
        
        return prices;
    }
    
    private double[] simulateAntitheticPaths(double spotPrice, double drift, double diffusion) {
        return simulateAntitheticPaths(spotPrice, drift, diffusion, random);
    }
    
    private double[] simulateAntitheticPaths(double spotPrice, double drift, double diffusion, Random random) {
        double[] prices = new double[timeSteps + 1];
        prices[0] = spotPrice;
        
        for (int t = 0; t < timeSteps; t++) {
            double z = -random.nextGaussian(); // Note the negative sign for antithetic
            prices[t + 1] = prices[t] * Math.exp(drift + diffusion * z);
        }
        
        return prices;
    }
    
    private double calculatePayoff(Option option, double finalPrice) {
        switch (option.getOptionType()) {
            case CALL:
                return Math.max(0, finalPrice - option.getStrikePrice());
            case PUT:
                return Math.max(0, option.getStrikePrice() - finalPrice);
            default:
                throw new IllegalArgumentException("Unsupported option type");
        }
    }
    
    /**
     * Run a sensitivity analysis by varying a parameter
     * 
     * @param option The option to analyze
     * @param valuationDate The valuation date
     * @param spotPrice Base spot price
     * @param volatility Base volatility
     * @param riskFreeRate Base risk-free rate
     * @param dividendYield Base dividend yield
     * @param parameterName Parameter to vary ("spotPrice", "volatility", "riskFreeRate", or "dividendYield")
     * @param range Range of parameter values to test
     * @return Map of parameter values to option prices
     */
    public Map<Double, Double> runSensitivityAnalysis(Option option, LocalDate valuationDate,
                                                     double spotPrice, double volatility,
                                                     double riskFreeRate, double dividendYield,
                                                     String parameterName, List<Double> range) {
        Map<Double, Double> results = new HashMap<>();
        
        for (double value : range) {
            double currentSpotPrice = parameterName.equals("spotPrice") ? value : spotPrice;
            double currentVolatility = parameterName.equals("volatility") ? value : volatility;
            double currentRiskFreeRate = parameterName.equals("riskFreeRate") ? value : riskFreeRate;
            double currentDividendYield = parameterName.equals("dividendYield") ? value : dividendYield;
            
            Map<String, Double> priceResult = priceOption(option, valuationDate, currentSpotPrice,
                                                         currentVolatility, currentRiskFreeRate, currentDividendYield);
            
            results.put(value, priceResult.get("price"));
        }
        
        return results;
    }
} 