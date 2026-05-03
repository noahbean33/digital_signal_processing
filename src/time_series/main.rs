use time_series_rust::{autocorrelation, stationarity, trend, seasonality, whitenoise, sma, wma, ses, arma};

fn main(){
    println!("=== Time Series Analysis Library Demo ===\n");
    
    let data: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    
    println!("1. Autocorrelation:");
    let lag: usize = 1;
    let result: f64 = autocorrelation::autocorrelation(&data, lag);
    println!("   Autocorrelation at lag {}: {}\n", lag, result);

    println!("2. Stationarity Check:");
    let window_size: usize = 2;
    let is_stationary = stationarity::check_stationarity(&data, window_size);
    println!("   Is series stationary: {}\n", is_stationary);

    println!("3. Trend Detection:");
    let trend_result = trend::detect_trend(&data);
    println!("   Trend: {}\n", trend_result);

    println!("4. Simple Moving Average:");
    let sma_result = sma::simple_moving_average(&data, 3);
    println!("   SMA (window=3): {:?}\n", sma_result);

    println!("5. Weighted Moving Average:");
    let weights = vec![0.1, 0.3, 0.6];
    let wma_result = wma::weighted_moving_average(&data, 3, &weights);
    println!("   WMA (window=3): {:?}\n", wma_result);

    println!("6. Exponential Smoothing:");
    let series = vec![10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
    let model = ses::ExponentialSmoothing{
        alpha: 0.5,
        beta: None,
    };
    let ses_result = model.single_exponential_smoothing(&series);
    println!("   SES result: {:?}\n", ses_result);

    println!("7. Seasonality Check:");
    let seasonal_data = [10.0, 15.0, 10.0, 15.0, 10.0, 15.0];
    let has_seasonality = seasonality::check_seasonality(&seasonal_data, 2);
    println!("   Has seasonality: {}\n", has_seasonality);

    println!("8. White Noise Check:");
    let noise_data = [0.01, -0.02, 0.03, -0.01, 0.02, -0.03, 0.01];
    let is_white_noise = whitenoise::check_white_noise(&noise_data, 1);
    println!("   Is white noise: {}\n", is_white_noise);

    println!("9. ARMA Model:");
    let ar_coeffs = vec![0.5, 0.3];
    let ma_coeffs = vec![0.2];
    let arma_result = arma::arma(&ar_coeffs, &ma_coeffs, 0.0, 10);
    println!("   ARMA series (first 5): {:?}\n", &arma_result[..5]);

    println!("=== Demo Complete ===");
}