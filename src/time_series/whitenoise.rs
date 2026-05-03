fn mean(data: &[f64]) -> f64{
    let sum: f64 = data.iter().sum();
    sum / data.len() as f64
}

fn variance(data: &[f64], mean: f64) -> f64{
    let variance: f64 = data.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / data.len() as f64;
    variance
}

fn autocorrelation(data: &[f64], lag: usize) -> f64{
    let n = data.len();
    if lag >= n{
        panic!("Lag must be less than the length of the data array")
    }

    let mean_val = mean(data);
    let mut numerator = 0.0;
    let mut denominator = 0.0;

    for elements in 0..(n - lag){
        numerator += (data[elements] - mean_val) * (data[elements + lag] - mean_val);
    }

    for elements in 0..n{
        denominator += (data[elements] - mean_val).powi(2);
    }

    numerator / denominator
}

pub fn check_white_noise(data: &[f64], lag: usize) -> bool{
    let mean_val = mean(data);
    if mean_val.abs() > 0.05{
        return false;
    }

    let var_val = variance(data, mean_val);
    if var_val < 0.05 {
        return false;
    }

    let autocorr = autocorrelation(data, lag);
    if autocorr.abs() > 0.05{
        return false;
    }

    return true;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_white_noise_detection() {
        // Data must satisfy: |mean| <= 0.05, variance >= 0.05, |autocorr(lag)| <= 0.05
        // Block-of-2 pattern [0.3,0.3,-0.3,-0.3,...] gives mean=0, var=0.09, low autocorr
        let mut data = Vec::new();
        for _ in 0..6 {
            data.extend_from_slice(&[0.3, 0.3, -0.3, -0.3]);
        }
        let result = check_white_noise(&data, 1);
        assert!(result);
    }

    #[test]
    fn test_non_white_noise_high_mean() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let result = check_white_noise(&data, 1);
        assert_eq!(result, false);
    }

    #[test]
    fn test_variance() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let m = mean(&data);
        let var = variance(&data, m);
        assert!((var - 2.0).abs() < 1e-10);
    }
}
