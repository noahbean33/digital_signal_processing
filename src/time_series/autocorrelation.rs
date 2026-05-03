fn mean(data: &[f64]) -> f64 {
    let sum: f64 = data.iter().sum();
    sum / data.len() as f64
}

pub fn autocorrelation(data: &[f64], lag: usize) -> f64 {
    let n = data.len();
    if lag >= n {
        panic!("lag must be less than the length of the data series")
    }

    let mean_val = mean(data);
    let mut numerator = 0.0;
    let mut denominator = 0.0;

    for elements in 0..(n - lag) {
        // lag = 1 and n = 5 -> elements from 0 to 3: it will consider the 5th index which will be 4
        // x[i] -> (x[0]+...+x[3] - mean) * (x[i + lag] -> x[0+lag] + ... x[2+lag] - mean )
        numerator += (data[elements] - mean_val) * (data[elements + lag] - mean_val)
    }

    for elements in 0..n {
        // (x[0 .. n] - mean_val)^2
        denominator += (data[elements] - mean_val).powi(2);
    }

    numerator / denominator
}

// 1...5 -> lag = 1 -> 1,2,3,,4,5 if you lag = 2 -> 1,3,5. Note: lag less than n

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correlation() {
        let data: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let lag = 1; // Defining the lag variable referenced in the image
        let result: f64 = autocorrelation(&data, lag);
        assert!((result - 0.4).abs() < 1e-10);
    }
}