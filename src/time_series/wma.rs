pub fn weighted_moving_average(data: &[f64], window_size: usize, weights: &[f64]) -> Vec<f64>{
    let mut result = Vec::new();
    if data.len() < window_size{
        return result;
    }

    for i in 0..=data.len() - window_size{
        let window = &data[i..i + window_size];
        let weighted_sum: f64 = window.iter().enumerate().map(|(j, &value)| value * weights[j]).sum();
        let sum_of_weights: f64 = weights.iter().sum();
        let weighted_avg = weighted_sum / sum_of_weights;
        result.push(weighted_avg);
    }

    return result;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_weighted_moving_average() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let weights = [0.1, 0.2, 0.7];
        let result = weighted_moving_average(&data, 3, &weights);
        assert_eq!(result.len(), 3);
        
        let expected_first = (1.0 * 0.1 + 2.0 * 0.2 + 3.0 * 0.7) / 1.0;
        assert!((result[0] - expected_first).abs() < 1e-10);
    }

    #[test]
    fn test_wma_equal_weights() {
        let data = [10.0, 20.0, 30.0];
        let weights = [1.0, 1.0, 1.0];
        let result = weighted_moving_average(&data, 3, &weights);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 20.0);
    }

    #[test]
    fn test_wma_empty_result() {
        let data = [1.0, 2.0];
        let weights = [0.5, 0.5, 0.5];
        let result = weighted_moving_average(&data, 3, &weights);
        assert_eq!(result.len(), 0);
    }
}
