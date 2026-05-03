pub fn simple_moving_average(data: &[f64], window_size: usize) -> Vec<f64>{
    let mut result = Vec::new();
    if data.len() < window_size{
        return result;
    }

    for i in 0..=data.len() - window_size{
        let window = &data[i..i + window_size];
        let sum: f64 = window.iter().sum();
        let avg = sum / window_size as f64;
        result.push(avg)
    }

    return result;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_moving_average() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let result = simple_moving_average(&data, 3);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], 2.0);
        assert_eq!(result[1], 3.0);
        assert_eq!(result[2], 4.0);
    }

    #[test]
    fn test_sma_empty_result() {
        let data = [1.0, 2.0];
        let result = simple_moving_average(&data, 5);
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_sma_single_window() {
        let data = [10.0, 20.0, 30.0];
        let result = simple_moving_average(&data, 3);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 20.0);
    }
}
