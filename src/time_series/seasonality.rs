fn mean(data: &[f64]) -> f64{
    let sum: f64 = data.iter().sum();
    sum / data.len() as f64
}

fn seasonal_difference(data: &[f64], period: usize) -> Vec<f64>{
    let mut difference = Vec::new();
    for elements in period..data.len(){
        difference.push(data[elements] - data[elements - period]);
    }
    return difference
}

fn std_dev(data: &[f64], mean: f64) -> f64{
    let std_dev: f64 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / data.len() as f64;
    std_dev.sqrt()
}

pub fn check_seasonality(data: &[f64], period: usize) -> bool{
    let differenced = seasonal_difference(data, period);

    let mean_val = mean(&differenced);
    let std_dev_val = std_dev(&differenced, mean_val);

    return std_dev_val > 0.05 * mean_val;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seasonal_data() {
        let data = [10.0, 15.0, 10.0, 15.0, 10.0, 15.0];
        let result = check_seasonality(&data, 2);
        assert_eq!(result, false);
    }

    #[test]
    fn test_non_seasonal_data() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let result = check_seasonality(&data, 2);
        assert_eq!(result, false);
    }

    #[test]
    fn test_seasonal_difference() {
        let data = [10.0, 20.0, 30.0, 40.0, 50.0];
        let diff = seasonal_difference(&data, 2);
        assert_eq!(diff.len(), 3);
        assert_eq!(diff[0], 20.0);
        assert_eq!(diff[1], 20.0);
    }
}
