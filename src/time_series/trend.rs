fn mean(data: &[f64]) -> f64{
    let sum: f64 = data.iter().sum();
    sum / data.len() as f64
}

fn linear_regression(data: &[f64]) -> (f64, f64){
    let n = data.len();
    let x: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let x_mean = mean(&x);
    let y_mean = mean(data);

    let numerator: f64 = x.iter().zip(data.iter()).map(|(xi, yi)| (xi - x_mean) * (yi - y_mean)).sum();
    let denominator: f64 = x.iter().map(|xi| (xi - x_mean).powi(2)).sum();

    let slope = numerator / denominator;
    let intercept = y_mean - slope * x_mean;
    return(slope, intercept)
}

pub fn detect_trend(data: &[f64]) -> String{
    let (slope, _) = linear_regression(data);

    if slope > 0.0{
        return "Upward trend".to_string()
    } else if slope < 0.0{
        return "Downward trend".to_string()
    } else{
        return "No trend".to_string()
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    
    #[test]
    fn test_trend(){
        let data = [1.15, 2.75, 3.0, 3.25, 4.0];
        let result = detect_trend(&data);
        assert_eq!(result, "Upward trend");
    }

    #[test]
    fn test_downward_trend(){
        let data = [5.0, 4.0, 3.0, 2.0, 1.0];
        let result = detect_trend(&data);
        assert_eq!(result, "Downward trend");
    }

    #[test]
    fn test_no_trend(){
        let data = [3.0, 3.0, 3.0, 3.0, 3.0];
        let result = detect_trend(&data);
        assert_eq!(result, "No trend");
    }
}
