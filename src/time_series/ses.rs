pub struct ExponentialSmoothing{
    pub alpha: f64,
    pub beta: Option<f64>,
}

impl ExponentialSmoothing{
    pub fn single_exponential_smoothing(&self, data: &[f64]) -> Vec<f64>{
        let alpha = self.alpha;
        let mut smoothed = Vec::new();
        let mut x = data[0];

        for &s in data{
            x = alpha * s + (1.0 - alpha) * x;
            smoothed.push(x);
        }

        return smoothed
    }

    pub fn holt_linear_model(&self, data: &[f64]) -> (Vec<f64>, Vec<f64>){
        let alpha = self.alpha;
        let beta = self.beta.expect("Beta param is required for Holt Linear Model");

        let mut smoothed = Vec::new();
        let mut trend_smoothed = Vec::new();
        let mut s = data[0];
        let mut t = data[1] - data[0];

        for &x in data{
            let prev_s = s;
            s = alpha * x + (1.0 - alpha) * (s + t);
            t = beta * (s - prev_s) + (1.0 - beta) * t;
            smoothed.push(s);
            trend_smoothed.push(t);
        }

        return (smoothed, trend_smoothed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_exponential_smoothing() {
        let data = vec![10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
        let model = ExponentialSmoothing{
            alpha: 0.5,
            beta: None,
        };
        let result = model.single_exponential_smoothing(&data);
        assert_eq!(result.len(), data.len());
        assert_eq!(result[0], 10.0);
        assert!(result[result.len() - 1] > result[0]);
    }

    #[test]
    fn test_holt_linear_model() {
        let data = vec![10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
        let model = ExponentialSmoothing{
            alpha: 0.5,
            beta: Some(0.3),
        };
        let (smoothed, trends) = model.holt_linear_model(&data);
        assert_eq!(smoothed.len(), data.len());
        assert_eq!(trends.len(), data.len());
    }

    #[test]
    #[should_panic(expected = "Beta param is required for Holt Linear Model")]
    fn test_holt_without_beta() {
        let data = vec![10.0, 12.0, 14.0];
        let model = ExponentialSmoothing{
            alpha: 0.5,
            beta: None,
        };
        model.holt_linear_model(&data);
    }
}
