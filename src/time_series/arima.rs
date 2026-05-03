pub struct Arima {
    p: usize,
    q: usize,
    d: usize,
    ar_coeffs: Vec<f64>,
    ma_coeffs: Vec<f64>,
    differenced_series: Vec<f64>,
}

impl Arima{
    pub fn new(p: usize, d: usize, q: usize, ar_coeffs: Vec<f64>, ma_coeffs: Vec<f64>) -> Self{
        Arima{
            p, 
            q,
            d,
            ar_coeffs,
            ma_coeffs,
            differenced_series: Vec::new(),
        }
    }

    pub fn difference(&mut self, series: &[f64]) -> Vec<f64>{
        let mut differenced = series.to_vec();
        for _ in 0..self.d{
            differenced = differenced.windows(2).map(|w| w[1] - w[0]).collect();
        }
        self.differenced_series = differenced.clone();
        return differenced
    }

    pub fn inverse_difference(&self, forecast: f64, original_series: &[f64]) -> f64{
        let last_original = original_series[original_series.len() - self.d];
        return last_original + forecast
    }

    pub fn predict(&self, series: &[f64]) -> f64{
        let diff_series = self.differenced_series.clone();
        let mut ar_term = 0.0;
        let mut ma_term = 0.0;
        let noise = 0.0;

        if self.p > 0{
            let p_values: Vec<f64> = diff_series.iter().rev().take(self.p).cloned().collect();
            for(i, &phi) in self.ar_coeffs.iter().enumerate(){
                if i < p_values.len(){
                    ar_term += phi * p_values[i];
                }
            }
        }

        if self.q > 0 {
            let q_values: Vec<f64> = diff_series.iter().rev().take(self.q).cloned().collect();
            for(i, &theta) in self.ma_coeffs.iter().enumerate(){
                if i < q_values.len(){
                    ma_term += theta * noise;
                }
            }
        }

        let forecast = ar_term + ma_term + noise;
        return self.inverse_difference(forecast, series);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arima_creation() {
        let arima = Arima::new(2, 1, 1, vec![0.5, 0.3], vec![0.2]);
        assert_eq!(arima.p, 2);
        assert_eq!(arima.d, 1);
        assert_eq!(arima.q, 1);
    }

    #[test]
    fn test_arima_differencing() {
        let mut arima = Arima::new(1, 1, 1, vec![0.5], vec![0.2]);
        let series = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let diff = arima.difference(&series);
        assert_eq!(diff.len(), 4);
        assert_eq!(diff[0], 1.0);
        assert_eq!(diff[1], 1.0);
    }

    #[test]
    fn test_arima_predict() {
        let mut arima = Arima::new(1, 1, 1, vec![0.5], vec![0.2]);
        let series = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        arima.difference(&series);
        let prediction = arima.predict(&series);
        assert!(!prediction.is_nan());
    }
}
