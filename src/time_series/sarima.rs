pub struct Sarima {
    p: usize,
    q: usize,
    d: usize,
    seasonal_p: usize,
    seasonal_q: usize,
    seasonal_d: usize,
    m: usize,
    ar_coeffs: Vec<f64>,
    ma_coeffs: Vec<f64>,
    seasonal_ar_coeffs: Vec<f64>,
    seasonal_ma_coeffs: Vec<f64>,
    differenced_series: Vec<f64>,
}

impl Sarima {
    pub fn new(
        p: usize, d: usize, q: usize,
        seasonal_p: usize, seasonal_q: usize, seasonal_d: usize,
        m: usize,
        ar_coeffs: Vec<f64>, ma_coeffs: Vec<f64>,
        seasonal_ar_coeffs: Vec<f64>, seasonal_ma_coeffs: Vec<f64>,
    ) -> Self {
        Sarima {
            p, d, q, seasonal_p, seasonal_q, seasonal_d, m,
            ar_coeffs, ma_coeffs, seasonal_ar_coeffs, seasonal_ma_coeffs,
            differenced_series: Vec::new(),
        }
    }

    pub fn difference(&mut self, series: &[f64]) -> Vec<f64>{
        let mut differenced = series.to_vec();

        for _ in 0..self.d{
            differenced = differenced.windows(2).map(|w| w[1] - w[0]).collect()
        }

        for _ in 0..self.seasonal_d {
            if differenced.len() > self.m {
                differenced = differenced.iter().skip(self.m).zip(differenced.iter()).map(|(a, b)| a - b).collect();
            }
        }

        self.differenced_series = differenced.clone();
        return differenced;
    }

    pub fn inverse_difference(&self, forecast: f64, original_series: &[f64]) -> f64{
        let last_original = original_series[original_series.len() - self.d];
        return last_original + forecast;
    }

    pub fn predict(&self, series: &[f64]) -> f64{
        let diff_series = self.differenced_series.clone();
        let mut ar_term = 0.0;
        let mut ma_term = 0.0;
        let mut seasonal_ar_term = 0.0;
        let mut seasonal_ma_term = 0.0;
        let noise = 0.0;  

        if self.p > 0{
            let p_values: Vec<f64> = diff_series.iter().rev().take(self.p).cloned().collect();
            for(i, &phi) in self.ar_coeffs.iter().enumerate(){
                if i < p_values.len(){
                    ar_term += phi * p_values[i]
                }
            }
        }

        if self.q > 0{
            let q_values: Vec<f64> = diff_series.iter().rev().take(self.q).cloned().collect();
            for(i, &theta) in self.ma_coeffs.iter().enumerate(){
                if i < q_values.len(){
                    ma_term += theta * noise;
                }
            }
        }

        if self.seasonal_p > 0 {
            let seasonal_p_values: Vec<f64> = diff_series.iter().rev().take(self.seasonal_p * self.m).step_by(self.m).cloned().collect();
            for (i, &phi_s) in self.seasonal_ar_coeffs.iter().enumerate() {
                if i < seasonal_p_values.len() {
                    seasonal_ar_term += phi_s * seasonal_p_values[i];
                }
            }
        }

        if self.seasonal_q > 0 {
            let seasonal_q_values: Vec<f64> = diff_series.iter().rev().take(self.seasonal_q * self.m).step_by(self.m).cloned().collect();
            for (i, &theta_s) in self.seasonal_ma_coeffs.iter().enumerate() {
                if i < seasonal_q_values.len() {
                    seasonal_ma_term += theta_s * noise;
                }
            }
        }

        let forecast = ar_term + ma_term + seasonal_ar_term + seasonal_ma_term;
        return self.inverse_difference(forecast, series);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sarima_creation() {
        let sarima = Sarima::new(
            1, 1, 1, 1, 1, 1, 12,
            vec![0.5], vec![0.3], vec![0.2], vec![0.1]
        );
        assert_eq!(sarima.p, 1);
        assert_eq!(sarima.m, 12);
        assert_eq!(sarima.seasonal_p, 1);
    }

    #[test]
    fn test_sarima_differencing() {
        let mut sarima = Sarima::new(
            1, 1, 1, 0, 0, 0, 4,
            vec![0.5], vec![0.3], vec![], vec![]
        );
        let series = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let diff = sarima.difference(&series);
        assert!(diff.len() > 0);
    }

    #[test]
    fn test_sarima_predict() {
        let mut sarima = Sarima::new(
            1, 1, 1, 1, 1, 1, 4,
            vec![0.5], vec![0.3], vec![0.2], vec![0.1]
        );
        let series = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        sarima.difference(&series);
        let prediction = sarima.predict(&series);
        assert!(!prediction.is_nan());
    }
}
