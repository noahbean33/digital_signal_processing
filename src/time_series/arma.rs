pub fn arma(ar_coeffs: &[f64], ma_coeffs: &[f64], c: f64, n: usize) -> Vec<f64>{
    let mut y_values = vec![0.0; n];
    let mut e_values = vec![0.0; n];

    for t in 0..n{
        e_values[t] = (t as f64).sin() * 0.1;
    }

    for t in 0..n{
        let mut ar_term = 0.0;
        for(i, &phi) in ar_coeffs.iter().enumerate(){
            if t >= i  + 1{
                ar_term += phi * y_values[t - i - 1];
            }
        }
        
        let mut ma_term = 0.0;
        for(j, &theta) in ma_coeffs.iter().enumerate(){
            if t >= j + 1{
                ma_term += theta * e_values[t - j - 1];
            }
        }

        y_values[t] = c + ar_term + ma_term + e_values[t];
    }
    return y_values;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arma_basic() {
        let ar_coeffs = vec![0.5, 0.3];
        let ma_coeffs = vec![0.2];
        let c = 0.0;
        let n = 10;
        let result = arma(&ar_coeffs, &ma_coeffs, c, n);
        assert_eq!(result.len(), n);
    }

    #[test]
    fn test_arma_with_constant() {
        let ar_coeffs = vec![0.7];
        let ma_coeffs = vec![0.3];
        let c = 2.0;
        let n = 5;
        let result = arma(&ar_coeffs, &ma_coeffs, c, n);
        assert_eq!(result.len(), n);
        assert!(result.iter().all(|&x| !x.is_nan()));
    }
}
