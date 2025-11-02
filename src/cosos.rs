use std::f64::consts::PI;
use num_complex::Complex64;
type Complex = Complex64;

#[derive(Debug, Clone)]
pub struct MieResult {
    pub scattering_cross_section: f64,
    pub extinction_cross_section: f64,
    pub scattering_efficiency: f64,
    pub extinction_efficiency: f64,
    pub absorption_cross_section: f64,
    pub absorption_efficiency: f64,
    pub scattering_amplitude_parallel: Complex,
    pub scattering_amplitude_perpendicular: Complex,
    pub scattering_intensity_parallel: f64,
    pub scattering_intensity_perpendicular: f64,
}

pub struct MieCalculator {
    pub wavelength: f64,
    pub particle_radius: f64,
    pub medium_refractive_index: f64,
    pub particle_refractive_index: Complex,
}

impl MieCalculator {
    pub fn new(
        wavelength: f64,
        particle_radius: f64,
        medium_refractive_index: f64,
        particle_refractive_index: Complex,
    ) -> Self {
        Self {
            wavelength,
            particle_radius,
            medium_refractive_index,
            particle_refractive_index,
        }
    }

    ///x = 2πr/λ
    pub fn size_parameter(&self) -> f64 {
        2.0 * PI * self.particle_radius / self.wavelength
    }

    /// refracción con respecto al medio m = n_particle / n_medium
    pub fn relative_refractive_index(&self) -> Complex {
        self.particle_refractive_index / self.medium_refractive_index
    }

    /// Riccati-Bessel ψ_n(z) = z * j_n(z)
    pub fn riccati_bessel_psi(&self, n: usize, z: f64) -> f64 {
        if n == 0 {
            z.sin()
        } else if n == 1 {
            z.sin() / z - z.cos()
        } else {
            let mut psi_n_minus_1 = z.sin();
            let mut psi_n = z.sin() / z - z.cos();
            
            for i in 2..=n {
                let i_f64 = i as f64;
                let psi_n_plus_1 = (2.0 * i_f64 - 1.0) / z * psi_n - psi_n_minus_1;
                psi_n_minus_1 = psi_n;
                psi_n = psi_n_plus_1;
            }
            psi_n
        }
    }

    /// Riccati-Bessel ξ_n(z) = z * h_n^(1)(z)
    pub fn riccati_bessel_xi(&self, n: usize, z: f64) -> Complex {
        if n == 0 {
            Complex::new(z.sin(), -z.cos())
        } else if n == 1 {
            Complex::new(z.sin() / z - z.cos(), -z.cos() / z - z.sin())
        } else {
            let mut psi_n_minus_1 = z.sin();
            let mut psi_n = z.sin() / z - z.cos();
            let mut chi_n_minus_1 = -z.cos();
            let mut chi_n = -z.cos() / z - z.sin();
            
            for i in 2..=n {
                let i_f64 = i as f64;
                let psi_n_plus_1 = (2.0 * i_f64 - 1.0) / z * psi_n - psi_n_minus_1;
                let chi_n_plus_1 = (2.0 * i_f64 - 1.0) / z * chi_n - chi_n_minus_1;
                psi_n_minus_1 = psi_n;
                psi_n = psi_n_plus_1;
                chi_n_minus_1 = chi_n;
                chi_n = chi_n_plus_1;
            }
            Complex::new(psi_n, -chi_n)
        }
    }

    /// Derivada de ψ
    pub fn riccati_bessel_psi_prime(&self, n: usize, z: f64) -> f64 {
        if n == 0 {
            z.cos()
        } else {
            let psi_n = self.riccati_bessel_psi(n, z);
            let psi_n_minus_1 = if n > 0 { self.riccati_bessel_psi(n - 1, z) } else { z.sin() };
            psi_n_minus_1 - (n as f64) / z * psi_n
        }
    }

    /// Derivada de ξ_n
    pub fn riccati_bessel_xi_prime(&self, n: usize, z: f64) -> Complex {
        let xi_n = self.riccati_bessel_xi(n, z);
        let xi_n_minus_1 = if n > 0 { self.riccati_bessel_xi(n - 1, z) } else { Complex::new(z.sin(), -z.cos()) };
        xi_n_minus_1 - (n as f64) / z * xi_n
    }

    /// Something something coeficientes
    pub fn mie_coefficients(&self, n_max: usize) -> (Vec<Complex>, Vec<Complex>) {
        let x = self.size_parameter();
        let m = self.relative_refractive_index();
        let mx = m * x;

        let mut a_n = Vec::with_capacity(n_max);
        let mut b_n = Vec::with_capacity(n_max);

        for n in 1..=n_max {
            // Calculate functions and their derivatives
            let psi_n_x = self.riccati_bessel_psi(n, x);
            let xi_n_x = self.riccati_bessel_xi(n, x);
            let psi_prime_n_x = self.riccati_bessel_psi_prime(n, x);
            let xi_prime_n_x = self.riccati_bessel_xi_prime(n, x);
            
            let psi_n_mx = self.riccati_bessel_psi(n, mx.re);
            let psi_prime_n_mx = self.riccati_bessel_psi_prime(n, mx.re);

            // Mie coefficients
            let numerator_a = m * psi_n_mx * psi_prime_n_x - psi_n_x * psi_prime_n_mx;
            let denominator_a = m * psi_n_mx * xi_prime_n_x - xi_n_x * psi_prime_n_mx;
            let a = numerator_a / denominator_a;
            
            let numerator_b = psi_n_mx * psi_prime_n_x - m * psi_n_x * psi_prime_n_mx;
            let denominator_b = psi_n_mx * xi_prime_n_x - m * xi_n_x * psi_prime_n_mx;
            let b = numerator_b / denominator_b;
            
            a_n.push(a);
            b_n.push(b);
        }

        (a_n, b_n)
    }

    /// Amplitudes de la dispersión como función del ángulo.
    pub fn scattering_amplitudes(&self, theta: f64, n_max: usize) -> (Complex, Complex) {
        let (a_n, b_n) = self.mie_coefficients(n_max);
        let mu = theta.cos();
        
        let mut s1 = Complex::new(0.0, 0.0);
        let mut s2 = Complex::new(0.0, 0.0);
        
        // Calculate Legendre polynomials and their derivatives
        let mut pi_n_minus_1 = 0.0;
        let mut pi_n = 1.0;
        
        for n in 1..=n_max {
            let n_f64 = n as f64;
            
            // Recurrence relations for Legendre functions
            let pi_n_plus_1 = ((2.0 * n_f64 + 1.0) * mu * pi_n - (n_f64 + 1.0) * pi_n_minus_1) / n_f64;
            let tau_n = (n_f64 + 1.0) * mu * pi_n_plus_1 - (n_f64 + 2.0) * pi_n;
            
            let factor = (2.0 * n_f64 + 1.0) / (n_f64 * (n_f64 + 1.0));
            
            s1 += factor * (a_n[n-1] * pi_n_plus_1 + b_n[n-1] * tau_n);
            s2 += factor * (a_n[n-1] * tau_n + b_n[n-1] * pi_n_plus_1);
            
            // Update for next iteration
            pi_n_minus_1 = pi_n;
            pi_n = pi_n_plus_1;
        }
        
        (s1, s2)
    }


    /// Calculate all Mie scattering properties
    pub fn calculate(&self, theta: Option<f64>) -> MieResult {
        let x = self.size_parameter();
        
        // Determine number of terms needed (Wiscombe criterion)
        let n_max = (x + 4.0 * x.powf(1.0/3.0) + 2.0).ceil() as usize;
        let n_max = n_max.max(1).min(1000); // Safety bounds
        
        let (a_n, b_n) = self.mie_coefficients(n_max);
        
        // Calculate cross sections
        let mut scattering_cross_section = 0.0;
        let mut extinction_cross_section = 0.0;
        
        for n in 1..=n_max {
            let n_f64 = n as f64;
            let term = (2.0 * n_f64 + 1.0) * (a_n[n-1].norm_sqr() + b_n[n-1].norm_sqr());
            scattering_cross_section += term;
            
            let extinction_term = (2.0 * n_f64 + 1.0) * (a_n[n-1].re + b_n[n-1].re);
            extinction_cross_section += extinction_term;
        }
        
        let geometric_cross_section = PI * self.particle_radius * self.particle_radius;
        
        scattering_cross_section *= 2.0 * PI / (x * x) * geometric_cross_section;
        extinction_cross_section *= 2.0 * PI / (x * x) * geometric_cross_section;
        
        let scattering_efficiency = scattering_cross_section / geometric_cross_section;
        let extinction_efficiency = extinction_cross_section / geometric_cross_section;
        let absorption_cross_section = extinction_cross_section - scattering_cross_section;
        let absorption_efficiency = absorption_cross_section / geometric_cross_section;
        
        // Calculate scattering amplitudes and intensities at given angle
        let (s1, s2, i1, i2) = if let Some(angle) = theta {
            let (s1, s2) = self.scattering_amplitudes(angle.to_radians(), n_max);
            let i1 = s1.norm_sqr();
            let i2 = s2.norm_sqr();
            (s1, s2, i1, i2)
        } else {
            (Complex::new(0.0, 0.0), Complex::new(0.0, 0.0), 0.0, 0.0)
        };
        
        MieResult {
            scattering_cross_section,
            extinction_cross_section,
            scattering_efficiency,
            extinction_efficiency,
            absorption_cross_section,
            absorption_efficiency,
            scattering_amplitude_parallel: s1,
            scattering_amplitude_perpendicular: s2,
            scattering_intensity_parallel: i1,
            scattering_intensity_perpendicular: i2,
        }
    }
}

// Info del agua
pub mod water_refractive_index {
    use super::Complex;
    pub fn get_refractive_index(wavelength_um: f64) -> Complex {
        // https://refractiveindex.info/?shelf=main&book=H2O&page=Hale
        match wavelength_um {
            w if w >= 0.2 && w < 0.3 => Complex::new(1.362, 3.3e-8),  // UV
            w if w >= 0.3 && w < 0.4 => Complex::new(1.343, 6.5e-9),
            w if w >= 0.4 && w < 0.5 => Complex::new(1.337, 1.0e-9),
            w if w >= 0.5 && w < 0.6 => Complex::new(1.333, 1.9e-9),  // Visible
            w if w >= 0.6 && w < 0.7 => Complex::new(1.331, 1.6e-8),
            w if w >= 0.7 && w < 0.8 => Complex::new(1.330, 1.5e-7),
            w if w >= 0.8 && w < 1.0 => Complex::new(1.328, 4.8e-7),
            w if w >= 1.0 && w < 2.0 => Complex::new(1.319, 1.0e-4),  // Near IR - some absorption
            _ => Complex::new(1.319, 1.0e-4)  // Default value
        }
    }
}
