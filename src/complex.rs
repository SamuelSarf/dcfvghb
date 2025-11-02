use std::f64::consts::PI;
use num_complex::Complex64;
use plotters::prelude::*;

type Complex = Complex64;

struct MieCalculator {
    wavelength: f64,
    particle_radius: f64,
    medium_refractive_index: f64,
    particle_refractive_index: Complex,
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

    fn size_parameter(&self) -> f64 {
        2.0 * PI * self.particle_radius / self.wavelength
    }

    fn relative_refractive_index(&self) -> Complex {
        self.particle_refractive_index / self.medium_refractive_index
    }

    fn riccati_bessel_psi(&self, n: usize, z: f64) -> f64 {
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

    fn riccati_bessel_xi(&self, n: usize, z: f64) -> Complex {
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

    fn riccati_bessel_psi_prime(&self, n: usize, z: f64) -> f64 {
        if n == 0 {
            z.cos()
        } else {
            let psi_n = self.riccati_bessel_psi(n, z);
            let psi_n_minus_1 = if n > 0 { self.riccati_bessel_psi(n - 1, z) } else { z.sin() };
            psi_n_minus_1 - (n as f64) / z * psi_n
        }
    }

    fn riccati_bessel_xi_prime(&self, n: usize, z: f64) -> Complex {
        let xi_n = self.riccati_bessel_xi(n, z);
        let xi_n_minus_1 = if n > 0 { self.riccati_bessel_xi(n - 1, z) } else { Complex::new(z.sin(), -z.cos()) };
        xi_n_minus_1 - (n as f64) / z * xi_n
    }

    fn mie_coefficients(&self, n_max: usize) -> (Vec<Complex>, Vec<Complex>) {
        let x = self.size_parameter();
        let m = self.relative_refractive_index();
        let mx = m * x;

        let mut a_n = Vec::with_capacity(n_max);
        let mut b_n = Vec::with_capacity(n_max);

        for n in 1..=n_max {
            let psi_n_x = self.riccati_bessel_psi(n, x);
            let xi_n_x = self.riccati_bessel_xi(n, x);
            let psi_prime_n_x = self.riccati_bessel_psi_prime(n, x);
            let xi_prime_n_x = self.riccati_bessel_xi_prime(n, x);
            
            let psi_n_mx = self.riccati_bessel_psi(n, mx.re);
            let psi_prime_n_mx = self.riccati_bessel_psi_prime(n, mx.re);

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

    pub fn scattering_amplitudes(&self, theta: f64) -> (Complex, Complex) {
        let x = self.size_parameter();
        let n_max = (x + 4.0 * x.powf(1.0/3.0) + 2.0).ceil() as usize;
        let n_max = n_max.max(1).min(1000);
        
        let (a_n, b_n) = self.mie_coefficients(n_max);
        let mu = theta.cos();
        
        let mut s1 = Complex::new(0.0, 0.0);
        let mut s2 = Complex::new(0.0, 0.0);
        
        let mut pi_n_minus_1 = 0.0;
        let mut pi_n = 1.0;
        
        for n in 1..=n_max {
            let n_f64 = n as f64;
            
            let pi_n_plus_1 = ((2.0 * n_f64 + 1.0) * mu * pi_n - (n_f64 + 1.0) * pi_n_minus_1) / n_f64;
            let tau_n = (n_f64 + 1.0) * mu * pi_n_plus_1 - (n_f64 + 2.0) * pi_n;
            
            let factor = (2.0 * n_f64 + 1.0) / (n_f64 * (n_f64 + 1.0));
            
            s1 += factor * (a_n[n-1] * pi_n_plus_1 + b_n[n-1] * tau_n);
            s2 += factor * (a_n[n-1] * tau_n + b_n[n-1] * pi_n_plus_1);
            
            pi_n_minus_1 = pi_n;
            pi_n = pi_n_plus_1;
        }
        
        (s1, s2)
    }
}

pub fn water_refractive_index(wavelength_um: f64) -> Complex {
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

fn complex_to_color(z: Complex, max_magnitude: f64) -> RGBColor {
    let magnitude = z.norm();
    let phase = z.arg();
    
    let hue = if phase < 0.0 { phase + 2.0 * PI } else { phase } / (2.0 * PI);
    
    let log_magnitude = (magnitude + 1.0).ln();
    let log_max = (max_magnitude + 1.0).ln();
    let brightness = (log_magnitude / log_max).min(1.0);
    
    hsv_to_rgb(hue, 1.0, brightness)
}

fn hsv_to_rgb(h: f64, s: f64, v: f64) -> RGBColor {
    let h = h * 6.0;
    let i = h.floor() as i32;
    let f = h - i as f64;
    let p = v * (1.0 - s);
    let q = v * (1.0 - s * f);
    let t = v * (1.0 - s * (1.0 - f));
    
    let (r, g, b) = match i % 6 {
        0 => (v, t, p),
        1 => (q, v, p),
        2 => (p, v, t),
        3 => (p, q, v),
        4 => (t, p, v),
        5 => (v, p, q),
        _ => (0.0, 0.0, 0.0),
    };
    
    RGBColor(
        (r * 255.0) as u8,
        (g * 255.0) as u8,
        (b * 255.0) as u8,
    )
}

pub fn comp(wavelength: f64, particle_radius: f64, medium_refractive_index: f64) -> Result<(), Box<dyn std::error::Error>> {
    let particle_refractive_index = water_refractive_index(wavelength);
    
    let calculator = MieCalculator::new(
        wavelength,
        particle_radius,
        medium_refractive_index,
        particle_refractive_index,
    );
    
    let root = BitMapBackend::new("c.png", (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;
    
    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Amplitud de S₁(θ)",
            ("sans-serif", 30).into_font(),
        )
        .margin(50)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..180.0, 0.0..1.0)?;
    
    chart
        .configure_mesh()
        .x_desc("Ángulo θ (°)")
        .y_desc("|S₁(θ)| (Normalizado)")
        .draw()?;
    
    let mut max_magnitude: f64 = 0.0;
    let num_points = 1000;
    
    for i in 0..=num_points {
        let theta_deg = (i as f64) * 180.0 / (num_points as f64);
        let theta_rad = theta_deg.to_radians();
        let (s1, _s2) = calculator.scattering_amplitudes(theta_rad);
        max_magnitude = max_magnitude.max(s1.norm());
    }
    
    println!("Magnitud máxima de S₁(θ): {:.6}", max_magnitude);
    
    let drawing_area = chart.plotting_area();
    let vertical_resolution = 50;
    
    for deg in 0..180 {
        let theta_deg = deg as f64;
        let theta_rad = theta_deg.to_radians();
        
        let (s1, s2) = calculator.scattering_amplitudes(theta_rad);
        let normalized_magnitude = s1.norm() / max_magnitude;
        
        for v_idx in 0..vertical_resolution {
            let y_pos = (v_idx as f64) / (vertical_resolution as f64) * normalized_magnitude;
            let color = complex_to_color(s2, max_magnitude);
            
            drawing_area.draw(&Circle::new(
                (theta_deg, y_pos),
                2,
                ShapeStyle::from(color).filled(),
            ))?;
        }
    }
    
    // Add text labels to the main drawing area (not the plotting area)
    root.draw(&Text::new(
        format!("Brillo = |S₂(θ)| (max = {:.3})", max_magnitude),
        (50, 30),
        ("sans-serif", 20).into_font(),
    ))?;
    
    root.draw(&Text::new(
        format!("Radio de partícula: {} µm, Longitud de onda: {} µm", particle_radius, wavelength),
        (50, 60),
        ("sans-serif", 16).into_font(),
    ))?;
    
    create_legend(&root)?;
    
    
    
    Ok(())
}

fn create_legend(root: &DrawingArea<BitMapBackend, plotters::coord::Shift>) -> Result<(), Box<dyn std::error::Error>> {
    let legend_width = 30;
    let legend_height = 150;
    let legend_x = 1000;
    let legend_y = 100;
    
    for h in 0..legend_height {
        let hue = h as f64 / legend_height as f64;
        let color = hsv_to_rgb(hue, 1.0, 1.0);
        
        root.draw(&Rectangle::new(
            [
                (legend_x, legend_y + h),
                (legend_x + legend_width, legend_y + h + 1)
            ],
            ShapeStyle::from(color).filled(),
        ))?;
    }
    
    root.draw(&Text::new(
        "Fase:".to_string(),
        (legend_x, legend_y - 20),
        ("sans-serif", 15).into_font(),
    ))?;
    
    root.draw(&Text::new(
        "0°".to_string(),
        (legend_x + legend_width + 5, legend_y),
        ("sans-serif", 12).into_font(),
    ))?;
    
    root.draw(&Text::new(
        "180°".to_string(),
        (legend_x + legend_width + 5, legend_y + legend_height),
        ("sans-serif", 12).into_font(),
    ))?;
    
    Ok(())
}

