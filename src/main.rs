use plotters::prelude::*;
mod cosos;
use cosos::{MieCalculator, water_refractive_index};
mod complex;
use complex::comp;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Verde:");
    let wavelength = 0.55; // microns
    let particle_radius = 0.5; // microns
    let medium_refractive_index = 1.0; // aire
    let particle_refractive_index = water_refractive_index::get_refractive_index(wavelength);
    let scattering_angle = 90.0;
    
    let calculator = MieCalculator::new(
        wavelength,
        particle_radius,
        medium_refractive_index,
        particle_refractive_index,
    );
    
    let result = calculator.calculate(Some(scattering_angle));
    
    println!("Longitud de onda: {:.3} µm", wavelength);
    println!("Radio: {:.3} µm", particle_radius);
    println!("Parametro de tamaño: {:.3}", calculator.size_parameter());
    println!("Indice de refracción relativo: {:.3} + {:.3e}i", 
             calculator.relative_refractive_index().re,
             calculator.relative_refractive_index().im);
    println!();
    println!("Crosección:");
    println!("  Dispersión: {:.6e} µm²", result.scattering_cross_section);
    println!("  Extinción: {:.6e} µm²", result.extinction_cross_section);
    println!("  Absorción: {:.6e} µm²", result.absorption_cross_section);
    println!();
    println!("Eficiencias:");
    println!("  Dispersión: {:.6}", result.scattering_efficiency);
    println!("  Extinción: {:.6}", result.extinction_efficiency);
    println!("  Absorción: {:.6}", result.absorption_efficiency);
    println!();
    println!("Dispersión a {}°:", scattering_angle);
    println!("  S1 (paralelo): {:.6e} + {:.6e}i", 
             result.scattering_amplitude_parallel.re,
             result.scattering_amplitude_parallel.im);
    println!("  S2 (perpendicular): {:.6e} + {:.6e}i", 
             result.scattering_amplitude_perpendicular.re,
             result.scattering_amplitude_perpendicular.im);
    println!("  Intensidad (paralelo): {:.6e}", result.scattering_intensity_parallel);
    println!("  Intensidad (perpendicular): {:.6e}", result.scattering_intensity_perpendicular);
    println!("\nIntensidad de dispersión por ángulo (for radius = 0.5 µm):");
    println!("Ángulo(°) | Intensidad");
    println!("--------|----------");
    for angle in [10, 30, 60, 90, 120, 150, 180].iter() {
        let result = calculator.calculate(Some(*angle as f64));
        println!("{:6}  | {:.4e}", angle, result.scattering_intensity_parallel);
    }
    let root = BitMapBackend::new("a.png", (1280, 960)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Intensidad de dispersión vs. ángulo para un radio de 1/2 µm", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f32..180f32, -0.1f32..1.5f32)?;

    chart.configure_mesh()
        .y_desc("Intensidad")
        .x_desc("Ángulo (°)")
        .draw()?;
    chart
        .draw_series(LineSeries::new(
            (0..=9000).map(|x| x as f32 / 50.0).map(|x| (x, (calculator.calculate(Some(x as f64)).scattering_intensity_perpendicular/2149.0) as f32)),
            &BLUE,
        ))?
    .label("Intensidad perpendicular")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));

    chart
        .draw_series(LineSeries::new(
            (0..=9000).map(|x| x as f32 / 50.0).map(|x| (x, (calculator.calculate(Some(x as f64)).scattering_intensity_parallel/2149.0) as f32)),
            &RED,
        ))?
    .label("Intensidad paralela")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));

        chart
        .configure_series_labels()
        .background_style(RGBColor(128, 128, 128))
        .draw()?;
    root.present()?;
    let _ = comp(wavelength, particle_radius, medium_refractive_index);
	Ok(())
}
