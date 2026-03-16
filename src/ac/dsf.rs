//! Dark Spectrum Fitting (DSF) atmospheric correction
//!
//! Implements the ACOLITE DSF algorithm:
//! 1. Extract dark spectrum (lowest percentile) per band
//! 2. For each aerosol model + band, invert LUT: dark_rhot → AOT
//! 3. Select minimum AOT across bands, best model by RMSD
//! 4. Apply correction: rhos = (rhot/tt_gas - romix) / (dutott + astot*(rhot/tt_gas - romix))

use crate::Result;
use ndarray::Array2;
use std::collections::HashMap;

#[cfg(feature = "full-io")]
use crate::ac::aerlut::AerosolLut;
#[cfg(feature = "full-io")]
use crate::ac::aerlut::{rsr_convolve_gauss, rsr_convolve_sensor, GenericAerosolLut};

/// Dark spectrum extraction method
#[derive(Debug, Clone)]
pub enum DarkSpectrumMethod {
    /// Nth percentile of valid pixels
    Percentile(f64),
    /// Linear regression intercept of N darkest pixels (ACOLITE default)
    Intercept(usize),
}

/// DSF AOT estimation mode (matches Python's dsf_aot_estimate setting)
#[derive(Debug, Clone)]
pub enum DsfMode {
    /// Single scene-wide dark spectrum → one AOT for the whole image
    Fixed,
    /// Per-tile dark spectrum → spatially varying AOT grid
    Tiled(usize, usize),
}

/// DSF configuration
#[derive(Debug, Clone)]
pub struct DsfConfig {
    pub dark_method: DarkSpectrumMethod,
    pub min_tgas_aot: f64,
    pub dsf_nbands: usize,
    pub dsf_nbands_fit: usize,
    pub aot_compute: AotCompute,
    pub wave_range: (f64, f64),
    /// If set, skip model voting and use only the LUT whose name ends with this suffix (e.g. "MOD1")
    pub fixed_model: Option<String>,
    /// AOT estimation mode: Fixed (whole-scene) or Tiled(rows, cols)
    pub mode: DsfMode,
}

/// How to compute the selected AOT from per-band AOTs
#[derive(Debug, Clone)]
pub enum AotCompute {
    /// Use the minimum AOT across all bands (ACOLITE default)
    Min,
    /// Use the mean of the N lowest AOT bands
    MeanNLowest,
}

impl Default for DsfConfig {
    fn default() -> Self {
        Self {
            dark_method: DarkSpectrumMethod::Intercept(200),
            min_tgas_aot: 0.85,
            dsf_nbands: 2,
            dsf_nbands_fit: 2,
            aot_compute: AotCompute::Min,
            wave_range: (400.0, 2500.0),
            fixed_model: None,
            mode: DsfMode::Tiled(200, 200),
        }
    }
}

/// Result of DSF AOT estimation
#[derive(Debug, Clone)]
pub struct DsfResult {
    pub aot: f64,
    pub model_idx: usize,
    pub model_name: String,
    pub rmsd: f64,
    pub band_aots: HashMap<String, f64>,
}

/// Extract dark spectrum per band using the configured method.
///
/// - `Percentile(p)`: Nth percentile of valid (finite, >0) pixels
/// - `Intercept(n)`: Linear regression y-intercept of N darkest valid pixels
///   (matches Python ACOLITE's `ac.shared.intercept`)
pub fn estimate_dark_spectrum(bands: &[Array2<f64>], method: &DarkSpectrumMethod) -> Vec<f64> {
    bands
        .iter()
        .map(|band| {
            let mut values: Vec<f64> = band
                .iter()
                .copied()
                .filter(|v| v.is_finite() && *v > 0.0)
                .collect();
            if values.is_empty() {
                return f64::NAN;
            }
            values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            match method {
                DarkSpectrumMethod::Percentile(pct) => {
                    let idx = ((values.len() as f64) * pct / 100.0) as usize;
                    values[idx.min(values.len() - 1)]
                }
                DarkSpectrumMethod::Intercept(npix) => {
                    let n = (*npix).min(values.len());
                    if n == 0 {
                        return 0.0;
                    }
                    // Linear regression: y = slope*x + intercept over indices 0..n
                    let xmean = (n - 1) as f64 / 2.0;
                    let ymean: f64 = values[..n].iter().sum::<f64>() / n as f64;
                    let mut sxx = 0.0_f64;
                    let mut sxy = 0.0_f64;
                    for (i, &y) in values[..n].iter().enumerate() {
                        let dx = i as f64 - xmean;
                        sxx += dx * dx;
                        sxy += dx * (y - ymean);
                    }
                    if sxx == 0.0 {
                        return ymean;
                    }
                    let slope = sxy / sxx;
                    ymean - slope * xmean // intercept
                }
            }
        })
        .collect()
}

/// Invert LUT to find AOT from observed dark reflectance for a single band.
/// Interpolates romix(tau) at given geometry, finds tau where romix == observed rhot.
#[cfg(feature = "full-io")]
fn invert_aot(
    lut: &AerosolLut,
    band: &str,
    rhot_dark: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> f64 {
    // Evaluate romix at each tau step
    let tau_steps = &lut.meta.tau;
    let romix_vals: Vec<f64> = tau_steps
        .iter()
        .map(|&tau| lut.romix(band, pressure, azi, thv, ths, tau))
        .collect();

    // Interpolate: find tau where romix == rhot_dark
    // romix is monotonically increasing with tau
    // Return NaN if out of range (matches Python np.interp with left=nan, right=nan)
    if rhot_dark <= romix_vals[0] {
        return f64::NAN;
    }
    if rhot_dark >= romix_vals[romix_vals.len() - 1] {
        return f64::NAN;
    }

    for i in 0..romix_vals.len() - 1 {
        if rhot_dark >= romix_vals[i] && rhot_dark <= romix_vals[i + 1] {
            let t = (rhot_dark - romix_vals[i]) / (romix_vals[i + 1] - romix_vals[i]);
            return tau_steps[i] * (1.0 - t) + tau_steps[i + 1] * t;
        }
    }
    f64::NAN
}

/// Run full DSF AOT estimation across multiple models.
///
/// For each model:
///   1. For each band, invert LUT to get per-band AOT from dark spectrum
///   2. Select AOT: min across bands (default) or mean of N lowest
///   3. Compute RMSD between observed dark spectrum and modeled path reflectance
///      for the N best-fitting bands at the selected AOT
///
/// `dark_spectrum` should be already gas-corrected (rhot / tt_gas).
///
/// Returns the best model + AOT.
#[cfg(feature = "full-io")]
pub fn optimize_aot(
    luts: &[AerosolLut],
    dark_spectrum: &[f64], // per-band dark rhot (already gas-corrected)
    band_names: &[String], // LUT band names (e.g. "1", "2", ...)
    wavelengths: &[f64],   // per-band wavelength in nm
    tt_gas: &[f64],        // per-band gas transmittance (for filtering only)
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
) -> DsfResult {
    let mut best = DsfResult {
        aot: 0.01,
        model_idx: 0,
        model_name: String::new(),
        rmsd: f64::MAX,
        band_aots: HashMap::new(),
    };

    for (li, lut) in luts.iter().enumerate() {
        // Filter bands by wavelength range and gas transmittance
        let mut band_aots: Vec<(usize, f64)> = Vec::new();
        for (bi, bname) in band_names.iter().enumerate() {
            if wavelengths[bi] < config.wave_range.0 || wavelengths[bi] > config.wave_range.1 {
                continue;
            }
            if tt_gas[bi] < config.min_tgas_aot {
                continue;
            }
            if !dark_spectrum[bi].is_finite() || dark_spectrum[bi] <= 0.0 {
                continue;
            }

            // dark_spectrum is already gas-corrected
            let rhot_dark = dark_spectrum[bi];

            let aot = invert_aot(lut, bname, rhot_dark, pressure, azi, thv, ths);
            if aot.is_finite() {
                band_aots.push((bi, aot));
            }
        }

        if band_aots.is_empty() {
            continue;
        }

        // Sort by AOT ascending
        band_aots.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        // Select AOT based on config
        let selected_aot = match config.aot_compute {
            AotCompute::Min => band_aots[0].1,
            AotCompute::MeanNLowest => {
                let n = config.dsf_nbands.min(band_aots.len());
                band_aots[..n].iter().map(|x| x.1).sum::<f64>() / n as f64
            }
        };

        // Compute RMSD for the N best-fitting bands (those with lowest AOT)
        let n_fit = config.dsf_nbands_fit.min(band_aots.len());
        let mut sum_sq = 0.0;
        let mut ba_map = HashMap::new();
        for &(bi, aot_bi) in &band_aots {
            ba_map.insert(band_names[bi].clone(), aot_bi);
        }
        for &(bi, _) in &band_aots[..n_fit] {
            let bname = &band_names[bi];
            let rhot_dark = dark_spectrum[bi]; // already gas-corrected
            let romix_model = lut.romix(bname, pressure, azi, thv, ths, selected_aot);
            sum_sq += (rhot_dark - romix_model).powi(2);
        }
        let rmsd = (sum_sq / n_fit as f64).sqrt();

        log::info!(
            "Model {} RMSD={:.6} AOT={:.4}",
            lut.name,
            rmsd,
            selected_aot
        );

        if rmsd < best.rmsd {
            best = DsfResult {
                aot: selected_aot,
                model_idx: li,
                model_name: lut.name.clone(),
                rmsd,
                band_aots: ba_map,
            };
        }
    }

    best
}

/// Apply atmospheric correction to a single band using LUT parameters.
///
/// Formula (matching Python ACOLITE):
///   rhot_noatm = (rhot / tt_gas) - romix
///   rhos = rhot_noatm / (dutott + astot * rhot_noatm)
#[cfg(feature = "full-io")]
pub fn dsf_correct_band(
    rhot: &Array2<f64>,
    lut: &AerosolLut,
    band: &str,
    tt_gas: f64,
    aot: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> Array2<f64> {
    let romix = lut.romix(band, pressure, azi, thv, ths, aot);
    let astot = lut.astot(band, pressure, azi, thv, ths, aot);
    let dutott = lut.dutott(band, pressure, azi, thv, ths, aot);

    rhot.mapv(|v| {
        if !v.is_finite() || v <= 0.0 {
            return f64::NAN;
        }
        let rhot_gc = v / tt_gas;
        let rhot_noatm = rhot_gc - romix;
        let denom = dutott + astot * rhot_noatm;
        if denom <= 0.0 {
            return f64::NAN;
        }
        rhot_noatm / denom
    })
}

/// Fixed (whole-scene) DSF: extract one dark spectrum from the full image,
/// invert AOT once, return a single scene-wide result.
/// Matches Python's `dsf_aot_estimate=fixed` behavior.
#[cfg(feature = "full-io")]
pub fn optimize_aot_fixed(
    luts: &[AerosolLut],
    toa_gc_bands: &[Array2<f64>],
    band_names: &[String],
    wavelengths: &[f64],
    tt_gas: &[f64],
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
) -> DsfResult {
    let dark_spectrum = estimate_dark_spectrum(toa_gc_bands, &config.dark_method);
    optimize_aot(
        luts,
        &dark_spectrum,
        band_names,
        wavelengths,
        tt_gas,
        pressure,
        azi,
        thv,
        ths,
        config,
    )
}

/// Result of tiled DSF estimation
#[cfg(feature = "full-io")]
#[derive(Debug, Clone)]
pub struct TiledDsfResult {
    /// Per-tile AOT grid (ni × nj)
    pub aot_grid: Vec<Vec<f64>>,
    /// Selected model index (single model, most-common across tiles)
    pub model_idx: usize,
    pub model_name: String,
    /// Tile dimensions used
    pub tile_rows: usize,
    pub tile_cols: usize,
    /// Number of tile rows/cols
    pub ni: usize,
    pub nj: usize,
}

/// Run tiled DSF: divide image into tiles, estimate AOT per tile, select most-common model.
///
/// Matches Python ACOLITE's `dsf_aot_estimate=tiled` with `dsf_spectrum_option=intercept`.
///
/// `toa_gc_bands`: gas-corrected TOA arrays (rhot / tt_gas), one per band.
/// Returns a `TiledDsfResult` with per-tile AOT on a grid.
#[cfg(feature = "full-io")]
pub fn optimize_aot_tiled(
    luts: &[AerosolLut],
    toa_gc_bands: &[Array2<f64>],
    band_names: &[String],
    wavelengths: &[f64],
    tt_gas: &[f64],
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
    tile_size: (usize, usize), // (tile_rows, tile_cols)
) -> TiledDsfResult {
    let (rows, cols) = toa_gc_bands[0].dim();
    let ni = (rows + tile_size.0 - 1) / tile_size.0;
    let nj = (cols + tile_size.1 - 1) / tile_size.1;
    let nbands = toa_gc_bands.len();

    // Per-tile, per-model: run DSF
    // model_votes[model_idx] = count of tiles where this model was best
    let mut model_votes = vec![0usize; luts.len()];
    // per_tile_results[ti][tj] = (best_model_idx, aot, rmsd) per model
    let mut tile_results: Vec<Vec<Vec<(usize, f64, f64)>>> = vec![vec![Vec::new(); nj]; ni];

    for ti in 0..ni {
        let r0 = ti * tile_size.0;
        let r1 = ((ti + 1) * tile_size.0).min(rows);
        for tj in 0..nj {
            let c0 = tj * tile_size.1;
            let c1 = ((tj + 1) * tile_size.1).min(cols);

            // Extract dark spectrum for this tile
            let tile_bands: Vec<Array2<f64>> = toa_gc_bands
                .iter()
                .map(|b| b.slice(ndarray::s![r0..r1, c0..c1]).to_owned())
                .collect();

            // Check minimum tile coverage
            let total_pixels = (r1 - r0) * (c1 - c0);
            let valid_count = tile_bands[0]
                .iter()
                .filter(|v| v.is_finite() && **v > 0.0)
                .count();
            if (valid_count as f64) < total_pixels as f64 * 0.1 {
                continue; // skip tiles with <10% valid data
            }

            let dark_spectrum = estimate_dark_spectrum(&tile_bands, &config.dark_method);

            // Run optimize_aot for this tile
            let result = optimize_aot(
                luts,
                &dark_spectrum,
                band_names,
                wavelengths,
                tt_gas,
                pressure,
                azi,
                thv,
                ths,
                config,
            );

            if result.rmsd < f64::MAX && result.aot.is_finite() {
                // Filter extreme AOTs (matching Python dsf_min/max_tile_aot)
                if result.aot >= 0.001 && result.aot <= 1.5 {
                    tile_results[ti][tj].push((result.model_idx, result.aot, result.rmsd));
                    model_votes[result.model_idx] += 1;
                }
            }
        }
    }

    // Select most common model (matching Python dsf_aot_most_common_model)
    // If fixed_model is set, find the matching LUT index instead of voting
    let best_model_idx = if let Some(ref model_name) = config.fixed_model {
        luts.iter()
            .position(|l| l.name.ends_with(model_name))
            .unwrap_or(0)
    } else {
        model_votes
            .iter()
            .enumerate()
            .max_by_key(|(_, &v)| v)
            .map(|(i, _)| i)
            .unwrap_or(0)
    };

    log::info!(
        "Model votes: {}",
        model_votes
            .iter()
            .enumerate()
            .map(|(i, v)| format!(
                "{}={}",
                luts.get(i).map(|l| l.name.as_str()).unwrap_or("?"),
                v
            ))
            .collect::<Vec<_>>()
            .join(", ")
    );
    log::info!(
        "Selected model: {} (idx={})",
        luts[best_model_idx].name,
        best_model_idx
    );

    // Now rebuild AOT grid using only the selected model
    // Re-run per-tile DSF with only the selected model
    let mut aot_grid = vec![vec![f64::NAN; nj]; ni];

    for ti in 0..ni {
        let r0 = ti * tile_size.0;
        let r1 = ((ti + 1) * tile_size.0).min(rows);
        for tj in 0..nj {
            let c0 = tj * tile_size.1;
            let c1 = ((tj + 1) * tile_size.1).min(cols);

            let tile_bands: Vec<Array2<f64>> = toa_gc_bands
                .iter()
                .map(|b| b.slice(ndarray::s![r0..r1, c0..c1]).to_owned())
                .collect();

            let total_pixels = (r1 - r0) * (c1 - c0);
            let valid_count = tile_bands[0]
                .iter()
                .filter(|v| v.is_finite() && **v > 0.0)
                .count();
            if (valid_count as f64) < total_pixels as f64 * 0.1 {
                continue;
            }

            let dark_spectrum = estimate_dark_spectrum(&tile_bands, &config.dark_method);

            // Invert AOT for each band using the selected model
            let mut band_aots: Vec<f64> = Vec::new();
            for (bi, bname) in band_names.iter().enumerate() {
                if wavelengths[bi] < config.wave_range.0 || wavelengths[bi] > config.wave_range.1 {
                    continue;
                }
                if tt_gas[bi] < config.min_tgas_aot {
                    continue;
                }
                if !dark_spectrum[bi].is_finite() || dark_spectrum[bi] <= 0.0 {
                    continue;
                }
                let aot = invert_aot(
                    &luts[best_model_idx],
                    bname,
                    dark_spectrum[bi],
                    pressure,
                    azi,
                    thv,
                    ths,
                );
                if aot.is_finite() && aot >= 0.001 && aot <= 1.5 {
                    band_aots.push(aot);
                }
            }
            if band_aots.is_empty() {
                continue;
            }
            band_aots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

            aot_grid[ti][tj] = match config.aot_compute {
                AotCompute::Min => band_aots[0],
                AotCompute::MeanNLowest => {
                    let n = config.dsf_nbands.min(band_aots.len());
                    band_aots[..n].iter().sum::<f64>() / n as f64
                }
            };
        }
    }

    // Fill NaN tiles with nearest valid value (matching Python's distance_transform_edt fill)
    fill_nan_nearest(&mut aot_grid, ni, nj);

    {
        let mut vals: Vec<f64> = aot_grid
            .iter()
            .flatten()
            .copied()
            .filter(|v| v.is_finite())
            .collect();
        if !vals.is_empty() {
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            let median = vals[vals.len() / 2];
            log::info!(
                "AOT grid ({}×{}): mean={:.4}, median={:.4}, n_valid={}",
                ni,
                nj,
                mean,
                median,
                vals.len()
            );
        }
    }

    TiledDsfResult {
        aot_grid,
        model_idx: best_model_idx,
        model_name: luts[best_model_idx].name.clone(),
        tile_rows: tile_size.0,
        tile_cols: tile_size.1,
        ni,
        nj,
    }
}

/// Fill NaN values in a 2D grid with nearest valid neighbor (BFS).
fn fill_nan_nearest(grid: &mut Vec<Vec<f64>>, ni: usize, nj: usize) {
    use std::collections::VecDeque;
    let mut queue = VecDeque::new();
    let mut filled = vec![vec![false; nj]; ni];

    // Seed with valid cells
    for i in 0..ni {
        for j in 0..nj {
            if grid[i][j].is_finite() {
                filled[i][j] = true;
                queue.push_back((i, j));
            }
        }
    }

    // BFS expansion
    while let Some((i, j)) = queue.pop_front() {
        let val = grid[i][j];
        for &(di, dj) in &[(0isize, 1isize), (0, -1), (1, 0), (-1, 0)] {
            let ni2 = i as isize + di;
            let nj2 = j as isize + dj;
            if ni2 >= 0 && ni2 < ni as isize && nj2 >= 0 && nj2 < nj as isize {
                let (ni2, nj2) = (ni2 as usize, nj2 as usize);
                if !filled[ni2][nj2] {
                    grid[ni2][nj2] = val;
                    filled[ni2][nj2] = true;
                    queue.push_back((ni2, nj2));
                }
            }
        }
    }
}

/// Interpolate a tile-grid value to a pixel position using nearest-neighbor.
/// Matches Python's `tiles_interp(..., method='nearest')`.
fn interp_tile_nearest(
    grid: &[Vec<f64>],
    ni: usize,
    nj: usize,
    row: usize,
    col: usize,
    tile_rows: usize,
    tile_cols: usize,
    _img_rows: usize,
    _img_cols: usize,
) -> f64 {
    // Map pixel to fractional tile coordinate
    let fi = row as f64 / tile_rows as f64;
    let fj = col as f64 / tile_cols as f64;
    // Nearest tile index
    let ti = (fi.round() as usize).min(ni - 1);
    let tj = (fj.round() as usize).min(nj - 1);
    grid[ti][tj]
}

/// Apply atmospheric correction with spatially varying (tiled) parameters.
///
/// For each pixel, looks up the tile AOT, interpolates LUT parameters, and applies correction.
#[cfg(feature = "full-io")]
pub fn dsf_correct_band_tiled(
    rhot: &Array2<f64>,
    lut: &AerosolLut,
    band: &str,
    tt_gas: f64,
    tiled: &TiledDsfResult,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> Array2<f64> {
    let (rows, cols) = rhot.dim();

    // Pre-compute LUT parameters per tile
    let mut romix_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];
    let mut astot_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];
    let mut dutott_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];

    for ti in 0..tiled.ni {
        for tj in 0..tiled.nj {
            let aot = tiled.aot_grid[ti][tj];
            romix_grid[ti][tj] = lut.romix(band, pressure, azi, thv, ths, aot);
            astot_grid[ti][tj] = lut.astot(band, pressure, azi, thv, ths, aot);
            dutott_grid[ti][tj] = lut.dutott(band, pressure, azi, thv, ths, aot);
        }
    }

    Array2::from_shape_fn((rows, cols), |(r, c)| {
        let v = rhot[(r, c)];
        if !v.is_finite() || v <= 0.0 {
            return f64::NAN;
        }

        let romix = interp_tile_nearest(
            &romix_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );
        let astot = interp_tile_nearest(
            &astot_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );
        let dutott = interp_tile_nearest(
            &dutott_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );

        let rhot_gc = v / tt_gas;
        let rhot_noatm = rhot_gc - romix;
        let denom = dutott + astot * rhot_noatm;
        if denom <= 0.0 {
            return f64::NAN;
        }
        rhot_noatm / denom
    })
}

// Keep the simple non-LUT versions for when full-io is not available

/// Invert AOT from observed dark reflectance using a generic (wavelength-indexed) LUT.
/// Convolves romix spectrum with Gaussian RSR for the given band, then interpolates.
#[cfg(feature = "full-io")]
fn invert_aot_generic(
    lut: &GenericAerosolLut,
    center_nm: f64,
    fwhm_nm: f64,
    rhot_dark: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> f64 {
    invert_aot_generic_rsr(lut, center_nm, fwhm_nm, None, rhot_dark, pressure, azi, thv, ths)
}

/// Invert AOT from dark reflectance, optionally using actual sensor RSR.
#[cfg(feature = "full-io")]
fn invert_aot_generic_rsr(
    lut: &GenericAerosolLut,
    center_nm: f64,
    fwhm_nm: f64,
    rsr: Option<(&[f64], &[f64])>, // (wave_nm, response)
    rhot_dark: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> f64 {
    let tau_steps = &lut.meta.tau;
    let romix_vals: Vec<f64> = tau_steps
        .iter()
        .map(|&tau| {
            let spectrum = lut.romix_spectrum(pressure, azi, thv, ths, tau);
            if let Some((rw, rr)) = rsr {
                rsr_convolve_sensor(&lut.wave_um, &spectrum, rw, rr)
            } else {
                rsr_convolve_gauss(&lut.wave_um, &spectrum, center_nm, fwhm_nm)
            }
        })
        .collect();

    if rhot_dark <= romix_vals[0] {
        return f64::NAN;
    }
    if rhot_dark >= romix_vals[romix_vals.len() - 1] {
        return f64::NAN;
    }

    for i in 0..romix_vals.len() - 1 {
        if rhot_dark >= romix_vals[i] && rhot_dark <= romix_vals[i + 1] {
            let t = (rhot_dark - romix_vals[i]) / (romix_vals[i + 1] - romix_vals[i]);
            return tau_steps[i] * (1.0 - t) + tau_steps[i + 1] * t;
        }
    }
    f64::NAN
}

/// Run DSF AOT estimation using generic (hyperspectral) LUTs.
///
/// Same algorithm as `optimize_aot` but uses wavelength-indexed LUTs with RSR convolution.
#[cfg(feature = "full-io")]
pub fn optimize_aot_generic(
    luts: &[GenericAerosolLut],
    dark_spectrum: &[f64], // per-band dark rhot (gas-corrected)
    wavelengths: &[f64],   // per-band center wavelength in nm
    bandwidths: &[f64],    // per-band FWHM in nm
    tt_gas: &[f64],
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
) -> DsfResult {
    let mut best = DsfResult {
        aot: 0.01,
        model_idx: 0,
        model_name: String::new(),
        rmsd: f64::MAX,
        band_aots: HashMap::new(),
    };

    for (li, lut) in luts.iter().enumerate() {
        let mut band_aots: Vec<(usize, f64)> = Vec::new();
        for bi in 0..wavelengths.len() {
            if wavelengths[bi] < config.wave_range.0 || wavelengths[bi] > config.wave_range.1 {
                continue;
            }
            if tt_gas[bi] < config.min_tgas_aot {
                continue;
            }
            if !dark_spectrum[bi].is_finite() || dark_spectrum[bi] <= 0.0 {
                continue;
            }

            let aot = invert_aot_generic(
                lut,
                wavelengths[bi],
                bandwidths[bi],
                dark_spectrum[bi],
                pressure,
                azi,
                thv,
                ths,
            );
            if aot.is_finite() {
                band_aots.push((bi, aot));
            }
        }

        if band_aots.is_empty() {
            continue;
        }
        band_aots.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let selected_aot = match config.aot_compute {
            AotCompute::Min => band_aots[0].1,
            AotCompute::MeanNLowest => {
                let n = config.dsf_nbands.min(band_aots.len());
                band_aots[..n].iter().map(|x| x.1).sum::<f64>() / n as f64
            }
        };

        let n_fit = config.dsf_nbands_fit.min(band_aots.len());
        let mut sum_sq = 0.0;
        let mut ba_map = HashMap::new();
        for &(bi, aot_bi) in &band_aots {
            ba_map.insert(format!("{}", bi), aot_bi);
        }
        for &(bi, _) in &band_aots[..n_fit] {
            let romix_spectrum = lut.romix_spectrum(pressure, azi, thv, ths, selected_aot);
            let romix_model = rsr_convolve_gauss(
                &lut.wave_um,
                &romix_spectrum,
                wavelengths[bi],
                bandwidths[bi],
            );
            sum_sq += (dark_spectrum[bi] - romix_model).powi(2);
        }
        let rmsd = (sum_sq / n_fit as f64).sqrt();

        log::info!(
            "Generic model {} RMSD={:.6} AOT={:.4}",
            lut.name,
            rmsd,
            selected_aot
        );

        if rmsd < best.rmsd {
            best = DsfResult {
                aot: selected_aot,
                model_idx: li,
                model_name: lut.name.clone(),
                rmsd,
                band_aots: ba_map,
            };
        }
    }

    best
}

/// Fixed (whole-scene) DSF using generic LUTs for hyperspectral sensors.
#[cfg(feature = "full-io")]
pub fn optimize_aot_fixed_generic(
    luts: &[GenericAerosolLut],
    toa_gc_bands: &[Array2<f64>],
    wavelengths: &[f64],
    bandwidths: &[f64],
    tt_gas: &[f64],
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
) -> DsfResult {
    let dark_spectrum = estimate_dark_spectrum(toa_gc_bands, &config.dark_method);
    optimize_aot_generic(
        luts,
        &dark_spectrum,
        wavelengths,
        bandwidths,
        tt_gas,
        pressure,
        azi,
        thv,
        ths,
        config,
    )
}

/// Fixed DSF using generic LUTs with actual sensor RSR for convolution.
#[cfg(feature = "full-io")]
pub fn optimize_aot_fixed_sensor_rsr(
    luts: &[GenericAerosolLut],
    toa_gc_bands: &[Array2<f64>],
    wavelengths: &[f64],
    bandwidths: &[f64],
    tt_gas: &[f64],
    band_rsr: &[Option<(&[f64], &[f64])>], // per-band RSR (wave_nm, response)
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
) -> DsfResult {
    let dark_spectrum = estimate_dark_spectrum(toa_gc_bands, &config.dark_method);

    let mut best = DsfResult {
        aot: 0.01, model_idx: 0, model_name: String::new(),
        rmsd: f64::MAX, band_aots: HashMap::new(),
    };

    for (li, lut) in luts.iter().enumerate() {
        let mut band_aots: Vec<(usize, f64)> = Vec::new();
        for bi in 0..wavelengths.len() {
            if wavelengths[bi] < config.wave_range.0 || wavelengths[bi] > config.wave_range.1 { continue; }
            if tt_gas[bi] < config.min_tgas_aot { continue; }
            if !dark_spectrum[bi].is_finite() || dark_spectrum[bi] <= 0.0 { continue; }

            let aot = invert_aot_generic_rsr(
                lut, wavelengths[bi], bandwidths[bi],
                band_rsr.get(bi).and_then(|r| *r),
                dark_spectrum[bi], pressure, azi, thv, ths,
            );
            if aot.is_finite() { band_aots.push((bi, aot)); }
        }
        if band_aots.is_empty() { continue; }
        band_aots.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

        let selected_aot = match config.aot_compute {
            AotCompute::Min => band_aots[0].1,
            AotCompute::MeanNLowest => {
                let n = config.dsf_nbands.min(band_aots.len());
                band_aots[..n].iter().map(|x| x.1).sum::<f64>() / n as f64
            }
        };

        let n_fit = config.dsf_nbands_fit.min(band_aots.len());
        let mut sum_sq = 0.0;
        let mut ba_map = HashMap::new();
        for &(bi, aot_bi) in &band_aots { ba_map.insert(format!("{}", bi), aot_bi); }
        for &(bi, _) in &band_aots[..n_fit] {
            let romix_spectrum = lut.romix_spectrum(pressure, azi, thv, ths, selected_aot);
            let romix_model = if let Some(Some((rw, rr))) = band_rsr.get(bi) {
                rsr_convolve_sensor(&lut.wave_um, &romix_spectrum, rw, rr)
            } else {
                rsr_convolve_gauss(&lut.wave_um, &romix_spectrum, wavelengths[bi], bandwidths[bi])
            };
            sum_sq += (dark_spectrum[bi] - romix_model).powi(2);
        }
        let rmsd = (sum_sq / n_fit as f64).sqrt();
        log::info!("Sensor-RSR model {} RMSD={:.6} AOT={:.4}", lut.name, rmsd, selected_aot);

        if rmsd < best.rmsd {
            best = DsfResult {
                aot: selected_aot, model_idx: li, model_name: lut.name.clone(),
                rmsd, band_aots: ba_map,
            };
        }
    }
    best
}

/// Tiled DSF using generic LUTs for hyperspectral sensors.
#[cfg(feature = "full-io")]
pub fn optimize_aot_tiled_generic(
    luts: &[GenericAerosolLut],
    toa_gc_bands: &[Array2<f64>],
    wavelengths: &[f64],
    bandwidths: &[f64],
    tt_gas: &[f64],
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
    config: &DsfConfig,
    tile_size: (usize, usize),
) -> TiledDsfResult {
    let (rows, cols) = toa_gc_bands[0].dim();
    let ni = (rows + tile_size.0 - 1) / tile_size.0;
    let nj = (cols + tile_size.1 - 1) / tile_size.1;

    // Phase 1: vote on best model across tiles
    let mut model_votes = vec![0usize; luts.len()];
    for ti in 0..ni {
        let r0 = ti * tile_size.0;
        let r1 = ((ti + 1) * tile_size.0).min(rows);
        for tj in 0..nj {
            let c0 = tj * tile_size.1;
            let c1 = ((tj + 1) * tile_size.1).min(cols);

            let tile_bands: Vec<Array2<f64>> = toa_gc_bands
                .iter()
                .map(|b| b.slice(ndarray::s![r0..r1, c0..c1]).to_owned())
                .collect();

            let total_pixels = (r1 - r0) * (c1 - c0);
            let valid_count = tile_bands[0]
                .iter()
                .filter(|v| v.is_finite() && **v > 0.0)
                .count();
            if (valid_count as f64) < total_pixels as f64 * 0.1 {
                continue;
            }

            let dark_spectrum = estimate_dark_spectrum(&tile_bands, &config.dark_method);
            let result = optimize_aot_generic(
                luts,
                &dark_spectrum,
                wavelengths,
                bandwidths,
                tt_gas,
                pressure,
                azi,
                thv,
                ths,
                config,
            );
            if result.rmsd < f64::MAX
                && result.aot.is_finite()
                && result.aot >= 0.001
                && result.aot <= 1.5
            {
                model_votes[result.model_idx] += 1;
            }
        }
    }

    let best_model_idx = if let Some(ref model_name) = config.fixed_model {
        luts.iter()
            .position(|l| l.name.ends_with(model_name))
            .unwrap_or(0)
    } else {
        model_votes
            .iter()
            .enumerate()
            .max_by_key(|(_, &v)| v)
            .map(|(i, _)| i)
            .unwrap_or(0)
    };

    log::info!(
        "Generic tiled: selected model {} (idx={})",
        luts[best_model_idx].name,
        best_model_idx
    );

    // Phase 2: rebuild AOT grid with selected model
    let mut aot_grid = vec![vec![f64::NAN; nj]; ni];
    let single_lut = std::slice::from_ref(&luts[best_model_idx]);

    for ti in 0..ni {
        let r0 = ti * tile_size.0;
        let r1 = ((ti + 1) * tile_size.0).min(rows);
        for tj in 0..nj {
            let c0 = tj * tile_size.1;
            let c1 = ((tj + 1) * tile_size.1).min(cols);

            let tile_bands: Vec<Array2<f64>> = toa_gc_bands
                .iter()
                .map(|b| b.slice(ndarray::s![r0..r1, c0..c1]).to_owned())
                .collect();

            let total_pixels = (r1 - r0) * (c1 - c0);
            let valid_count = tile_bands[0]
                .iter()
                .filter(|v| v.is_finite() && **v > 0.0)
                .count();
            if (valid_count as f64) < total_pixels as f64 * 0.1 {
                continue;
            }

            let dark_spectrum = estimate_dark_spectrum(&tile_bands, &config.dark_method);
            let result = optimize_aot_generic(
                single_lut,
                &dark_spectrum,
                wavelengths,
                bandwidths,
                tt_gas,
                pressure,
                azi,
                thv,
                ths,
                config,
            );
            if result.aot.is_finite() && result.aot >= 0.001 && result.aot <= 1.5 {
                aot_grid[ti][tj] = result.aot;
            }
        }
    }

    fill_nan_nearest(&mut aot_grid, ni, nj);

    TiledDsfResult {
        aot_grid,
        model_idx: best_model_idx,
        model_name: luts[best_model_idx].name.clone(),
        tile_rows: tile_size.0,
        tile_cols: tile_size.1,
        ni,
        nj,
    }
}

/// Apply atmospheric correction to a single band using a generic LUT.
/// Convolves LUT parameters with Gaussian RSR for the band.
#[cfg(feature = "full-io")]
pub fn dsf_correct_band_generic(
    rhot: &Array2<f64>,
    lut: &GenericAerosolLut,
    center_nm: f64,
    fwhm_nm: f64,
    tt_gas: f64,
    aot: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> Array2<f64> {
    dsf_correct_band_generic_rsr(rhot, lut, center_nm, fwhm_nm, None, tt_gas, aot, pressure, azi, thv, ths)
}

/// Apply atmospheric correction using a generic LUT, optionally with actual sensor RSR.
#[cfg(feature = "full-io")]
pub fn dsf_correct_band_generic_rsr(
    rhot: &Array2<f64>,
    lut: &GenericAerosolLut,
    center_nm: f64,
    fwhm_nm: f64,
    rsr: Option<(&[f64], &[f64])>,
    tt_gas: f64,
    aot: f64,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> Array2<f64> {
    let (romix_spec, astot_spec, dutott_spec) = lut.params_spectrum(pressure, azi, thv, ths, aot);
    let (romix, astot, dutott) = if let Some((rw, rr)) = rsr {
        (rsr_convolve_sensor(&lut.wave_um, &romix_spec, rw, rr),
         rsr_convolve_sensor(&lut.wave_um, &astot_spec, rw, rr),
         rsr_convolve_sensor(&lut.wave_um, &dutott_spec, rw, rr))
    } else {
        (rsr_convolve_gauss(&lut.wave_um, &romix_spec, center_nm, fwhm_nm),
         rsr_convolve_gauss(&lut.wave_um, &astot_spec, center_nm, fwhm_nm),
         rsr_convolve_gauss(&lut.wave_um, &dutott_spec, center_nm, fwhm_nm))
    };

    rhot.mapv(|v| {
        if !v.is_finite() || v <= 0.0 { return f64::NAN; }
        let rhot_gc = v / tt_gas;
        let rhot_noatm = rhot_gc - romix;
        let denom = dutott + astot * rhot_noatm;
        if denom <= 0.0 { return f64::NAN; }
        rhot_noatm / denom
    })
}

/// Apply atmospheric correction with tiled AOT using a generic LUT.
#[cfg(feature = "full-io")]
pub fn dsf_correct_band_tiled_generic(
    rhot: &Array2<f64>,
    lut: &GenericAerosolLut,
    center_nm: f64,
    fwhm_nm: f64,
    tt_gas: f64,
    tiled: &TiledDsfResult,
    pressure: f64,
    azi: f64,
    thv: f64,
    ths: f64,
) -> Array2<f64> {
    let (rows, cols) = rhot.dim();

    // Pre-compute LUT parameters per tile
    let mut romix_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];
    let mut astot_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];
    let mut dutott_grid = vec![vec![0.0f64; tiled.nj]; tiled.ni];

    for ti in 0..tiled.ni {
        for tj in 0..tiled.nj {
            let aot = tiled.aot_grid[ti][tj];
            let (romix_spec, astot_spec, dutott_spec) =
                lut.params_spectrum(pressure, azi, thv, ths, aot);
            romix_grid[ti][tj] = rsr_convolve_gauss(&lut.wave_um, &romix_spec, center_nm, fwhm_nm);
            astot_grid[ti][tj] = rsr_convolve_gauss(&lut.wave_um, &astot_spec, center_nm, fwhm_nm);
            dutott_grid[ti][tj] =
                rsr_convolve_gauss(&lut.wave_um, &dutott_spec, center_nm, fwhm_nm);
        }
    }

    Array2::from_shape_fn((rows, cols), |(r, c)| {
        let v = rhot[(r, c)];
        if !v.is_finite() || v <= 0.0 {
            return f64::NAN;
        }

        let romix = interp_tile_nearest(
            &romix_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );
        let astot = interp_tile_nearest(
            &astot_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );
        let dutott = interp_tile_nearest(
            &dutott_grid,
            tiled.ni,
            tiled.nj,
            r,
            c,
            tiled.tile_rows,
            tiled.tile_cols,
            rows,
            cols,
        );

        let rhot_gc = v / tt_gas;
        let rhot_noatm = rhot_gc - romix;
        let denom = dutott + astot * rhot_noatm;
        if denom <= 0.0 {
            return f64::NAN;
        }
        rhot_noatm / denom
    })
}

/// Simplified AOT estimation (no LUT)
pub fn optimize_aot_simple(
    dark_spectrum: &[f64],
    _wavelengths: &[f64],
    _sun_zenith: f64,
    _view_zenith: f64,
) -> f64 {
    let mean_dark = dark_spectrum.iter().sum::<f64>() / dark_spectrum.len() as f64;
    (mean_dark * 10.0).min(1.0).max(0.01)
}

/// Simplified aerosol correction (no LUT)
pub fn dsf_correction_simple(
    toa: &Array2<f64>,
    aot: f64,
    wavelength: f64,
    _sun_zenith: f64,
    _view_zenith: f64,
) -> Array2<f64> {
    let path_reflectance = aot * 0.1 * (550.0 / wavelength).powf(1.3);
    toa.mapv(|v| (v - path_reflectance).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::arr2;

    #[test]
    fn test_estimate_dark_spectrum() {
        let band1 = arr2(&[[0.1, 0.2], [0.3, 0.4]]);
        let band2 = arr2(&[[0.05, 0.15], [0.25, 0.35]]);
        let dark = estimate_dark_spectrum(&[band1, band2], &DarkSpectrumMethod::Percentile(25.0));
        assert_eq!(dark.len(), 2);
        assert!(dark[0] > 0.0);
        assert!(dark[1] > 0.0);
    }
}
