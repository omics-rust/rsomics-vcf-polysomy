//! Copy-number estimation from BAF histograms.
//!
//! Algorithm mirrors bcftools/polysomy.c (MIT).  For each chromosome:
//!   1. Smooth the histogram with a moving average (default window=7).
//!   2. Find irr = minimum in [0, n/2), iaa = minimum in [n/2, n).
//!   3. Compute srr/sra/saa integrals and classify CN=1/-1 or proceed.
//!   4. Fit CN2/CN3/CN4 Gaussian mixture models against normalised smoothed data.
//!   5. Select the best model with a CN-penalty tiebreaker.
//!   6. Write DIST, FIT, and CN records to dist.dat.

use std::io::Write;

use crate::fitting::{Gaussian, fit_mc, sprint_func};
use crate::histogram::{BafHistogram, NBINS};
use rsomics_common::Result;

/// Tuning parameters for polysomy analysis.
#[derive(Clone, Debug)]
pub struct PolysomyArgs {
    pub sample: Option<String>,
    pub fit_th: f64,
    pub cn_penalty: f64,
    pub peak_symmetry: f64,
    pub min_peak_size: f64,
    pub min_fraction: f64,
    pub include_aa: bool,
    pub smooth: i32,
}

impl Default for PolysomyArgs {
    fn default() -> Self {
        Self {
            sample: None,
            fit_th: 3.3,
            cn_penalty: 0.7,
            peak_symmetry: 0.5,
            min_peak_size: 0.1,
            min_fraction: 0.1,
            include_aa: false,
            smooth: 3,
        }
    }
}

/// Result of per-chromosome copy-number fitting.
#[derive(Debug, Clone)]
pub struct CnEstimate {
    pub chrom: String,
    pub cn: f64,
    pub fit: f64,
}

/// Fit models and write DIST / FIT / CN lines for one chromosome.
pub fn fit_and_write(w: &mut impl Write, hist: &BafHistogram, args: &PolysomyArgs) -> Result<()> {
    let n = NBINS;
    let win = (args.smooth.unsigned_abs() as usize) * 2 + 1;

    let smoothed = hist.smooth(args.smooth);
    let xs: Vec<f64> = (0..n).map(BafHistogram::xval).collect();

    // Valley positions: irr separates RR from RA, iaa separates RA from AA.
    let irr = (0..n / 2)
        .min_by(|&a, &b| smoothed[a].partial_cmp(&smoothed[b]).unwrap())
        .unwrap_or(0);
    let iaa = (n / 2..n)
        .min_by(|&a, &b| smoothed[a].partial_cmp(&smoothed[b]).unwrap())
        .unwrap_or(n - 1);

    // bcftools shifts valley indices by win/2 to account for smoothing lag.
    let irr_adj = (irr + win / 2).min(n - 1);
    let iaa_adj = (iaa + win / 2).min(n - 1);

    let irr_adj = irr_adj.min(iaa_adj.saturating_sub(1));
    let iaa_adj = iaa_adj.max(irr_adj + 1).min(n - 1);

    let srr: f64 = smoothed[0..=irr_adj].iter().sum();
    let sra: f64 = if irr_adj < iaa_adj {
        smoothed[irr_adj + 1..=iaa_adj].iter().sum()
    } else {
        0.0
    };
    let saa: f64 = if iaa_adj + 1 < n {
        smoothed[iaa_adj + 1..n].iter().sum()
    } else {
        0.0
    };

    let max_rr = smoothed[0..=irr_adj].iter().cloned().fold(0.0f64, f64::max);
    let max_ra = if irr_adj < iaa_adj {
        smoothed[irr_adj + 1..=iaa_adj]
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
    } else {
        0.0
    };
    let max_aa = if iaa_adj + 1 < n {
        smoothed[iaa_adj + 1..n]
            .iter()
            .cloned()
            .fold(0.0f64, f64::max)
    } else {
        0.0
    };

    let mut ys_norm = smoothed.clone();
    if max_rr > 0.0 {
        for v in ys_norm[0..=irr_adj].iter_mut() {
            *v /= max_rr;
        }
    }
    if max_ra > 0.0 && irr_adj < iaa_adj {
        for v in ys_norm[irr_adj + 1..=iaa_adj].iter_mut() {
            *v /= max_ra;
        }
    }
    if max_aa > 0.0 && iaa_adj + 1 < n {
        for v in ys_norm[iaa_adj + 1..n].iter_mut() {
            *v /= max_aa;
        }
    }

    for (i, &yn) in ys_norm.iter().enumerate() {
        writeln!(
            w,
            "DIST\t{}\t{}\t{}",
            hist.chrom,
            format_f(xs[i]),
            format_f(yn)
        )?;
    }

    // CN=1: no RA peak at all, or RA negligible while AA dominates.
    // CN=-1: RA small or AA large but not both conditions (bcftools polysomy.c).
    if sra == 0.0 || (srr > 0.0 && sra / srr < args.min_fraction && saa / sra > 1.0) {
        writeln!(w, "CN\t{}\t{:.2}\t{}", hist.chrom, 1.0f64, format_e(0.0))?;
        return Ok(());
    }
    if srr > 0.0 && (sra / srr < args.min_fraction || saa / sra > 1.0) {
        writeln!(w, "CN\t{}\t{:.2}\t{}", hist.chrom, -1.0f64, format_e(0.0))?;
        return Ok(());
    }

    // Fit over the normalised RA region [irr_adj, iaa_adj], as bcftools does.
    let ra_xs = xs[irr_adj..=iaa_adj].to_vec();
    let ra_ys = ys_norm[irr_adj..=iaa_adj].to_vec();

    let mut rng: u64 = 0x12345678_9abcdef0
        ^ hist
            .chrom
            .bytes()
            .fold(0u64, |acc, b| acc.wrapping_mul(131).wrapping_add(b as u64));

    // CN2: single Gaussian at 0.5
    let cn2_init = vec![Gaussian {
        params: [1.0, 0.5, 0.05],
        bounds: Some([0.0, 10.0, 0.45, 0.55, 0.01, 0.3]),
    }];
    let (cn2_gs, cn2_fit) = fit_mc(&cn2_init, &ra_xs, &ra_ys, 50, &mut rng);
    writeln!(
        w,
        "FIT\t{}\t{}\t{}\t{}\t{}",
        hist.chrom,
        format_e(cn2_fit),
        irr_adj,
        iaa_adj,
        sprint_func(&cn2_gs)
    )?;

    // CN3: two symmetric Gaussians around 0.5, separation dx ≤ 0.5/3
    let cn3_init = vec![
        Gaussian {
            params: [0.7, 0.33, 0.05],
            bounds: Some([0.0, 10.0, 0.01, 0.49, 0.01, 0.3]),
        },
        Gaussian {
            params: [0.7, 0.67, 0.05],
            bounds: Some([0.0, 10.0, 0.51, 0.99, 0.01, 0.3]),
        },
    ];
    let (mut cn3_gs, cn3_fit) = fit_mc(&cn3_init, &ra_xs, &ra_ys, 50, &mut rng);
    let cn3_dx = (0.5 - cn3_gs[0].params[1]).abs();
    if cn3_dx > 0.5 / 3.0 {
        cn3_gs[0].params[1] = 0.5 - 0.5 / 3.0;
        cn3_gs[1].params[1] = 0.5 + 0.5 / 3.0;
    }
    writeln!(
        w,
        "FIT\t{}\t{}\t{}\t{}\t{}",
        hist.chrom,
        format_e(cn3_fit),
        irr_adj,
        iaa_adj,
        sprint_func(&cn3_gs)
    )?;

    // CN4: three Gaussians (0.5 ± dx plus centre), dx ≤ 0.25
    let cn4_init = vec![
        Gaussian {
            params: [0.5, 0.25, 0.05],
            bounds: Some([0.0, 10.0, 0.01, 0.49, 0.01, 0.3]),
        },
        Gaussian {
            params: [1.0, 0.5, 0.05],
            bounds: Some([0.0, 10.0, 0.45, 0.55, 0.01, 0.3]),
        },
        Gaussian {
            params: [0.5, 0.75, 0.05],
            bounds: Some([0.0, 10.0, 0.51, 0.99, 0.01, 0.3]),
        },
    ];
    let (cn4_gs, cn4_fit) = fit_mc(&cn4_init, &ra_xs, &ra_ys, 50, &mut rng);
    writeln!(
        w,
        "FIT\t{}\t{}\t{}\t{}\t{}",
        hist.chrom,
        format_e(cn4_fit),
        irr_adj,
        iaa_adj,
        sprint_func(&cn4_gs)
    )?;

    let cn3_a0 = cn3_gs[0].params[0].powi(2);
    let cn3_a1 = cn3_gs[1].params[0].powi(2);
    let cn3_sym_ok = if cn3_a0.max(cn3_a1) > 0.0 {
        cn3_a0.min(cn3_a1) / cn3_a0.max(cn3_a1) >= args.peak_symmetry
    } else {
        false
    };

    let cn4_a: Vec<f64> = cn4_gs.iter().map(|g| g.params[0].powi(2)).collect();
    let cn4_min = cn4_a.iter().cloned().fold(f64::INFINITY, f64::min);
    let cn4_max = cn4_a.iter().cloned().fold(0.0f64, f64::max);
    let cn4_peak_ok = cn4_min >= args.min_peak_size;
    let cn4_sym_ok = if cn4_max > 0.0 {
        cn4_min / cn4_max >= args.peak_symmetry
    } else {
        false
    };
    let cn4_asym = (cn4_gs[0].params[1] - cn4_gs[2].params[1]).abs();
    let cn4_asym_ok = cn4_asym <= 0.2;

    let cn2_ok = cn2_fit <= args.fit_th;
    let cn3_ok = cn3_fit <= args.fit_th && cn3_sym_ok;
    let cn4_ok = cn4_fit <= args.fit_th && cn4_peak_ok && cn4_sym_ok && cn4_asym_ok;

    // Accept higher CN only if GOF improves by at least (1 - cn_penalty) factor.
    let (best_cn, best_fit) = if cn2_ok {
        let mut winner_cn = 2.0f64;
        let mut winner_fit = cn2_fit;

        if cn3_ok && cn3_fit < (1.0 - args.cn_penalty) * cn2_fit {
            // CN fraction from peak centre: (1 - 2·b) / b
            let centre = cn3_gs[0].params[1];
            let frac = ((1.0 - 2.0 * centre) / centre).clamp(0.0, 1.0);
            winner_cn = 2.0 + frac;
            winner_fit = cn3_fit;
        }
        if cn4_ok && cn4_fit < (1.0 - args.cn_penalty) * winner_fit {
            let frac = (cn4_gs[2].params[1] - cn4_gs[0].params[1]).clamp(0.0, 1.0);
            winner_cn = 3.0 + frac;
            winner_fit = cn4_fit;
        }

        (winner_cn, winner_fit)
    } else {
        (-1.0, f64::MAX)
    };

    let fit_out = if best_fit == f64::MAX { 0.0 } else { best_fit };
    writeln!(
        w,
        "CN\t{}\t{:.2}\t{}",
        hist.chrom,
        best_cn,
        format_f(fit_out)
    )?;

    Ok(())
}

/// Format a float with %f-like precision (6 decimal places).
fn format_f(v: f64) -> String {
    format!("{v:.6}")
}

/// Format with C-style %e: 6 decimal places, sign-padded 2-digit exponent.
fn format_e(v: f64) -> String {
    if v == 0.0 {
        return "0.000000e+00".to_string();
    }
    let neg = v < 0.0;
    let abs_v = v.abs();
    let exp = abs_v.log10().floor() as i32;
    let mantissa = abs_v / 10.0f64.powi(exp);
    let sign = if neg { "-" } else { "" };
    format!("{sign}{mantissa:.6}e{:+03}", exp)
}
