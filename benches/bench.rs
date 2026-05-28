use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile::tempdir;

fn bench_vcf_polysomy(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-polysomy");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest.join("tests/golden/diploid_ra_500.vcf");
    let dir = tempdir().unwrap();
    c.bench_function("rsomics-vcf-polysomy golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([vcf.to_str().unwrap(), "-o", dir.path().to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_polysomy);
criterion_main!(benches);
