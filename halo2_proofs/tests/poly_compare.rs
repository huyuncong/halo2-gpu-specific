use ff::Field;
use halo2_proofs::poly::EvaluationDomain;
use halo2_proofs::poly::{Coeff, Polynomial};
use pairing::bn256::Fr;
use rand::thread_rng;
use std::time::Instant;

type CoeffPoly = Polynomial<Fr, Coeff>;
// Generates a CoeffPoly with random Fr elements.
fn generate_random_poly(size: usize) -> CoeffPoly {
    let mut rng = thread_rng();
    let values: Vec<Fr> = (0..size).map(|_| Fr::random(&mut rng)).collect();
    Polynomial::new(values)
}

// Elementwise multiplication for ExtendedLagrangeCoeff polynomials.

#[test]
fn test() {
    let k = 22;
    let size = 1 << k;

    // Generate random polynomials a, b, c, d.
    let a = generate_random_poly(size);
    let b = generate_random_poly(size);
    let c = generate_random_poly(size);
    let d = generate_random_poly(size);

    // Clone polynomials for timing to avoid clone overhead inside measurement.
    let a1 = a.clone();
    let b1 = b.clone();
    let c1 = c.clone();
    let d1 = d.clone();

    let a2 = a.clone();
    let b2 = b.clone();
    let c2 = c.clone();
    let d2 = d.clone();

    // --- Method 1 ---
    let domain1 = EvaluationDomain::<Fr>::new(4, k);
    let start1 = Instant::now();

    let prod_ext = domain1.coeff_to_extended(a1)
        * domain1.coeff_to_extended(b1)
        * domain1.coeff_to_extended(c1)
        * domain1.coeff_to_extended(d1);
    let _h_poly_method1 = domain1.extended_to_coeff(domain1.divide_by_vanishing_poly(prod_ext));
    let duration1 = start1.elapsed();

    let domain2 = EvaluationDomain::<Fr>::new(3, k);
    let domain3 = EvaluationDomain::<Fr>::new(3, k + 1);
    let start2 = Instant::now();

    let poly_ab =
        domain2.extended_to_coeff(domain2.coeff_to_extended(a2) * domain2.coeff_to_extended(b2));
    let poly_ab = domain3.coeff_from_vec(poly_ab);
    let poly_cd =
        domain2.extended_to_coeff(domain2.coeff_to_extended(c2) * domain2.coeff_to_extended(d2));
    let poly_cd = domain3.coeff_from_vec(poly_cd);
    let combined_coeff = domain3
        .extended_to_coeff(domain3.coeff_to_extended(poly_ab) * domain3.coeff_to_extended(poly_cd));
    let mut final_coeffs = combined_coeff;
    long_division(&mut final_coeffs, domain3.n as usize);
    let duration2 = start2.elapsed();

    println!("Method 1 took: {:?}", duration1);
    println!("Method 2 took: {:?}", duration2);
}

fn long_division(poly: &mut Vec<Fr>, n: usize) {
    let m = poly.len();
    for i in (n..m).rev() {
        let factor = poly[i];
        poly[i - n] = poly[i - n] + factor;
        poly[i] = Fr::zero();
    }
    poly.truncate(m - n);
}
