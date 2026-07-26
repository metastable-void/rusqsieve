use super::*;

pub(super) fn extract(ctx: &Context, columns: &[Column]) -> Result<Natural, EngineError> {
    let matrix_cols: Vec<Vec<u32>> = columns
        .iter()
        .map(|c| {
            let mut v = Vec::new();
            if c.sign {
                v.push(0)
            }
            for &(i, e) in &c.powers {
                if e & 1 != 0 {
                    v.push(i + 1)
                }
            }
            v
        })
        .collect();
    let matrix = SparseBinaryMatrix::from_columns(ctx.base.len() + 1, &matrix_cols)
        .map_err(|_| EngineError::InvalidDependency)?;
    let dependencies = matrix
        .filtered_dependencies()
        .map_err(|_| EngineError::ResourceLimit)?;
    for dep in dependencies.iter() {
        if !matrix.verify_dependency(dep) {
            return Err(EngineError::InvalidDependency);
        }
        let mut x = Natural::ONE;
        let mut y = Natural::ONE;
        let mut sums = vec![0u32; ctx.base.len()];
        for (j, c) in columns.iter().enumerate() {
            if (dep[j / 64] >> (j % 64)) & 1 == 0 {
                continue;
            }
            x = x.mul_mod(&c.root, &ctx.n);
            for &lp in &c.extra_sqrt {
                y = y.mul_mod(&Natural::from_u64(lp), &ctx.n);
            }
            for &(i, e) in &c.powers {
                sums[i as usize] += e
            }
        }
        for (e, &s) in ctx.base.iter().zip(&sums) {
            for _ in 0..s / 2 {
                y = y.mul_mod(&Natural::from_u64(e.prime as u64), &ctx.n)
            }
        }
        let d = if x >= y {
            x.wrapping_sub(&y)
        } else {
            y.wrapping_sub(&x)
        };
        let g = d.gcd(&ctx.n);
        if !g.is_one() && g != ctx.n {
            return Ok(g);
        }
        let g = x.add_mod(&y, &ctx.n).gcd(&ctx.n);
        if !g.is_one() && g != ctx.n {
            return Ok(g);
        }
    }
    Err(EngineError::NoFactor)
}
