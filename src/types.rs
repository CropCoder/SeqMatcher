use std::collections::HashMap;

use crate::bio;

#[derive(Debug, Clone)]
pub struct Primer {
    pub id: String,
    pub f: String,
    pub r: String,
    pub rc_f: String,
    pub rc_r: String,
}

impl Primer {
    pub fn new(id: String, f: String, r: String) -> Self {
        let rc_f = bio::reverse_complement(&f);
        let rc_r = bio::reverse_complement(&r);
        Self {
            id,
            f: f.to_uppercase(),
            r: r.to_uppercase(),
            rc_f,
            rc_r,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Variant {
    pub raw: String,
    pub rc: String,
}

impl Variant {
    pub fn new(seq: &str) -> Self {
        let raw = seq.to_uppercase();
        let rc = bio::reverse_complement(&raw);
        Self { raw, rc }
    }
}

#[derive(Default, Clone)]
pub struct ThreadResult {
    pub primer_counts: HashMap<String, usize>,
    pub variant_counts: HashMap<String, HashMap<usize, usize>>,
}

impl ThreadResult {
    pub fn merge(&mut self, other: &ThreadResult) {
        for (k, v) in &other.primer_counts {
            *self.primer_counts.entry(k.clone()).or_default() += v;
        }
        for (p_id, vars) in &other.variant_counts {
            let entry = self.variant_counts.entry(p_id.clone()).or_default();
            for (v_idx, count) in vars {
                *entry.entry(*v_idx).or_default() += count;
            }
        }
    }
}
