use crate::types::Primer;

const COMPLEMENT_TABLE: [u8; 128] = {
    let mut table = [0u8; 128];
    // Initialize identity mapping
    let mut i = 0;
    while i < 128 {
        table[i] = i as u8;
        i += 1;
    }
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'R' as usize] = b'Y';
    table[b'Y' as usize] = b'R';
    table[b'K' as usize] = b'M';
    table[b'M' as usize] = b'K';
    table[b'N' as usize] = b'N';
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'r' as usize] = b'y';
    table[b'y' as usize] = b'r';
    table[b'k' as usize] = b'm';
    table[b'm' as usize] = b'k';
    table[b'n' as usize] = b'n';
    table
};

pub fn reverse_complement(seq: &str) -> String {
    let bytes = seq.as_bytes();
    let mut result = Vec::with_capacity(bytes.len());
    for &b in bytes.iter().rev() {
        if (b as usize) < 128 {
            result.push(COMPLEMENT_TABLE[b as usize]);
        } else {
            result.push(b);
        }
    }
    String::from_utf8(result).unwrap_or_default()
}

pub fn check_primer_match(seq: &str, primer: &Primer) -> bool {
    (seq.starts_with(&primer.f) && seq.ends_with(&primer.rc_r))
        || (seq.starts_with(&primer.r) && seq.ends_with(&primer.rc_f))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AATT"), "AATT");
        assert_eq!(reverse_complement("RYKM"), "KMRY");
    }

    #[test]
    fn test_primer_match() {
        let primer = Primer::new("test".into(), "ATCG".into(), "GCTA".into());
        // seq starts with f=ATCG and ends with rc_r=TAGC
        assert!(check_primer_match("ATCGNNNNTAGC", &primer));
        // no match
        assert!(!check_primer_match("GGGGNNNNCCCC", &primer));
    }
}
