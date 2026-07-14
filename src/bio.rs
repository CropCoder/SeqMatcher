//! Biological sequence algorithms: reverse complement computation
//! and primer match detection for DNA/RNA sequences.
//!
//! Supports all standard IUPAC nucleotide codes including degenerate bases.

use crate::types::Primer;

/// Compile-time reverse complement lookup table (128 bytes).
///
/// Maps each ASCII byte to its complement nucleotide:
/// - Standard bases: A↔T, C↔G
/// - IUPAC degenerate: R↔Y, K↔M, S↔S, W↔W, B↔V, D↔H, H↔D, V↔B
/// - N (any base) → N (self-complementary)
/// - All other bytes → identity (pass-through)
///
/// Both upper and lower case are supported. Non-ASCII bytes (>127)
/// are handled by the caller (passed through unchanged).
const COMPLEMENT_TABLE: [u8; 128] = {
    let mut table = [0u8; 128];
    // Initialize identity mapping for all bytes
    let mut i = 0;
    while i < 128 {
        table[i] = i as u8;
        i += 1;
    }
    // Standard bases
    table[b'A' as usize] = b'T';
    table[b'T' as usize] = b'A';
    table[b'C' as usize] = b'G';
    table[b'G' as usize] = b'C';
    // IUPAC degenerate bases (purine/pyrimidine families)
    table[b'R' as usize] = b'Y'; // A|G → T|C
    table[b'Y' as usize] = b'R'; // T|C → A|G
    table[b'K' as usize] = b'M'; // G|T → C|A
    table[b'M' as usize] = b'K'; // C|A → G|T
    table[b'S' as usize] = b'S'; // G|C → C|G (self-complementary)
    table[b'W' as usize] = b'W'; // A|T → T|A (self-complementary)
    table[b'B' as usize] = b'V'; // C|G|T → G|C|A
    table[b'D' as usize] = b'H'; // A|G|T → T|C|A
    table[b'H' as usize] = b'D'; // A|C|T → T|G|A
    table[b'V' as usize] = b'B'; // A|C|G → T|G|C
    table[b'N' as usize] = b'N'; // A|C|G|T → any (self-complementary)
                                 // Lowercase versions
    table[b'a' as usize] = b't';
    table[b't' as usize] = b'a';
    table[b'c' as usize] = b'g';
    table[b'g' as usize] = b'c';
    table[b'r' as usize] = b'y';
    table[b'y' as usize] = b'r';
    table[b'k' as usize] = b'm';
    table[b'm' as usize] = b'k';
    table[b's' as usize] = b's';
    table[b'w' as usize] = b'w';
    table[b'b' as usize] = b'v';
    table[b'd' as usize] = b'h';
    table[b'h' as usize] = b'd';
    table[b'v' as usize] = b'b';
    table[b'n' as usize] = b'n';
    table
};

/// Compute the reverse complement of a DNA sequence.
///
/// Iterates bytes in reverse order, mapping each nucleotide to its
/// complement via a compile-time lookup table. Non-ASCII bytes pass
/// through unchanged. Returns an empty string on invalid UTF-8
/// (should not occur for valid DNA input).
///
/// # Performance
///
/// O(n) time and memory — one reverse pass, one allocation.
/// The LUT lookup is a single indexed load with no branching.
///
/// # Examples
///
/// ```
/// use seq_matcher::bio::reverse_complement;
/// assert_eq!(reverse_complement("ATCG"), "CGAT");
/// assert_eq!(reverse_complement("SSWW"), "WWSS"); // S↔S, W↔W, order reversed
/// assert_eq!(reverse_complement("BDHV"), "BDHV"); // B↔V, D↔H, full round-trip
/// ```
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

/// Check whether a sequence matches a primer pair using anchored prefix/suffix matching.
///
/// A sequence matches if either:
/// - It starts with the forward primer and ends with the reverse-complement of the reverse primer
/// - It starts with the reverse primer and ends with the reverse-complement of the forward primer
///
/// This models PCR amplification where primers must bind at the exact termini of the amplicon.
/// Uses `starts_with`/`ends_with` for exact terminal matching — no internal or partial binding.
pub fn check_primer_match(seq: &str, primer: &Primer) -> bool {
    (seq.starts_with(&primer.f) && seq.ends_with(&primer.rc_r))
        || (seq.starts_with(&primer.r) && seq.ends_with(&primer.rc_f))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement_standard() {
        assert_eq!(reverse_complement("ATCG"), "CGAT");
        assert_eq!(reverse_complement("AATT"), "AATT");
        assert_eq!(reverse_complement("GGCC"), "GGCC");
    }

    #[test]
    fn test_reverse_complement_iupac_degenerate() {
        // Purine/pyrimidine families
        assert_eq!(reverse_complement("R"), "Y"); // A|G → T|C
        assert_eq!(reverse_complement("K"), "M"); // G|T → C|A

        // Self-complementary degenerate bases
        assert_eq!(reverse_complement("S"), "S"); // G|C → C|G
        assert_eq!(reverse_complement("W"), "W"); // A|T → T|A

        // Three-base degenerate families
        assert_eq!(reverse_complement("B"), "V"); // C|G|T → G|C|A
        assert_eq!(reverse_complement("D"), "H"); // A|G|T → T|C|A
        assert_eq!(reverse_complement("H"), "D"); // A|C|T → T|G|A
        assert_eq!(reverse_complement("V"), "B"); // A|C|G → T|G|C

        // Round-trip: complement twice = identity
        let deg = "RYKMSWBDHVN";
        assert_eq!(reverse_complement(&reverse_complement(deg)), deg);

        // Lowercase
        assert_eq!(reverse_complement("atcg"), "cgat");
        assert_eq!(reverse_complement("b"), "v");
    }

    #[test]
    fn test_reverse_complement_mixed_case() {
        assert_eq!(reverse_complement("AtCg"), "cGaT");
        assert_eq!(reverse_complement("aTcG"), "CgAt");
    }

    #[test]
    fn test_reverse_complement_empty() {
        assert_eq!(reverse_complement(""), "");
    }

    #[test]
    fn test_reverse_complement_non_iupac() {
        // Non-IUPAC characters pass through (identity)
        assert_eq!(reverse_complement("ATXG"), "CXAT");
        assert_eq!(reverse_complement("123"), "321");
    }

    #[test]
    fn test_primer_match_forward_first() {
        let primer = Primer::new("test".into(), "ATCG".into(), "GCTA".into());
        // seq starts with f=ATCG and ends with rc_r=TAGC
        assert!(check_primer_match("ATCGNNNNTAGC", &primer));
    }

    #[test]
    fn test_primer_match_reverse_first() {
        let primer = Primer::new("test".into(), "ATCG".into(), "GCTA".into());
        // f=ATCG, r=GCTA → rc_f=CGAT, rc_r=TAGC
        // seq starts with r=GCTA and ends with rc_f=CGAT
        assert!(check_primer_match("GCTANNNNCGAT", &primer));
    }

    #[test]
    fn test_primer_match_no_match() {
        let primer = Primer::new("test".into(), "ATCG".into(), "GCTA".into());
        // neither end matches
        assert!(!check_primer_match("GGGGNNNNCCCC", &primer));
    }

    #[test]
    fn test_primer_match_partial() {
        // partial prefix only — must match BOTH ends
        let primer = Primer::new("test".into(), "ATCG".into(), "GCTA".into());
        assert!(!check_primer_match("ATCGNNNNCCCC", &primer));
    }
}
