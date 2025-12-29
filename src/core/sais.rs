use crate::error::CoreError;
use libsais::context::Context;
use libsais::{LibsaisError, SuffixArrayConstruction, ThreadCount};
use std::borrow::Cow;
use std::ops::RangeInclusive;

pub fn build_sa_lcp(
    text: &[u8],
    parallel: bool,
    _println: impl Fn(Cow<'_, str>),
    eprintln: impl Fn(Cow<'_, str>),
) -> Result<(Vec<u32>, Vec<u32>), CoreError> {
    let map_libsais_err = |step: &str, err: LibsaisError| {
        eprintln(
            format!("parallel {step} construction failed, falling back to single-threaded: {err}")
                .into(),
        );
    };

    let mt = || -> Result<_, _> {
        let thread_count = ThreadCount::openmp_default();
        let mut ctx = Context::<u8, i32, _>::new_multi_threaded(thread_count);

        SuffixArrayConstruction::for_text(text)
            .in_owned_buffer32()
            .multi_threaded(thread_count)
            .with_context(&mut ctx)
            .run()
            .map_err(|err| map_libsais_err("suffix array", err))?
            .plcp_construction()
            .multi_threaded(thread_count)
            .run()
            .map_err(|err| map_libsais_err("PLCP", err))?
            .lcp_construction()
            .multi_threaded(thread_count)
            .run()
            .map_err(|err| map_libsais_err("LCP", err))
    };

    let st = || -> Result<_, _> {
        SuffixArrayConstruction::for_text(text)
            .in_owned_buffer32()
            .single_threaded()
            .run()?
            .plcp_construction()
            .single_threaded()
            .run()?
            .lcp_construction()
            .single_threaded()
            .run()
    };

    let (sa, lcp, _, _) = if parallel {
        match mt() {
            Ok(sa) => sa.into_parts(),
            Err(()) => st()
                .map_err(|err| CoreError::Other {
                    err: err.into(),
                    msg: Some("single-threaded fallback failed too".to_owned()),
                })?
                .into_parts(),
        }
    } else {
        st().map_err(|err| CoreError::Other {
            err: err.into(),
            msg: Some("suffix array, LCP construction failed in single-threaded mode".to_owned()),
        })?
        .into_parts()
    };

    let sa: Vec<_> = sa.into_iter().map(|v| v as u32).collect();
    let lcp = lcp.into_iter().map(|v| v as u32).collect();

    Ok((sa, lcp))
}

fn binary_search<T, F: Fn(&T) -> bool>(x: &[T], cmp: F) -> usize {
    let mut start = 0;
    let mut m = x.len();
    while m > 0 {
        let mid = start + m / 2;
        let this = &x[mid];
        if cmp(this) {
            start = mid + 1;
            m -= m / 2 + 1;
        } else {
            m /= 2;
        }
    }
    start
}

pub fn find_kmer_range(
    text: &[u8],
    sa: &[u32],
    lcp: &[u32],
    kmer: &[u8],
    contains_n: fn(&[u8]) -> bool,
) -> Option<RangeInclusive<usize>> {
    if contains_n(kmer) {
        return None;
    }

    let lower = binary_search(sa, |&i| &text[i as usize..] < kmer);
    if lower >= sa.len() || !&text[sa[lower] as usize..].starts_with(kmer) {
        return None;
    }

    let count = lcp[(lower + 1).min(lcp.len())..]
        .iter()
        .take_while(|&&t| t as usize >= kmer.len())
        .count();

    Some(lower..=lower + count)
}

#[test]
fn bs_prefix_lower_bound_cases() {
    let cases: &[(&[&str], &str, usize)] = &[
        (&["b", "c", "d"], "a", 0),
        (&["b", "c", "d", "e"], "a", 0),
        (&["b", "d"], "c", 1),
        (&["b", "c", "d"], "e", 3),
        (&["b", "b", "c", "d"], "e", 4),
        (&["b", "d", "e"], "c", 1),
        (&["b", "bb", "d", "e"], "c", 2),
        (&["a", "b", "c", "d"], "a", 0),
        (&["a", "b", "b", "c", "d"], "a", 0),
        (&["a", "a", "b", "c", "d"], "a", 0),
        (&["a", "aa", "aa", "b", "c", "d"], "aa", 1),
        (&["a", "aa", "aa", "b", "c", "d"], "d", 5),
        (&["a", "aa", "aa", "b", "c", "d", "d"], "d", 5),
    ];

    for &(v, needle, expect) in cases {
        for i in 0..v.len().saturating_sub(1) {
            assert!(
                v[i] <= v[i + 1],
                "input slice is not sorted: v=[{v:?}], w[{i}] </= w[{}]",
                i + 1
            );
        }

        let idx = binary_search(v, |&s| s < needle);
        assert_eq!(idx, expect, "v={:?} needle={}", v, needle);
    }
}
