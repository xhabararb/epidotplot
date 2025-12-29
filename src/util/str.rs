use std::borrow::Cow;

pub fn truncate_ellipsis(s: &str, max_len: usize, count_ellipsis: bool) -> Cow<'_, str> {
    const ELLIPSIS: &str = "...";
    if s.chars().count() <= max_len {
        return Cow::Borrowed(s);
    }

    let use_ellipsis = max_len > ELLIPSIS.len();
    let n = if count_ellipsis && use_ellipsis {
        max_len - ELLIPSIS.len()
    } else {
        max_len
    };

    s.char_indices().nth(n).map_or_else(
        || Cow::Owned(s.to_owned()),
        |(i, _)| {
            if use_ellipsis {
                Cow::Owned(format!("{}{}", &s[..i], ELLIPSIS))
            } else {
                Cow::Borrowed(&s[..i])
            }
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ascii_truncation_with_ellipsis() {
        assert_eq!(truncate_ellipsis("abcdef", 6, true), "abcdef");
        assert_eq!(truncate_ellipsis("abcdef", 6, false), "abcdef");
        assert_eq!(truncate_ellipsis("abcdef", 3, true), "abc");
        assert_eq!(truncate_ellipsis("abcdef", 2, true), "ab");
    }
}
