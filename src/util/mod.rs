use num_format::{Locale, ToFormattedString};
use std::any::type_name;
use std::time::Duration;

pub mod color;
pub mod str;

#[inline]
pub fn to_num_pretty(num: &impl ToFormattedString) -> String {
    num.to_formatted_string(&Locale::en)
}

#[inline]
pub fn timeit<T, F: FnOnce() -> T>(f: F) -> (T, Duration) {
    use std::time::Instant;
    let start = Instant::now();
    let result = f();
    (result, start.elapsed())
}

#[inline]
pub fn timeit_pretty<T, F: FnOnce() -> T>(f: F) -> (T, Duration, String) {
    let (res, d) = timeit(f);
    (res, d, format_approx(&d))
}

#[inline]
pub fn timeit_print<T, F: FnOnce() -> T>(label: &str, f: F) -> (T, Duration) {
    let (res, d, s) = timeit_pretty(f);
    println!("{label}: {s}",);
    (res, d)
}

#[allow(unused)]
#[inline]
pub fn print_type_of<T>(_: &T) {
    println!("{}", type_name::<T>());
}

fn format_approx(duration: &Duration) -> String {
    let mut whole = duration.as_nanos();
    if whole < 1 {
        return "< 1 ns".to_string();
    }

    let mut rem = 0;
    // todo: use u as micro?
    let units = ["n", "μ", "m", ""];
    let mut unit = units[0];
    for u in &units[1..] {
        if whole >= 1_000 {
            rem = whole % 1_000;
            whole /= 1_000;
            unit = u;
        } else {
            break;
        }
    }

    let frac: Option<String> = (rem != 0).then_some(
        format!("{rem:03}")
            .trim_end_matches('0')
            .chars()
            .take(2)
            .collect(),
    );

    let whole = to_num_pretty(&whole);
    frac.map_or_else(
        || format!("{whole} {unit}s",),
        |f| format!("{whole}.{f} {unit}s",),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn exact_ns() {
        assert_eq!(format_approx(&Duration::from_nanos(0)), "< 1 ns");
        assert_eq!(format_approx(&Duration::from_nanos(1)), "1 ns");
        assert_eq!(format_approx(&Duration::from_nanos(15)), "15 ns");
        assert_eq!(format_approx(&Duration::from_nanos(999)), "999 ns");
    }

    #[test]
    fn micros_and_millis() {
        assert_eq!(format_approx(&Duration::from_nanos(1_000)), "1 μs");
        assert_eq!(format_approx(&Duration::from_nanos(1_000_000)), "1 ms");
        assert_eq!(format_approx(&Duration::from_nanos(1_234_567)), "1.23 ms");
    }

    #[test]
    fn pad_decimals() {
        assert_eq!(format_approx(&Duration::from_nanos(1_010)), "1.01 μs");
        assert_eq!(format_approx(&Duration::from_nanos(1_100)), "1.1 μs");
    }

    #[test]
    fn drop_last_decimal() {
        assert_eq!(format_approx(&Duration::from_nanos(1_234)), "1.23 μs");
        assert_eq!(format_approx(&Duration::from_nanos(12_345)), "12.34 μs");
    }

    #[test]
    fn seconds() {
        assert_eq!(
            format_approx(&Duration::from_nanos(12_345_678_123)),
            "12.34 s"
        );
        assert_eq!(
            format_approx(&Duration::from_nanos(12_345_678_123_456)),
            "12,345.67 s"
        );

        assert_eq!(
            format_approx(&Duration::from_nanos(1_345_678_123_456_123)),
            "1,345,678.12 s"
        );
        assert_eq!(
            format_approx(&Duration::from_nanos(12_345_678_123_456_123)),
            "12,345,678.12 s"
        );
        assert_eq!(
            format_approx(&Duration::from_nanos(121_345_678_123_456_123)),
            "121,345,678.12 s"
        );
    }
}
