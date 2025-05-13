//! Iterator functionality.

use std::mem::replace;

/// Wraps a source iterator, transforming it into an iterator that yields pairs of type `(V, u64)` such
/// that each pair corresponds to a grouping of the contiguous items from the source iterator that have the
/// same value, where the pair's first component is the value and the pair's second component is the count of
/// items from the source iterator in the grouping.
struct IterWithCounts<I, V> {
    source: I,
    prev_value: Option<V>,
}

impl<I, V> IterWithCounts<I, V> {
    fn new(source: I) -> Self {
        Self {
            source,
            prev_value: None,
        }
    }
}

impl<I, V> Iterator for IterWithCounts<I, V>
where
    V: PartialEq,
    I: Iterator<Item = V>,
{
    type Item = (V, u64);

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 1;
        loop {
            let curr_value = self.source.next();
            match (curr_value, &self.prev_value) {
                (Some(v), Some(prev)) if v == *prev => count += 1,
                (Some(v), Some(_)) => {
                    let ret_v = replace(&mut self.prev_value, Some(v));
                    return Some((ret_v.unwrap_or_else(|| panic!("can't fail")), count));
                }
                (Some(v), None) => {
                    self.prev_value = Some(v);
                }
                (None, Some(_)) => {
                    let ret_v = self.prev_value.take();
                    return Some((ret_v.unwrap_or_else(|| panic!("can't fail")), count));
                }
                (None, None) => return None,
            }
        }
    }
}

/// Transforms an iterator of values into an iterator of value-count pairs.
///
/// The pairs of type `(V, u64)` are such that:
/// - Each pair corresponds to a grouping of the contiguous items from `source` that have the
///   same value.
/// - The pair's first component is the value from `source`.
/// - The pair's second component is the count of items from `source` in the grouping.
pub fn iter_with_counts<V: PartialEq>(
    source: impl Iterator<Item = V>,
) -> impl Iterator<Item = (V, u64)> {
    IterWithCounts::new(source)
}

#[cfg(test)]
mod iter_test {
    use super::*;

    #[test]
    fn iter_test() {
        let dat = [1., 3., 10., 10., 10., 9., 9., 10., 20.];
        let dat_c = iter_with_counts(dat.into_iter()).collect::<Vec<_>>();
        let exp_dat_c = vec![(1., 1), (3., 1), (10., 3), (9., 2), (10., 1), (20., 1)];
        assert_eq!(exp_dat_c, dat_c);
    }
}
