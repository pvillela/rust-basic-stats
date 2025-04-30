pub struct IterWithCounts<I, V> {
    source: I,
    prev_value: Option<V>,
}

impl<I, V> IterWithCounts<I, V> {
    pub fn new(source: I) -> Self {
        Self {
            source,
            prev_value: None,
        }
    }
}

impl<I, V> Iterator for IterWithCounts<I, V>
where
    V: PartialEq + Clone,
    I: Iterator<Item = V>,
{
    type Item = (V, u64);

    fn next(&mut self) -> Option<Self::Item> {
        let mut count = 1;
        loop {
            let curr_value = self.source.next();
            let prev_value = self.prev_value.clone();
            match (curr_value, prev_value) {
                (Some(v), Some(prev)) if v == prev => count += 1,
                (Some(v), Some(prev)) => {
                    self.prev_value = Some(v);
                    return Some((prev, count));
                }
                (Some(v), None) => {
                    self.prev_value = Some(v);
                }
                (None, Some(prev)) => {
                    self.prev_value = None;
                    return Some((prev, count));
                }
                (None, None) => return None,
            }
        }
    }
}

#[cfg(test)]
mod iter_test {
    use super::*;

    #[test]
    fn iter_test() {
        let dat = [1., 3., 9., 9., 10., 10., 10., 10., 20.];
        let dat_c = IterWithCounts::new(dat.into_iter()).collect::<Vec<_>>();
        let exp_dat_c = vec![(1., 1), (3., 1), (9., 2), (10., 4), (20., 1)];
        assert_eq!(exp_dat_c, dat_c);
    }
}
