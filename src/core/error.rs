//! Errors returned by functions in this library.

use std::{
    borrow::Cow,
    error::Error,
    fmt::{Debug, Display},
};

/// Alias for results in this library.
pub(crate) type StatsResult<V> = Result<V, StatsError>;

/// Used for errors from functions in this library.
#[derive(Debug)]
pub struct StatsError {
    msg: Cow<'static, str>,
}

impl StatsError {
    /// Instantiates Self. `msg` can be either a `String` or a `&'static str`.
    pub fn new(msg: impl Into<Cow<'static, str>>) -> Self {
        Self { msg: msg.into() }
    }

    pub fn msg(&self) -> &Cow<'static, str> {
        &self.msg
    }
}

impl Display for StatsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self, f)
    }
}

impl Error for StatsError {}

#[allow(unused)]
pub(crate) trait AsStatsResult<V> {
    fn stats_result(self, msg: &'static str) -> StatsResult<V>;
}

impl<V, E> AsStatsResult<V> for Result<V, E>
where
    E: Error,
{
    fn stats_result(self: Result<V, E>, msg: &'static str) -> StatsResult<V> {
        self.map_err(|_| StatsError::new(msg))
    }
}

#[cfg(test)]
mod test {
    use crate::core::StatsError;

    #[test]
    fn test_stats_error() {
        let e1 = StatsError::new("foo");
        let e2 = StatsError::new("bar".to_string());

        let msg1 = e1.msg();
        assert_eq!(msg1, "foo");

        let msg2 = e2.msg();
        assert_eq!(msg2, "bar");
    }
}
