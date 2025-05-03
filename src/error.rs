//! Errors returned by functions in this library.

use std::{
    error::Error,
    fmt::{Debug, Display},
};

/// Alias for results in this library.
pub(crate) type StatsResult<V> = Result<V, StatsError>;

/// Signals an error in a function in this library.
#[derive(Debug)]
pub struct StatsError(pub &'static str);

impl Display for StatsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self, f)
    }
}

impl Error for StatsError {}

pub(crate) trait AsStatsResult<V> {
    fn stats_result(self, msg: &'static str) -> StatsResult<V>;
}

impl<V, E> AsStatsResult<V> for Result<V, E>
where
    E: Error,
{
    fn stats_result(self: Result<V, E>, msg: &'static str) -> StatsResult<V> {
        self.map_err(|_| StatsError(msg))
    }
}
