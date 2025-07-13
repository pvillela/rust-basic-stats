use super::{StatsError, StatsResult};

pub fn check_alpha_in_open_0_1(alpha: f64) -> StatsResult<()> {
    if 0.0 < alpha && alpha < 1.0 {
        return Ok(());
    }
    Err(StatsError::new("arg `alpha` must be in interval (0, 1)"))
}
