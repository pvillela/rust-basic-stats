use super::{StatsError, StatsResult};

pub(crate) fn check_p0_in_open_0_1(p0: f64) -> StatsResult<()> {
    if (0.0..1.0).contains(&p0) && p0 > 0. {
        return Ok(());
    }
    Err(StatsError("arg `p0` must be in interval (0, 1)"))
}

pub(crate) fn check_alpha_in_closed_0_1(alpha: f64) -> StatsResult<()> {
    if (0.0..=1.0).contains(&alpha) {
        return Ok(());
    }
    Err(StatsError("arg `alpha` must be in interval [0, 1]"))
}

pub(crate) fn check_alpha_in_open_0_1(alpha: f64) -> StatsResult<()> {
    if (0.0..1.0).contains(&alpha) && alpha > 0. {
        return Ok(());
    }
    Err(StatsError("arg `alpha` must be in interval (0, 1)"))
}
