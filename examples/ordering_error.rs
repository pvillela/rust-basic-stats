use basic_stats::error::OrderingError;

/// Simplistic function that adds the products of ranks and values for a data set represented as
/// an iterator sorted in ascending order.
///
/// # Errors:
/// - Returns [`OrderingError`] if `data_set` is not in ascending order.
fn rank_value_sum_prod(data_set: impl Iterator<Item = f64>) -> Result<f64, OrderingError> {
    let mut sum_prod = 0.;
    let mut last_value = f64::MIN;

    for (i, value) in data_set.enumerate() {
        if value < last_value {
            return Err(OrderingError);
        }
        sum_prod += (i + 1) as f64 * value;
        last_value = value;
    }

    Ok(sum_prod)
}

fn main() {
    {
        let data_set = [-1.5, 2.0, 3.0, 5.25].into_iter();
        let good = rank_value_sum_prod(data_set);
        println!("good={good:?}");
        // good=Ok(32.5)
    }

    {
        let data_set = [-1.5, 3.0, 2.0, 5.25].into_iter();
        let bad = rank_value_sum_prod(data_set);
        println!("bad={bad:?}");
        // bad=Err(OrderingError)
    }
}
