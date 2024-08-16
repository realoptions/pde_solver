## Option Price PDE

This is a basic library for solving 1-dimensional PDEs of the form `f_t = alpha(x, t) f_x + sigma(x, t) f_xx + rho(x, t) f`.  It uses finite differences and implicit solutions using the Thomas algorithm.  It is tested on the "basic" Black Scholes equation (with constant coefficients via the transformation `x = log(S)`) and with a transformation that maps the infinite Black Scholes domain to a finite domain (via the transformation `x = S/(S+K)`)


## Example

```rust
extern crate black_scholes;
extern crate option_price_pde;
let strike = 4.5;
let rate = 0.01;
let sigma = 0.3;
let maturity = 2.0;
let t_lim: solve_second_order_pde::Boundary = solve_second_order_pde::Boundary {
    min: 0.0,
    max: maturity,
};
let x_lim = solve_second_order_pde::Boundary { min: 0.0, max: 1.0 };
let size_t = 700;
let size_x = 255;
let result = solve_second_order_pde::solve_second_order_pde(
    |_t, x| -rate * (1.0 - x),                                        // coefficient on f(x)
    |_t, x| rate * (1.0 - x) * x,                                     // coefficient on f'(x)
    |_t, x: f64| (1.0 - x) * (1.0 - x) * x * x * sigma * sigma * 0.5, // coefficient on f''(x)
    |_t, x| (2.0 * x - 1.0),                                          // upper boundary
    |_t, _x| 0.0,                                                     // lower boundary
    |x| {                                                             // time boundary (t=0)
        let v = 2.0 * x - 1.0;
        if v > 0.0 {
            v
        } else {
            0.0
        }
    },
    700,
    &t_lim,
    255,
    &x_lim,
)
.unwrap();
let (_t, x, v) = result
    .get_result_at_x_and_t(size_t - 1, ((size_x - 1) / 2) as usize)
    .unwrap();
let s = (x * strike) / (1.0 - x);
let approx_bs = (s + strike) * v;
assert!((approx_bs - black_scholes::call(s, strike, rate, sigma, maturity)).abs() < 0.001);
```