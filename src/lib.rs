//v is updated in place
//a is length n-1, b length n, c length n-1, v length n
//overwrite is appropriate if the matrix coefficients change (eg, if they depend on x), since they have to be recreated anyway
fn thomas_algorithm_overwrite(a: &[f64], b: &[f64], c: &mut [f64], v: &mut [f64]) {
    let mut solution = vec![0.0; v.len()];
    c[0] = c[0] / b[0];
    v[0] = v[0] / b[0];

    for i in 1..c.len() {
        let den = b[i] - a[i - 1] * c[i - 1];
        c[i] = c[i] / den;
        v[i] = (v[i] - a[i - 1] * v[i - 1]) / den;
    }
    let last_index = v.len() - 1;
    v[last_index] = (v[last_index] - a[last_index - 1] * v[last_index - 1])
        / (b[last_index] - a[last_index - 1] * c[last_index - 1]);

    solution[last_index] = v[last_index];

    for i in (0..last_index).rev() {
        v[i] = v[i] - c[i] * v[i + 1];
    }
}

pub struct Surface {
    pub surface: Vec<Vec<f64>>,
    pub size_t: usize,
    pub size_x: usize,
    pub dx: f64,
    pub dt: f64,
    pub min_t: f64,
    pub min_x: f64,
}
impl Surface {
    pub fn init(size_t: usize, size_x: usize, dt: f64, dx: f64, min_t: f64, min_x: f64) -> Self {
        Surface {
            surface: vec![],
            size_t,
            size_x,
            dx,
            dt,
            min_t,
            min_x,
        }
    }
    //note that this forces a move
    pub fn set_next_x(&mut self, x_at_t: Vec<f64>) {
        self.surface.push(x_at_t);
    }
    pub fn get_last_x(&self) -> Option<&[f64]> {
        self.get_x_at_t(self.surface.len() - 1)
    }
    pub fn get_x_at_t(&self, t_index: usize) -> Option<&[f64]> {
        let curr_length = self.surface.len();
        if t_index < curr_length {
            Some(&self.surface[t_index])
        } else {
            None
        }
    }
    pub fn get_result_at_x_and_t(&self, t_index: usize, x_index: usize) -> Option<(f64, f64, f64)> {
        let x = self.get_x_at_t(t_index)?;
        let curr_length = x.len();
        if x_index < curr_length {
            Some((
                self.min_t + self.dt * (t_index as f64 + 1.0),
                self.min_x + self.dx * (x_index as f64 + 1.0),
                x[x_index],
            ))
        } else {
            None
        }
    }
}
pub struct Boundary {
    min: f64,
    max: f64,
}
pub fn solve_second_order_pde(
    fn_fn: impl Fn(f64, f64) -> f64, //for black scholes, this is a constant (-r)
    fn_first_deriv: impl Fn(f64, f64) -> f64, //for black scholes, this is r*S (or, in logs, a constant r-sigma*sigma/2)
    fn_second_deriv: impl Fn(f64, f64) -> f64, //for black scholes, this is sigma*sigma*S*S/2
    boundary_up: impl Fn(f64, f64) -> f64,    //for black scholes, this is s-k
    boundary_down: impl Fn(f64, f64) -> f64,  //for black scholes, this is 0
    init_fn: impl Fn(f64) -> f64,             //for black scholes, this is (s-k)^+
    size_t: usize,
    t_lim: &Boundary,
    size_x: usize,
    x_lim: &Boundary,
) -> Option<Surface> {
    let dx = (x_lim.max - x_lim.min) / (size_x as f64 + 1.0); //+ two so that boundary is "outside" of mesh
    let dt = (t_lim.max - t_lim.min) / (size_t as f64); //+ one so that boundary is "outside" of mesh
    let mut surface = Surface::init(size_t, size_x, dt, dx, t_lim.min, x_lim.min);
    let mut construct_starting_surface = vec![0.0; size_x]; //does not include boundaries
    for i in 0..size_x {
        construct_starting_surface[i] = init_fn(x_lim.min + dx * (i as f64 + 1.0));
        //start at 1 since 0 is boundary
    }

    let fxminus =
        |t: f64, x: f64| (dt / dx) * (0.5 * fn_first_deriv(t, x) - fn_second_deriv(t, x) / dx);
    let fxplus =
        |t: f64, x: f64| -(dt / dx) * (0.5 * fn_first_deriv(t, x) + fn_second_deriv(t, x) / dx);
    let fx =
        |t: f64, x: f64| 1.0 + 2.0 * fn_second_deriv(t, x) * (dt / (dx * dx)) - fn_fn(t, x) * dt;

    let put_boundarys = |v: &mut [f64], t: f64| {
        v[0] = v[0] - fxminus(t, x_lim.min + dx) * boundary_down(t, x_lim.min); //a at x+1 for the "first" row in mesh
        v[size_x - 1] = v[size_x - 1] - fxplus(t, x_lim.max - dx) * boundary_up(t, x_lim.max);
        //c at x-1 for the "last" row in mesh
    };
    surface.set_next_x(construct_starting_surface);
    //t=0 is the start of the surface, continuing with t=mint+dt
    for j in 1..size_t {
        let t = t_lim.min + dt * (j as f64);
        let a: Vec<f64> =
            (1..size_x) //start at 1 since 0 is boundary
                .map(|index| fxminus(t, x_lim.min + dx * (index as f64 + 1.0)))
                .collect(); //add 1 since the first "a" is actually the second row

        let b: Vec<f64> = (1..size_x + 1) //start at 1 since 0 is boundary, end at size_x since size_x+1 is boundary
            .map(|index| fx(t, x_lim.min + dx * (index as f64)))
            .collect();
        let mut c: Vec<f64> = (1..size_x) //start at 1 since 0 is boundary
            .map(|index| fxplus(t, x_lim.min + dx * (index as f64)))
            .collect();

        let mut v = surface.get_last_x()?.to_vec().clone();
        put_boundarys(&mut v, t);
        thomas_algorithm_overwrite(&a, &b, &mut c, &mut v);
        surface.set_next_x(v);
    }
    Some(surface)
}

#[cfg(test)]
mod tests {

    use super::*;
    extern crate black_scholes;

    #[test]
    fn thomas_algorithm() {
        let a = vec![2.0, -4.0];
        let b = vec![1.0, 2.0, 1.0];
        let mut c = vec![3.0, -1.0];
        let mut d = vec![6.0, 4.0, 0.0];
        /*let result = */
        thomas_algorithm_overwrite(&a, &b, &mut c, &mut d);
        for res in d.iter() {
            println!("res {}", res);
        }
        assert_eq!(d[0], 3.0);
        assert_eq!(d[1], 1.0);
        assert_eq!(d[2], 4.0);
    }

    #[test]
    fn test_pde_on_black_scholes() {
        let stock = 5.0;
        let strike = 4.5;
        let rate = 0.01;
        let sigma = 0.3;
        let maturity = 2.0;
        let t_lim = Boundary {
            min: 0.0,
            max: maturity,
        };
        let x_lim = Boundary {
            min: -10.0,
            max: 10.0,
        };
        let size_t = 700;
        let size_x = 255;
        let result = solve_second_order_pde(
            |_t, _x| -rate,
            |_t, _x| rate - sigma * sigma * 0.5,
            |_t, _x: f64| sigma * sigma * 0.5,
            |_t, x| (stock * x.exp() - strike),
            |_t, _x| 0.0,
            |x| {
                let v = stock * x.exp() - strike;
                if v > 0.0 {
                    v
                } else {
                    0.0
                }
            },
            size_t,
            &t_lim,
            size_x,
            &x_lim,
        )
        .unwrap();
        assert_eq!(result.surface.len(), size_t);
        let (_t, x, approx_bs) = result
            .get_result_at_x_and_t(size_t - 1, ((size_x - 1) / 2) as usize)
            .unwrap();
        assert!(
            (approx_bs - black_scholes::call(stock * x.exp(), strike, rate, sigma, maturity)).abs()
                < 0.01
        );
    }
    #[test]
    fn test_pde_on_full_range_bs() {
        let strike = 4.5;
        let rate = 0.01;
        let sigma = 0.3;
        let maturity = 2.0;
        let t_lim: Boundary = Boundary {
            min: 0.0,
            max: maturity,
        };
        let x_lim = Boundary { min: 0.0, max: 1.0 };
        let size_t = 700;
        let size_x = 255;
        let result = solve_second_order_pde(
            |_t, x| -rate * (1.0 - x),
            |_t, x| rate * (1.0 - x) * x,
            |_t, x: f64| (1.0 - x) * (1.0 - x) * x * x * sigma * sigma * 0.5,
            |_t, x| (2.0 * x - 1.0),
            |_t, _x| 0.0,
            |x| {
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
        //note that this is more accurate than the approximation of the boundary in the previous test
    }
}
