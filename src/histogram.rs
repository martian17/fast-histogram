use plotly::{Plot, Scatter};
use plotly::common::{Mode, Line, color::Rgb, Visible};
use std::iter::Cycle;
use std::path::Path;
use crate::utils::serve_once;

pub struct HistogramData<'a> {
    pub values: Vec<f64>,
    pub timespan: u64, // picoseconds
    pub out_dir: &'a Path,
}

// Rust version of Plotly doesn't provide `express` subpackage, which provides color schemes
// So I yoinked the colors from there and using it here
// original import in python: px.colors.qualitative.Plotly
// color scheme location: https://github.com/plotly/plotly.py/blob/main/_plotly_utils/colors/qualitative.py
fn get_quantative_color_scheme() -> Cycle<impl Iterator<Item = Rgb> + Clone> {
    fn rgb_hex(hex: &str) -> Rgb {
        let r = u8::from_str_radix(&hex[1..3], 16).unwrap();
        let g = u8::from_str_radix(&hex[3..5], 16).unwrap();
        let b = u8::from_str_radix(&hex[5..7], 16).unwrap();
        Rgb::new(r, g, b)
    }
    let color_scheme: [Rgb; _] = [
        rgb_hex("#636EFA"),
        rgb_hex("#EF553B"),
        rgb_hex("#00CC96"),
        rgb_hex("#AB63FA"),
        rgb_hex("#FFA15A"),
        rgb_hex("#19D3F3"),
        rgb_hex("#FF6692"),
        rgb_hex("#B6E880"),
        rgb_hex("#FF97FF"),
        rgb_hex("#FECB52"),
    ];
    return color_scheme.into_iter().cycle();
}

pub fn display_histogram (data: HistogramData) {
    let mut color_cycle = get_quantative_color_scheme();
    let mut plot = Plot::new();
    let ch_ref = 8;
    let ch_sync = 3;
    let dist = 0.0;
    let trace = Scatter::new(data.values.iter().enumerate().map(|(i, _)| i as f64).collect(), data.values)
        .mode(Mode::Lines)
        .name(format!("Ch: {ch_ref}, Dist: {dist}, Sync: {ch_sync}"))
        .visible(Visible::True)
        .line(Line::new().color(color_cycle.next().unwrap()));
        //.line(color_cycle.next());
    plot.add_trace(trace);
    let out_path_buffer = data.out_dir.with_file_name("./out.html");
    let out_path = out_path_buffer.as_path();
    println!("out path: {:?}", out_path);
    plot.write_html(out_path);
    use std::env;
    let path = env::current_dir().unwrap();
    println!("The current directory is {}", path.display());
    serve_once(out_path).unwrap();
}
