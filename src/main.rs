use std::env;
use std::fs::File;
use std::path::Path;

mod commands;
mod parquet_iterator;
mod types;
mod histogram;
mod utils;
mod stats;
mod signal_processing;

use crate::stats::display_parquet_stats;

fn main () {
    let args: Vec<String> = env::args().collect();
    println!("reading: {:?}", args[2]);
    let path = Path::new(&args[2]);
    let file = File::open(path).unwrap();

    display_parquet_stats(&file);

    if args[1] == "--stats".to_string() {
        std::process::exit(1);
    }

    if args[1] == "--histogram".to_string() {
        println!("mode: histogram");
        commands::histogram::handle_histogram(path);
    }

    if args[1] == "--uniform".to_string() {
        println!("mode: uniform");
        commands::experimental_fft::handle_uniform(path);
    }

    if args[1] == "--pulsed".to_string() {
        println!("mode: pulsed");
        commands::experimental_fft::handle_pulsed(path);
    }
}
