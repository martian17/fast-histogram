use std::fs::File;
use std::path::Path;
use owo_colors::OwoColorize;



use crate::parquet_iterator::{
    parquet_to_time_tag_iter,
    log_parquet_stat,
    parquet_stat,
};
// use crate::types::NormalizedTimeTag;
use crate::histogram::{HistogramData, display_histogram};
use crate::signal_processing::{fft_convolve, apply_pdf, PulseTracker};


pub fn handle_pulsed(path: &Path) {
    let file = &File::open(path).unwrap();
    log_parquet_stat(file);
    let stat = parquet_stat(file);
    let channel0_id = stat.channels[0];
    let channel1_id = stat.channels[1];
    println!("Using channel {} and {}", channel0_id.yellow(), channel1_id.yellow());

    let bucket_size: usize = 1000;

    let mut bucket0: Vec<f64> = vec![0.0; bucket_size];
    let mut bucket1: Vec<f64> = vec![0.0; bucket_size];

    let _kernel_pdf_function = |x: f64| -> f64 {
        use std::f64::consts::PI;
        let sigma = 20.0;
        let coefficient = 1.0 / (sigma * (2.0 * PI).sqrt());
        let exponent = -(x*x) / (2.0 * sigma*sigma);
        return coefficient * exponent.exp();
    };


    let mut tracker0 = PulseTracker::new();
    let mut tracker1 = PulseTracker::new();
    for tag in parquet_to_time_tag_iter(file) {
        if tag.channel_id == channel0_id {
            tracker0.update(tag.time_tag_ps as f64);
        } else if tag.channel_id == channel1_id {
            tracker1.update(tag.time_tag_ps as f64);
        }
    }
    println!("tracker result: {:?}", tracker0.yellow());

    // let sample_interval = (tracker0.interval + tracker1.interval)/2.0;
    // let sample_offset = tracker0.phase;

    for tag in parquet_to_time_tag_iter(file) {
        let sample_interval = tracker0.interval * 1.0;
        let r = (tag.time_tag_ps as f64 - tracker0.phase).rem_euclid(sample_interval)/sample_interval;
        let mut index = (r * bucket_size as f64).floor() as usize;
        if index == bucket_size {
            index = bucket_size - 1;
        }
        if tag.channel_id == channel0_id {
            bucket0[index] += 1.0;
        } else if tag.channel_id == channel1_id {
            bucket1[index] += 1.0;
        } else {
            panic!("what???");
        }
    }
    // apply_pdf(&mut bucket0, kernel_pdf_function);
    // apply_pdf(&mut bucket1, kernel_pdf_function);

    let result = fft_convolve(bucket0.clone(), bucket1.clone());

    display_histogram(HistogramData{
        values: result,
        timespan: 1000,// hard code it for now
        out_dir: path,
    });
}

pub fn handle_uniform(path: &Path) {
    let file = &File::open(path).unwrap();
    log_parquet_stat(file);
    let stat = parquet_stat(&file);
    let channel0_id = stat.channels[0];
    let channel1_id = stat.channels[1];
    println!("Using channel {} and {}", channel0_id.yellow(), channel1_id.yellow());

    let bucket_size: usize = 123_456;//100_000;

    let mut bucket0: Vec<f64> = vec![0.0; bucket_size];
    let mut bucket1: Vec<f64> = vec![0.0; bucket_size];

    let kernel_pdf_function = |x: f64| -> f64 {
        use std::f64::consts::PI;
        let sigma = 2000.0;
        let coefficient = 1.0 / (sigma * (2.0 * PI).sqrt());
        let exponent = -(x*x) / (2.0 * sigma*sigma);
        return coefficient * exponent.exp();
    };

    for tag in parquet_to_time_tag_iter(file) {
        let index = (tag.time_tag_ps as usize).rem_euclid(bucket_size);
        if tag.channel_id == channel0_id {
            bucket0[index] += 1.0;
        } else if tag.channel_id == channel1_id {
            bucket1[index] += 1.0;
        }
    }
    apply_pdf(&mut bucket0, kernel_pdf_function);
    apply_pdf(&mut bucket1, kernel_pdf_function);

    let result = fft_convolve(bucket0.clone(), bucket1.clone());

    display_histogram(HistogramData{
        values: result,
        timespan: bucket_size as u64,
        out_dir: path,
    });
}
