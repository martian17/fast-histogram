use std::env;
use std::fs::File;
use std::ops::Range;
use rustfft::{FftPlanner, num_complex::Complex};

use macroquad::prelude::*;


mod parquet_iterator;
mod types;

use crate::parquet_iterator::{
    parquet_to_time_tag_iter,
    log_parquet_stat,
    parquet_stat,
};
use crate::types::NormalizedTimeTag;

fn window_conf() -> Conf {
    Conf {
        window_title: "Fast-Histogram Rust".to_owned(),
        window_width: 1024,
        window_height: 700,
        window_resizable: true,
        ..Default::default()
    }
}

struct MinMaxInfo<T> {
    max: T,
    min: T,
    max_idx: usize,
    min_idx: usize,
}

fn max_min_and_index<T>(values: &[T]) -> MinMaxInfo<T>
where
    T: PartialOrd + Copy,
{
    let mut min = values[0];
    let mut max = values[0];
    let mut min_idx = 0;
    let mut max_idx = 0;

    for (i, &value) in values.iter().enumerate().skip(1) {
        if value > max {
            max = value;
            max_idx = i;
        }
        if value < min {
            min = value;
            min_idx = i;
        }
    }

    MinMaxInfo {
        max,
        min,
        max_idx,
        min_idx,
    }
}


async fn display_histogram (values: Vec<f64>) {
    // log_peak(&values);
    let min_max = max_min_and_index(&values);

    let len = values.len();
    let bin_size = len/1000;

    let view_size = values.len()/bin_size;
    let mut coalesced: Vec<f64> = vec![0.0; view_size];
    for i in 0..values.len() {
        let idx = i / bin_size;
        if idx >= view_size {
            break;
        }
        coalesced[idx] += values[i];
    }
    let view_min_max = max_min_and_index(&coalesced);
    println!("max: {}", view_min_max.max);
    println!("min: {}", view_min_max.min);
    loop {
        let min = view_min_max.min;
        let max = view_min_max.max;
        clear_background(WHITE);
        let width = screen_width() as usize;
        let height = screen_height() as usize;
        for i in 0..(view_size-1) {
            let idx0 = (i+view_size+view_size/2)%view_size;
            let idx1 = (i+1+view_size+view_size/2)%view_size;
            let y0 = (height as f64) - (coalesced[idx0] - min)/(max - min) * (height as f64);
            let y1 = (height as f64) - (coalesced[idx1] - min)/(max - min) * (height as f64);
            draw_line(i as f32, y0 as f32, (i+1) as f32, y1 as f32, 1.0, BLACK);
        }
        next_frame().await;
    }
}

fn normalize_complex(c: Complex<f64>) -> Complex<f64> {
    let norm = c.norm();
    if norm < 0.01 {
        return Complex {re: 0.0, im: 0.0};
    }
    return c / norm;
}

fn fft_convolve(buf1: Vec<f64>, buf2: Vec<f64>) -> Vec<f64> {
    let len = buf1.len();
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(len);
    let ifft = planner.plan_fft_inverse(len);

    let mut complex_buf1: Vec<Complex<f64>> = buf1.iter().map(|value| Complex{re: value.clone(), im: 0.0}).collect();
    let mut complex_buf2: Vec<Complex<f64>> = buf2.iter().map(|value| Complex{re: value.clone(), im: 0.0}).collect();

    fft.process(&mut complex_buf1);
    fft.process(&mut complex_buf2);

    // calculate dot between the two
    let mut complex_buf3 = vec![Complex{re: 0.0, im: 0.0}; len];
    for i in 0..len {
        complex_buf3[i] = (complex_buf1[i] * complex_buf2[i].conj()) / (len as f64);
    }

    // calculate ifft on complex_buf3 to get the convolution
    ifft.process(&mut complex_buf3);

    let mut reals: Vec<f64> = complex_buf3.iter().map(|complex| complex.re).collect();
    // let mut imags: Vec<f64> = complex_buf2.iter().map(|complex| complex.im).collect();
    return reals;
}

fn index_center_offset(idx: usize, len: usize) -> i64 {
    let idx = idx as i64;
    let len = len as i64;
    // the second half will be negative
    return (idx + len/2) % len - len/2;
}

fn apply_pdf<F>(buffer: &mut Vec<f64>, pdf: F)
where
    F: Fn(f64) -> f64,
{
    let len = buffer.len();
    let mut pdf_buffer = vec![0.0; len];

    for i in 0..len {
        pdf_buffer[i] = pdf(index_center_offset(i, len) as f64);
    }
    // TODO: see if there is any copy with .to_vec()
    let result = fft_convolve(buffer.to_vec(), pdf_buffer);
    for i in 0..len {
        buffer[i] = result[i];
    }
}

// filtered tracker to track the pulse phase and interval
#[derive(Debug)]
struct PulseTracker {
    interval: f64, // picoseconds
    phase: f64,    // picoseconds
    initialized: bool,
}

impl PulseTracker {
    fn update(&mut self, t: f64){
        if !self.initialized {
            self.phase = t;
            self.initialized = true;
        }else{
            let intervals = (t-self.phase)/self.interval;
            let rounded_intervals = intervals.round();
            let derived_interval = (t-self.phase)/rounded_intervals;
            let derived_phase = t - rounded_intervals * self.interval;

            // perform gained adjustment
            self.interval = (self.interval * 0.1 + derived_interval * 0.9);
            self.phase = (self.phase * 0.1 + derived_phase * 0.9);
        }
    }
    fn new() -> Self {
        return Self {
            interval: 1000.0, // picoseconds
            phase: 0.0,       // picoseconds
            initialized: false,
        };
    }
}

fn display_parquet_stats(args: Vec<String>){
    if args.len() < 3 {
        std::process::exit(1);
    }
    let file = File::open(args[2].clone()).unwrap();

    use std::collections::{HashSet, HashMap};

    struct ChannelStat{
        t_min: u64,// = f64::MAX
        t_max: u64,
        count: u64,
    };

    let mut channel_ids = HashSet::<u16>::new();
    let mut channel_pulse_trackers = HashMap::<u16, PulseTracker>::new();
    let mut channel_stats = HashMap::<u16, ChannelStat>::new();
    for tag in parquet_to_time_tag_iter(file.try_clone().unwrap()) {
        if !channel_ids.contains(&tag.channel_id) {
            channel_ids.insert(tag.channel_id);
            channel_pulse_trackers.insert(tag.channel_id, PulseTracker::new());
            channel_stats.insert(tag.channel_id, ChannelStat{
                t_min: u64::MAX,
                t_max: u64::MIN,
                count: 0,
            });
        }

        channel_pulse_trackers.get_mut(&tag.channel_id).unwrap().update(tag.time_tag_ps as f64);
        let stat = channel_stats.get_mut(&tag.channel_id).unwrap();
        stat.count += 1;
        if tag.time_tag_ps < stat.t_min {
            stat.t_min = tag.time_tag_ps;
        }
        if tag.time_tag_ps > stat.t_max {
            stat.t_max = tag.time_tag_ps;
        }
    }
    let mut channels = channel_ids.into_iter().collect::<Vec<u16>>();
    channels.sort();
    use ansi_term::Colour::{Red, Yellow, Blue};
    use ansi_term::Style;
    println!("{}", Style::new().bold().paint(format!("\nAll channels scanned: {:?}", channels)));
    for id in channels {
        let pulse_tracker = channel_pulse_trackers.get(&id).unwrap();
        let stat = channel_stats.get(&id).unwrap();
        println!("{}, {} Counts, Start: {}, End: {}", 
            Style::new().bold().paint(format!("Channel {}", id)), 
            Yellow.paint(format!("{}", stat.count)), 
            Yellow.paint(format!("{}ps", stat.t_min)), 
            Yellow.paint(format!("{}ps", stat.t_max)));
        println!("Derived pulse interval: {}", Yellow.paint(format!("{}ps", pulse_tracker.interval)));
        println!("Derived pulse phase: {}", Yellow.paint(format!("{}ps", pulse_tracker.phase % pulse_tracker.interval)));
    }
}


#[macroquad::main(window_conf)]
async fn main () {
    let args: Vec<String> = env::args().collect();

    if args[1] == "--stats".to_string() {
        println!("mode: stats");
        display_parquet_stats(args.clone());
        std::process::exit(1);
    }

    if args[1] == "--histogram".to_string() {
        println!("mode: histogram");
        handle_histogram(args.clone()).await;
    }

    if args[1] == "--uniform".to_string() {
        println!("mode: uniform");
        handle_uniform(args.clone()).await;
    }

    if args[1] == "--pulsed".to_string() {
        println!("mode: pulsed");
        handle_pulsed(args.clone()).await;
    }
}

async fn handle_uniform(args: Vec<String>) {
    let file = File::open(args[2].clone()).unwrap();
    println!("reading: {:?}", args[2]);
    log_parquet_stat(&file);
    let stat = parquet_stat(&file);
    let channel0_id = stat.channels[0];
    let channel1_id = stat.channels[1];
    println!("Using channel {} and {}", channel0_id, channel1_id);

    let bucket_size: usize = 123_456;//100_000;

    let mut bucket0: Vec<f64> = vec![0.0; bucket_size];
    let mut bucket1: Vec<f64> = vec![0.0; bucket_size];

    let mut kernel_pdf_function = |x: f64| -> f64 {
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

    display_histogram(result).await;
}


async fn handle_pulsed(args: Vec<String>) {
    let file = File::open(args[2].clone()).unwrap();
    println!("reading: {:?}", args[2]);
    log_parquet_stat(&file);
    let stat = parquet_stat(&file);
    let channel0_id = stat.channels[0];
    let channel1_id = stat.channels[1];
    println!("Using channel {} and {}", channel0_id, channel1_id);

    let bucket_size: usize = 1000;

    let mut bucket0: Vec<f64> = vec![0.0; bucket_size];
    let mut bucket1: Vec<f64> = vec![0.0; bucket_size];

    let mut kernel_pdf_function = |x: f64| -> f64 {
        use std::f64::consts::PI;
        let sigma = 20.0;
        let coefficient = 1.0 / (sigma * (2.0 * PI).sqrt());
        let exponent = -(x*x) / (2.0 * sigma*sigma);
        return coefficient * exponent.exp();
    };


    let mut tracker0 = PulseTracker::new();
    let mut tracker1 = PulseTracker::new();
    for tag in parquet_to_time_tag_iter(file.try_clone().unwrap()) {
        if tag.channel_id == channel0_id {
            tracker0.update(tag.time_tag_ps as f64);
        } else if tag.channel_id == channel1_id {
            tracker1.update(tag.time_tag_ps as f64);
        }
    }
    println!("tracker result: {:?}", tracker0);

    let sample_interval = (tracker0.interval + tracker1.interval)/2.0;
    let sample_offset = tracker0.phase;

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

    display_histogram(result).await;
}


fn log_peak (values: &Vec<i32>, bin_size: usize) {
    use ansi_term::Colour::{Red, Yellow, Blue};
    let min_max = max_min_and_index(&values);
    let len: i32 = values.len() as i32;
    let max_idx: i32 = min_max.max_idx as i32;
    let peak_idx = (max_idx + len/2).rem_euclid(len) - len/2;
    //+ len)%len - len/2;
    println!("{}", Yellow.bold().paint(format!("Peak found at {}ps", peak_idx * bin_size as i32)));
}


async fn handle_histogram(args: Vec<String>) {
    let file = File::open(args[2].clone()).unwrap();
    println!("reading: {:?}", args[2]);
    log_parquet_stat(&file);
    let stat = parquet_stat(&file);
    let channel1_id = stat.channels[1];
    let channel2_id = stat.channels[0];
    println!("Using channel {} and {}", channel1_id, channel2_id);

    let bucket_size: usize = 200_0;// 5 x 200n = 1 micro seconds of range
    let bin_size: usize = 5;// 5 picoseconds
    let mut bucket: Vec<i32> = vec![0; bucket_size];
    
    let mut syn_start: u64 = 0;
    
    let mut stored_tags: Vec<NormalizedTimeTag> = Vec::new();
    for tag in parquet_to_time_tag_iter(file.try_clone().unwrap()) {
        stored_tags.push(tag.clone());
        if tag.channel_id == channel1_id {
            syn_start = tag.time_tag_ps;
        } else if tag.channel_id == channel2_id {
            if syn_start == 0 {
                continue;
            }
            let diff = (tag.time_tag_ps as i64 - syn_start as i64)/(bin_size as i64);
            if diff < (bucket_size/2) as i64 {
                bucket[diff as usize] += 1;
            }
            // channel2.push(tag.time_tag_ps);
        }
    }
    syn_start = u64::MAX;
    for tag in stored_tags.iter().rev() {
        if tag.channel_id == channel1_id {
            syn_start = tag.time_tag_ps;
        } else if tag.channel_id == channel2_id {
            if syn_start == u64::MAX {
                continue;
            }
            let mut diff = (syn_start as i64 - tag.time_tag_ps as i64)/(bin_size as i64);// tag is always smaller here
            if diff < (bucket_size/2) as i64 {
                bucket[((bucket_size - 1) as i64 - diff) as usize] += 1;
            }
            // channel2.push(tag.time_tag_ps);
        }
    }
    display_parquet_stats(args.clone());
    log_peak(&bucket, bin_size);
    
    display_histogram(bucket.iter().map(|v| *v as f64).collect()).await;
}


// async fn handle_coincidence(args: Vec<String>) {
//     let file = File::open(args[2].clone()).unwrap();
//     println!("reading: {:?}", args[2]);
//     log_parquet_stat(&file);
//     let stat = parquet_stat(&file);
//     let channel1_id = stat.channels[0];
//     let channel2_id = stat.channels[1];
//     println!("Using channel {} and {}", channel1_id, channel2_id);
// 
//     let bucket_size: usize = 200_0;// 5 x 200n = 1 micro seconds of range
//     let bin_size: usize = 5;// 5 picoseconds
//     let mut bucket: Vec<i32> = vec![0; bucket_size];
//     
//     let mut syn_start: u64 = 0;
//     
//     let mut stored_tags: Vec<NormalizedTimeTag> = Vec::new();
//     for tag in parquet_to_time_tag_iter(file.try_clone().unwrap()) {
//         stored_tags.push(tag.clone());
//         if tag.channel_id == channel1_id {
//             syn_start = tag.time_tag_ps;
//         } else if tag.channel_id == channel2_id {
//             if syn_start == 0 {
//                 continue;
//             }
//             let diff = (tag.time_tag_ps as i64 - syn_start as i64)/(bin_size as i64);
//             if diff < (bucket_size/2) as i64 {
//                 bucket[diff as usize] += 1;
//             }
//             // channel2.push(tag.time_tag_ps);
//         }
//     }
//     syn_start = u64::MAX;
//     for tag in stored_tags.iter().rev() {
//         if tag.channel_id == channel1_id {
//             syn_start = tag.time_tag_ps;
//         } else if tag.channel_id == channel2_id {
//             if syn_start == u64::MAX {
//                 continue;
//             }
//             let mut diff = (syn_start as i64 - tag.time_tag_ps as i64)/(bin_size as i64);// tag is always smaller here
//             if diff < (bucket_size/2) as i64 {
//                 bucket[((bucket_size - 1) as i64 - diff) as usize] += 1;
//             }
//             // channel2.push(tag.time_tag_ps);
//         }
//     }
//     log_peak(&bucket, bin_size);
//     
//     display_histogram(bucket.iter().map(|v| *v as f64).collect()).await;
// }
