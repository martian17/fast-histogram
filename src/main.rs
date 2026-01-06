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

struct MinMaxInfo {
    max: f64,
    min: f64,
    max_idx: usize,
    min_idx: usize,
}

fn max_min_and_index (values: &Vec<f64>) -> MinMaxInfo {
    let mut min: f64 = f64::INFINITY;
    let mut max: f64 = -f64::INFINITY;
    let mut min_idx: usize = 0;
    let mut max_idx: usize = 0;
    for i in 0..values.len() {
        if values[i] > max {
            max = values[i];
            max_idx = i;
        }
        if values[i] < min {
            min = values[i];
            min_idx = i;
        }
    }
    return MinMaxInfo {
        max: max,
        min: min,
        max_idx: max_idx,
        min_idx: min_idx,
    }
}

fn log_peak (values: &Vec<f64>) {
    let min_max = max_min_and_index(&values);
    let len: i32 = values.len() as i32;
    let max_idx: i32 = min_max.max_idx as i32;
    let peak_location = (max_idx + len/2 + len)%len - len/2;
    println!("peak at {} ps", peak_location);
}

async fn display_histogram (values: Vec<f64>, bin_size: usize) {
    log_peak(&values);
    let min_max = max_min_and_index(&values);

    let view_size = values.len()/bin_size;
    let mut coalesced: Vec<f64> = vec![0.0; view_size];
    for i in 0..values.len() {
        let idx = i / bin_size;
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
            let y0 = (height as f64) - (coalesced[i+0] - min)/(max - min) * (height as f64);
            let y1 = (height as f64) - (coalesced[i+1] - min)/(max - min) * (height as f64);
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

    // normalize doesn't work
    // // normalize
    // for i in 0..bucket_size {
    //     complex_buf0[i] = normalize_complex(complex_buf0[i]);
    //     complex_buf1[i] = normalize_complex(complex_buf1[i]);
    //     println!("{}", complex_buf0[i]);
    //     // complex_buf1[i] = complex_buf1[i] / complex_buf1[i].norm();
    // }

    //let min_max = max_min_and_index(complex_buf0.iter().map(|complex| =>));
    //let min_max = max_min_and_index(&values);

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


#[macroquad::main(window_conf)]
async fn main () {
    let args: Vec<String> = env::args().collect();
    // println!("args: {:?}", args);
    println!("reading: {:?}", args[1]);

    let file = File::open(args[1].clone()).unwrap();
    log_parquet_stat(&file);
    let stat = parquet_stat(&file);
    let channel0_id = stat.channels[0];
    let channel1_id = stat.channels[1];
    println!("Using channel {} and {}", channel0_id, channel1_id);

    let bucket_size: usize = 10000;

    let mut bucket0: Vec<f64> = vec![0.0; bucket_size];
    let mut bucket1: Vec<f64> = vec![0.0; bucket_size];

    // let mut kernel_range: Range<i64> = -160..160;
    // let mut kernel_pdf_function = |x: f64| -> f64 {
    //     if x < -15.0 || x > 15.0 {
    //         return 0.0;
    //     }
    //     return 1.0/30.0;
    // };
    let mut kernel_pdf_function = |x: f64| -> f64 {
        use std::f64::consts::PI;
        let sigma = 200.0;
        let coefficient = 1.0 / (sigma * (2.0 * PI).sqrt());
        let exponent = -(x*x) / (2.0 * sigma*sigma);
        return coefficient * exponent.exp();
    };

    // let mut cnt: i32 = 0;
    for tag in parquet_to_time_tag_iter(file) {
        // cnt += 1;
        // if cnt < 1000 {
        //     println!("tag: {:?}", tag);
        // }
        let index = tag.time_tag_ps as usize % bucket_size;
        if tag.channel_id == channel0_id {
            bucket0[index] += 1.0;
        } else if tag.channel_id == channel1_id {
            bucket1[index] += 1.0;
        }
    }
    apply_pdf(&mut bucket0, kernel_pdf_function);
    apply_pdf(&mut bucket1, kernel_pdf_function);

    let result = fft_convolve(bucket0.clone(), bucket1.clone());

    display_histogram(result, 4).await;
    // display_histogram(bucket0, 4).await;

}
