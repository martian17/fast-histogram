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
    if(norm < 0.01){
        return Complex {re: 0.0, im: 0.0};
    }
    return c / norm;
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

    
    let mut kernel_range: Range<i64> = -15..15;
    // let mut kernel_range: Range<i64> = -160..160;
    let mut kernel_pdf_function = |_x: f64| {
        // let mut x: f64 = x as f64;
        // may want to try out normal distribution in the future but flat for now
        let start = kernel_range.start as f64;
        let end = kernel_range.end as f64;
        return 1.0/(end - start);
    };

    // let mut cnt: i32 = 0;
    for tag in parquet_to_time_tag_iter(file) {
        // cnt += 1;
        // if cnt < 1000 {
        //     println!("tag: {:?}", tag);
        // }
        // let index = tag.time_tag_ps % bucket_size;
        if tag.channel_id == channel0_id {
            for offset in kernel_range.clone() {
                let index = ((bucket_size as i64 + tag.time_tag_ps as i64 + offset) as usize) % bucket_size;
                bucket0[index] += kernel_pdf_function(offset as f64);
            }
        } else if tag.channel_id == channel1_id {
            for offset in kernel_range.clone() {
                let index = ((bucket_size as i64 + tag.time_tag_ps as i64 + offset) as usize) % bucket_size;
                bucket1[index] += kernel_pdf_function(offset as f64);
            }
        }
    }

    // println!("dasdfsadf");
    // for i in 0..100 {
    //     println!("{}", bucket0[i]);
    // }
    {
        let mut planner = FftPlanner::<f64>::new();
        let fft = planner.plan_fft_forward(bucket_size);
        let ifft = planner.plan_fft_inverse(bucket_size);

        let mut buff0: Vec<Complex<f64>> = bucket0.iter().map(|value| Complex{re: value.clone(), im: 0.0}).collect();
        let mut buff1: Vec<Complex<f64>> = bucket1.iter().map(|value| Complex{re: value.clone(), im: 0.0}).collect();

        fft.process(&mut buff0);
        fft.process(&mut buff1);

        // normalize doesn't work
        // // normalize
        // for i in 0..bucket_size {
        //     buff0[i] = normalize_complex(buff0[i]);
        //     buff1[i] = normalize_complex(buff1[i]);
        //     println!("{}", buff0[i]);
        //     // buff1[i] = buff1[i] / buff1[i].norm();
        // }

        //let min_max = max_min_and_index(buff0.iter().map(|complex| =>));
        //let min_max = max_min_and_index(&values);

        // calculate dot between the two
        let mut buff2 = vec![Complex{re: 0.0, im: 0.0}; bucket_size];
        for i in 0..bucket_size {
            buff2[i] = (buff0[i] * buff1[i].conj()) / (bucket_size as f64);
        }

        // calculate ifft on buff2 to get the convolution
        ifft.process(&mut buff2);

        let mut reals: Vec<f64> = buff2.iter().map(|complex| complex.re).collect();
        let mut imags: Vec<f64> = buff2.iter().map(|complex| complex.im).collect();

        // display_histogram(imags, 4);
        // display_histogram(bucket0, 4).await;
        // display_histogram(bucket1, 4).await;
        display_histogram(reals, 4).await;
    }
    //display_histogram(bucket0, 4);
    // let histogram_bin_size = 4;// only when displaying, non-critical

}
