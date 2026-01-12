use std::fs::File;
use std::path::Path;
use owo_colors::OwoColorize;

use crate::parquet_iterator::{
    parquet_to_time_tag_iter,
    parquet_stat,
};
use crate::types::NormalizedTimeTag;
use crate::histogram::{HistogramData, display_histogram};
use crate::utils::log_peak;


// pub fn handle_histogram(args: Vec<String>) {
pub fn handle_histogram(path: &Path) {
    let file = &File::open(path).unwrap();

    let stat = parquet_stat(file);
    let channel1_id = stat.channels[1];
    let channel2_id = stat.channels[0];
    println!("Using channel {} and {}", channel1_id.yellow(), channel2_id.yellow());

    let tags = parquet_to_time_tag_iter(file);

    let bucket_size: usize = 200_0;// 5 x 200n = 1 micro seconds of range
    let bin_size: usize = 5;// 5 picoseconds
    let mut bucket: Vec<i32> = vec![0; bucket_size];
    
    let mut syn_start: u64 = 0;
    
    let mut stored_tags: Vec<NormalizedTimeTag> = Vec::new();
    for tag in tags {
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
            let diff = (syn_start as i64 - tag.time_tag_ps as i64)/(bin_size as i64);// tag is always smaller here
            if diff < (bucket_size/2) as i64 {
                bucket[((bucket_size - 1) as i64 - diff) as usize] += 1;
            }
        }
    }
    log_peak(&bucket, bin_size);
    
    display_histogram(HistogramData{
        values: bucket.iter().map(|v| *v as f64).collect(),
        timespan: (bucket_size * bin_size) as u64,
        out_dir: path,
    });
}
