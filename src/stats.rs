use std::fs::File;
use owo_colors::OwoColorize;



use crate::parquet_iterator::{
    parquet_to_time_tag_iter,
};
// use crate::types::NormalizedTimeTag;
use crate::signal_processing::PulseTracker;


pub fn display_parquet_stats(file: &File){
    use std::collections::{HashSet, HashMap};

    struct ChannelStat{
        t_min: u64,// = f64::MAX
        t_max: u64,
        count: u64,
    }

    let mut channel_ids = HashSet::<u16>::new();
    let mut channel_pulse_trackers = HashMap::<u16, PulseTracker>::new();
    let mut channel_stats = HashMap::<u16, ChannelStat>::new();
    for tag in parquet_to_time_tag_iter(file) {
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
    println!("{}", format!("\nAll channels scanned: {:?}", channels).bold());
    for id in channels {
        let pulse_tracker = channel_pulse_trackers.get(&id).unwrap();
        let stat = channel_stats.get(&id).unwrap();
        println!("{}, {} Counts, Start: {}, End: {}", 
            format!("Channel {}", id).bold(), 
            stat.count.yellow(), 
            format!("{}ps", stat.t_min).yellow(), 
            format!("{}ps", stat.t_max).yellow());
        println!("Derived pulse interval: {}", format!("{}ps", pulse_tracker.interval).yellow());
        println!("Derived pulse phase: {}", format!("{}ps", pulse_tracker.phase % pulse_tracker.interval).yellow());
    }
}
