use std::env;
use std::fs::File;


mod parquet_iterator;
mod types;

use crate::parquet_iterator::parquet_to_time_tag_iter;
use crate::types::NormalizedTimeTag;


fn main () {
    let args: Vec<String> = env::args().collect();
    // println!("args: {:?}", args);
    println!("reading: {:?}", args[1]);


    let file = File::open(args[1].clone()).unwrap();
    let mut cnt: i32 = 0;
    for tag in parquet_to_time_tag_iter(file) {
        cnt += 1;
        if cnt > 100 {
            break;
        }
        println!("tag: {:?}", tag);
    }
}
