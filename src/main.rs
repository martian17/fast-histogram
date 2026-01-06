use parquet::arrow::arrow_reader::{
    ParquetRecordBatchReaderBuilder,
    ParquetRecordBatchReader
};
use arrow::array::{UInt16Array, UInt64Array};
use arrow::array::RecordBatchReader;
use arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use arrow::record_batch::RecordBatch;
use std::fs::File;
use std::path::Path;
use std::env;

#[derive(Debug)]
pub struct NormalizedTimeTag {
    pub channel_id: u16,

    /// The time tag, in picoseconds, counting up from the start of the measurement.
    pub time_tag_ps: u64,
}

pub struct TimeTagBatchIterator {
    col_channel: Vec<u16>,
    col_time_tag: Vec<u64>,
    index: usize,
    size: usize,
}

impl Iterator for TimeTagBatchIterator {
    type Item = NormalizedTimeTag;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.size {
            return None;
        }
        let result = Some(NormalizedTimeTag {
            channel_id: self.col_channel[self.index],
            time_tag_ps: self.col_time_tag[self.index],
        });
        self.index += 1;
        return result;
    }
}

fn parquet_to_time_tag_iter(file: File) -> impl Iterator<Item = NormalizedTimeTag> {
    let mut builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let mut reader: ParquetRecordBatchReader = builder.build().unwrap();
    return reader.flat_map(|batch_result| {
        let batch = batch_result.expect("Expected batch");
        let col_channel_ref = batch.column_by_name("channel").expect("Expected channel ArrayRef");
        let col_channel: &UInt16Array = col_channel_ref.as_any().downcast_ref().unwrap();
        let col_channel_vec: Vec<u16> = col_channel.iter().map(|option| match option {
            Some(value) => value,
            None => panic!("None found in channel column field"),
        }).collect();
        let col_time_tag_ref = batch.column_by_name("time_tag").expect("Expected time_tag ArrayRef");
        let col_time_tag: &UInt64Array = col_time_tag_ref.as_any().downcast_ref().unwrap();
        let col_time_tag_vec: Vec<u64> = col_time_tag.iter().map(|option| match option {
            Some(value) => value,
            None => panic!("None found in channel column field"),
        }).collect();
        return TimeTagBatchIterator {
            col_channel: col_channel_vec,
            col_time_tag: col_time_tag_vec,
            index: 0,
            size: batch.num_rows(),
        }
    })
}

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
