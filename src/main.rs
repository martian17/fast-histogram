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

pub struct TimeTagIter {
    reader: ParquetRecordBatchReader,
    col_channel: Vec<u16>,
    col_time_tag: Vec<u64>,
    index: usize,
}

impl TimeTagIter {
    fn next_columns(reader: &mut ParquetRecordBatchReader) -> Option<(Vec<u16>, Vec<u64>)> {
        match reader.next() {
            Some(batch_result) => {
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
                return Some((col_channel_vec, col_time_tag_vec));
            }
            None => {
                return None;
            }
        }
    }
}

impl Iterator for TimeTagIter {
    type Item = NormalizedTimeTag;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.col_channel.len() {
            match TimeTagIter::next_columns(&mut self.reader) {
                Some((col_channel, col_time_tag)) => {
                    self.col_channel = col_channel;
                    self.col_time_tag = col_time_tag;
                    self.index = 0;
                },
                None => {
                    // end of index
                    return None;
                },
            }
        }
        let result = Some(NormalizedTimeTag {
            channel_id: self.col_channel[self.index],
            time_tag_ps: self.col_time_tag[self.index],
        });
        self.index += 1;
        return result;
    }
}

pub fn parquet_to_time_tag_iter(file: File) -> TimeTagIter {
    let mut builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let mut reader: ParquetRecordBatchReader = builder.build().unwrap();
    match TimeTagIter::next_columns(&mut reader) {
        Some((col_channel, col_time_tag)) => {
            return TimeTagIter {
                reader: reader,
                col_channel: col_channel,
                col_time_tag: col_time_tag,
                index: 0,
            };
        },
        None => {
            // this iterator will always return None
            return TimeTagIter {
                reader: reader,
                col_channel: Vec::new(),
                col_time_tag: Vec::new(),
                index: 0,
            };
        },
    }
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


    //let file_reader = SerializedFileReader::new(file).unwrap();
    
    // let mut builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    // let mut reader: ParquetRecordBatchReader = builder.build().unwrap();
    // let schema: &SchemaRef = &reader.schema();
    // println!("schema: {:?}", schema);
    // for batch_option in reader {
    //     let batch = batch_option.expect("Expected batch");
    //     let col_channel_ref = batch.column_by_name("channel").expect("Expected channel ArrayRef");
    //     let col_channel: &UInt16Array = col_channel_ref.as_any().downcast_ref().unwrap();
    //     let col_time_tag_ref = batch.column_by_name("time_tag").expect("Expected time_tag ArrayRef");
    //     let col_time_tag: &UInt64Array = col_time_tag_ref.as_any().downcast_ref().unwrap();
    //     for i in 0..batch.num_rows() {
    //         let channel: u16 = col_channel.value(i);
    //         let time_tag: u64 = col_time_tag.value(i);
    //         let tag = NormalizedTimeTag {
    //             channel_id: channel,
    //             time_tag_ps: time_tag,
    //         };
    //         println!("tag: {:?}", tag);
    //         break;
    //     }
    // }
}
