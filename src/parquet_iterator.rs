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

use crate::types::NormalizedTimeTag;

struct TimeTagBatchIterator {
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

pub struct ParquetStat {
    pub channels: Vec<u16>,
    pub schema: SchemaRef,
}

pub fn parquet_stat(file: &File) ->  ParquetStat {
    use std::collections::HashSet;
    let mut builder = ParquetRecordBatchReaderBuilder::try_new(file.try_clone().unwrap()).unwrap();
    let schema = builder.schema();

    let mut channel_ids = HashSet::<u16>::new();
    let mut cnt = 0;
    let max_cnt = 100;
    for tag in parquet_to_time_tag_iter(file.try_clone().unwrap()) {
        channel_ids.insert(tag.channel_id);
        cnt += 1;
        if cnt >= max_cnt {
            break;
        }
    }
    let mut channels = channel_ids.into_iter().collect::<Vec<u16>>();
    channels.sort();
    return ParquetStat {
        channels: channels,
        schema: schema.clone(),
    };
}

pub fn log_parquet_stat(file: &File) {
    let stat = parquet_stat(file);
    println!("parquet schema {:?}", stat.schema);
    println!("channel_id found in the first 100 tags: {:?}", stat.channels);
}

pub fn parquet_to_time_tag_iter(file: File) -> impl Iterator<Item = NormalizedTimeTag> {
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
