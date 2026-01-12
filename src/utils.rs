use owo_colors::OwoColorize;

// TODO: some fields are dead, but clean it up later
#[allow(dead_code)]
pub struct MinMaxInfo<T> {
    pub max: T,
    pub min: T,
    pub max_idx: usize,
    pub min_idx: usize,
}

pub fn min_max_and_index<T>(values: &[T]) -> MinMaxInfo<T>
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

pub fn log_peak(values: &Vec<i32>, bin_size: usize) {
    let min_max = min_max_and_index(&values);
    let len: i32 = values.len() as i32;
    let max_idx: i32 = min_max.max_idx as i32;
    let peak_idx = (max_idx + len/2).rem_euclid(len) - len/2;
    //+ len)%len - len/2;
    println!("{}", format!("Peak found at {}ps", peak_idx * bin_size as i32).yellow().bold());
}


use std::path::Path;
use anyhow;
use open;
use tiny_http;
use std::fs::File;
use std::str::FromStr;

pub fn serve_once(path: &Path) -> anyhow::Result<()> {
    // TODO: use try operator "?"
    let server = tiny_http::Server::http("127.0.0.1:7878").unwrap();
    let header = tiny_http::Header::from_str("Cache-Control: no-cache, no-store, must-revalidate").unwrap();
    let response = tiny_http::Response::from_file(File::open(path).unwrap()).with_header(header);

    open::that("http://127.0.0.1:7878").unwrap();

    // usually this goes inside a loop{}, but as we're only serving once, it's not necessary
    let request = server.recv().expect("Serve once failed");
    let _ = request.respond(response);

    return anyhow::Ok(());
}


// // TODO: This part is likely unstable. Use spawn to make it more robust
// // Got the code from: https://doc.rust-lang.org/book/ch21-01-single-threaded.html
// use std::{
//     fs,
//     io::{BufReader, prelude::*},
//     net::{TcpListener},
// };
// use std::path::Path;
// use anyhow;
// use open;
// pub fn serve_once<P: AsRef<Path>>(path: P) -> anyhow::Result<()> {
//     let contents = fs::read_to_string(path).unwrap();
//     let listener = TcpListener::bind("127.0.0.1:7878").unwrap();
// 
//     open::that("http://127.0.0.1:7878")?;
// 
//     let mut stream = listener.incoming().next().unwrap().unwrap();
// 
//     let buf_reader = BufReader::new(&stream);
//     let _http_request: Vec<_> = buf_reader
//         .lines()
//         .map(|result| result.unwrap())
//         .take_while(|line| !line.is_empty())
//         .collect();
// 
//     let status_line = "HTTP/1.1 200 OK";
//     let length = contents.len();
// 
//     let response =
//         format!("{status_line}\r\nContent-Length: {length}\r\n\r\n{contents}");
// 
//     stream.write_all(response.as_bytes()).unwrap();
//     return anyhow::Ok(());
// }
