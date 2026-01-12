use rustfft::{FftPlanner, num_complex::Complex};

pub fn fft_convolve(buf1: Vec<f64>, buf2: Vec<f64>) -> Vec<f64> {
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

    let reals: Vec<f64> = complex_buf3.iter().map(|complex| complex.re).collect();
    // let mut imags: Vec<f64> = complex_buf2.iter().map(|complex| complex.im).collect();
    return reals;
}

fn index_center_offset(idx: usize, len: usize) -> i64 {
    let idx = idx as i64;
    let len = len as i64;
    // the second half will be negative
    return (idx + len/2) % len - len/2;
}

pub fn apply_pdf<F>(buffer: &mut Vec<f64>, pdf: F)
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
pub struct PulseTracker {
    pub interval: f64, // picoseconds
    pub phase: f64,    // picoseconds
    pub initialized: bool,
}

impl PulseTracker {
    pub fn update(&mut self, t: f64){
        if !self.initialized {
            self.phase = t;
            self.initialized = true;
        }else{
            let intervals = (t-self.phase)/self.interval;
            let rounded_intervals = intervals.round();
            let derived_interval = (t-self.phase)/rounded_intervals;
            let derived_phase = t - rounded_intervals * self.interval;

            // perform gained adjustment
            self.interval = self.interval * 0.1 + derived_interval * 0.9;
            self.phase = self.phase * 0.1 + derived_phase * 0.9;
        }
    }
    pub fn new() -> Self {
        return Self {
            interval: 1000.0, // picoseconds
            phase: 0.0,       // picoseconds
            initialized: false,
        };
    }
}
