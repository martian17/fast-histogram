# fast-histogram

Usage for pulsed entangled photon sources (such as SPDC, phase information only): 
```bash
cargo run --features pulsed -- xxx.parquet
```

Usage for uniform entangled photon sources (fast broad range coincidence detection): 
```bash
cargo run --features pulsed -- xxx.parquet
```

Note: Parquet files must contain "channel" and "time_tags" columns
