//! Data acquisition — search and download from remote sources

pub mod cmr;
pub mod stac;
pub mod download;

pub use cmr::{search_cmr, CmrGranule, search_pace_l1b, search_pace_scene};
pub use stac::StacClient;
pub use download::download_file;
