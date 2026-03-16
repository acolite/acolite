//! Data acquisition — search and download from remote sources

pub mod cmr;
pub mod download;
pub mod stac;

pub use cmr::{search_cmr, search_pace_l1b, search_pace_scene, CmrGranule};
pub use download::download_file;
pub use stac::StacClient;
